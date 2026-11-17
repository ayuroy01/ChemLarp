// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
// ChemLarp for macOS: the web app (WebAssembly engine and UI, bundled in Contents/Resources/web) in a native
// window. A small HTTP server bound to 127.0.0.1 serves the bundled files, so the page runs exactly as in a
// browser (workers, WebAssembly, file reading); nothing is fetched from the network.
//
//   ChemLarp.app/Contents/MacOS/ChemLarp --selftest   checks the server and exits (0 on success)
//   ChemLarp.app/Contents/MacOS/ChemLarp --smoke      opens the window, runs an engine call through the web
//                                                     app's worker (WebAssembly) and exits (0 on success)
import AppKit
import Network
import WebKit

// MARK: - Local server

final class LocalServer {
  let root: URL
  private var listener: NWListener?
  private let queue = DispatchQueue(label: "chemlarp.server", attributes: .concurrent)

  init(root: URL) { self.root = root.standardizedFileURL }

  /// Starts listening on 127.0.0.1 at a port the system chooses; `ready` gets the port, or nil on failure.
  func start(ready: @escaping (UInt16?) -> Void) {
    let params = NWParameters.tcp
    params.requiredLocalEndpoint = NWEndpoint.hostPort(host: "127.0.0.1", port: .any)
    guard let l = try? NWListener(using: params) else { return ready(nil) }
    listener = l
    var reported = false
    l.stateUpdateHandler = { state in
      switch state {
      case .ready:
        if !reported { reported = true; ready(l.port?.rawValue) }
      case .failed:
        if !reported { reported = true; ready(nil) }
      default: break
      }
    }
    l.newConnectionHandler = { [weak self] c in
      c.start(queue: self?.queue ?? .global())
      self?.receive(c, Data())
    }
    l.start(queue: queue)
  }

  private func receive(_ c: NWConnection, _ buffer: Data) {
    c.receive(minimumIncompleteLength: 1, maximumLength: 65536) { [weak self] data, _, done, error in
      guard let self else { return c.cancel() }
      var buf = buffer
      if let data { buf.append(data) }
      if let end = buf.range(of: Data("\r\n\r\n".utf8)) {
        self.respond(c, String(decoding: buf.subdata(in: 0..<end.lowerBound), as: UTF8.self))
      } else if done || error != nil || buf.count > 65536 {
        c.cancel()
      } else {
        self.receive(c, buf)
      }
    }
  }

  /// The bundled file for a request path, or nil (unknown file, or a path outside the bundle).
  func file(for target: String) -> URL? {
    var path = String(target.prefix { $0 != "?" && $0 != "#" })
    path = path.removingPercentEncoding ?? path
    if path.hasSuffix("/") { path += "index.html" }
    let url = root.appendingPathComponent(String(path.drop { $0 == "/" })).standardizedFileURL
    guard url.path.hasPrefix(root.path + "/") else { return nil }
    var isDir: ObjCBool = false
    guard FileManager.default.fileExists(atPath: url.path, isDirectory: &isDir), !isDir.boolValue else { return nil }
    return url
  }

  private func respond(_ c: NWConnection, _ head: String) {
    let parts = (head.components(separatedBy: "\r\n").first ?? "").split(separator: " ")
    let method = parts.first.map(String.init) ?? ""
    var status = "404 Not Found", type = "text/plain; charset=utf-8", body = Data("not found\n".utf8)
    if method != "GET" && method != "HEAD" {
      status = "405 Method Not Allowed"
      body = Data()
    } else if parts.count >= 2, let url = file(for: String(parts[1])), let data = try? Data(contentsOf: url) {
      status = "200 OK"
      type = LocalServer.mime(url.pathExtension)
      body = data
    }
    var out = Data(("HTTP/1.1 \(status)\r\nContent-Type: \(type)\r\nContent-Length: \(body.count)\r\n" +
                    "Cache-Control: no-cache\r\nX-Content-Type-Options: nosniff\r\nConnection: close\r\n\r\n").utf8)
    if method != "HEAD" { out.append(body) }
    c.send(content: out, completion: .contentProcessed { _ in c.cancel() })
  }

  static func mime(_ ext: String) -> String {
    switch ext.lowercased() {
    case "html": return "text/html; charset=utf-8"
    case "js", "mjs": return "text/javascript; charset=utf-8"
    case "css": return "text/css; charset=utf-8"
    case "wasm": return "application/wasm"
    case "json", "webmanifest": return "application/json"
    case "svg": return "image/svg+xml"
    case "png": return "image/png"
    default: return "text/plain; charset=utf-8"
    }
  }
}

// MARK: - Self-test (used by CI): serve the bundle and fetch from it

func selftest() -> Never {
  let root = Bundle.main.resourceURL!.appendingPathComponent("web")
  let server = LocalServer(root: root)
  let sem = DispatchSemaphore(value: 0)
  var port: UInt16?
  server.start { p in port = p; sem.signal() }
  _ = sem.wait(timeout: .now() + 10)
  guard let port else { print("selftest: the server did not start"); exit(1) }
  let base = "http://127.0.0.1:\(port)"
  let checks: [(String, Int, String?)] = [
    ("/", 200, "text/html"), ("/chemlarp.wasm", 200, "application/wasm"), ("/engine-worker.js", 200, "text/javascript"),
    ("/%2e%2e/Info.plist", 404, nil), ("/no-such-file", 404, nil),
  ]
  var ok = true
  for (path, want, type) in checks {
    var got = 0, ctype = ""
    let s = DispatchSemaphore(value: 0)
    URLSession.shared.dataTask(with: URL(string: base + path)!) { _, r, _ in
      if let h = r as? HTTPURLResponse { got = h.statusCode; ctype = h.value(forHTTPHeaderField: "Content-Type") ?? "" }
      s.signal()
    }.resume()
    _ = s.wait(timeout: .now() + 10)
    let pass = got == want && (type == nil || ctype.hasPrefix(type!))
    ok = ok && pass
    print("selftest: \(pass ? "ok " : "BAD") \(path) -> \(got) \(ctype)")
  }
  print(ok ? "selftest passed" : "selftest FAILED")
  exit(ok ? 0 : 1)
}

// MARK: - App

final class AppDelegate: NSObject, NSApplicationDelegate, WKUIDelegate, WKNavigationDelegate, WKDownloadDelegate {
  var window: NSWindow!
  var web: WKWebView!
  var server: LocalServer!
  var origin = ""
  let smoke = CommandLine.arguments.contains("--smoke")

  func applicationDidFinishLaunching(_ notification: Notification) {
    buildMenu()
    server = LocalServer(root: Bundle.main.resourceURL!.appendingPathComponent("web"))
    server.start { port in
      DispatchQueue.main.async {
        guard let port else {
          let a = NSAlert()
          a.messageText = "ChemLarp could not start its local server."
          a.runModal()
          NSApp.terminate(nil)
          return
        }
        self.open(port)
      }
    }
  }

  func open(_ port: UInt16) {
    origin = "http://127.0.0.1:\(port)"
    let config = WKWebViewConfiguration()
    web = WKWebView(frame: .zero, configuration: config)
    web.uiDelegate = self
    web.navigationDelegate = self
    window = NSWindow(contentRect: NSRect(x: 0, y: 0, width: 1280, height: 820),
                      styleMask: [.titled, .closable, .miniaturizable, .resizable], backing: .buffered, defer: false)
    window.title = "ChemLarp"
    window.contentView = web
    window.minSize = NSSize(width: 480, height: 400)
    window.center()
    window.setFrameAutosaveName("ChemLarpMain")
    window.makeKeyAndOrderFront(nil)
    web.load(URLRequest(url: URL(string: origin + "/")!))
    NSApp.activate(ignoringOtherApps: true)
  }

  func applicationShouldTerminateAfterLastWindowClosed(_ sender: NSApplication) -> Bool { true }

  // --smoke: once the page has loaded, start the engine worker and canonicalize benzene through it.
  func webView(_ webView: WKWebView, didFinish navigation: WKNavigation!) {
    guard smoke else { return }
    let js = """
      const w = new Worker('engine-worker.js');
      const r = await new Promise((resolve, reject) => {
        w.onmessage = (e) => resolve(e.data);
        w.onerror = (e) => reject(String(e.message || e));
        setTimeout(() => reject('no answer from the engine within 60 s'), 60000);
        w.postMessage({ id: 1, op: 'canonical', args: { smiles: 'C1=CC=CC=C1' } });
      });
      return JSON.stringify({ title: document.title, ok: r.ok, canonical: r.result && r.result.canonical, ms: r.ms });
      """
    webView.callAsyncJavaScript(js, arguments: [:], in: nil, in: .page) { result in
      switch result {
      case .success(let value):
        let text = value as? String ?? "\(String(describing: value))"
        let pass = text.contains("\"canonical\":\"c1ccccc1\"")
        print("smoke: \(text)")
        print(pass ? "smoke passed" : "smoke FAILED")
        exit(pass ? 0 : 1)
      case .failure(let error):
        print("smoke FAILED: \(error)")
        exit(1)
      }
    }
  }

  func webView(_ webView: WKWebView, didFail navigation: WKNavigation!, withError error: Error) {
    if smoke { print("smoke FAILED: \(error)"); exit(1) }
  }

  func webView(_ webView: WKWebView, didFailProvisionalNavigation navigation: WKNavigation!, withError error: Error) {
    if smoke { print("smoke FAILED: \(error)"); exit(1) }
  }

  // Links outside the app open in the default browser.
  func webView(_ webView: WKWebView, decidePolicyFor action: WKNavigationAction, preferences: WKWebpagePreferences,
               decisionHandler: @escaping (WKNavigationActionPolicy, WKWebpagePreferences) -> Void) {
    if action.shouldPerformDownload { return decisionHandler(.download, preferences) }
    if let url = action.request.url, let scheme = url.scheme, ["http", "https"].contains(scheme),
       !url.absoluteString.hasPrefix(origin) {
      NSWorkspace.shared.open(url)
      return decisionHandler(.cancel, preferences)
    }
    decisionHandler(.allow, preferences)
  }

  func webView(_ webView: WKWebView, decidePolicyFor response: WKNavigationResponse,
               decisionHandler: @escaping (WKNavigationResponsePolicy) -> Void) {
    decisionHandler(response.canShowMIMEType ? .allow : .download)
  }

  func webView(_ webView: WKWebView, createWebViewWith configuration: WKWebViewConfiguration, for action: WKNavigationAction,
               windowFeatures: WKWindowFeatures) -> WKWebView? {
    if let url = action.request.url { NSWorkspace.shared.open(url) }
    return nil
  }

  // <input type="file">: the open panel.
  func webView(_ webView: WKWebView, runOpenPanelWith parameters: WKOpenPanelParameters, initiatedByFrame frame: WKFrameInfo,
               completionHandler: @escaping ([URL]?) -> Void) {
    let panel = NSOpenPanel()
    panel.allowsMultipleSelection = parameters.allowsMultipleSelection
    panel.canChooseDirectories = false
    panel.beginSheetModal(for: window) { r in completionHandler(r == .OK ? panel.urls : nil) }
  }

  func webView(_ webView: WKWebView, runJavaScriptAlertPanelWithMessage message: String, initiatedByFrame frame: WKFrameInfo,
               completionHandler: @escaping () -> Void) {
    let a = NSAlert()
    a.messageText = message
    a.beginSheetModal(for: window) { _ in completionHandler() }
  }

  func webView(_ webView: WKWebView, runJavaScriptConfirmPanelWithMessage message: String, initiatedByFrame frame: WKFrameInfo,
               completionHandler: @escaping (Bool) -> Void) {
    let a = NSAlert()
    a.messageText = message
    a.addButton(withTitle: "OK")
    a.addButton(withTitle: "Cancel")
    a.beginSheetModal(for: window) { r in completionHandler(r == .alertFirstButtonReturn) }
  }

  // Exports (SVG, PNG, MOL, ...): a save panel.
  func webView(_ webView: WKWebView, navigationAction: WKNavigationAction, didBecome download: WKDownload) { download.delegate = self }
  func webView(_ webView: WKWebView, navigationResponse: WKNavigationResponse, didBecome download: WKDownload) { download.delegate = self }

  func download(_ download: WKDownload, decideDestinationUsing response: URLResponse, suggestedFilename: String,
                completionHandler: @escaping (URL?) -> Void) {
    let panel = NSSavePanel()
    panel.nameFieldStringValue = suggestedFilename
    panel.beginSheetModal(for: window) { r in
      guard r == .OK, let url = panel.url else { return completionHandler(nil) }
      try? FileManager.default.removeItem(at: url)  // the panel has already asked whether to replace it
      completionHandler(url)
    }
  }

  @objc func reload(_ sender: Any?) { web?.reload() }

  func buildMenu() {
    let main = NSMenu()
    func item(_ title: String, _ action: Selector?, _ key: String, _ mods: NSEvent.ModifierFlags = .command) -> NSMenuItem {
      let i = NSMenuItem(title: title, action: action, keyEquivalent: key)
      i.keyEquivalentModifierMask = mods
      return i
    }
    let app = NSMenu()
    app.addItem(item("About ChemLarp", #selector(NSApplication.orderFrontStandardAboutPanel(_:)), ""))
    app.addItem(.separator())
    app.addItem(item("Hide ChemLarp", #selector(NSApplication.hide(_:)), "h"))
    app.addItem(item("Quit ChemLarp", #selector(NSApplication.terminate(_:)), "q"))
    let edit = NSMenu(title: "Edit")
    edit.addItem(item("Undo", Selector(("undo:")), "z"))
    edit.addItem(item("Redo", Selector(("redo:")), "z", [.command, .shift]))
    edit.addItem(.separator())
    edit.addItem(item("Cut", #selector(NSText.cut(_:)), "x"))
    edit.addItem(item("Copy", #selector(NSText.copy(_:)), "c"))
    edit.addItem(item("Paste", #selector(NSText.paste(_:)), "v"))
    edit.addItem(item("Select All", #selector(NSText.selectAll(_:)), "a"))
    let view = NSMenu(title: "View")
    let reload = item("Reload", #selector(reload(_:)), "r")
    reload.target = self
    view.addItem(reload)
    view.addItem(item("Enter Full Screen", #selector(NSWindow.toggleFullScreen(_:)), "f", [.command, .control]))
    let win = NSMenu(title: "Window")
    win.addItem(item("Minimize", #selector(NSWindow.performMiniaturize(_:)), "m"))
    win.addItem(item("Zoom", #selector(NSWindow.performZoom(_:)), ""))
    for (title, menu) in [("ChemLarp", app), ("Edit", edit), ("View", view), ("Window", win)] {
      let top = NSMenuItem(title: title, action: nil, keyEquivalent: "")
      top.submenu = menu
      main.addItem(top)
    }
    NSApp.mainMenu = main
    NSApp.windowsMenu = win
  }
}

if CommandLine.arguments.contains("--selftest") { selftest() }
let application = NSApplication.shared
let delegate = AppDelegate()
application.delegate = delegate
application.setActivationPolicy(.regular)
application.run()
