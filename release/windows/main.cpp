// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
// ChemLarp for Windows: the web app (WebAssembly engine and UI, in the "web" folder next to this program) in a
// native window through WebView2. The folder is mapped to https://chemlarp.local, so nothing is fetched from
// the network. WebView2 is part of Windows 11 and of current Windows 10; without it the program says where to
// get it.
//
//   ChemLarp.exe --smoke <file>   loads the app, runs an engine call through its worker (WebAssembly), writes
//                                 the result to <file> and exits with 0 on success (used by CI)
#include <windows.h>
#include <shellapi.h>
#include <shlobj.h>
#include <wrl.h>

#include <fstream>
#include <string>

#include "WebView2.h"

using Microsoft::WRL::Callback;
using Microsoft::WRL::ComPtr;

namespace {

HWND g_window = nullptr;
ComPtr<ICoreWebView2Controller> g_controller;
ComPtr<ICoreWebView2> g_view;
std::wstring g_smoke_file;  // non-empty in --smoke mode

std::wstring exe_dir() {
  wchar_t path[MAX_PATH];
  const DWORD n = GetModuleFileNameW(nullptr, path, MAX_PATH);
  std::wstring p(path, n);
  return p.substr(0, p.find_last_of(L"\\/"));
}

std::wstring user_data_dir() {
  PWSTR local = nullptr;
  std::wstring dir;
  if (SUCCEEDED(SHGetKnownFolderPath(FOLDERID_LocalAppData, 0, nullptr, &local))) dir = std::wstring(local) + L"\\ChemLarp\\WebView2";
  CoTaskMemFree(local);
  return dir;
}

void finish_smoke(const std::wstring& text, bool ok) {
  std::wofstream out(g_smoke_file);
  out << text << L"\n" << (ok ? L"smoke passed" : L"smoke FAILED") << L"\n";
  out.close();
  ExitProcess(ok ? 0 : 1);
}

void fail(const wchar_t* what, HRESULT hr) {
  if (!g_smoke_file.empty()) finish_smoke(std::wstring(what) + L" (HRESULT " + std::to_wstring(hr) + L")", false);
  std::wstring msg = std::wstring(what) +
                     L"\n\nChemLarp needs the Microsoft Edge WebView2 Runtime, which comes with Windows 11 and current "
                     L"Windows 10. It can be installed from https://developer.microsoft.com/microsoft-edge/webview2/";
  MessageBoxW(g_window, msg.c_str(), L"ChemLarp", MB_ICONERROR);
  PostQuitMessage(1);
}

const wchar_t* kSmokeScript = LR"JS(
(async () => {
  try {
    const w = new Worker('engine-worker.js');
    const r = await new Promise((resolve, reject) => {
      w.onmessage = (e) => resolve(e.data);
      w.onerror = (e) => reject(String(e.message || e));
      setTimeout(() => reject('no answer from the engine within 60 s'), 60000);
      w.postMessage({ id: 1, op: 'canonical', args: { smiles: 'C1=CC=CC=C1' } });
    });
    chrome.webview.postMessage(JSON.stringify({ title: document.title, ok: r.ok, canonical: r.result && r.result.canonical, ms: r.ms }));
  } catch (err) {
    chrome.webview.postMessage(JSON.stringify({ error: String(err) }));
  }
})();
)JS";

void create_webview() {
  const std::wstring web = exe_dir() + L"\\web";
  CreateCoreWebView2EnvironmentWithOptions(
      nullptr, user_data_dir().c_str(), nullptr,
      Callback<ICoreWebView2CreateCoreWebView2EnvironmentCompletedHandler>(
          [web](HRESULT hr, ICoreWebView2Environment* env) -> HRESULT {
            if (FAILED(hr) || !env) {
              fail(L"The WebView2 runtime could not be started.", hr);
              return S_OK;
            }
            env->CreateCoreWebView2Controller(
                g_window, Callback<ICoreWebView2CreateCoreWebView2ControllerCompletedHandler>(
                              [web](HRESULT hr2, ICoreWebView2Controller* controller) -> HRESULT {
                                if (FAILED(hr2) || !controller) {
                                  fail(L"The WebView2 window could not be created.", hr2);
                                  return S_OK;
                                }
                                g_controller = controller;
                                g_controller->get_CoreWebView2(&g_view);
                                RECT bounds;
                                GetClientRect(g_window, &bounds);
                                g_controller->put_Bounds(bounds);

                                ComPtr<ICoreWebView2_3> view3;
                                if (FAILED(g_view.As(&view3))) {
                                  fail(L"This WebView2 runtime is too old (folder mapping is missing).", E_NOINTERFACE);
                                  return S_OK;
                                }
                                view3->SetVirtualHostNameToFolderMapping(L"chemlarp.local", web.c_str(),
                                                                         COREWEBVIEW2_HOST_RESOURCE_ACCESS_KIND_ALLOW);
                                // Links to other sites and window.open go to the default browser.
                                g_view->add_NewWindowRequested(
                                    Callback<ICoreWebView2NewWindowRequestedEventHandler>(
                                        [](ICoreWebView2*, ICoreWebView2NewWindowRequestedEventArgs* args) -> HRESULT {
                                          LPWSTR uri = nullptr;
                                          args->get_Uri(&uri);
                                          if (uri) ShellExecuteW(nullptr, L"open", uri, nullptr, nullptr, SW_SHOWNORMAL);
                                          CoTaskMemFree(uri);
                                          args->put_Handled(TRUE);
                                          return S_OK;
                                        })
                                        .Get(),
                                    nullptr);
                                g_view->add_NavigationStarting(
                                    Callback<ICoreWebView2NavigationStartingEventHandler>(
                                        [](ICoreWebView2*, ICoreWebView2NavigationStartingEventArgs* args) -> HRESULT {
                                          LPWSTR uri = nullptr;
                                          args->get_Uri(&uri);
                                          const std::wstring u = uri ? uri : L"";
                                          CoTaskMemFree(uri);
                                          if (u.rfind(L"https://chemlarp.local/", 0) != 0 && u.rfind(L"http", 0) == 0) {
                                            ShellExecuteW(nullptr, L"open", u.c_str(), nullptr, nullptr, SW_SHOWNORMAL);
                                            args->put_Cancel(TRUE);
                                          }
                                          return S_OK;
                                        })
                                        .Get(),
                                    nullptr);
                                if (!g_smoke_file.empty()) {
                                  g_view->add_WebMessageReceived(
                                      Callback<ICoreWebView2WebMessageReceivedEventHandler>(
                                          [](ICoreWebView2*, ICoreWebView2WebMessageReceivedEventArgs* args) -> HRESULT {
                                            LPWSTR msg = nullptr;
                                            args->TryGetWebMessageAsString(&msg);
                                            const std::wstring text = msg ? msg : L"";
                                            CoTaskMemFree(msg);
                                            finish_smoke(text, text.find(L"\"canonical\":\"c1ccccc1\"") != std::wstring::npos);
                                            return S_OK;
                                          })
                                          .Get(),
                                      nullptr);
                                  g_view->add_NavigationCompleted(
                                      Callback<ICoreWebView2NavigationCompletedEventHandler>(
                                          [](ICoreWebView2* sender, ICoreWebView2NavigationCompletedEventArgs* args) -> HRESULT {
                                            BOOL ok = FALSE;
                                            args->get_IsSuccess(&ok);
                                            if (!ok) finish_smoke(L"the page did not load", false);
                                            sender->ExecuteScript(kSmokeScript, nullptr);
                                            return S_OK;
                                          })
                                          .Get(),
                                      nullptr);
                                }
                                g_view->Navigate(L"https://chemlarp.local/index.html");
                                return S_OK;
                              })
                              .Get());
            return S_OK;
          })
          .Get());
}

LRESULT CALLBACK window_proc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
  switch (msg) {
    case WM_SIZE:
      if (g_controller) {
        RECT bounds;
        GetClientRect(hwnd, &bounds);
        g_controller->put_Bounds(bounds);
      }
      return 0;
    case WM_DESTROY:
      PostQuitMessage(0);
      return 0;
    default:
      return DefWindowProcW(hwnd, msg, wp, lp);
  }
}

}  // namespace

int WINAPI wWinMain(HINSTANCE instance, HINSTANCE, PWSTR, int show) {
  int argc = 0;
  LPWSTR* argv = CommandLineToArgvW(GetCommandLineW(), &argc);
  for (int i = 1; i + 1 < argc; ++i)
    if (std::wstring(argv[i]) == L"--smoke") g_smoke_file = argv[i + 1];
  LocalFree(argv);

  SetProcessDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);
  CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED);
  WNDCLASSEXW wc{};
  wc.cbSize = sizeof(wc);
  wc.lpfnWndProc = window_proc;
  wc.hInstance = instance;
  wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
  wc.hIcon = LoadIconW(instance, L"APPICON");
  wc.hIconSm = wc.hIcon;
  wc.lpszClassName = L"ChemLarpWindow";
  RegisterClassExW(&wc);
  g_window = CreateWindowExW(0, wc.lpszClassName, L"ChemLarp", WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT, 1280, 820,
                             nullptr, nullptr, instance, nullptr);
  ShowWindow(g_window, show);
  UpdateWindow(g_window);
  if (!g_smoke_file.empty()) SetTimer(g_window, 1, 120000, [](HWND, UINT, UINT_PTR, DWORD) { finish_smoke(L"timed out", false); });
  create_webview();
  MSG m;
  while (GetMessageW(&m, nullptr, 0, 0)) {
    TranslateMessage(&m);
    DispatchMessageW(&m);
  }
  CoUninitialize();
  return static_cast<int>(m.wParam);
}
