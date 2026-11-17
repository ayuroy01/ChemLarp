// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
package io.github.ayuroy01.chemlarp;

import android.content.res.AssetManager;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.InetAddress;
import java.net.ServerSocket;
import java.net.Socket;
import java.net.URLDecoder;
import java.nio.charset.StandardCharsets;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 * Serves the bundled web app (assets/web) on 127.0.0.1 at a port the system chooses, so the WebView loads it
 * exactly as a browser would: workers, WebAssembly and file reading all work, and nothing comes from the network.
 */
final class LocalServer {
    private final AssetManager assets;
    private final String root;
    private ServerSocket socket;
    private final ExecutorService pool = Executors.newCachedThreadPool();

    LocalServer(AssetManager assets, String root) {
        this.assets = assets;
        this.root = root;
    }

    /** Starts listening; returns the port. */
    int start() throws IOException {
        socket = new ServerSocket(0, 50, InetAddress.getByName("127.0.0.1"));
        Thread t = new Thread(() -> {
            while (!socket.isClosed()) {
                try {
                    Socket c = socket.accept();
                    pool.execute(() -> serve(c));
                } catch (IOException e) {
                    return;  // closed
                }
            }
        }, "chemlarp-server");
        t.setDaemon(true);
        t.start();
        return socket.getLocalPort();
    }

    void stop() {
        try {
            if (socket != null) socket.close();
        } catch (IOException ignored) {
        }
        pool.shutdownNow();
    }

    /** The asset for a request path, or null (unknown, or a path that leaves the app's folder). */
    String assetPath(String target) {
        String path = target;
        int q = path.indexOf('?');
        if (q >= 0) path = path.substring(0, q);
        int h = path.indexOf('#');
        if (h >= 0) path = path.substring(0, h);
        try {
            path = URLDecoder.decode(path.replace("+", "%2B"), "UTF-8");
        } catch (Exception e) {
            return null;
        }
        if (path.endsWith("/")) path += "index.html";
        StringBuilder out = new StringBuilder(root);
        for (String part : path.split("/")) {
            if (part.isEmpty() || part.equals(".")) continue;
            if (part.equals("..") || part.contains("\\")) return null;
            out.append('/').append(part);
        }
        return out.toString();
    }

    private void serve(Socket c) {
        try (Socket s = c; InputStream in = s.getInputStream(); OutputStream out = s.getOutputStream()) {
            ByteArrayOutputStream head = new ByteArrayOutputStream();
            int prev3 = 0, prev2 = 0, prev1 = 0, b;
            while ((b = in.read()) != -1 && head.size() < 65536) {
                head.write(b);
                if (prev3 == '\r' && prev2 == '\n' && prev1 == '\r' && b == '\n') break;
                prev3 = prev2;
                prev2 = prev1;
                prev1 = b;
            }
            String first = head.toString("ISO-8859-1").split("\r\n", 2)[0];
            String[] parts = first.split(" ");
            String method = parts.length > 0 ? parts[0] : "";
            String status = "404 Not Found", type = "text/plain; charset=utf-8";
            byte[] body = "not found\n".getBytes(StandardCharsets.UTF_8);
            if (!method.equals("GET") && !method.equals("HEAD")) {
                status = "405 Method Not Allowed";
                body = new byte[0];
            } else if (parts.length >= 2) {
                String asset = assetPath(parts[1]);
                if (asset != null) {
                    try (InputStream a = assets.open(asset)) {
                        body = readAll(a);
                        status = "200 OK";
                        type = mime(asset);
                    } catch (IOException ignored) {
                        // not found
                    }
                }
            }
            String headers = "HTTP/1.1 " + status + "\r\nContent-Type: " + type + "\r\nContent-Length: " + body.length +
                    "\r\nCache-Control: no-cache\r\nX-Content-Type-Options: nosniff\r\nConnection: close\r\n\r\n";
            out.write(headers.getBytes(StandardCharsets.ISO_8859_1));
            if (!method.equals("HEAD")) out.write(body);
            out.flush();
        } catch (IOException ignored) {
            // the client went away
        }
    }

    private static byte[] readAll(InputStream in) throws IOException {
        ByteArrayOutputStream out = new ByteArrayOutputStream();
        byte[] buf = new byte[65536];
        int n;
        while ((n = in.read(buf)) > 0) out.write(buf, 0, n);
        return out.toByteArray();
    }

    static String mime(String path) {
        String p = path.toLowerCase();
        if (p.endsWith(".html")) return "text/html; charset=utf-8";
        if (p.endsWith(".js") || p.endsWith(".mjs")) return "text/javascript; charset=utf-8";
        if (p.endsWith(".css")) return "text/css; charset=utf-8";
        if (p.endsWith(".wasm")) return "application/wasm";
        if (p.endsWith(".json") || p.endsWith(".webmanifest")) return "application/json";
        if (p.endsWith(".svg")) return "image/svg+xml";
        if (p.endsWith(".png")) return "image/png";
        return "text/plain; charset=utf-8";
    }
}
