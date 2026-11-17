// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
package io.github.ayuroy01.chemlarp;

import android.app.Activity;
import android.app.AlertDialog;
import android.content.ActivityNotFoundException;
import android.content.Intent;
import android.net.Uri;
import android.os.Bundle;
import android.util.Log;
import android.webkit.JavascriptInterface;
import android.webkit.ValueCallback;
import android.webkit.WebChromeClient;
import android.webkit.WebResourceRequest;
import android.webkit.WebSettings;
import android.webkit.WebView;
import android.webkit.WebViewClient;

import java.io.IOException;

/**
 * ChemLarp for Android: the web app (WebAssembly engine and UI, bundled in assets/web) in a WebView, served by
 * LocalServer on 127.0.0.1. Started with the boolean extra "smoke", it runs one engine call through the app's
 * worker and logs the result under the tag ChemLarpSmoke (used by CI on an emulator).
 */
public class MainActivity extends Activity {
    private static final int PICK_FILE = 1;
    private static final String SMOKE_JS =
            "(async () => { try {" +
            " const w = new Worker('engine-worker.js');" +
            " const r = await new Promise((resolve, reject) => {" +
            "   w.onmessage = (e) => resolve(e.data); w.onerror = (e) => reject(String(e.message || e));" +
            "   setTimeout(() => reject('no answer from the engine within 60 s'), 60000);" +
            "   w.postMessage({ id: 1, op: 'canonical', args: { smiles: 'C1=CC=CC=C1' } }); });" +
            " ChemLarpSmoke.report(JSON.stringify({ title: document.title, ok: r.ok, canonical: r.result && r.result.canonical, ms: r.ms }));" +
            "} catch (err) { ChemLarpSmoke.report(JSON.stringify({ error: String(err) })); } })();";

    private WebView web;
    private LocalServer server;
    private ValueCallback<Uri[]> fileCallback;

    @Override
    protected void onCreate(Bundle state) {
        super.onCreate(state);
        final boolean smoke = getIntent().getBooleanExtra("smoke", false);
        server = new LocalServer(getAssets(), "web");
        final int port;
        try {
            port = server.start();
        } catch (IOException e) {
            new AlertDialog.Builder(this).setMessage("ChemLarp could not start its local server.")
                    .setPositiveButton(android.R.string.ok, (d, w) -> finish()).show();
            return;
        }
        final String origin = "http://127.0.0.1:" + port;
        web = new WebView(this);
        WebSettings s = web.getSettings();
        s.setJavaScriptEnabled(true);
        s.setDomStorageEnabled(true);
        s.setAllowFileAccess(false);
        s.setAllowContentAccess(false);
        web.setWebViewClient(new WebViewClient() {
            @Override
            public boolean shouldOverrideUrlLoading(WebView view, WebResourceRequest request) {
                Uri u = request.getUrl();
                if (u.toString().startsWith(origin)) return false;
                try {
                    startActivity(new Intent(Intent.ACTION_VIEW, u));  // other sites: the browser
                } catch (ActivityNotFoundException ignored) {
                }
                return true;
            }

            @Override
            public void onPageFinished(WebView view, String url) {
                if (smoke) view.evaluateJavascript(SMOKE_JS, null);
            }
        });
        web.setWebChromeClient(new WebChromeClient() {
            @Override
            public boolean onShowFileChooser(WebView view, ValueCallback<Uri[]> callback, FileChooserParams params) {
                if (fileCallback != null) fileCallback.onReceiveValue(null);
                fileCallback = callback;
                try {
                    startActivityForResult(params.createIntent(), PICK_FILE);
                } catch (ActivityNotFoundException e) {
                    fileCallback = null;
                    return false;
                }
                return true;
            }
        });
        if (smoke) {
            web.addJavascriptInterface(new Object() {
                @JavascriptInterface
                public void report(String json) {
                    Log.i("ChemLarpSmoke", json);
                    Log.i("ChemLarpSmoke", json.contains("\"canonical\":\"c1ccccc1\"") ? "smoke passed" : "smoke FAILED");
                }
            }, "ChemLarpSmoke");
        }
        setContentView(web);
        web.loadUrl(origin + "/index.html");
    }

    @Override
    protected void onActivityResult(int request, int result, Intent data) {
        if (request == PICK_FILE && fileCallback != null) {
            fileCallback.onReceiveValue(WebChromeClient.FileChooserParams.parseResult(result, data));
            fileCallback = null;
            return;
        }
        super.onActivityResult(request, result, data);
    }

    @Override
    @SuppressWarnings("deprecation")
    public void onBackPressed() {
        if (web != null && web.canGoBack()) web.goBack();
        else super.onBackPressed();
    }

    @Override
    protected void onDestroy() {
        if (server != null) server.stop();
        if (web != null) web.destroy();
        super.onDestroy();
    }
}
