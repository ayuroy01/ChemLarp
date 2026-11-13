// Offline support: cache the app shell, the WebAssembly engine and example data.
// The build stamps a new cache name on every deploy, so updates replace old files.
const CACHE = 'chemlarp-__BUILD_ID__';
const FILES = [
  './', 'index.html', 'style.css', 'app.js', 'viewer.js', 'editor.js', 'chart.js', 'engine-worker.js',
  'chemlarp.js', 'chemlarp.wasm', 'manifest.webmanifest',
  'icons/icon.svg', 'icons/icon-192.png', 'icons/icon-512.png', 'icons/maskable-512.png',
  'data/roadmap.json', 'data/1crn.pdb', 'data/mechanisms/h2o2.ocm', 'data/mechanisms/h2o2.yaml',
];

self.addEventListener('install', (event) => {
  // cache: 'reload' skips the browser's HTTP cache (GitHub Pages sends max-age=600),
  // otherwise a new worker could store the previous deploy's files.
  event.waitUntil(
    caches.open(CACHE)
      .then((c) => c.addAll(FILES.map((f) => new Request(f, { cache: 'reload' }))))
      .then(() => self.skipWaiting()),
  );
});

self.addEventListener('activate', (event) => {
  event.waitUntil(
    caches.keys()
      .then((keys) => Promise.all(keys.filter((k) => (k.startsWith('chemlarp-') || k.startsWith('omnichem-')) && k !== CACHE).map((k) => caches.delete(k))))
      .then(() => self.clients.claim()),
  );
});

// Network first for data that changes (roadmap), cache first for everything else.
self.addEventListener('fetch', (event) => {
  const url = new URL(event.request.url);
  if (event.request.method !== 'GET' || url.origin !== location.origin) return;
  if (url.pathname.endsWith('roadmap.json')) {
    event.respondWith(fetch(event.request, { cache: 'no-store' }).then((r) => {
      const copy = r.clone();
      caches.open(CACHE).then((c) => c.put(event.request, copy));
      return r;
    }).catch(() => caches.match(event.request)));
    return;
  }
  event.respondWith(caches.match(event.request).then((hit) => hit || fetch(event.request)));
});
