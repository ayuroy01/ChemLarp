#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
# macOS app and disk image: release/macos/main.swift around the built web app.
#
#   scripts/build_macos_app.sh <web dir> [out dir]
#
# <web dir> is the output of scripts/build_web.sh (index.html, chemlarp.js, chemlarp.wasm, ...). Output:
# <out>/ChemLarp.app (universal: arm64 and x86_64, macOS 12+) and <out>/ChemLarp-<version>-macos.dmg. The app is
# signed ad hoc only (no Developer ID, not notarized): Gatekeeper asks before the first launch of a downloaded
# copy (right-click > Open, or System Settings > Privacy & Security > Open Anyway).
set -euo pipefail
cd "$(dirname "$0")/.."
WEB="${1:?usage: scripts/build_macos_app.sh <web dir> [out dir]}"
OUT="${2:-build/macos}"
test -f "$WEB/index.html" && test -f "$WEB/chemlarp.wasm" || { echo "$WEB is not a built web app" >&2; exit 1; }
VERSION=$(python3 scripts/check_packages.py --version-only)
APP="$OUT/ChemLarp.app"
rm -rf "$OUT"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources" "$OUT/obj"

for arch in arm64 x86_64; do
  swiftc -O -target "$arch-apple-macos12.0" release/macos/main.swift -o "$OUT/obj/ChemLarp-$arch"
done
lipo -create "$OUT/obj/ChemLarp-arm64" "$OUT/obj/ChemLarp-x86_64" -output "$APP/Contents/MacOS/ChemLarp"

cp -R "$WEB" "$APP/Contents/Resources/web"

# Icon from the web app's 512 px icon.
ICONSET="$OUT/obj/AppIcon.iconset"
mkdir -p "$ICONSET"
for s in 16 32 128 256 512; do
  sips -z "$s" "$s" web/icons/icon-512.png --out "$ICONSET/icon_${s}x${s}.png" >/dev/null
  d=$((s * 2)); [ "$d" -le 512 ] && sips -z "$d" "$d" web/icons/icon-512.png --out "$ICONSET/icon_${s}x${s}@2x.png" >/dev/null
done
iconutil -c icns "$ICONSET" -o "$APP/Contents/Resources/AppIcon.icns"

cat > "$APP/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleIdentifier</key><string>io.github.ayuroy01.chemlarp</string>
  <key>CFBundleName</key><string>ChemLarp</string>
  <key>CFBundleDisplayName</key><string>ChemLarp</string>
  <key>CFBundleExecutable</key><string>ChemLarp</string>
  <key>CFBundleIconFile</key><string>AppIcon</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleShortVersionString</key><string>$VERSION</string>
  <key>CFBundleVersion</key><string>$VERSION</string>
  <key>LSMinimumSystemVersion</key><string>12.0</string>
  <key>LSApplicationCategoryType</key><string>public.app-category.education</string>
  <key>NSHighResolutionCapable</key><true/>
  <key>NSHumanReadableCopyright</key><string>Copyright © 2026 ayuroy01. GNU AGPL-3.0-only</string>
  <key>NSAppTransportSecurity</key><dict><key>NSAllowsLocalNetworking</key><true/></dict>
</dict>
</plist>
PLIST

codesign --force --deep --sign - "$APP"
"$APP/Contents/MacOS/ChemLarp" --selftest

STAGE="$OUT/obj/dmg"
mkdir -p "$STAGE"
cp -R "$APP" "$STAGE/"
ln -s /Applications "$STAGE/Applications"
hdiutil create -volname "ChemLarp" -srcfolder "$STAGE" -ov -format UDZO "$OUT/ChemLarp-$VERSION-macos.dmg" >/dev/null
rm -rf "$OUT/obj"
echo "built $APP and $OUT/ChemLarp-$VERSION-macos.dmg"
