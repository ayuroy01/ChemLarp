#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
# Android app (APK): release/android around the built web app.
#
#   scripts/build_android.sh <web dir> [out dir]
#
# Needs JDK 17+, Gradle 8.9+ on PATH and the Android SDK (ANDROID_HOME, with platform 35 and build-tools).
# The APK is signed with the keystore in CHEMLARP_KEYSTORE (password CHEMLARP_KEYSTORE_PASS, alias chemlarp)
# when set; otherwise with a throwaway key made here, which is enough to install it by hand ("unknown
# sources") but not for Google Play, and an update signed with a different key will not install over it.
set -euo pipefail
cd "$(dirname "$0")/.."
WEB="${1:?usage: scripts/build_android.sh <web dir> [out dir]}"
OUT="${2:-build/android}"
test -f "$WEB/index.html" && test -f "$WEB/chemlarp.wasm" || { echo "$WEB is not a built web app" >&2; exit 1; }
VERSION=$(python3 scripts/check_packages.py --version-only)
# versionCode: 1.2.3 -> 10203 (must increase with every release on the Play Store).
CODE=$(python3 -c "import sys; a=[int(x) for x in sys.argv[1].split('.')]; print(a[0]*10000+a[1]*100+a[2])" "$VERSION")
ASSETS=release/android/app/src/main/assets
rm -rf "$ASSETS" && mkdir -p "$ASSETS" && cp -R "$WEB" "$ASSETS/web"
gradle -p release/android --no-daemon -q assembleRelease -PchemlarpVersion="$VERSION" -PchemlarpVersionCode="$CODE"
UNSIGNED=release/android/app/build/outputs/apk/release/app-release-unsigned.apk
BT=$(ls -d "$ANDROID_HOME"/build-tools/* | sort -V | tail -1)
mkdir -p "$OUT"
APK="$OUT/ChemLarp-$VERSION-android.apk"
"$BT/zipalign" -f -p 4 "$UNSIGNED" "$OUT/aligned.apk"
KS="${CHEMLARP_KEYSTORE:-}"
PASS="${CHEMLARP_KEYSTORE_PASS:-}"
if [ -z "$KS" ]; then
  KS="$OUT/throwaway.jks"
  PASS=throwaway
  rm -f "$KS"
  keytool -genkeypair -keystore "$KS" -storepass "$PASS" -keypass "$PASS" -alias chemlarp -keyalg RSA -keysize 2048 \
    -validity 10000 -dname "CN=ChemLarp test build" >/dev/null 2>&1
  echo "signed with a throwaway key (set CHEMLARP_KEYSTORE for a release key)"
fi
"$BT/apksigner" sign --ks "$KS" --ks-pass "pass:$PASS" --ks-key-alias chemlarp --out "$APK" "$OUT/aligned.apk"
"$BT/apksigner" verify "$APK"
rm -f "$OUT/aligned.apk" "$OUT"/*.idsig
echo "built $APK"
