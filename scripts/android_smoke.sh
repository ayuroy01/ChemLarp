#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (C) 2026 ayuroy01. See LICENSE and NOTICE.
# Installs the APK on a running emulator or device, starts the app with the "smoke" extra and waits for its
# engine call to report in logcat (tag ChemLarpSmoke). Exit 0 on "smoke passed".
set -euo pipefail
APK="$(ls $1 | head -1)"
adb install -r "$APK"
adb logcat -c
adb shell am start -n io.github.ayuroy01.chemlarp/.MainActivity --ez smoke true
for i in $(seq 1 150); do
  out=$(adb logcat -d -s ChemLarpSmoke:I || true)
  if echo "$out" | grep -q "smoke passed"; then echo "$out"; exit 0; fi
  if echo "$out" | grep -q "smoke FAILED"; then echo "$out"; exit 1; fi
  sleep 2
done
echo "no result from the app within 300 s"; adb logcat -d | tail -50; exit 1
