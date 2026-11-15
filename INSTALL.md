# Installing ChemLarp: a step-by-step guide

This guide is written for chemists, students and teachers, not programmers. You do not need to know anything
about GitHub, code or the command line. Pick **one** option.

| Option | Best for | Time | Needs |
|---|---|---|---|
| [1. Web browser](#option-1-use-it-in-your-web-browser-no-installation) | Everyone; try this first | 1 minute | A browser |
| [2. Desktop or phone app](#option-2-install-the-app) | Using ChemLarp in its own window, offline | 5 minutes | Windows 10/11, macOS 12+, or Android 8+ |
| [3. Python](#option-3-python-for-scripts-and-notebooks) | Scripts, Jupyter notebooks, many molecules at once | 10 minutes | Python 3.10 or newer |

Options 1 and 2 give you exactly the same program. Everything runs on your own device: no account is needed and
your molecules are never uploaded anywhere.

---

## Option 1: use it in your web browser (no installation)

1. Open **Google Chrome**, **Microsoft Edge**, **Safari** or **Firefox**.
2. Go to **https://ayuroy01.github.io/ChemLarp/**
3. Wait a few seconds the first time. Below the ChemLarp logo, top left, it says *"Starting engine…"* and then
   *"Engine 0.1.0 · offline ready"*.
4. That's it. Type a molecule in the **Compound** box at the top, for example `CCO` (ethanol) or the name of
   one of the samples, and press **Open**. Use the menu on the left (on a phone, the bar at the bottom) to switch
   tools.

**Keep it like an app (optional, recommended).** It then opens from your dock, Start menu or home screen and works
without internet:

- **Windows, Mac or Linux with Chrome or Edge:** click **Install ChemLarp** at the bottom of the left menu, or the
  small install icon (a screen with a down arrow) at the right end of the address bar, then **Install**.
- **Android (Chrome):** tap the **⋮** menu (top right) → **Install app** or **Add to Home screen**.
- **iPhone or iPad (Safari):** tap the **Share** button (a square with an arrow pointing up) → scroll down →
  **Add to Home Screen** → **Add**.

> **Tip:** if a calculation seems slow, close other heavy browser tabs. Quantum chemistry on large molecules can
> take a minute or two; it runs in the background, so the page stays usable.

---

## Option 2: install the app

All downloads are on the **Releases page**: **https://github.com/ayuroy01/ChemLarp/releases/latest**

On that page, scroll down to the section called **Assets** (click **Assets** if it is folded), and click the
file for your device. You can ignore all the other files.

| Your device | Click this file |
|---|---|
| Windows PC | the one ending in **`-windows-x64-setup.exe`** |
| Mac (any, Apple silicon or Intel) | the one ending in **`-macos.dmg`** |
| Android phone or tablet | the one ending in **`-android.apk`** |
| iPhone or iPad | none: use [Option 1](#option-1-use-it-in-your-web-browser-no-installation) and *Add to Home Screen* |
| Linux | use [Option 1](#option-1-use-it-in-your-web-browser-no-installation), or [Option 3](#option-3-python-for-scripts-and-notebooks) |

**About security warnings.** ChemLarp is made by an independent developer. Its apps are not yet "signed" with a
paid Apple or Microsoft developer certificate, so your computer will warn you the first time. That is normal for
free, open-source programs distributed this way, and the steps below show how to continue. You can read every
line of the source code in this repository, and the apps are built automatically from it on GitHub's servers.

### Windows

1. Download the file ending in **`-windows-x64-setup.exe`** from the Releases page.
2. Open your **Downloads** folder and double-click that file.
3. If a blue window says **"Windows protected your PC"**, click **More info**, then **Run anyway**.
4. The installer opens. If it asks who to install for, choose **Install for me only** (no administrator rights
   needed). Click **Next** and then **Install**.
5. Click **Finish**. ChemLarp is now in your **Start menu**: press the Windows key and type `ChemLarp`.

To uninstall later: **Settings → Apps → Installed apps → ChemLarp → Uninstall**.

> If the app window stays blank, your PC is missing Microsoft Edge WebView2 (included in Windows 11 and in up-to-date
> Windows 10). Install the free "Evergreen Bootstrapper" from
> https://developer.microsoft.com/microsoft-edge/webview2/ and open ChemLarp again.

### Mac

1. Download the file ending in **`-macos.dmg`** from the Releases page.
2. Open your **Downloads** folder and double-click the `.dmg` file. A window opens with the ChemLarp icon.
3. **Drag the ChemLarp icon onto the Applications folder** shown in the same window.
4. Open your **Applications** folder and double-click **ChemLarp**.
5. The first time, macOS says it *"cannot verify"* the app or that it is *"from an unidentified developer"*.
   Click **Done** (or **OK**). Do **not** click *Move to Trash*. Then:
   - Open **System Settings** (Apple menu  → System Settings) → **Privacy & Security**.
   - Scroll down to the message about **ChemLarp** and click **Open Anyway**.
   - Enter your Mac password if asked, then click **Open Anyway** once more.

   On older macOS (12 Monterey to 14 Sonoma) you can instead **right-click** (or Control-click) ChemLarp in
   Applications, choose **Open**, then click **Open** in the dialog.
6. From now on ChemLarp opens normally. To keep it in the Dock, right-click its Dock icon → **Options → Keep in Dock**.

To uninstall: drag ChemLarp from Applications to the Trash.

### Android

1. On your phone, open the Releases page in Chrome and tap the file ending in **`-android.apk`**.
2. When the download finishes, tap **Open** in the notification (or find the file in the **Files** app, under
   **Downloads**).
3. Android may say that installing apps from this source is not allowed. Tap **Settings**, switch on **Allow
   from this source**, then go back.
4. Tap **Install**, then **Open**. If Google Play Protect shows a warning about an unknown app, tap **More
   details → Install anyway**.

To uninstall: press and hold the ChemLarp icon → **Uninstall** (or **App info → Uninstall**).

---

## Option 3: Python (for scripts and notebooks)

Use this when you want to process many molecules at once, or use ChemLarp from Jupyter.

**Step 1: check that you have Python 3.10 or newer.**

- **Windows:** open the Start menu, type `cmd`, press Enter, and in the black window type `python --version`
  and press Enter. If it says "not found" or opens the Microsoft Store, install Python from
  https://www.python.org/downloads/ and **tick "Add python.exe to PATH"** on the first screen of the installer.
- **Mac:** open **Terminal** (Applications → Utilities → Terminal) and type `python3 --version`. If it is
  missing or older than 3.10, install it from https://www.python.org/downloads/.

**Step 2: download the right file.** On the [Releases page](https://github.com/ayuroy01/ChemLarp/releases/latest),
under **Assets**, download the file starting with `chemlarp-` and ending in `.whl` whose name matches your
computer:

| Your computer | The `.whl` name contains |
|---|---|
| Windows (64-bit) | `win_amd64` |
| Mac with Apple silicon (M1, M2, M3, M4…) | `macosx` and `arm64` |
| Mac with an Intel processor | `macosx` and `x86_64` |
| Linux | `manylinux` and `x86_64` (or `aarch64` on ARM) |

(Not sure which Mac you have? Apple menu  → **About This Mac**: "Chip: Apple M…" means Apple silicon.)

**Step 3: install it.** In the same Terminal or command window, type `pip install ` (with a space at the end),
then drag the downloaded `.whl` file onto the window so its path appears, and press Enter. On a Mac use
`pip3` instead of `pip`.

**Step 4: try it.**

```bash
python -c "import chemlarp; print(chemlarp.iupac_name('CC(=O)Oc1ccccc1C(=O)O'))"
```

It should print `2-(acetyloxy)benzoic acid`. The install also adds a `chemlarp` command: type `chemlarp` to
list everything it can do, for example `chemlarp name CCO`.

In Python or Jupyter:

```python
import chemlarp
chemlarp.inchi("CCO")["key"]                 # InChIKey of ethanol
chemlarp.qc("O", "b3lyp", "6-31G*")["energy"]  # DFT energy of water, in hartree
chemlarp.embed("CCO")["sdf"]                 # 3D structure as an SDF block
```

---

## Troubleshooting

| Problem | What to do |
|---|---|
| The web page stays on "Starting engine…" | Use an up-to-date Chrome, Edge, Safari (iOS 15 or later) or Firefox. Private/incognito windows and some school or company networks block parts the app needs; try a normal window or another network. |
| "Invalid structure" or a red error under the compound box | The text is read as a SMILES string. Check the message: it points to the atom with the problem (for example a carbon with five bonds). You can also draw the molecule in the **Draw** tool instead. |
| A calculation takes a long time | Quantum chemistry and molecular dynamics depend on molecule size. Try a smaller basis set (6-31G or STO-3G) or a shorter simulation first. |
| The Releases page has no files | Packages are built automatically for each version; if none is published yet, use the web app (Option 1), which is always up to date. |
| Something else, or a result looks wrong | Open an issue at https://github.com/ayuroy01/ChemLarp/issues with what you typed and what you expected. |

## Building from source (developers)

See [the developer section of the README](README.md#for-developers) and [docs/PACKAGING.md](docs/PACKAGING.md).
