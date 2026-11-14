; Inno Setup script for ChemLarp (Windows x64). Built by the packages workflow:
;   iscc /DAppVersion=0.1.0 /DStage=<folder with ChemLarp.exe, web\, bin\chemlarp.exe> /DOutDir=<dir> chemlarp.iss
; Installs for the current user (no administrator rights) under %LOCALAPPDATA%\Programs\ChemLarp, with a Start
; menu entry and an uninstaller. The installer is not code-signed, so SmartScreen may warn before it runs.
#ifndef AppVersion
  #define AppVersion "0.0.0"
#endif

[Setup]
AppId={{7C2B5E1A-4F3D-4C8B-9A6E-2D1F0B3C5A71}
AppName=ChemLarp
AppVersion={#AppVersion}
AppVerName=ChemLarp {#AppVersion}
AppPublisher=ayuroy01
AppPublisherURL=https://ayuroy01.github.io/ChemLarp/
AppSupportURL=https://github.com/ayuroy01/ChemLarp/issues
DefaultDirName={autopf}\ChemLarp
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
DisableProgramGroupPage=yes
OutputDir={#OutDir}
OutputBaseFilename=ChemLarp-{#AppVersion}-windows-x64-setup
Compression=lzma2
SolidCompression=yes
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
SetupIconFile=chemlarp.ico
UninstallDisplayIcon={app}\ChemLarp.exe
WizardStyle=modern
LicenseFile=..\..\LICENSE

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "{#Stage}\*"; DestDir: "{app}"; Flags: recursesubdirs ignoreversion

[Icons]
Name: "{autoprograms}\ChemLarp"; Filename: "{app}\ChemLarp.exe"
Name: "{autodesktop}\ChemLarp"; Filename: "{app}\ChemLarp.exe"; Tasks: desktopicon

[Run]
Filename: "{app}\ChemLarp.exe"; Description: "{cm:LaunchProgram,ChemLarp}"; Flags: nowait postinstall skipifsilent
