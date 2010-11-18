; example2.nsi
;
; This script is based on example1.nsi, but it remember the directory, 
; has uninstall support and (optionally) installs start menu shortcuts.
;
; It will install example2.nsi into a directory that the user selects,

;--------------------------------

!include EnvVarUpdate.nsh
!include "MUI.nsh"

; The name of the installer
Name "FRANz for Windows"

; The file to write
OutFile "FRANz.exe"

; The default installation directory
InstallDir $PROGRAMFILES\FRANz

; Registry key to check for directory (so if you install again, it will 
; overwrite the old one automatically)
InstallDirRegKey HKLM "Software\FRANz" "Install_Dir"

; Request application privileges for Windows Vista
RequestExecutionLevel admin

;--------------------------------

; Pages


  !insertmacro MUI_PAGE_WELCOME
  !insertmacro MUI_PAGE_LICENSE "COPYING"
  !insertmacro MUI_PAGE_COMPONENTS
  !insertmacro MUI_PAGE_DIRECTORY
  !insertmacro MUI_PAGE_INSTFILES
  !insertmacro MUI_PAGE_FINISH

  !insertmacro MUI_UNPAGE_CONFIRM
  !insertmacro MUI_UNPAGE_INSTFILES

;--------------------------------

!insertmacro MUI_LANGUAGE "English"

; The stuff to install
Section "FRANz (required)" SecFRANz 

  SectionIn RO
  
  ; Set output path to the installation directory.
  SetOutPath $INSTDIR
  
  ; Put file there
  File "src\FRANz.exe"
  
  ; Write the installation path into the registry
  WriteRegStr HKLM SOFTWARE\FRANz "Install_Dir" "$INSTDIR"
  
  ; Write the uninstall keys for Windows
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\FRANz" "DisplayName" "FRANz"
  WriteRegStr HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\FRANz" "UninstallString" '"$INSTDIR\uninstall.exe"'
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\FRANz" "NoModify" 1
  WriteRegDWORD HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\FRANz" "NoRepair" 1
  WriteUninstaller "Uninstall.exe"
  
SectionEnd

; Optional section (can be disabled by the user)
Section "Start Menu Shortcuts"

  CreateDirectory "$SMPROGRAMS\FRANz"
  CreateShortCut "$SMPROGRAMS\FRANz\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0
  CreateShortCut "$SMPROGRAMS\FRANz\FRANz.lnk" "$INSTDIR\FRANz.exe" "" "$INSTDIR\FRANz.exe" 0
SectionEnd

; Optional section (can be disabled by the user)
Section "Examples"
  CreateDirectory "$INSTDIR\examples"
  CreateDirectory "$INSTDIR\examples\simpsons"
  CreateDirectory "$INSTDIR\examples\penaeus_monodon"
  CreateDirectory "$INSTDIR\examples\penaeus_monodon\input"
  SetOutPath "$INSTDIR\examples\simpsons"
  FILE "examples\simpsons\simpsons.dat"
  FILE "examples\simpsons\simpsons.mothers"
  FILE "examples\simpsons\futurama.dat"
  SetOutPath "$INSTDIR\examples\penaeus_monodon\input"
  FILE "examples\penaeus_monodon\input\penaeus_monodon.dat"
  FILE "examples\penaeus_monodon\input\penaeus_monodon.pedigree"
  SetOutPath "$INSTDIR\examples"
  FILE "examples\README"
  SetOutPath $INSTDIR
SectionEnd


Section "Extras (requires Perl)" SecExtras
  CreateDirectory "$INSTDIR\extras"
  CreateDirectory "$INSTDIR\extras\input"
  SetOutPath "$INSTDIR\extras\input"
  FILE "extras\input\csv.pl"
  FILE "extras\input\microsat.pl"
  FILE "extras\input\migrate.pl"
  SetOutPath $INSTDIR
SectionEnd

Section "Add to path" SecPath
${EnvVarUpdate} $0 "PATH" "A" "HKLM" $INSTDIR
SectionEnd
 
;--------------------------------
;Descriptions

  ;Language strings
  LangString DESC_SecFRANz ${LANG_ENGLISH} "The FRANz program and an uninstaller."
  LangString DESC_SecPath ${LANG_ENGLISH} "By adding the FRANz install dir to the path, you can type FRANz instead of C:\Program Files \FRANz\FRANz in the Command Prompt."
  LangString DESC_SecExtras ${LANG_ENGLISH} "Will install some Perl scripts for input file conversion. Requires Perl (available under http://www.activestate.com/activeperl/)"

  ;Assign language strings to sections
  !insertmacro MUI_FUNCTION_DESCRIPTION_BEGIN
    !insertmacro MUI_DESCRIPTION_TEXT ${SecFRANz} $(DESC_SecFRANz)
    !insertmacro MUI_DESCRIPTION_TEXT ${SecPath} $(DESC_SecPath)
    !insertmacro MUI_DESCRIPTION_TEXT ${SecExtras} $(DESC_SecExtras)
  !insertmacro MUI_FUNCTION_DESCRIPTION_END

;--------------------------------


Section "Uninstall"
  
  ; Remove registry keys
  DeleteRegKey HKLM "Software\Microsoft\Windows\CurrentVersion\Uninstall\FRANz"
  DeleteRegKey HKLM SOFTWARE\FRANz
  
  ; remove path
  ${un.EnvVarUpdate} $0 "PATH" "R" "HKLM" $INSTDIR 

  ; Remove files and uninstaller
  Delete $INSTDIR\FRANz.exe
  Delete $INSTDIR\Uninstall.exe
  Delete "$INSTDIR\examples\README"
  Delete "$INSTDIR\examples\simpsons\*.*"
  Delete "$INSTDIR\examples\penaeus_monodon\input\*.*"
  Delete "$INSTDIR\extras\input\*.*"

  ; Remove shortcuts, if any
  Delete "$SMPROGRAMS\FRANz\*.*"

  ; Remove directories used
  RMDir "$INSTDIR\extras\input"
  RMDir "$INSTDIR\extras"
  RMDir "$INSTDIR\examples\simpsons"
  RMDir "$INSTDIR\examples\penaeus_monodon\input"
  RMDir "$INSTDIR\examples\penaeus_monodon"
  RMDir "$INSTDIR\examples"
  RMDir "$SMPROGRAMS\FRANz"
  RMDir "$INSTDIR"

SectionEnd
