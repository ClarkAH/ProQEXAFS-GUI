@echo=off
Set Shortcut=ProQEXAFS-GUI.lnk
echo set WshShell = WScript.CreateObject("WScript.Shell")>DecodeShortCut.vbs
echo set Lnk = WshShell.CreateShortcut(WScript.Arguments.Unnamed(0))>>DecodeShortCut.vbs
echo wscript.Echo Lnk.TargetPath>>DecodeShortCut.vbs
set vbscript=cscript //nologo DecodeShortCut.vbs
For /f "delims=" %%T in ( ' %vbscript% "%Shortcut%" ' ) do set target=%%T
del DecodeShortCut.vbs
Echo Shortcut %shortcut%
Echo Target   "%target%"

For %%A in ("%target%") do set ptt=%%~dpA

echo F| xcopy "%cdir%\ProQEXAFS-GUI\env-install\ProQEXAFS.bat" "%target%"
echo D| xcopy /E "%cdir%\ProQEXAFS-GUI" "%ptt%\main-pkg"