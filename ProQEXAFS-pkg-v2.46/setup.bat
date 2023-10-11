Set cdir=%cd%
CALL "%cdir%\ProQEXAFS-GUI\env-install\ProQEXAFS-2.46-Windows-x86_64.exe"
cd %userprofile%\desktop
CALL "%cdir%\ProQEXAFS-GUI\env-install\get-link.bat"
cd "%cdir%"