@echo off
setlocal

set "SETVARS=%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat"
if not exist "%SETVARS%" set "SETVARS=%ProgramFiles%\Intel\oneAPI\setvars.bat"

if not exist "%SETVARS%" (
    echo ERROR: Intel oneAPI setvars.bat was not found.
    exit /b 1
)

call "%SETVARS%" intel64 >nul
if errorlevel 1 exit /b %errorlevel%

powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0scripts\release.ps1" %*
exit /b %errorlevel%
