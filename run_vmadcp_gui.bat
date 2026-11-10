@echo off
REM Sruth VMADCP Processor launcher. Optional: drag a cruise folder onto this file.
cd /d "%~dp0"
set PY=py
where py >nul 2>nul || set PY=python
%PY% -c "import numpy, scipy, matplotlib" 2>nul || %PY% -m pip install -r requirements.txt
%PY% sruth_vmadcp_gui.py %1
if errorlevel 1 pause
