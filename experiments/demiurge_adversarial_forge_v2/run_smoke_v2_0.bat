@echo off
setlocal EnableExtensions
cd /d "%~dp0"
chcp 65001 >nul
set PYTHONUTF8=1
set PYTHONUNBUFFERED=1
where py >nul 2>&1
if not errorlevel 1 (set "PY=py -3") else (set "PY=python")
%PY% -m pip install --disable-pip-version-check -r requirements_v2_0.txt
if errorlevel 1 goto end
%PY% -u self_test_v2_0.py
if errorlevel 1 goto end
%PY% -u demiurge_forge_v2.py --quiet --verify-expected
if errorlevel 1 goto end
%PY% -u download_sky_v2.py
if errorlevel 1 goto end
%PY% -u janus_cosmos_v2_0.py --smoke --nulls 32 --cal-nulls 16
:end
pause
