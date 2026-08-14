@echo off
setlocal EnableExtensions
cd /d "%~dp0"
chcp 65001 >nul
set PYTHONUTF8=1
set PYTHONUNBUFFERED=1
title JANUS COSMOS v2.0.2 - DEMIURGE ADVERSARIAL FORGE

echo ============================================================
echo JANUS COSMOS v2.0.2 - DEMIURGE ADVERSARIAL COSMOS CHECK
echo ============================================================
echo Working directory: %CD%
echo.

where py >nul 2>&1
if not errorlevel 1 (
  set "PY=py -3"
) else (
  where python >nul 2>&1
  if errorlevel 1 (
    echo [ERROR] Python 3 not found.
    pause
    exit /b 1
  )
  set "PY=python"
)

echo [1/5] Checking / installing dependencies...
%PY% -m pip install --disable-pip-version-check -r requirements_v2_0.txt
if errorlevel 1 goto fail

echo.
echo [2/5] Offline scientific + integrity self-test...
%PY% -u self_test_v2_0.py
if errorlevel 1 goto fail

echo.
echo [3/5] Re-forging detector and verifying deterministic freeze hash...
%PY% -u demiurge_forge_v2.py --quiet --verify-expected
if errorlevel 1 goto forgefail

echo.
echo [4/5] Downloading 22 fixed astronomy FITS products...
echo       Orion: DSS2 Red/Blue + 2MASS J/K
echo       Blind controls: 4 fields x 4 bands
echo       NGC1425: HST WFPC2 F555W + F814W
%PY% -u download_sky_v2.py
if errorlevel 1 goto downloadfail

echo.
echo [5/5] Starting powered frozen-detector sky evaluation...
echo       768 test nulls/model, 96 calibration nulls, 3 seed chunks.
echo       Checkpoints are saved. If interrupted, run this BAT again.
echo.
%PY% -u janus_cosmos_v2_0.py
set "RC=%errorlevel%"

echo.
echo ============================================================
echo FINISHED, EXIT CODE = %RC%
echo Report : results_v2_0\janus-cosmos-v2.0-report.json
echo Summary: results_v2_0\SUMMARY_v2.0.txt
echo Events : results_v2_0\janus-cosmos-v2.0-events.jsonl
echo Log    : results_v2_0\terminal.log
echo ============================================================
pause
exit /b %RC%

:forgefail
echo.
echo [ERROR] The deterministic forge did not reproduce the packaged frozen hash.
echo The sky run is blocked to prevent target-informed detector drift.
pause
exit /b 3

:downloadfail
echo.
echo [ERROR] One or more astronomy FITS downloads failed.
echo Re-run this BAT. Valid completed downloads are cached.
echo See external_data\download_errors_v2_0.json
pause
exit /b 2

:fail
echo.
echo [ERROR] Dependency or self-test failure. Sky analysis was NOT started.
pause
exit /b 1
