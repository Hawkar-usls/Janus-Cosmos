@echo off
setlocal EnableExtensions
cd /d "%~dp0"
chcp 65001 >nul
set PYTHONUTF8=1
set PYTHONUNBUFFERED=1
set JANUS_COSMOS_WORKERS=10
set OMP_NUM_THREADS=1
set OPENBLAS_NUM_THREADS=1
set MKL_NUM_THREADS=1
set NUMEXPR_NUM_THREADS=1
title JANUS COSMOS v2.1.1 - PARALLEL SPECIFICITY REPAIR

echo ============================================================
echo JANUS COSMOS v2.1.1 - PARALLEL DETECTOR SPECIFICITY REPAIR
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
echo [2/5] Verifying negative certificate and frozen v2.1 protocol...
%PY% -u self_test_v2_1.py
if errorlevel 1 goto fail

echo.
echo [3/5] Re-forging unchanged parent detector identity...
%PY% -u demiurge_forge_v2.py --quiet --verify-expected
if errorlevel 1 goto forgefail

echo.
echo [4/5] Downloading 168 frozen source products, up to 10 concurrently...
echo       Orion: 4 target bands + 20 fresh fields x 4 matching bands
echo       HST: NGC1425 WF3 science/weight pairs
echo       HST controls: 20 SGAL fields x 2 bands x science/weight
echo       Completed files are cached; interrupted downloads can be resumed.
%PY% -u download_sky_v2_1.py --workers 10
if errorlevel 1 goto downloadfail

echo.
echo [5/5] Starting real-field detector-specificity evaluation...
echo       Synthetic nulls remain diagnostics; admission uses real controls.
echo       10 independent fields run concurrently; frozen result order is restored.
echo       Existing v2.1.0 model checkpoints are reused when hashes match.
echo       Live progress, ETA and atomic partial reports are enabled.
echo.
%PY% -u janus_cosmos_v2_1.py --workers 10
set "RC=%errorlevel%"

echo.
echo ============================================================
echo FINISHED, EXIT CODE = %RC%
echo Report : results_v2_1\janus-cosmos-v2.1-report.json
echo Summary: results_v2_1\SUMMARY_v2.1.txt
echo Events : results_v2_1\janus-cosmos-v2.1-events.jsonl
echo Log    : results_v2_1\terminal_v2.1.log
echo Progress: results_v2_1\progress_v2.1.json
echo ============================================================
pause
exit /b %RC%

:forgefail
echo.
echo [ERROR] The unchanged parent detector did not reproduce its portable identity.
echo v2.1 is blocked to prevent target-informed detector drift.
pause
exit /b 3

:downloadfail
echo.
echo [ERROR] One or more frozen FITS downloads failed.
echo Re-run this BAT. Valid completed downloads are cached.
echo See external_data\download_errors_v2_1.json
pause
exit /b 2

:fail
echo.
echo [ERROR] Dependency, certificate, protocol, or self-test failure.
echo Sky analysis was NOT started.
pause
exit /b 1
