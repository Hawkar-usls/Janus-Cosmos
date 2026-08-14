@echo off
setlocal EnableExtensions
cd /d "%~dp0"
chcp 65001 >nul
set PYTHONUTF8=1
set PYTHONUNBUFFERED=1

if not exist "results\canonical_v1" mkdir "results\canonical_v1"
where py >nul 2>&1
if %errorlevel%==0 goto use_py
where python >nul 2>&1
if %errorlevel%==0 goto use_python

echo [ERROR] Python 3 not found.
pause
exit /b 1

:use_py
py -3 -m pip install --disable-pip-version-check -r requirements.txt
if errorlevel 1 goto dependency_error
powershell.exe -NoProfile -ExecutionPolicy Bypass -Command "& { py -3 -u -m janus_cosmos.pipeline --manifest data/hst_blind_corpus.json --output-dir results/canonical_v1 --cache-dir .cache/janus_cosmos --nulls 1024 --seeds 20260810,20260811,20260812 2^>^&1 | Tee-Object -FilePath 'results/canonical_v1/terminal.log'; exit $LASTEXITCODE }"
goto done

:use_python
python -m pip install --disable-pip-version-check -r requirements.txt
if errorlevel 1 goto dependency_error
powershell.exe -NoProfile -ExecutionPolicy Bypass -Command "& { python -u -m janus_cosmos.pipeline --manifest data/hst_blind_corpus.json --output-dir results/canonical_v1 --cache-dir .cache/janus_cosmos --nulls 1024 --seeds 20260810,20260811,20260812 2^>^&1 | Tee-Object -FilePath 'results/canonical_v1/terminal.log'; exit $LASTEXITCODE }"
goto done

:dependency_error
echo [ERROR] Dependency installation failed.
pause
exit /b 1

:done
set RC=%errorlevel%
echo.
echo ============================================================
echo JANUS COSMOS FINISHED, EXIT CODE = %RC%
echo Receipt: results\canonical_v1\janus-cosmos-receipt.json
echo Events : results\canonical_v1\janus-cosmos-events.jsonl
echo Log    : results\canonical_v1\terminal.log
echo ============================================================
pause
exit /b %RC%
