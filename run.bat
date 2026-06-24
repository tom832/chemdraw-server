@echo off
REM ============================================================================
REM ChemDraw Server launcher (multi-worker)
REM
REM   WORKERS      Number of uvicorn worker processes. Each worker is a full OS
REM                process with its own GIL and its own ChemScript/.NET instance,
REM                so CPU-bound throughput scales ~linearly with cores. Set it to
REM                roughly the CPU core count, leaving 1-2 cores for the OS / IO.
REM                More workers = more RAM (each loads ChemScript + RDKit).
REM   HOST / PORT  Listen address. Defaults to 0.0.0.0:1145 (see main_server.py).
REM
REM   PYTHONUTF8=1 Forces UTF-8 for stdout/stderr. ChemScript prints a copyright
REM                banner containing non-GBK characters (c) on import; without
REM                this, Chinese (GBK) Windows consoles raise UnicodeEncodeError
REM                in every (worker) process at startup.
REM
REM Logs:
REM   log\app_<pid>.log          Per-worker app log (loguru, WARNING+,
REM                              auto-rotated 20MB / 7 days / zip).
REM   log\native.log             Native ChemScript/.NET stdout+stderr. OVERWRITE
REM                              mode (single ">"): reset every start, never grows
REM                              unbounded. Replace with "> NUL 2>&1" to discard.
REM   log\prometheus_multiproc\  Per-worker Prometheus shards, aggregated by the
REM                              /metrics endpoint (cleaned on each start).
REM ============================================================================

cd /d "%~dp0"

set PYTHONUTF8=1
set PYTHONIOENCODING=utf-8
set WORKERS=4
REM set PORT=1145
REM set HOST=0.0.0.0

if not exist "log" mkdir "log"

uv run main_server.py > "log\native.log" 2>&1
