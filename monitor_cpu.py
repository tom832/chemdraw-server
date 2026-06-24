# monitor_cpu.py
"""
Diagnostic: measure how many CPU cores the server process actually uses under
load. This tells us whether the throughput ceiling is a single-core/GIL limit
or simply running out of physical cores.

It launches loadtest.py as a child process and samples the CPU% of:
  - the uvicorn server process (the one LISTENING on the target port)
  - the load-test client process itself (to show same-machine contention)

Interpretation (Windows / psutil): 100% == one fully-busy logical core.
The max possible is 100% * logical_core_count.

Usage:
    uv run monitor_cpu.py --endpoint name2smiles -c 100 -d 15
"""

import argparse
import subprocess
import sys
import time

import psutil


def find_listen_pid(port: int):
    """Parse `netstat -ano` to find the PID LISTENING on the given port."""
    out = subprocess.run(["netstat", "-ano"], capture_output=True, text=True).stdout
    for line in out.splitlines():
        if f":{port} " in line and "LISTENING" in line.upper():
            parts = line.split()
            try:
                return int(parts[-1])
            except ValueError:
                continue
    return None


def stats(xs):
    return (min(xs), sum(xs) / len(xs), max(xs)) if xs else (0.0, 0.0, 0.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--port", type=int, default=1145)
    ap.add_argument("--endpoint", default="name2smiles")
    ap.add_argument("-c", "--concurrency", type=int, default=100)
    ap.add_argument("-d", "--duration", type=int, default=15)
    args = ap.parse_args()

    n_logical = psutil.cpu_count(logical=True)
    n_physical = psutil.cpu_count(logical=False)

    pid = find_listen_pid(args.port)
    if not pid:
        print(f"ERROR: no process LISTENING on :{args.port}. Is the server up?")
        sys.exit(1)
    server = psutil.Process(pid)
    print(f"Server PID         : {pid}")
    print(f"Logical cores      : {n_logical}  (100% = 1 core, "
          f"max possible = {n_logical * 100}%)")
    print(f"Physical cores     : {n_physical}")
    print(f"Load profile       : endpoint={args.endpoint} c={args.concurrency} "
          f"d={args.duration}s\n")

    client = subprocess.Popen(
        [sys.executable, "loadtest.py", "--endpoint", args.endpoint,
         "-c", str(args.concurrency), "-d", str(args.duration)],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
    )
    client_proc = psutil.Process(client.pid)

    # Prime cpu_percent (first call always returns 0.0).
    server.cpu_percent(None)
    client_proc.cpu_percent(None)
    psutil.cpu_percent(None)
    time.sleep(2.0)  # let the test ramp up

    srv_samples, cli_samples, sys_samples = [], [], []
    while client.poll() is None:
        srv_samples.append(server.cpu_percent(interval=0.5))
        try:
            cli_samples.append(client_proc.cpu_percent(None))
        except psutil.Error:
            pass
        sys_samples.append(psutil.cpu_percent(None))

    out, _ = client.communicate()

    s_min, s_avg, s_max = stats(srv_samples)
    c_min, c_avg, c_max = stats(cli_samples)
    _, sys_avg, sys_max = stats(sys_samples)

    print("---- load test result (tail) ----")
    tail = [ln for ln in out.strip().splitlines()
            if ln.strip()][-13:]
    print("\n".join(tail))

    print("\n---- CPU usage DURING load ----")
    print(f"Server process : avg={s_avg:6.0f}%  max={s_max:6.0f}%  "
          f"-> ~{s_avg / 100:.1f} cores avg, ~{s_max / 100:.1f} cores peak")
    print(f"Client process : avg={c_avg:6.0f}%  max={c_max:6.0f}%  "
          f"-> ~{c_avg / 100:.1f} cores avg")
    print(f"System total   : avg={sys_avg:6.1f}%  max={sys_max:6.1f}%  "
          f"(of all {n_logical} logical cores)")

    print("\n---- verdict ----")
    if s_avg < 130:
        print("Server uses ~1 core => single-core / GIL-bound. "
              "More threads won't help; multiprocessing would.")
    else:
        print(f"Server spreads across ~{s_avg / 100:.1f} cores => threads ARE "
              "running in parallel (GIL largely released during ChemScript "
              "calls). The ceiling is core count + same-machine contention, "
              "NOT the GIL.")


if __name__ == "__main__":
    main()
