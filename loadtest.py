# loadtest.py
"""
Async load tester for the ChemDraw server (REST API path).

It spins up N concurrent virtual clients (one persistent connection each via a
shared httpx connection pool) and hammers a chosen endpoint for a fixed
duration, then reports throughput, latency percentiles, status-code breakdown
and a sample of error bodies (handy for spotting ChemScript license failures).

Designed to answer: "how many concurrent connections / RPS can this hold?"

Usage (run inside the project venv):
    uv run loadtest.py --endpoint health -c 1000 -d 15
    uv run loadtest.py --endpoint name2smiles -c 1000 -d 20
    uv run loadtest.py --endpoint name2smiles -c 1000 -d 30 --think 1.0   # 1000 idle-ish conns

The script auto-detects whether the API lives under "/api" (direct uvicorn)
or "/chemdraw/api" (behind a reverse proxy), and reads the API key from
config.settings unless --api-key is given.
"""

import argparse
import asyncio
import random
import time
from collections import Counter

import httpx

# ---- Test payload pools (kept small & valid for ChemScript) ------------------
NAMES = [
    "aspirin", "benzene", "caffeine", "glucose", "ethanol", "toluene",
    "acetone", "phenol", "aniline", "pyridine", "naphthalene", "methane",
]
SMILES = [
    "CCO", "c1ccccc1", "CC(=O)O", "CCN", "C1CCCCC1", "O", "CO", "CCCC",
    "c1ccncc1", "CC(=O)Oc1ccccc1C(=O)O",
]

# (method, path, payload_factory) per endpoint. path is relative to the prefix.
ENDPOINTS = {
    "health":      ("GET",  "/api/health",      lambda: None),
    "name2smiles": ("POST", "/api/name2smiles", lambda: {"name": random.choice(NAMES)}),
    "smiles2name": ("POST", "/api/smiles2name", lambda: {"smiles": random.choice(SMILES)}),
    "mol_compare": ("POST", "/api/mol_compare", lambda: {
        "mol1": random.choice(SMILES), "mol2": random.choice(SMILES)
    }),
}


class Stats:
    def __init__(self):
        self.latencies = []          # ms, all completed requests
        self.status = Counter()      # http status code -> count
        self.errors = Counter()      # exception class name -> count
        self.sample_errors = []      # up to N (status, body[:200]) samples
        self.total = 0
        self.ok = 0
        self.failed = 0
        self.inflight = 0
        self.max_inflight = 0


async def detect_prefix(client: httpx.AsyncClient, base: str) -> str:
    """Return the URL prefix ('' or '/chemdraw') that serves /api/health."""
    for prefix in ("", "/chemdraw"):
        try:
            r = await client.get(f"{base}{prefix}/api/health", timeout=5.0)
            if r.status_code == 200:
                return prefix
        except Exception:
            pass
    raise RuntimeError(
        f"Could not reach /api/health on {base} (tried '' and '/chemdraw'). "
        "Is the server running?"
    )


async def worker(client, base, prefix, endpoint, headers, think, stop_at, stats):
    method, path, payload_factory = ENDPOINTS[endpoint]
    url = f"{base}{prefix}{path}"
    while time.perf_counter() < stop_at:
        body = payload_factory()
        stats.inflight += 1
        if stats.inflight > stats.max_inflight:
            stats.max_inflight = stats.inflight
        t0 = time.perf_counter()
        try:
            r = await client.request(method, url, json=body, headers=headers)
            dt = (time.perf_counter() - t0) * 1000.0
            stats.latencies.append(dt)
            stats.total += 1
            stats.status[r.status_code] += 1
            if 200 <= r.status_code < 300:
                stats.ok += 1
            else:
                stats.failed += 1
                if len(stats.sample_errors) < 5:
                    stats.sample_errors.append((r.status_code, r.text[:200]))
        except Exception as e:  # noqa: BLE001 - we want to bucket every failure
            dt = (time.perf_counter() - t0) * 1000.0
            stats.latencies.append(dt)
            stats.total += 1
            stats.failed += 1
            stats.errors[type(e).__name__] += 1
        finally:
            stats.inflight -= 1
        if think > 0:
            await asyncio.sleep(think)


async def progress_printer(stats: Stats, stop_at: float):
    last_total = 0
    last_t = time.perf_counter()
    try:
        while True:
            await asyncio.sleep(2.0)
            now = time.perf_counter()
            cur = stats.total
            rps = (cur - last_total) / (now - last_t)
            remaining = max(0.0, stop_at - now)
            print(f"  [t-{remaining:5.1f}s] total={cur:>8}  ok={stats.ok:>8}  "
                  f"failed={stats.failed:>6}  ~{rps:8.1f} req/s  inflight={stats.inflight}")
            last_total, last_t = cur, now
    except asyncio.CancelledError:
        pass


def pct(sorted_ms, p):
    if not sorted_ms:
        return 0.0
    k = int(round((len(sorted_ms) - 1) * p))
    return sorted_ms[k]


def report(stats: Stats, cfg, elapsed: float):
    lat = sorted(stats.latencies)
    print("\n==== Load Test Report ====")
    print(f"Endpoint     : {cfg.endpoint}")
    print(f"Target       : {cfg.base}{cfg.prefix}{ENDPOINTS[cfg.endpoint][1]}")
    print(f"Concurrency  : {cfg.concurrency}")
    print(f"Think time   : {cfg.think}s")
    print(f"Duration     : requested {cfg.duration}s, actual {elapsed:.1f}s")
    print(f"Max in-flight: {stats.max_inflight}")
    print("-" * 40)
    print(f"Requests     : total={stats.total}  ok={stats.ok}  failed={stats.failed}")
    if elapsed > 0:
        print(f"Throughput   : {stats.total / elapsed:.1f} req/s")
    if stats.total:
        print(f"Error rate   : {100.0 * stats.failed / stats.total:.2f}%")
    print(f"Status codes : {dict(stats.status)}")
    if stats.errors:
        print(f"Exceptions   : {dict(stats.errors)}")
    print("-" * 40)
    if lat:
        mean = sum(lat) / len(lat)
        print("Latency (ms) :")
        print(f"  min={lat[0]:.1f}  mean={mean:.1f}  p50={pct(lat,0.50):.1f}  "
              f"p90={pct(lat,0.90):.1f}  p95={pct(lat,0.95):.1f}  "
              f"p99={pct(lat,0.99):.1f}  max={lat[-1]:.1f}")
    if stats.sample_errors:
        print("-" * 40)
        print("Sample error responses (first 5):")
        for code, text in stats.sample_errors:
            print(f"  [{code}] {text!r}")
    print("=" * 40)


async def run(cfg):
    limits = httpx.Limits(
        max_connections=cfg.concurrency + 10,
        max_keepalive_connections=cfg.concurrency + 10,
    )
    timeout = httpx.Timeout(cfg.timeout)
    headers = {"Authorization": f"Bearer {cfg.api_key}"} if cfg.api_key else {}

    async with httpx.AsyncClient(limits=limits, timeout=timeout) as client:
        cfg.prefix = await detect_prefix(client, cfg.base)
        print(f"Detected API prefix: '{cfg.prefix or '(none)'}' "
              f"-> {cfg.base}{cfg.prefix}/api/...")

        # Pre-flight one real request to surface auth / license problems early.
        if cfg.endpoint != "health":
            method, path, factory = ENDPOINTS[cfg.endpoint]
            pre = await client.request(method, f"{cfg.base}{cfg.prefix}{path}",
                                       json=factory(), headers=headers)
            print(f"Pre-flight {cfg.endpoint}: HTTP {pre.status_code} "
                  f"{pre.text[:160]!r}")
            if pre.status_code == 401:
                print("\nAUTH FAILED (401). Set --api-key or check .env API_KEY. Aborting.")
                return

        stats = Stats()
        print(f"\nStarting: {cfg.concurrency} workers x {cfg.duration}s "
              f"(think={cfg.think}s) ...")
        start = time.perf_counter()
        stop_at = start + cfg.duration
        progress = asyncio.create_task(progress_printer(stats, stop_at))
        workers = [
            asyncio.create_task(
                worker(client, cfg.base, cfg.prefix, cfg.endpoint,
                       headers, cfg.think, stop_at, stats)
            )
            for _ in range(cfg.concurrency)
        ]
        await asyncio.gather(*workers)
        elapsed = time.perf_counter() - start
        progress.cancel()
        report(stats, cfg, elapsed)


def parse_args():
    p = argparse.ArgumentParser(description="Async load tester for ChemDraw server")
    p.add_argument("--base", default="http://127.0.0.1:1145",
                   help="Base URL (default: http://127.0.0.1:1145)")
    p.add_argument("--endpoint", default="health", choices=list(ENDPOINTS),
                   help="Endpoint to hit (default: health)")
    p.add_argument("-c", "--concurrency", type=int, default=100,
                   help="Number of concurrent virtual clients (default: 100)")
    p.add_argument("-d", "--duration", type=float, default=15.0,
                   help="Test duration in seconds (default: 15)")
    p.add_argument("--think", type=float, default=0.0,
                   help="Sleep (s) after each request to simulate idle "
                        "connections (default: 0 = max load)")
    p.add_argument("--timeout", type=float, default=30.0,
                   help="Per-request timeout in seconds (default: 30)")
    p.add_argument("--api-key", default=None,
                   help="Bearer API key (default: read from config.settings)")
    cfg = p.parse_args()

    if cfg.api_key is None:
        try:
            from config import settings
            cfg.api_key = settings.API_KEY
        except Exception:
            cfg.api_key = None
    cfg.prefix = ""
    return cfg


if __name__ == "__main__":
    cfg = parse_args()
    try:
        asyncio.run(run(cfg))
    except KeyboardInterrupt:
        print("\nInterrupted.")
