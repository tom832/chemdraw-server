# main_server.py

import os
import multiprocessing

# Prometheus multiprocess mode MUST be configured before prometheus_client is
# imported (it is pulled in by the Instrumentator below). In multi-worker mode
# every worker records metrics into this directory and the /metrics endpoint
# aggregates them. The directory must already exist or the Instrumentator
# constructor raises ValueError.
os.environ.setdefault(
    "PROMETHEUS_MULTIPROC_DIR", os.path.join("log", "prometheus_multiproc")
)
_PROM_DIR = os.environ["PROMETHEUS_MULTIPROC_DIR"]
os.makedirs(_PROM_DIR, exist_ok=True)

# Master process only: clear stale metric shards from previous runs BEFORE any
# metric object is created. Worker processes are spawned (parent_process() is
# not None) and must NOT wipe each other's shards.
if multiprocessing.parent_process() is None:
    import glob

    for _stale in glob.glob(os.path.join(_PROM_DIR, "*")):
        try:
            os.remove(_stale)
        except OSError:
            pass

import functools
import logging

import uvicorn
from prometheus_fastapi_instrumentator import Instrumentator
from fastapi.middleware.cors import CORSMiddleware
from fastapi import FastAPI
from loguru_logging_intercept import setup_loguru_logging_intercept
from uvicorn.supervisors import Multiprocess
from config import settings
from starlette.middleware.base import BaseHTTPMiddleware
from fastapi import Request, Response
from loguru import logger

# Import shared modules
from mcp_server import mcp_server
from api_server import api_app


# 从生成的 MCP 服务中获取可挂载的 ASGI 应用
mcp_asgi_app = mcp_server.http_app()

DOCS_ACCESS_TOKEN = settings.DOCS_ACCESS_TOKEN

class DocsProtectionMiddleware(BaseHTTPMiddleware):
    """保护文档访问的中间件"""

    def __init__(self, app, docs_token: str = None):
        super().__init__(app)
        self.docs_token = docs_token

    async def dispatch(self, request: Request, call_next):
        # 检查是否是访问docs相关路径
        if request.url.path in ["/chemdraw/docs", "/chemdraw/redoc"] or request.url.path.startswith("/chemdraw/openapi") or request.url.path.startswith("/chemdraw/api/openapi"):
            # 如果设置了docs token，则需要验证
            if self.docs_token:
                # 从查询参数中获取token
                token = request.query_params.get("token")

                # 对于openapi.json请求，检查Referer头部是否包含正确的token
                if request.url.path.startswith("/chemdraw/openapi") or request.url.path.startswith("/chemdraw/api/openapi"):
                    referer = request.headers.get("referer", "")
                    if f"token={self.docs_token}" in referer:
                        # 如果Referer包含正确的token，允许访问
                        pass
                    elif token != self.docs_token:
                        return Response(
                            content="Unauthorized access to docs. Please provide token as query parameter or ensure referer contains the token.",
                            status_code=401,
                            media_type="text/plain; charset=utf-8"
                        )
                elif token != self.docs_token:
                    # 对于docs和redoc页面，直接检查token参数
                    return Response(
                        content="Unauthorized access to docs",
                        status_code=401,
                        media_type="text/plain; charset=utf-8"
                    )

        # 继续处理请求
        response = await call_next(request)
        return response

root_app = FastAPI(
    title="ChemDraw Unified Server",
    lifespan=mcp_asgi_app.lifespan,
    docs_url="/docs",
    root_path="/chemdraw"
)

if DOCS_ACCESS_TOKEN:
    root_app.add_middleware(DocsProtectionMiddleware, docs_token=DOCS_ACCESS_TOKEN)

# 对根应用进行统一监控
instrumentator = Instrumentator(
    excluded_handlers=["/metrics", "/health", "/mcp"]
).instrument(root_app).expose(root_app)

# 挂载两个子应用
root_app.mount("/api", api_app)
root_app.mount("/", mcp_asgi_app)


# 在根应用上添加 CORS
root_app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

def _get_log_level(config: uvicorn.Config) -> int:
    if isinstance(config.log_level, str):
        return logging.getLevelName(config.log_level)
    return config.log_level or logging.INFO


def _filter_metrics(record):
    """Drop log lines that mention /metrics to keep the logs clean."""
    return "/metrics" not in str(record.get("message", ""))


def _server_target(config: uvicorn.Config, sockets=None):
    """Run a uvicorn server. Executed in EACH worker process (and in the
    single-process path).

    This MUST be a module-level function (not a closure) so it can be pickled
    and shipped to worker processes under the Windows 'spawn' start method.
    """
    # loguru interception must be (re)configured inside every worker process.
    setup_loguru_logging_intercept(
        level=_get_log_level(config),
        modules=("uvicorn.error", "uvicorn.asgi", "uvicorn.access"),
    )
    logger.remove()
    # Per-worker file: avoids multi-process rotation/compression races that
    # would otherwise occur when several workers write one shared log file.
    logger.add(
        f"log/app_{os.getpid()}.log",
        filter=_filter_metrics,
        level="WARNING",
        rotation="20 MB",
        retention="7 days",
        compression="zip",
        enqueue=True,
        format="{time:YYYY-MM-DD HH:mm:ss.SSS} | {level} | {name}:{function}:{line} - {message}",
    )
    server = uvicorn.Server(config=config)
    server.run(sockets=sockets)


def run_uvicorn_with_metrics_filter(config: uvicorn.Config):
    """Start uvicorn, wiring loguru interception into each worker.

    For workers > 1 the listening socket is bound once in the parent and handed
    to the Multiprocess supervisor. The per-worker target is a picklable
    functools.partial so it survives the Windows 'spawn' start method.
    """
    if config.workers and config.workers > 1:
        sock = config.bind_socket()
        target = functools.partial(_server_target, config)
        supervisor = Multiprocess(config, target=target, sockets=[sock])
        supervisor.run()
    else:
        _server_target(config)


def main():
    # Worker count: WORKERS env var overrides; default = CPU count.
    workers = int(os.environ.get("WORKERS", "0")) or (os.cpu_count() or 1)
    host = os.environ.get("HOST", "0.0.0.0")
    port = int(os.environ.get("PORT", "1145"))
    run_uvicorn_with_metrics_filter(
        uvicorn.Config(
            "main_server:root_app",
            host=host,
            port=port,
            access_log=False,  # uvicorn's own access log stays disabled
            workers=workers,
        )
    )

if __name__ == "__main__":
    
    main()