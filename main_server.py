# main_server.py

import uvicorn
from prometheus_fastapi_instrumentator import Instrumentator
from fastapi.middleware.cors import CORSMiddleware
from fastapi import FastAPI
from loguru_logging_intercept import setup_loguru_logging_intercept
from uvicorn.supervisors import Multiprocess, ChangeReload
import logging
from config import settings
from starlette.middleware.base import BaseHTTPMiddleware
from fastapi import Request, Response
from loguru import logger
import sys
import time

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

def run_uvicorn_with_metrics_filter(config: uvicorn.Config, force_exit=False):
    """自定义的 uvicorn 运行函数，在拦截设置后添加过滤器"""
    
    def _get_log_level(config: uvicorn.Config) -> int:
        if isinstance(config.log_level, str):
            return logging.getLevelName(config.log_level)
        else:
            return config.log_level or logging.INFO
    
    def _get_supervisor_type(config: uvicorn.Config):
        if config.should_reload:
            return ChangeReload
        if config.workers > 1:
            return Multiprocess
        return None
    
    def _run_server_with_intercept(**kwargs):
        # 设置 loguru 拦截
        setup_loguru_logging_intercept(
            level=_get_log_level(config), 
            modules=("uvicorn.error", "uvicorn.asgi", "uvicorn.access")
        )
        
        # 在拦截设置后，添加过滤器来过滤掉 /metrics 日志
        def filter_metrics(record):
            """过滤掉包含 /metrics 的日志消息"""
            message = record.get("message", "")
            return "/metrics" not in str(message)
        
        # 移除所有现有处理器
        logger.remove()
        # 重新添加带过滤器的处理器
        logger.add(
            sys.stderr,
            filter=filter_metrics,
            level="INFO",
            format="{time:YYYY-MM-DD HH:mm:ss.SSS} | {level} | {name}:{function}:{line} - {message}",
        )
        
        server = uvicorn.Server(config=config)
        server.force_exit = force_exit
        server.run(**kwargs)
    
    supervisor_type = _get_supervisor_type(config)
    if supervisor_type:
        sock = config.bind_socket()
        supervisor = supervisor_type(
            config, target=_run_server_with_intercept, sockets=[sock]
        )
        supervisor.run()
    else:
        _run_server_with_intercept()

def main():
    run_uvicorn_with_metrics_filter(
        uvicorn.Config(
            "main_server:root_app",
            host="0.0.0.0",
            port=1145,
            access_log=False,  # 禁用 uvicorn 自带的 access 日志
            # reload=True,
            # workers=1,
        )
    )

if __name__ == "__main__":
    
    main()