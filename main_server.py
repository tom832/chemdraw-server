# main_server.py

import uvicorn
from prometheus_fastapi_instrumentator import Instrumentator
from fastapi.middleware.cors import CORSMiddleware
from fastapi import FastAPI
from uvicorn_loguru_integration import run_uvicorn_loguru

# Import shared modules
from mcp_server import mcp_server
from api_server import api_app


# 从生成的 MCP 服务中获取可挂载的 ASGI 应用
mcp_asgi_app = mcp_server.http_app()


root_app = FastAPI(
    title="ChemDraw Unified Server",
    lifespan=mcp_asgi_app.lifespan,
    docs_url="/docs",
    root_path="/chemdraw"
)

# 对根应用进行统一监控
instrumentator = Instrumentator(
    excluded_handlers=["/metrics", "/health"]
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

def main():
    run_uvicorn_loguru(
        uvicorn.Config(
            "main_server:root_app",
            host="0.0.0.0",
            port=1145,
            # reload=True,
            # workers=1,
        )
    )

if __name__ == "__main__":
    
    main()