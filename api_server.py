from fastapi import FastAPI, Depends, APIRouter, Request, Response
from fastapi.responses import JSONResponse
from starlette.middleware.base import BaseHTTPMiddleware
from config import settings
from dependencies import verify_api_key
from schemas import (
    NameInput, SmilesInput, SmilesResponse,
    NameResponse, MoleculeResponse,
    MolCompareInput, MolCompareResponse,
    BatchMolCompareInput, BatchMolCompareResponse,
)
from services import (
    name_to_smiles, smiles_to_name, smiles_to_rdkit_mol,
    compare_molecules, batch_compare_molecules,
    ChemServiceError,
)


DOCS_ACCESS_TOKEN = settings.DOCS_ACCESS_TOKEN

class DocsProtectionMiddleware(BaseHTTPMiddleware):
    """保护文档访问的中间件"""
    
    def __init__(self, app, docs_token: str = None):
        super().__init__(app)
        self.docs_token = docs_token
    
    async def dispatch(self, request: Request, call_next):
        # 检查是否是访问docs相关路径
        if request.url.path in ["/chemdraw/api/docs", "/chemdraw/api/redoc"] or request.url.path.startswith("/chemdraw/api/openapi"):
            # 如果设置了docs token，则需要验证
            if self.docs_token:
                # 从查询参数中获取token
                token = request.query_params.get("token")
                
                # 对于openapi.json请求，检查Referer头部是否包含正确的token
                if request.url.path.startswith("/chemdraw/api/openapi"):
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

api_app = FastAPI(
    title=f"Chemdraw Tools API",
    docs_url=None if not DOCS_ACCESS_TOKEN else "/docs",
    redoc_url=None if not DOCS_ACCESS_TOKEN else "/redoc",
)

if DOCS_ACCESS_TOKEN:
    api_app.add_middleware(DocsProtectionMiddleware, docs_token=DOCS_ACCESS_TOKEN)


@api_app.exception_handler(ChemServiceError)
async def chem_service_error_handler(request: Request, exc: ChemServiceError):
    """ChemScript 无法命名/解析等业务失败属于输入问题，返回 422 而非 500。

    这样可避免把每个无法处理的分子记成 500 + 海量 diagnose 堆栈，
    既止住日志膨胀，也不向调用方暴露内部堆栈。
    """
    return JSONResponse(status_code=422, content={"detail": str(exc)})


@api_app.get("/", tags=["Public"])
def root():
    return {"message": f"Welcome to Chemdraw Tools API!"}

@api_app.get("/health", tags=["Public"])
def health_check_api():
    """Provides a public health check for the API service."""
    return {"status": "ok", "service": "API Server"}

api_mcp_router = APIRouter(tags=["MCP-Compatible"], dependencies=[Depends(verify_api_key)])

@api_mcp_router.post("/name2smiles", response_model=SmilesResponse)
def name_to_smiles_api(data: NameInput):
    return SmilesResponse(smiles=name_to_smiles(data.name))

@api_mcp_router.post("/smiles2name", response_model=NameResponse)
def smiles_to_name_api(data: SmilesInput):
    return NameResponse(name=smiles_to_name(data.smiles))

@api_mcp_router.post("/mol_compare", response_model=MolCompareResponse)
def mol_compare_api(data: MolCompareInput):
    result = compare_molecules(data.mol1, data.mol2)
    return MolCompareResponse(**result)

@api_mcp_router.post("/batch_mol_compare", response_model=BatchMolCompareResponse)
def batch_mol_compare_api(data: BatchMolCompareInput):
    """
    批量分子比较接口
    
    可以一次性比较多对分子，每对分子包含：
    - mol1: 第一个分子的表示（SMILES、化学名称等）
    - mol2: 第二个分子的表示（SMILES、化学名称等）
    - pair_id: 可选的标识符，用于标识该分子对
    
    返回结果包含每对分子的比较结果，即使某些对比较失败也会继续处理其他对。
    """
    pairs = [pair.model_dump() for pair in data.pairs]
    result = batch_compare_molecules(pairs)
    return BatchMolCompareResponse(**result)

api_router = APIRouter(tags=["ChemDraw Tools"], dependencies=[Depends(verify_api_key)])

@api_router.post("/smiles2rdkit", response_model=MoleculeResponse)
def smiles_to_rdkit_api(data: SmilesInput):
    return MoleculeResponse(molecule=smiles_to_rdkit_mol(data.smiles))

api_app.include_router(api_mcp_router)
api_app.include_router(api_router)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("api_server:api_app", host="0.0.0.0", port=1145, workers=1)
