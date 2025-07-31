# mcp_server.py
from fastapi import FastAPI
from fastmcp import FastMCP
from config import settings
from api_server import api_mcp_router

api_mcp_app = FastAPI(
    title="Chemdraw Tools API compatible with MCP",
)

api_mcp_app.include_router(api_mcp_router, prefix="/chemdraw")

mcp_server = FastMCP.from_fastapi(
    app=api_mcp_app,
    # mcp_names={
    #     "name_to_smiles_api": "name2smiles", 
    #     "smiles_to_name_api": "smiles2name",
    # },
    httpx_client_kwargs={
        "headers": {
            "Authorization": f"Bearer {settings.API_KEY}",
        }
    }
)


# # mcp_server.py
# from fastapi import Depends, HTTPException
# from fastmcp import FastMCP
# from schemas import NameInput, SmilesInput, SmilesResponse, NameResponse
# from services import name_to_smiles, smiles_to_name
# from config import settings
# from dependencies import verify_api_key


# # 创建 FastMCP 应用
# mcp_server = FastMCP(
#     "ChemDraw Server",
#     version=settings.PROJECT_VERSION,
# )

# # 定义工具函数
# @mcp_server.tool(name="name2smiles")  # 明确指定工具名称
# def name_to_smiles_tool(data: NameInput, _: bool = Depends(verify_api_key)) -> SmilesResponse:
#     """
#     Converts a chemical name to a SMILES string.
#     """
#     if _:
#         smiles_result = name_to_smiles(data.name)
#         return SmilesResponse(smiles=smiles_result)
#     else:
#         raise HTTPException(status_code=401, detail="Unauthorized")

# @mcp_server.tool(name="smiles2name")  # 明确指定工具名称
# def smiles_to_name_tool(data: SmilesInput, _: bool = Depends(verify_api_key)) -> NameResponse:
#     """
#     Converts a SMILES string to a chemical name.
#     """
#     if _:
#         name_result = smiles_to_name(data.smiles)
#         return NameResponse(name=name_result)
#     else:
#         raise HTTPException(status_code=401, detail="Unauthorized")

if __name__ == "__main__":
    mcp_server.run(transport="http", host="0.0.0.0", port=1147)
