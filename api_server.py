from fastapi import FastAPI, Depends, APIRouter

from dependencies import verify_api_key
from schemas import (
    NameInput, SmilesInput, SmilesResponse, 
    NameResponse, MoleculeResponse
)
from services import (
    name_to_smiles, smiles_to_name, smiles_to_rdkit_mol, 
)


api_app = FastAPI(
    title=f"Chemdraw Tools API",
)

@api_app.get("/", tags=["Public"])
def root():
    return {"message": f"Welcome to Chemdraw Tools API!"}

@api_app.get("/health", tags=["Public"])
def health_check_api():
    """Provides a public health check for the API service."""
    return {"status": "ok", "service": "API Server"}

api_mcp_router = APIRouter(tags=["MCP-Compatible"], dependencies=[Depends(verify_api_key)])

@api_mcp_router.post("/name2smiles", response_model=SmilesResponse)
async def name_to_smiles_api(data: NameInput):
    return SmilesResponse(smiles=name_to_smiles(data.name))

@api_mcp_router.post("/smiles2name", response_model=NameResponse)
async def smiles_to_name_api(data: SmilesInput):
    return NameResponse(name=smiles_to_name(data.smiles))

api_router = APIRouter(tags=["ChemDraw Tools"], dependencies=[Depends(verify_api_key)])

@api_router.post("/smiles2rdkit", response_model=MoleculeResponse)
async def smiles_to_rdkit_api(data: SmilesInput):
    return MoleculeResponse(molecule=smiles_to_rdkit_mol(data.smiles))

api_app.include_router(api_mcp_router)
api_app.include_router(api_router)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run("api_server:api_app", host="0.0.0.0", port=1145, workers=1)
