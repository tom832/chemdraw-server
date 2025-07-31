# schemas.py
from pydantic import BaseModel

class NameInput(BaseModel):
    name: str

class SmilesInput(BaseModel):
    smiles: str

class SmilesResponse(BaseModel):
    smiles: str

class NameResponse(BaseModel):
    name: str

class MoleculeResponse(BaseModel):
    smiles: str
    mol_pickle: str # base64 encoded pickle string