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

class MolCompareInput(BaseModel):
    mol1: str
    mol2: str

class MolCompareResponse(BaseModel):
    exact_match: bool
    tanimoto: float
    warning: str | None = None

class MolPair(BaseModel):
    """单个分子对"""
    mol1: str
    mol2: str
    pair_id: str | None = None  # 可选的标识符，用于标识每对分子

class BatchMolCompareInput(BaseModel):
    """批量分子比较输入"""
    pairs: list[MolPair]

class BatchMolCompareResult(BaseModel):
    """单个分子对的比较结果"""
    pair_id: str | None = None
    mol1: str
    mol2: str
    exact_match: bool
    tanimoto: float
    warning: str | None = None  # 警告信息（如手性不匹配或互变异构体）
    error: str | None = None  # 如果比较失败，记录错误信息

class BatchMolCompareResponse(BaseModel):
    """批量分子比较响应"""
    total: int  # 总对数
    success: int  # 成功比较的对数
    failed: int  # 失败的对数
    results: list[BatchMolCompareResult]

class MoleculeResponse(BaseModel):
    smiles: str
    mol_pickle: str # base64 encoded pickle string