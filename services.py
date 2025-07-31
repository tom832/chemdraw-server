# services.py
import pickle
import base64
import ChemScript as cs
from rdkit import Chem

# 定义自定义异常，以便在 API 层捕获
class ChemServiceError(Exception):
    pass

def name_to_smiles(name: str) -> str:
    """使用 ChemScript 将化学名称转换为 SMILES"""
    try:
        return cs.StructureData.LoadData(name, "name").Smiles
    except Exception as e:
        raise ChemServiceError(f"ChemScript 名称转SMILES失败: {e}")

def smiles_to_name(smiles: str) -> str:
    """使用 ChemScript 将 SMILES 转换为化学名称"""
    try:
        return cs.StructureData.LoadData(smiles, "smiles").ChemicalName()
    except Exception as e:
        raise ChemServiceError(f"ChemScript SMILES转名称失败: {e}")

def smiles_to_rdkit_mol(smiles: str) -> Chem.Mol:
    """将 SMILES 转换为 RDKit 分子对象"""
    try:
        if not smiles:
            raise ValueError("SMILES 字符串不能为空")
            
        csmol = cs.StructureData.LoadData(smiles, "smiles")
        rdmol = Chem.RWMol()
        atom_map = {}
        
        for idx, atom in enumerate(csmol.Atoms):
            atom_map[atom.Name] = idx
            rdatom = Chem.Atom(atom.Element)
            rdatom.SetFormalCharge(atom.FormalCharge)
            rdmol.AddAtom(rdatom)
        
        bond_map = {
            1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
            3: Chem.BondType.TRIPLE, 6: Chem.BondType.AROMATIC
        }
        
        for bond in csmol.Bonds:
            atom1_idx = atom_map[bond.Atom1.Name]
            atom2_idx = atom_map[bond.Atom2.Name]
            bond_type = bond_map.get(bond.Order.BondCode, Chem.BondType.SINGLE)
            rdmol.AddBond(atom1_idx, atom2_idx, bond_type)
        
        mol = rdmol.GetMol()
        conf = Chem.Conformer(mol.GetNumAtoms())
        for atom in csmol.Atoms:
            pos = atom.GetCartesian()
            conf.SetAtomPosition(atom_map[atom.Name], (pos.x, pos.y, pos.z))
        mol.AddConformer(conf)

        Chem.SanitizeMol(mol)
        Chem.Kekulize(mol)
        return mol
    except Exception as e:
        raise ChemServiceError(f"RDKit 转换失败: {e}")

def serialize_molecule(mol: Chem.Mol) -> str:
    """将 RDKit 分子对象序列化为 base64 字符串"""
    return base64.b64encode(pickle.dumps(mol)).decode('utf-8')