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
        # 尝试加载化学结构
        structure = cs.StructureData.LoadData(name, "name")

        # 检查是否成功加载
        if structure is None:
            raise ChemServiceError(f"无法识别化学名称: '{name}'。请检查名称格式是否正确，或尝试使用SMILES格式。")

        # 检查是否有 Smiles 属性
        if not hasattr(structure, 'Smiles'):
            raise ChemServiceError(f"加载的结构对象没有 Smiles 属性: '{name}'")

        smiles = structure.Smiles
        if not smiles:
            raise ChemServiceError(f"转换结果为空SMILES: '{name}'")

        return smiles

    except ChemServiceError:
        return ''
        
    except Exception as e:
        raise ChemServiceError(f"ChemScript 名称转SMILES失败: {e}")

def smiles_to_name(smiles: str) -> str:
    """使用 ChemScript 将 SMILES 转换为化学名称"""
    try:
        if not smiles:
            raise ChemServiceError("SMILES 字符串不能为空")

        # 尝试加载化学结构
        structure = cs.StructureData.LoadData(smiles, "smiles")

        # 检查是否成功加载
        if structure is None:
            raise ChemServiceError(f"无法解析SMILES: '{smiles}'。请检查SMILES格式是否正确。")

        # 检查是否有 ChemicalName 方法
        if not hasattr(structure, 'ChemicalName'):
            raise ChemServiceError(f"加载的结构对象没有 ChemicalName 方法: '{smiles}'")

        chemical_name = structure.ChemicalName()
        if not chemical_name:
            raise ChemServiceError(f"转换结果为空化学名称: '{smiles}'")

        return chemical_name

    except ChemServiceError:
        # 重新抛出我们自定义的异常
        raise
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

def compare_molecules(mol1_repr: str, mol2_repr: str) -> dict:
    """
    对比两个分子表示是否相同，并计算 Tanimoto 分数。
    要求：使用 ChemDraw LoadData，不指定类型，直接让 ChemScript 自行解析。
    相等判断基于 ChemScript 的 Inchi，Tanimoto 由 ChemScript 内置方法计算。
    """
    try:
        if not mol1_repr or not mol2_repr:
            raise ChemServiceError("分子表示文本不能为空")

        def _load_structure(repr_text: str, label: str):
            structure = cs.StructureData.LoadData(repr_text)
            if structure is None:
                raise ChemServiceError(f"无法解析第{label}个分子表示")
            if not hasattr(structure, "Inchi"):
                raise ChemServiceError(f"第{label}个分子缺少 Inchi 属性")
            return structure

        struct1 = _load_structure(mol1_repr, "一")
        struct2 = _load_structure(mol2_repr, "二")

        inchi1 = getattr(struct1, "Inchi", None)
        inchi2 = getattr(struct2, "Inchi", None)
        if not inchi1 or not inchi2:
            raise ChemServiceError("任一分子的 Inchi 为空，无法比较")

        is_equal = inchi1 == inchi2

        if not hasattr(struct1, "Tanimoto"):
            raise ChemServiceError("ChemScript 结构对象缺少 Tanimoto 方法")
        tanimoto = struct1.Tanimoto(struct2)
        try:
            tanimoto = round(float(tanimoto), 4)
        except Exception:
            raise ChemServiceError("Tanimoto 结果无法转换为浮点数")

        return {
            "is_equal": is_equal,
            "tanimoto": tanimoto,
        }

    except ChemServiceError:
        raise
    except Exception as e:
        raise ChemServiceError(f"分子比较失败: {e}")