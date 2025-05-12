# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 10:35:10 2022

@author: yk

"""
import pandas as pd
import rdkit.Chem as rkc
from rdkit.Chem import rdchem
from rdkit.Chem import RWMol
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import PandasTools
from rdkit.Chem import Descriptors
import re

#mol = rkc.MolFromSmiles('[*:0]c1ccc2ccnc(NC3CCNCC3)c2c1')

ATTACHMENT_POINT_TOKEN = "*"
ATTACHMENT_POINT_NUM_REGEXP = r"\[{}:(\d+)\]".format(re.escape(ATTACHMENT_POINT_TOKEN))
ATTACHMENT_POINT_REGEXP = r"(?:{0}|\[{0}[^\]]*\])".format(re.escape(ATTACHMENT_POINT_TOKEN))
ATTACHMENT_POINT_NO_BRACKETS_REGEXP = r"(?<!\[){}".format(re.escape(ATTACHMENT_POINT_TOKEN))
def remove_attachment_point_numbers(mol_or_smi):
    """
    Removes the numbers for the attachment points throughout the molecule.
    :param mol_or_smi: SMILES string or mol object to convert.
    :return : A converted SMILES string.
    """
    if isinstance(mol_or_smi, rkc.Mol):
        for atom in mol_or_smi.GetAtoms():
            atom.ClearProp("molAtomMapNumber")
        return mol_or_smi
    return re.sub(ATTACHMENT_POINT_NUM_REGEXP, "[{}]".format(ATTACHMENT_POINT_TOKEN), mol_or_smi)

def replace_attachment(mol):
    if mol is not None:
        mw = rkc.RWMol(mol)
        atta_idx = []
        for i in range(0,mw.GetNumAtoms()):
            if mw.GetAtomWithIdx(i).GetSymbol() == '*':
                atta_idx.append(i)
        for i in reversed(atta_idx):
            mw.GetAtomWithIdx(i).SetAtomicNum(1)
        writer.write(mw)
        return mw
    else:
        return '[xe]'                

suppl= rkc.SDMolSupplier(r'generated_molecules_nik.sdf')
writer = rkc.SDWriter(r'generated_molecules_nik_add_H.sdf')
for mol in suppl:   
#    remove_attachment_point_numbers(mol)
    replace_attachment(mol)
writer.close()
print("输出完成，请到输出路径中查看")

'''
writer.SetProps(['LOGP', 'MW'])
for i, mol in enumerate(mols):
    mw = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    mol.SetProp('MW', '%.2f' %(mw))
    mol.SetProp('LOGP', '%.2f' %(logp))
    mol.SetProp('_Name', 'No_%s' %(i))
#获取分子中的原子数目

atom_num = mol.GetNumAtoms()
#获取分子中的键数目
bond_num = mol.GetNumBonds()
nei_atom = []
nei_bond = []
#获取分子中的原子的相邻原子的数目以及序号
for i in range(atom_num):
  nei_atom.append([(i.GetSmarts(),i.GetIdx()) for i in mol.GetAtomWithIdx(i).GetNeighbors()])
#获取分子中键的相关类型以及键的特征，以及成键原子序号
for i in range(bond_num):
  nei_bond.append((mol.GetBondWithIdx(i).GetBondType().name,mol.GetBondWithIdx(i).GetBeginAtomIdx(),mol.GetBondWithIdx(i).GetEndAtomIdx()))


m2 = Chem.AddHs(m)
print("m Smiles:",Chem.MolToSmiles(m))
print("m2 Smiles:",Chem.MolToSmiles(m2))
print(Draw.ShowMol(m))
print(Chem.rdchem.SubstanceGroup.GetAttachPoints(m))
mol=Chem.SDMolSupplier(generated_molecule1.sdf')
print(type(mol))
print(Draw.ShowMol(mol[0]))

df1= pd.read_csv(')
PandasTools.AddMoleculeColumnToFrame(df1,'smiles','mol',includeFingerprints=True)
df1['MW'] = df1['mol'].apply(Descriptors.MolWt)  
import rdkit
from rdkit import Chem
#导入一个分子
smi = 'c1ccccc1'
#rdkit读取
mol = Chem.MolFromSmiles(smi)
#获取分子中的原子数目
atom_num = mol.GetNumAtoms()
#获取分子中的键数目
bond_num = mol.GetNumBonds()
nei_atom = []
nei_bond = []
#获取分子中的原子的相邻原子的数目以及序号
for i in range(atom_num):
  nei_atom.append([(i.GetSmarts(),i.GetIdx()) for i in mol.GetAtomWithIdx(i).GetNeighbors()])
#获取分子中键的相关类型以及键的特征，以及成键原子序号
for i in range(bond_num):
  nei_bond.append((mol.GetBondWithIdx(i).GetBondType().name,mol.GetBondWithIdx(i).GetBeginAtomIdx(),mol.GetBondWithIdx(i).GetEndAtomIdx()))
#我们就获得了分子中的原子以及成键信息，后续操作待定
In [90]: nei_atom
Out[90]:
[[('C', 1), ('C', 5)],
[('C', 0), ('C', 2)],
[('C', 1), ('C', 3)],
[('C', 2), ('C', 4)],
[('C', 3), ('C', 5)],
[('C', 4), ('C', 0)]]


In [91]: nei_bond
Out[91]:
[('SINGLE', 0, 1),
('DOUBLE', 1, 2),
('SINGLE', 2, 3),
('DOUBLE', 3, 4),
('SINGLE', 4, 5),
('SINGLE', 5, 0)]

xw:
def delete_attachment(smi):
    m = Chem.MolFromSmiles(smi)
    if m is not None:
        mw = Chem.RWMol(m)
        atta_idx = []
        for i in range(0,mw.GetNumAtoms()):
            if mw.GetAtomWithIdx(i).GetSymbol() == '*':
                atta_idx.append(i)
        for i in reversed(atta_idx):
            mw.RemoveAtom(i)
        return Chem.MolToSmiles(mw)
    else:
        return '[Xe]'
'''

