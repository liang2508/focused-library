# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 15:07:52 2024

@author: HY
"""
from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd

# 读取CSV文件
data = pd.read_csv('./selected_pdbfile_by_keyaas_overlap_descending_tyr.csv')

# 创建RDKit的分子对象并添加到SDF中

output_sdf_file = './similar_pdb_ligands_tyr.sdf'
smiles_list = data['ligand_smiles']
sdf_writer = Chem.SDWriter(output_sdf_file)

for i, smiles in enumerate(smiles_list):
    # 创建RDKit的分子对象
    molecule = Chem.MolFromSmiles(smiles)
    
    # 如果分子对象创建成功，则将其添加到SDWriter中，并设置分子的Name属性
    if molecule is not None:
    # 设置分子的名称为该行数据的第一列和第二列的字符组合
        molecule_name = str(data.iloc[i, 2]) + '_' + str(data.iloc[i,4])
        molecule.SetProp('_Name', molecule_name)
    
    # 将修改后的分子写入SDF
    sdf_writer.write(molecule)

sdf_writer.close()

