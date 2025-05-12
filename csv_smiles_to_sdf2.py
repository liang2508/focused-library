# -*- coding: utf-8 -*-
"""
Created on Sat Sep 30 21:11:38 2023

@author: HY
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, SDWriter

# 读取CSV文件
#data = pd.read_csv('output_deduplicated_deoriginal_pdb_ligands_lib_4_methods.csv')
data = pd.read_csv('sampled_nik_mol_gen_by_decorator.csv')
# 指定SDF文件保存的路径和文件名
output_file = 'sampled_nik_mol_gen_by_decorator.sdf'
#output_file = 'pdb_ligands_lib_4_methods.sdf'
# 创建SDWriter对象，用于保存所有分子的SDF数据
writer = SDWriter(output_file)

# 提取SMILES列和Name列数据
#smiles_list = data['ligand_smiles']
#name_list = data['name']
#smiles_list = data['new_smiles']
smiles_list = data['SMILES']
name = 'decorator'
#name = 'compared_pocket_similarity'
#name = 'pocket_similarity'
#name = 'keyaas_similarity'
#name = 'pocket_seq_similarity'
'''
用于pdb_ligand命名
# 遍历每个SMILES和Name，将其转换为分子对象并添加到SDWriter中
for i, (smiles, name) in enumerate(zip(smiles_list, name_list)):
    # 创建RDKit的分子对象
    molecule = Chem.MolFromSmiles(smiles)
    
    # 如果分子对象创建成功，则将其添加到SDWriter中，并设置分子的Name属性
    if molecule is not None:
        molecule.SetProp('_Name', name + '_' + str(i+1))
        writer.write(molecule)
        print(f"Converted SMILES to SDF: {name}_{i+1}")
    else:
        print(f"Failed to convert SMILES: {name}")
'''
# 遍历每个SMILES和Name，将其转换为分子对象并添加到SDWriter中
for i, smiles in enumerate(smiles_list):
    # 创建RDKit的分子对象
    molecule = Chem.MolFromSmiles(smiles)
    
    # 如果分子对象创建成功，则将其添加到SDWriter中，并设置分子的Name属性
    if molecule is not None:
        molecule.SetProp('_Name', name + '_' + str(i+1))
        writer.write(molecule)
        print(f"Converted SMILES to SDF: {name}_{i+1}")
    else:
        print(f"Failed to convert SMILES: {name}")
# 关闭SDWriter
writer.close()

print("Conversion completed.")

