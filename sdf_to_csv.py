# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 15:54:18 2023

@author: HY
"""

from rdkit import Chem
from rdkit.Chem import PandasTools
import pandas as pd
# 读取SDF文件
sdf_file_path = 'generated_molecules_nik_add_H.sdf'
supplier = Chem.SDMolSupplier(sdf_file_path)

# 提取SMILES并存储到列表中
smiles_list = [Chem.MolToSmiles(mol) for mol in supplier if mol is not None]

# 创建包含 SMILES 列的 DataFrame
df = pd.DataFrame({'SMILES': smiles_list})

# 将 DataFrame 保存为 CSV 文件

output_csv_path = 'generated_molecules_nik_add_H.csv'
df.to_csv(output_csv_path, index=False)
