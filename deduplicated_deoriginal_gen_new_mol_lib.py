# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 16:48:41 2023

@author: HY
"""

import pandas as pd

# 读取第一个CSV文件
df1 = pd.read_csv('new_molecules_file_by_keyaas_overlap_filter_10000.csv', usecols=['new_smiles'])
df1["method"]='keyaas'

# 读取第二个CSV文件
df2 = pd.read_csv('new_molecules_file_by_pocket_seq_similarity_filter_10000.csv', usecols=['new_smiles'])
df2["method"]='pocket_seq_sim'
# 读取第三个CSV文件
df3 = pd.read_csv('new_molecules_file_by_pocket_similarity_filter_10000.csv', usecols=['new_smiles'])
df3["method"]='pocket_sim'
# 读取第四个CSV文件
df4 = pd.read_csv('new_molecules_file_by_compared_pocket_similarity_refined_filter_10000.csv', usecols=['new_smiles'])
df4["method"]='compared_pocket_sim'

# 读取第五个CSV文件
df5 = pd.read_csv('random_sample_gen_by_decorator_10000.csv', usecols=['new_smiles'])
df5["method"]='decorator'

# 合并多个DataFrame
merged_df = pd.concat([df1, df2,df3,df4,df5], axis=0)

# 去重
deduplicated_df = merged_df.drop_duplicates(subset=['new_smiles'])

#去原靶点数据
#deoriginal_df = deduplicated_df[deduplicated_df['Uniprot_accession']!='Q07889']

# 读取第三个CSV文件
other_df = pd.read_csv('SOS1_pdb_molecule.csv')

# 去除存在于另一个CSV文件中的数据
final_df = deduplicated_df[~deduplicated_df['new_smiles'].isin(other_df['ligand_smiles'])].dropna()

# 输出到新的CSV文件
final_df.to_csv('output_deduplicated_deoriginal_gen_mol_4_methods.csv', index=True)