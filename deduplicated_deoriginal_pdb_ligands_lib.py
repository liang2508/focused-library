# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 16:48:41 2023

@author: HY
"""

import pandas as pd

# 读取第一个CSV文件
df1 = pd.read_csv('selected_pdbfile_by_keyaas_overlap_descending.csv', usecols=['Uniprot_accession', 'pdb_id', 'ligand_name', 'ligand_smiles'])
df1 = df1.head(1000)
# 读取第二个CSV文件
df2 = pd.read_csv('selected_pdbfile_by_pocket_seq_similarity_1000.csv', usecols=['Uniprot_accession', 'pdb_id', 'ligand_name', 'ligand_smiles'])
# 读取第三个CSV文件
df3 = pd.read_csv('selected_pdbfile_by_pocket_similarity.csv', usecols=['Uniprot_accession', 'pdb_id', 'ligand_name', 'ligand_smiles'])

# 读取第四个CSV文件
df4 = pd.read_csv('selected_pdbfile_by_compared_pocket_similarity_refined_1000.csv', usecols=['Uniprot_accession', 'pdb_id', 'ligand_name', 'ligand_smiles'])

# 合并两个DataFrame
merged_df = pd.concat([df1, df2,df3,df4], axis=0)

# 去重
deduplicated_df = merged_df.drop_duplicates(subset=['ligand_smiles'])

#去原靶点数据
deoriginal_df = deduplicated_df[deduplicated_df['Uniprot_accession']!='Q07889']

# 读取第三个CSV文件
other_df = pd.read_csv('SOS1_pdb_molecule.csv')

# 去除存在于另一个CSV文件中的数据
final_df = deoriginal_df[~deoriginal_df['ligand_smiles'].isin(other_df['ligand_smiles'])].dropna()

#添加靶点+pdb名+配体名作为化合物名方便后续进行反向研究
final_df['name'] = final_df['Uniprot_accession']+ '_' +final_df['pdb_id'].str.strip('"') + '_' + final_df['ligand_name'].str.strip('"')
# 输出到新的CSV文件
final_df.to_csv('output_deduplicated_deoriginal_pdb_ligands_lib_4_methods.csv', index=False)