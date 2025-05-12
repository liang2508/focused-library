# -*- coding: utf-8 -*-
"""
Created on Fri Sep 15 20:21:44 2023

@author: HY
"""

import csv
from rdkit import Chem

# 读取 CSV 文件，并获取指定列的数据
def read_csv(filename, column_index):
    data = []
    with open(filename, 'r', newline='') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if len(row) > column_index:
                data.append(row[column_index])
    return data

# 标准化 SMILES
def normalize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is not None:
        return Chem.MolToSmiles(mol,canonical=True)
    return None

# 统计重复的小分子 SMILES
def find_duplicate_smiles(smiles_list1, smiles_list2):
    normalized_smiles1_list=[]
    for smile1 in smiles_list1:
        normalized_smiles1 = normalize_smiles(smile1)
        if normalized_smiles1 is not None:
            normalized_smiles1_list.append(normalized_smiles1)
    duplicate_smiles = []
    duplicate_count = 0    
    smiles_set1 = set(normalized_smiles1_list)
    for smiles2 in smiles_list2:
        normalized_smiles2 = normalize_smiles(smiles2)
        if normalized_smiles2 is not None and normalized_smiles2 in smiles_set1:
            duplicate_smiles.append(normalized_smiles2)
            duplicate_count += 1
    return duplicate_smiles, duplicate_count

# 文件 1 的路径和列索引
file1_path = 'output_PDB_all_processed_select_ligandsplit_4_5_4.csv'
file1_column_index = 6

# 文件 2 的路径和列索引
file2_path = 'new_molecules_file_by_compared_pocket_similarity_refined_filter_10000_nik.csv'
file2_column_index = 3

# 读取文件 1 和文件 2 的小分子 SMILES 列数据
smiles_list1 = read_csv(file1_path, file1_column_index)
smiles_list2 = read_csv(file2_path, file2_column_index)

# 判断重复的小分子 SMILES 并统计个数
duplicate_smiles, duplicate_count = find_duplicate_smiles(smiles_list1, smiles_list2)

# 输出结果
print("重复的小分子 SMILES：")
for smiles in duplicate_smiles:
    print(smiles)
print("重复个数：", duplicate_count)