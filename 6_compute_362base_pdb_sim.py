# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 19:48:53 2023

@author: HY
"""


import os

from Bio.PDB import PDBParser, NeighborSearch, Selection, Polypeptide,standard_aa_names,PDBIO,Select,Vector

from rdkit import Chem
from rdkit.Chem import AllChem,Descriptors
from rdkit.Chem import Crippen
import pandas as pd

from rdkit.Chem import rdFMCS, MolFromSmiles, MolToSmiles, rdmolops, rdchem, Recap
import numpy as np

import random
from rdkit import DataStructs
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import blosum62
from sklearn.metrics.pairwise import cosine_similarity
from rdkit.Chem import MolSurf
from rdkit.Chem import MACCSkeys
import ast


# 文件夹路径
folder_path = 'D:/llinfile/dataset/PDBbind_v2020_refined/refined-set/selected_base_pdb_pocket'  # 替换为实际的文件夹路径

# 列出文件夹中的所有文件
file_list = os.listdir(folder_path)

# 存储处理后的结果的列表
result_vector = []


        
def pocket_similarity(row,protein2_pdb):
    # 计算蛋白质口袋的形状描述符（Morgan指纹）
    pdb_id = row['pdb_id'].strip('"')
    ligand_name = row['ligand_name'].strip('"')
    protein1_pdb = os.path.join("D:/llinfile/dataset/PDB/PDB_homo_processed_select_ligandsplit_4_5_3",row['Uniprot_accession'],f"{pdb_id}_pocket_{ligand_name}.pdb")
    print(protein1_pdb)
    mol1 = Chem.MolFromPDBFile(protein1_pdb)
    mol2 = Chem.MolFromPDBFile(protein2_pdb)
    fp_radius = 2
    fp_nBits = 2048
    try:
        pocket_fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, fp_radius, nBits=fp_nBits)
        pocket_fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, fp_radius, nBits=fp_nBits)
        
        # 计算口袋形状特征
        #pocket_atoms1 = [atom for atom in protein1.get_atoms()]
        #pocket_atoms2 = [atom for atom in protein2.get_atoms()]
        
        pocket_macc1 = MACCSkeys.GenMACCSKeys(mol1)
        pocket_macc2 = MACCSkeys.GenMACCSKeys(mol2)
     
        #similarity_score = cosine_similarity(features1.reshape(1, -1), features2.reshape(1, -1))[0][0]
        # 计算形状相似性
        similarity_1 = DataStructs.DiceSimilarity(pocket_fp1 , pocket_fp2)
        similarity_2 = DataStructs.DiceSimilarity(pocket_macc1, pocket_macc2)
        similarity_score = (similarity_1 + similarity_2)/2
        print("结合位点相似性分数:", similarity_score,similarity_1,similarity_2)
        return similarity_score
    except :
        return -100  # 或者根据需求返回其他默认值
'''
def pocket_similarity(protein1_pdb,protein2_pdb):
    # 计算蛋白质口袋的形状描述符（Morgan指纹）
    #pdb_id = row['pdb_id'].strip('"')
    #ligand_name = row['ligand_name'].strip('"')
    #protein1_pdb = os.path.join("D:/llinfile/dataset/PDB/PDB_homo_processed_select_ligandsplit_4_5_3",row['Uniprot_accession'],f"{pdb_id}_pocket_{ligand_name}.pdb")
    #print(protein1_pdb)
    mol1 = Chem.MolFromPDBFile(protein1_pdb)
    mol2 = Chem.MolFromPDBFile(protein2_pdb)
    fp_radius = 2
    fp_nBits = 2048
    try:
        pocket_fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, fp_radius, nBits=fp_nBits)
        pocket_fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, fp_radius, nBits=fp_nBits)
        
        # 计算口袋形状特征
        #pocket_atoms1 = [atom for atom in protein1.get_atoms()]
        #pocket_atoms2 = [atom for atom in protein2.get_atoms()]
        
        pocket_macc1 = MACCSkeys.GenMACCSKeys(mol1)
        pocket_macc2 = MACCSkeys.GenMACCSKeys(mol2)
     
        #similarity_score = cosine_similarity(features1.reshape(1, -1), features2.reshape(1, -1))[0][0]
        # 计算形状相似性
        similarity_1 = DataStructs.DiceSimilarity(pocket_fp1 , pocket_fp2)
        similarity_2 = DataStructs.DiceSimilarity(pocket_macc1, pocket_macc2)
        similarity_score = (similarity_1 + similarity_2)/2
        #print("结合位点相似性分数:", similarity_score,similarity_1,similarity_2)
        return similarity_score
    except :
        print(protein1_pdb)
        return -100  # 或者根据需求返回其他默认值

# 逐一处理文件
for file_name in file_list:
    file_path = os.path.join(folder_path, file_name)
    if os.path.isfile(file_path):  # 确保处理的是文件而非文件夹
        result = pocket_similarity(file_path,src_pdb)  # 调用处理函数处理文件
        result_vector.append(result)  # 将处理结果添加到结果向量中
print(result_vector)'''

 # 读取原始csv文件
df = pd.read_csv('./output_PDB_homo_processed_select_ligandsplit_4_5_4.csv')

##########################################################################################################################
# 按条件筛选行,选取与需计算的靶点口袋氨基酸序列局部相似度前1000的PDB
# 计算口袋相似性
# 逐一处理文件
for file_name in file_list:
    file_path = os.path.join(folder_path, file_name)
    if os.path.isfile(file_path):
        column_name = str(file_name)[:4]+"_sim"  # 将列表值的前4位转换为字符串作为新列的列名
        df[column_name] = df.apply(pocket_similarity, args=( file_path,), axis=1, result_type='expand')
        
df.to_csv('output_PDB_homo_processed_select_ligandsplit_4_5_4_compute_362base_sim_refined.csv', index=None)
