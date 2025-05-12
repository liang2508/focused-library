# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 10:10:15 2023

@author: HY
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# 读取多个CSV文件，每个文件的SMILES所在的列数不同
csv_files = [
    
    #{'file': './new_molecules_file_by_pocket_seq_similarity_filter_10000.csv', 'smiles_column': 3},
    #{'file': './new_molecules_file_by_keyaas_overlap_filter_10000.csv', 'smiles_column':3},
    #{'file': './new_molecules_file_by_pocket_similarity_filter_10000.csv', 'smiles_column': 3},
    #{'file': './new_molecules_file_by_compared_pocket_similarity_refined_filter_10000.csv', 'smiles_column': 3},
    #{'file': './random_sample_gen_by_decorator_10000.csv', 'smiles_column': 1},
    #{'file': './SOS1_pdb_molecule.csv', 'smiles_column': 6},
    {'file': './output_deduplicated_deoriginal_gen_mol_5_methods.csv', 'smiles_column': 1},
    {'file': './external_dataset/random_asinex_50000.csv', 'smiles_column': 0},
    {'file': './external_dataset/random_chemdiv_50000.csv', 'smiles_column': 0},
    # Add more files as needed
]


# 循环使用颜色列表
colors = ['#845EC2', '#FF9671', '#F9F871']
#colors = ['#4B4453','#845EC2','#D65DB1','#FF6F91', '#FF9671','#FFC75F', '#F9F871']

# 绘制化学空间分布对比图
#labels = [f'File {i+1}' for i in range(len(csv_files))]  # 文件标签
labels = ["gen_molecules","random_asinex","random_chemdiv"]
#labels = ['pocket_seq_similarity', 'keyaas_similarity','pocket_similarity','compared_pocket_similarity','decorator', 'SOS1_pdb_molecule',"random_asinex","random_chemdiv"]
#labels = ['pocket_seq_similarity', 'keyaas_similarity','pocket_similarity','compared_pocket_similarity','decorator', "random_asinex","random_chemdiv"]
# 存储所有SMILES字符串和文件编号的列表
all_smiles = []
file_labels = []

# 遍历每个CSV文件，提取SMILES字符串
for i, csv_file in enumerate(csv_files):
    data = pd.read_csv(csv_file['file'])
    smiles_column = csv_file['smiles_column']  # 第3列的索引是2，根据你的文件结构调整
    smiles_list = list(data.iloc[:, smiles_column])
    all_smiles.extend(smiles_list)
    file_labels.extend([i] * len(smiles_list))  # 为每个文件添加对应的编号

# 将SMILES字符串转换为RDKit的Mol对象
molecules = [Chem.MolFromSmiles(smiles) for smiles in all_smiles if Chem.MolFromSmiles(smiles) is not None]

# 计算分子的分子指纹（Morgan指纹）
fingerprints = [AllChem.GetMorganFingerprintAsBitVect(mol, 2) for mol in molecules]

# 使用PCA进行降维
pca = PCA(n_components=2)
reduced_fingerprints = pca.fit_transform(fingerprints)

# 绘制化学空间分布对比图，每个文件用不同颜色表示
plt.figure(figsize=(10, 6))
for i in range(len(csv_files)):
    indices = [idx for idx, label in enumerate(file_labels) if label == i]
    plt.scatter(
        reduced_fingerprints[indices, 0],
        reduced_fingerprints[indices, 1],
        alpha=0.7,
        label=labels[i],
        color=colors[i]  # 使用循环颜色列表
    )

plt.title('Chemical Space Distribution Comparison')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend()
plt.show()


