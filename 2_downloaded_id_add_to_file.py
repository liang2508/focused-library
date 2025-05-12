# -*- coding: utf-8 -*-
"""
Created on Thu May  4 16:34:24 2023

@author: HY
"""

import os
import csv

# 读取CSV文件
with open('uniprot_pdb_list_rattus_protein_downloaded.csv', 'r') as file:
    reader = csv.reader(file)
    rows = [row for row in reader]

# 给表格添加新列downloaded_pdb_list
#rows[0].append('download_num')
rows[0].append('downloaded_pdb_list')  # 表头添加新列
for row in rows[1:]:
#    row.append('')  # 初始化新列为空
    row.append('')
# 遍历文件夹
folder_path = './PDB_all/rattus'
for folder_name in os.listdir(folder_path):
    for row in rows[1:]:
        if folder_name == row[1]:
            # 遍历文件夹中所有文件，将文件名前四位写入列表
            file_path = os.path.join(folder_path, folder_name)
            file_list = []
            for file_name in os.listdir(file_path):
                file_list.append(file_name[:4])
            # 将列表写入对应的行列B中
            row[5] = len(file_list)
            row[6] = ','.join(file_list)
            

# 将修改后的表格写入CSV文件
with open('uniprot_pdb_list_rattus_protein_downloaded_new.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(rows)