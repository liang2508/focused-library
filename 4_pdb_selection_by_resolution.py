# -*- coding: utf-8 -*-
"""
Created on Thu May 18 10:51:22 2023

@author: HY
"""

import shutil
import os
from Bio.PDB.PDBParser import PDBParser

# 定义函数，递归读取文件夹中的文件，并对每个文件进行操作，最后复制到指定目录
def process_and_copy_files(src_folder, dst_folder):
    # 获取源文件夹的上级文件夹名称
    count_to_resolution = 0 
    for item in os.listdir(src_folder):
        # 构建源文件路径和目标文件路径
        src_path = os.path.join(src_folder, item)
        parent_folder_name = os.path.basename(os.path.dirname(src_path))
        #print(parent_folder_name+"********")
        dst_path = os.path.join(dst_folder, parent_folder_name+"/")#, item
        #print(dst_path)
        # 如果是文件夹，则递归调用该函数
         
        if os.path.isdir(src_path):
            process_and_copy_files(src_path, dst_folder)
        # 如果是文件，则进行操作后复制到目标文件夹中
        elif os.path.isfile(src_path):
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure(src_path,src_path)
            
            # 读取 PDB 文件中的生物信息和分辨率
            #organism = structure.header.get('organism', '?')
            resolution = structure.header.get('resolution', 0.0)
            
            #-------------------------------------------------------------------
            # 根据分辨率是否小于等于 3A进行筛选
            if resolution is not None:
                if  resolution <= select_criteria['resolution']: 
                    print("["+src_path+"]"+"分辨率符合要求:"+str(resolution)+"埃",file = f)
                    count_to_resolution += 1
                
                    # 根据分辨率是否小于等于 3A进行筛选
                    os.makedirs(dst_path, exist_ok=True)
                    shutil.copy2(src_path, dst_path)
                    print(f'select {src_path}',file = f)
                    
                else:
                    print("["+src_path+"]"+"分辨率不符合要求:"+str(resolution)+"埃",file = f)
            else:
                print("["+src_path+"]"+"无分辨率信息",file = f)
            
        
    print("分辨率符合标准的PDB数量有："+str(count_to_resolution)+"个",file = f)     
                
# 设置筛选条件
select_criteria = {
    #'organism': 'Homo sapiens',
    'resolution': 4.5,
}
# 定义源文件夹路径和目标文件夹路径
src_folder = './PDB_processed_1/other'
dst_folder = './PDB_processed_resolution_4_5_2/other'

# 调用函数进行处理和复制
with open('./log_PDB_other_processed_resolution_4_5_2.txt', 'w') as f:
    process_and_copy_files(src_folder, dst_folder)
