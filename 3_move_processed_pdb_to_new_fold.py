# -*- coding: utf-8 -*-
"""
Created on Thu May 11 14:04:26 2023

@author: HY
"""

import os
import shutil

source_dir = "D:/llinfile/dataset/PDB/PDB_all/other"#
target_dir = "D:/llinfile/dataset/PDB/PDB_processed_1/other"

for root, dirs, files in os.walk(source_dir):
    for dir_name in dirs:
        sub_dir = os.path.join(root, dir_name)
        target_sub_dir = os.path.join(target_dir, dir_name)
        if not os.path.exists(target_sub_dir):
            os.makedirs(target_sub_dir)
        for file_name in os.listdir(sub_dir):
            if file_name.endswith("_1.pdb"):
                source_file = os.path.join(sub_dir, file_name)
                target_file = os.path.join(target_sub_dir, file_name)
                shutil.move(source_file, target_file)
                # 将文件名从“xxxx_1.pdb”改为“xxxx.pdb”
                new_file_name = file_name.replace("_1.pdb", ".pdb")
                os.rename(os.path.join(target_sub_dir, file_name), os.path.join(target_sub_dir, new_file_name))