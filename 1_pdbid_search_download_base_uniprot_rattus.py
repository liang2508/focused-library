# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 14:42:37 2023

@author: HY
"""
#use pdbe
import requests
import os

import urllib.request

import csv
import pandas as pd

def get_pdb_ids(uniprot_id): #, resolution_threshold=3.0
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/best_structures/{uniprot_id}?pretty=True"
    
    response = requests.get(url)
    data = response.json()
    pdb_ids = []
    if uniprot_id not in data:
        print(f"No PDB structures found for Uniprot ID: {uniprot_id}")
    else:
        for pdb_entry in data[uniprot_id]:
            #if pdb_entry["resolution"] <= resolution_threshold and (pdb_entry["pdb_id"] not in pdb_ids):
            if (pdb_entry["pdb_id"] not in pdb_ids):
                pdb_ids.append(pdb_entry["pdb_id"])

    return pdb_ids


def download_pdb(uniprot,pdb_ids):
    count = 0
    for pdbid in pdb_ids:
        # 构建PDB文件下载链接
        #url = f'https://www.ebi.ac.uk/pdbe/entry-files/download/pdb{pdbid}.ent'
        url = f'https://files.rcsb.org/download/{pdbid}.pdb'
    
        # 指定PDB文件保存路径# 指定下载目录
        download_dir = f'./PDB_all/rattus/{uniprot}'
        filename = f'{pdbid}.pdb'
        os.makedirs(download_dir, exist_ok=True)
        filepath = os.path.join(download_dir, filename)
    
        try:
            # 下载PDB文件
            urllib.request.urlretrieve(url, filepath)
    
            print(f'{pdbid} 下载完成')
            count += 1
        except urllib.error.HTTPError as e:
            # 如果无法下载，输出相应的消息
            print(f'{pdbid} 无法下载: {e}')
    return count

with open('./log_rattus_download.txt', 'w') as f:
    df = pd.read_csv("./chembl31_single_protein_list_deduplicated_rattus.csv", usecols=["UniProt_Accessions","sequence"])
    uniprot_ids = list(df["UniProt_Accessions"])
    df["pdb_id"]=''
    df["pdb_num"]=''
    df["download_num"]=''
    
    
    for uniprot_id in uniprot_ids:
        pdb_ids = get_pdb_ids(uniprot_id)
        df["pdb_id"][df.loc[df['UniProt_Accessions'] == uniprot_id].index.tolist()]= str(pdb_ids)
        df["pdb_num"][df.loc[df['UniProt_Accessions'] == uniprot_id].index.tolist()]= len(pdb_ids)
        print(f"PDB IDs for Uniprot ID '{uniprot_id}'  :",file = f)
        for pdb_id in pdb_ids:
            print(pdb_id,file = f)
            
        downloadnum = download_pdb(uniprot_id,pdb_ids)
        df["download_num"][df.loc[df['UniProt_Accessions'] == uniprot_id].index.tolist()]= downloadnum
        
    df.to_csv("./uniprot_pdb_list_rattus_protein.csv")
    
        


