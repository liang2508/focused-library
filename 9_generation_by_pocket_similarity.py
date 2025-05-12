# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:35:31 2023

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
#import csv
#import math
#from sklearn.cluster import AgglomerativeClustering
#from sklearn.cluster import KMeans
#import matplotlib.pyplot as plt
#from sklearn.preprocessing import MinMaxScaler
#import matplotlib.colors as colors




# 定义源文件
src_pdb = './pdb_files/4g3e.pdb'



def pocket_similarity(row,protein2_pdb):
    # 计算蛋白质口袋的形状描述符（Morgan指纹）
    pdb_id = row['pdb_id'].strip('"')
    ligand_name = row['ligand_name'].strip('"')
    protein1_pdb = os.path.join("D:/llinfile/dataset/PDB/PDB_processed_select_ligandsplit_4_5_3",row['Uniprot_accession'],f"{pdb_id}_pocket_{ligand_name}.pdb")
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
# 计算疏水性特征（使用电荷）
def calculate_hydrophobicity(molecule):
    mol = Chem.AddHs(molecule)
    AllChem.ComputeGasteigerCharges(mol)
    charges = [atom.GetDoubleProp('_GasteigerCharge') for atom in mol.GetAtoms()]
    hydrophobicity = sum(charges)
    return hydrophobicity

pocket_hydrophobicity1 = calculate_hydrophobicity(mol1)
pocket_hydrophobicity2 = calculate_hydrophobicity(mol2)

pocket_hydrophobicity1 = MolSurf.pyLabuteASA(mol1)
pocket_hydrophobicity2 = MolSurf.pyLabuteASA(mol2
                                             
                                             # 将计算得到的特征拼接成向量
features1 = np.concatenate([
    np.array(pocket_fp1),
    np.array([pocket_hydrophobicity1]),
])

features2 = np.concatenate([
    np.array(pocket_fp2),
    np.array([pocket_hydrophobicity2]),
])
'''


   
#局部序列比对
def seq_similarity_local(seq1, seq2):
    try:
        alignment = pairwise2.align.localms(seq1, seq2, 2, -1, -.5, -.1)
        score = alignment[0][2]
        start = alignment[0][3]
        end = alignment[0][4]
        return score,start,end
    except IndexError:
        return -100  # 或者根据需求返回其他默认值

#全局序列比对
def seq_similarity_global(seq1, seq2):
    try:
        alignment = pairwise2.align.globalms(seq1, seq2, 2, -1, -.5, -.1)
        score = alignment[0][2]
        start = alignment[0][3]
        end = alignment[0][4]
        return score,start,end
    except IndexError:
        return -100

def calculate_similarity(row, sequence):
    similarities = seq_similarity_local(row['pocket_sequence'], sequence)
    return pd.Series(similarities, index=['pocket_seq_similarity', 'pocket_seq_similarity_start', 'pocket_seq_similarity_end'])

def calculate_global_similarity(row, sequence):
    similarities = seq_similarity_global(row['pocket_sequence'], sequence)
    return pd.Series(similarities, index=['pocket_seq_similarity', 'pocket_seq_similarity_start', 'pocket_seq_similarity_end'])

#-------------------------------------------------------------------
# 根据配体进行筛选
# 找到所有配体
excluded_residues = ["1PE","2HT","2PE","7PE","ACT","ACY","AKG","BCT","BMA","BME","BOG","BU3","BUD",
                     "CAC","CIT","CME","CO3","DMS","DTT","DTV","EDO","EPE","FES","FMT","GBL","GOL",
                     "GSH","HEC","HED","HEM","IMD","IOD","IPA","MAN","MES","MG8","MLI","MO6","MPD",
                     "MYR","NAG","NCO","NH3","NO3","OCT","OGA","OPG","P2U","PG4","PGE","PGO","PHO",
                     "PLP","PO4","POP","PSE","PSU","PTL","SGM","SO4","SPD","SPM","SRT","TAM","TAR",
                     "TFA","TLA","TPP","TRS"]

#,"ATP","GTP","UMP","CMP","AMP","GMP","IMP","ADP","ATP","dAMP","dADP","dATP"保留核苷酸及其衍生物
#氨基酸缩写
aa_dict = {
            "ALA": "A",
            "ARG": "R",
            "ASN": "N",
            "ASP": "D",
            "CYS": "C",
            "GLU": "E",
            "GLN": "Q",
            "GLY": "G",
            "HIS": "H",
            "ILE": "I",
            "LEU": "L",
            "LYS": "K",
            "MET": "M",
            "PHE": "F",
            "PRO": "P",
            "SER": "S",
            "THR": "T",
            "TRP": "W",
            "TYR": "Y",
            "VAL": "V"
        }




# 定义 Pi-Pi 堆积的氨基酸类型
PI_RESIDUES = set(['PHE', 'TYR', 'TRP', 'HIS'])

# 定义哪些氨基酸是芳香的，以及它们的环中的原子
aromatic_residues = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TRP': ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
    'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2']  # HIS的咪唑环是一个五元环
}

# Define interaction types and parameters
interaction_types = {
    #"hbond": {"donor_atoms": ["NE", "NH1", "NH2", "ND1", "OD1", "SD"], "acceptor_atoms": ["O", "OD1", "OD2", "OE1", "OE2", "NE2", "ND1", "SG","N"]},
    "hbond": {"donor_atoms": ["N",  "O", "S"], "acceptor_atoms": ["O", "N"],"cutoff_angle":45.0},
    "ionic": {"positive_residues": ["ARG", "HIS", "LYS"], "negative_residues": ["ASP", "GLU"], "cutoff_distance": 4.0, "cutoff_angle": 120.0},
    "hydrophobic": {"residues": ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"], "cutoff_distance": 4.5},
    "pi_stacking": {"aromatic_residues": ["PHE", "TYR", "TRP", "HIS"], "cutoff_distance": 6.0, "cutoff_angle": 30.0},
    "hydrophobe":{"residues": ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TYR", "TRP"],"ligands": ["C","S","F","CL","BR","I"], "cutoff_distance": 4.5},
    "pi_cation": {"aromatic_residues": ["PHE", "TYR", "TRP", "HIS"], "positive_residues": ["ARG", "LYS"], "cutoff_distance": 4.0},
    #"salt_bridge": {"positive_residues": ["ARG", "HIS", "LYS"], "negative_residues": ["ASP", "GLU"], "cutoff_distance": 6.0},
    "metal_coordination": {"metal_ions": ["CA","CD","CO", "FE", "CU", "MG","MN","NI","ZN" ],"acceptor_atoms": ["O", "N"], "cutoff_distance": 2.8},
}


def compute_centroid(coords):
    return np.mean(coords, axis=0)


# 定义函数，输入3个原子的坐标，返回该平面的法向量
def calculate_normal(p1, p2, p3):
    v1 = p2 - p1
    v2 = p3 - p1
    normal = np.cross(v1, v2)
    return normal / np.linalg.norm(normal)

# 定义函数，输入芳香环上的原子坐标，返回该环的平面法向量
def calculate_plane_normal(atom_coords):
    p1, p2, p3 = atom_coords[:3]
    normal = calculate_normal(p1, p2, p3)
    if np.allclose(normal, [0, 0, 0]):
        # 如果计算得到的法向量为零向量，说明三个原子共线，选择不同的三个原子重新计算
        p1, p2, p3 = atom_coords[1:4]
        normal = calculate_normal(p1, p2, p3)
    return normal

# 定义函数，输入两个芳香环的原子坐标，返回它们的平面夹角（单位为弧度）
def calculate_plane_angle(coords1, coords2):
    normal1 = calculate_plane_normal(coords1)
    normal2 = calculate_plane_normal(coords2)
    dot_product = np.dot(normal1, normal2)
    angle_rad = np.arccos(dot_product / np.linalg.norm(normal1) / np.linalg.norm(normal2))
    return angle_rad


class ResidueSelector(Select):
    def __init__(self, residue_ids):
        self.residue_ids = residue_ids

    def accept_residue(self, residue):
        if residue.get_id()[1] in self.residue_ids:
            return True
        else:
            return False

def listdir(path,item):  #传入存储的list    
    file_path = os.path.join(path, item)
    if os.path.isfile(file_path): 
        file_path = file_path
        #print(file_path)
        return file_path
    else:
        for file in os.listdir(file_path):
            return listdir(file_path, file)
        
# 创建一个函数，计算两个列表的重合度
def overlap_ratio(list1, list2):
    set1, set2 = set(list1), set(list2)
    overlap = set1 & set2
    return len(overlap)/((len(set1)+len(set2))/2) #重合度计算法2：同时考虑对比双方的关键氨基酸列表长度
    #return len(overlap) / len(set2) #重合度计算法1：只计算重合的氨基酸数量在研究的靶点的关键氨基酸列表占比
    
def breakmol(df):
    # 基于recap规则打碎分子
    fragments = []
    for i, row in df.iterrows():
        mol = Chem.MolFromSmiles(row[1])
        if mol:
            for fragment in AllChem.RecapDecompose(mol):
                smiles = Chem.MolToSmiles(fragment)
                fragments.append([row[1], smiles])



similarity_way = 0
    
parser = PDBParser(QUIET=True)
structure = parser.get_structure("target",src_pdb)
# 获取所有的氨基酸残基
residues = list(structure.get_residues())
ligands = []
selected_ligands =[]
selected_ligands_name =[]
selected_ligands_id =[]
selected_chain = []
resolution = structure.header.get('resolution', 0.0)
for model in structure:
    for chain in model:
        for residue in chain:
            if residue.id[0].startswith('H_'):
                # Exclude water molecules and other non-ligand residues
                #print(residue.id[0])
                ligands.append(residue)
#if len(ligands) == 0:
    #print("无配体")
if len(ligands) != 0:
    for ligand in ligands:
        #判断原子数是否符合范围
        if len(ligand) >= 5 and len(ligand) <= 130:
            #print("原子数符合范围："+str(len(ligand)))
            #判断是否为非氨基酸残基（非多肽）
            if Polypeptide.is_aa(ligand, standard=False) is not True:#standard = true 只检查20个标准氨基酸
                #print("非氨基酸残基")
                #判断是否为可参考配体
                if ligand.resname not in excluded_residues:
                    #提取配体
                    ligand_extract = ''
                    #重复配体只取第一个并排除命名不规范的ligand，如7wv3
                    if ligand.resname not in selected_ligands_name and len(ligand.resname)== 3:
                        with open(src_pdb, 'r') as f:
                            for line in f:
                                #print(" "+ligand.get_parent().id+("{:>4}".format(str(ligand.get_id()[1]))+"    "))
                                if line.startswith('HETATM') and (ligand.resname in line ) and (" "+ligand.get_parent().id+("{:>4}".format(str(ligand.get_id()[1]))+"    ") in line) :
                                    
                                    ligand_extract += line
                                    #print(line)
                                elif line.startswith('CONECT') :
                                    ligand_extract += line
                                elif line.startswith('END'):
                                    # 如果 PDB 文件中包含多个模型，分别提取小配体
                                    break 
                            ligand_file_path = os.path.splitext(src_pdb)[0] + f'_ligand_{ligand.resname}.pdb'
                            # 确保文件夹路径存在
                            os.makedirs(os.path.dirname(ligand_file_path), exist_ok=True)
                            with open(ligand_file_path, 'w') as f:
                                f.write(ligand_extract)
                            #print(f"Small molecule ligand has been extracted and saved in {ligand_file_path}.")
                            # 尝试读取PDB文件，如果失败则删除该文件
                            try:
                                # Convert the ligand PDB file into an RDKit molecule
                                ligand_mol = Chem.MolFromPDBFile(ligand_file_path)
                                
                                # Add hydrogen atoms to the ligand molecule
                                #ligand_mol = Chem.AddHs(ligand_mol)
                                
                                # Kekulize the ligand molecule to restore the double bonds in aromatic rings
                                #rdmolops.Kekulize(ligand_mol, clearAromaticFlags=False)
                                
                                # Generate the canonical SMILES string for the ligand
                                
                                ligand_smiles = Chem.MolToSmiles(ligand_mol, isomericSmiles=True, canonical=True)
                            except:
                                #print("无法正常读取PDB文件，删除文件:", ligand_file_path)
                                
                                os.remove(ligand_file_path)
                                continue
                            
                                
                                
                                
                            #print(f"Canonical SMILES: {ligand_smiles}")
                        selected_ligands_name.append(ligand.resname)
                        selected_ligands.append(ligand)
                        selected_ligands_id.append(ligand.get_id()[1])
                        selected_chain.append(" "+ligand.get_parent().id+" ")
                        # 获取配体的原子列表
                        ligand_atoms = list(ligand.get_atoms())
                        #print(ligand_atoms)
                        #print(list(atom.element for atom in ligand.get_atoms()))
                        # 输出配体信息
                        #print(f"Ligand name: {ligand.resname}")
                        #print(f"Ligand id: {ligand.get_id()[1]}")
                        #print(f"Ligand chain ID: {ligand.get_parent().id}")
                        
                    
                    
                        #计算口袋序列
                        ns = NeighborSearch(list(structure.get_atoms()))
                        nearby_residues = set()
                        for atom in ligand_atoms:
                            nearby_residues.update(ns.search(atom.coord, 6.0, level='R'))
                    
                        # 过滤掉非氨基酸残基，如水分子和配体
                        binding_pocket_residues = [r for r in nearby_residues if Polypeptide.is_aa(r)]
                        binding_pocket_residues_ids = [r.get_id()[1] for r in nearby_residues if Polypeptide.is_aa(r)]
                        selector = ResidueSelector(binding_pocket_residues_ids)
                            

                        # 将氨基酸序号和名称存储到一个字典中
                        
                        if binding_pocket_residues is not None:
                            # 将结合口袋保存为新的PDB文件
                            io = PDBIO()
                            io.set_structure(structure)
                            io.save(f"{os.path.splitext(src_pdb)[0]}_pocket_{ligand.resname}.pdb", selector)
                            pdb_pocket_filepath = f"{os.path.splitext(src_pdb)[0]}_pocket_{ligand.resname}.pdb"
                            bindsite_aas = []
                            start = 100000
                            end = 0
                               
                            res_dict = {}
                            for residue in binding_pocket_residues: 
                                residue_name = residue.get_resname()
                                residue_id = residue.get_id()[1]
                                chain_id = residue.get_full_id()[2]
                                #print(f"{residue_name} {chain_id} {residue_id}")
                                res_dict[residue.id[1]] = residue.get_resname()
                                # 获取结合域中编号最小和最大的氨基酸编号
                                if residue_id < start:
                                    start = residue_id
                                if residue_id > end:
                                    end = residue_id
                                    
                            # 按氨基酸序号排序并输出为缩写序列
                            sorted_residues = sorted(res_dict.items(), key=lambda x: x[0])
                            sequence = ""
                            for res in sorted_residues:
                                if res[1] in aa_dict:
                                    sequence += aa_dict[res[1]]
                                else:
                                    sequence +='X'
                            #print(sequence)
                            
                            # 将结合域中编号最小的氨基酸至编号最大的氨基酸整段序列的每个氨基酸和名称存储到字典中
                            res_dict_comp = {}
                            for residue in residues:
                                if residue.get_id()[1] >= start and residue.get_id()[1] <= end:
                                    res_dict_comp[residue.id[1]] = residue.get_resname()
                            
                            # 按氨基酸序号排序并输出为缩写序列
                            sorted_residues_comp = sorted(res_dict_comp.items(), key=lambda x: x[0])
                            sequence_comp = ""
                            for res in sorted_residues_comp:
                                if res[1] in aa_dict:
                                    sequence_comp += aa_dict[res[1]]
                                else:
                                    sequence_comp +='X'
                            #print(sequence_comp)
                            
                            # 判断可能有相互作用的关键氨基酸并输出
                            key_aas = []
                            # 计算配体中每个环的所有原子坐标
                            ligand_rings_coords = []
                            for ring in ligand_mol.GetRingInfo().AtomRings():
                                x, y, z = 0, 0, 0
                                ligand_ring_coords=[]
                                for atom_idx in ring:
                                    #print(atom_idx)
                                    pos = ligand_mol.GetConformer().GetAtomPosition(atom_idx)
                                    ligand_ring_coords.append((pos.x,pos.y,pos.z))
                                ligand_rings_coords.append(ligand_ring_coords)
                             
                            for residue in binding_pocket_residues: 
                                #hbond
                                for protein_atom in residue.get_atoms():
                                    #print(protein_atom.element)
                                    #print(protein_atom.get_name())
                                    if protein_atom.element in interaction_types["hbond"]["donor_atoms"]:
                                        for ligand_atom in ligand.get_atoms():
                                            #print(ligand_atom)
                                            if ligand_atom.element in interaction_types["hbond"]["acceptor_atoms"]:
                                                if protein_atom - ligand_atom <= 3.5:  # Cutoff distance for hydrogen bond
                                                    #print("Hbond")
                                                    angle = protein_atom.get_vector().angle(ligand_atom.get_vector())
                                                    if angle <= interaction_types["hbond"]["cutoff_angle"]:
                                                        if residue not in key_aas:
                                                            key_aas.append(residue)
                                                            #print(f"hbond -- {residue.resname} {residue.id[1]} ")
                                #hydrophobe
                                if residue.resname in interaction_types["hydrophobe"]["residues"]:
                                    for protein_atom in residue.get_atoms():
                                        for ligand_atom in ligand.get_atoms():
                                            #print(ligand_atom)
                                            if ligand_atom.element in interaction_types["hydrophobe"]["ligands"]:
                                                r = protein_atom - ligand_atom
                                                if r <= interaction_types["hydrophobe"]["cutoff_distance"]:
                                                    if residue not in key_aas:
                                                        key_aas.append(residue)
                                                        #print(f"hydrophobe -- {residue.resname} {residue.id[1]} ")
                                #pi-pi stacking
                                if residue.resname in interaction_types["pi_stacking"]["aromatic_residues"]:
                                    for ligand_ring_coords in ligand_rings_coords:
                        
                                        #计算小配体中芳香环中心原子的坐标
                                        ligand_aromatic_ring_center = compute_centroid(ligand_ring_coords)
                                        
                                        # 计算氨基酸中芳香环的中心原子坐标
                                        protein_aromatic_ring_centers =[]
                                        protein_aromatic_atoms=[]
                                        for atom in residue.get_atoms():
                                            if atom.get_name() in aromatic_residues[residue.resname]:
                                                protein_aromatic_atoms.append(atom)
                                        
                                        #print(protein_aromatic_atoms)
                                        protein_coords = [atom.get_coord() for atom in protein_aromatic_atoms]
                                        if len(protein_coords) > 0:
                                            protein_aromatic_ring_centers= compute_centroid(protein_coords)
                                            #print(protein_aromatic_ring_centers)
                                            #计算两中心原子间的距离
                                            distance = (Vector(protein_aromatic_ring_centers) - Vector(ligand_aromatic_ring_center)).norm()
                                            
                                            if distance <= interaction_types["pi_stacking"]["cutoff_distance"]:
                                                #print(distance)
                                                
                                                #print(str(residue.resname)+str(residue.get_id()[1])+"距离符合条件")
                                            
                                                # 定义两个三维坐标列表，其中分别包含氨基酸和配体芳香环上的所有原子坐标
                                                coord1 = np.array(protein_coords)
                                                coord2 = np.array(ligand_ring_coords)
                                                #计算两环平面夹角
                                                angle= calculate_plane_angle(coord1, coord2)
                                                # 将弧度转换为角度
                                                angle_degrees = np.degrees(angle)
                                                if angle_degrees > 90:
                                                    angle_degrees = 360-angle_degrees
                                                
                                                # 输出两个向量之间的夹角
                                                
                                                #print(angle)
                                                #print(angle_degrees)
                                                if angle_degrees <= interaction_types["pi_stacking"]["cutoff_angle"]:
                                                    #print(f"Pi-Pi stacking between {residue.resname} {residue.id[1]} and ligand aromatic ring")               
                                                    if residue not in key_aas:
                                                        key_aas.append(residue)
                                # metal
                                
                                for ligand_atom in ligand.get_atoms():
                                    #print(ligand_atom)
                                    if ligand_atom.element in interaction_types["metal_coordination"]["metal_ions"]:
                                        for protein_atom in residue.get_atoms():
                                            if protein_atom.element in interaction_types["metal_coordination"]["acceptor_atoms"]:
                                                r = protein_atom - ligand_atom
                                                if r <= interaction_types["metal_coordination"]["cutoff_distance"]:
                                                    if residue not in key_aas:
                                                        key_aas.append(residue)
                                                        #print(f"metal -- {residue.resname} {residue.id[1]} ")
                            #print(key_aas)
                            key_aas_list=[]
                            key_aas_with_id = []
                            for  key_aa in key_aas:
                                key_aas_list.append(key_aa.resname)
                                aas_id = str(key_aa.resname)+str(key_aa.id[1])
                                key_aas_with_id.append(aas_id)
                            #[str(target),str(full_sequence),'"'+str(pdb_name)+'"',str(resolution),'"'+str(ligand.resname)+'"',str(ligand.get_parent().id),str(ligand_smiles),str(sequence),str(sequence_comp),str(key_aas_list),str(key_aas_with_id)]
                            
                            
                            # 读取原始csv文件
                            df = pd.read_csv('./output_PDB_all_processed_select_ligandsplit_4_5_4_pocketseq_notnull.csv')
    
##########################################################################################################################
                            # 按条件筛选行,选取与需计算的靶点口袋氨基酸序列局部相似度前1000的PDB
                            # 计算口袋相似性
                            
                            df[['pocket_similarity']] = df.apply(pocket_similarity, args=(pdb_pocket_filepath,), axis=1, result_type='expand')
                            #计算口袋氨基酸序列的全局相似性
                            #df[['pocket_seq_similarity', 'pocket_seq_similarity_start', 'pocket_seq_similarity_end']] = df.apply(calculate_global_similarity, args=(sequence,), axis=1, result_type='expand')
                            
                            # 计算口袋完整序列的局部相似性
                            #df['pocket_comp_seq_similarity'] = df['pocket_comp_sequence'].apply(lambda x: seq_similarity_local(x, sequence_comp))
                            #df[['pocket_comp_seq_similarity', 'pocket_comp_seq_similarity_start', 'pocket_comp_seq_similarity_end']] = df.apply(calculate_similarity, args=(sequence_comp,), axis=1, result_type='expand')
                            # 计算口袋完整序列的全局相似性
                            #df[['pocket_comp_seq_similarity', 'pocket_comp_seq_similarity_start', 'pocket_comp_seq_similarity_end']] = df.apply(calculate_global_similarity, args=(sequence_comp,), axis=1, result_type='expand')
                            
                            # 按照相似性降序排列
                            seq_similarity_df = df.sort_values(by='pocket_similarity', ascending=False)
                            # 提取前1000行
                            top_1000_rows = seq_similarity_df.head(1000)
                            
                            # 将前1000行转换为字典
                            #top_1000_dict = top_1000_rows.to_dict(orient='records')

                            # 选取第2、4列并另存为csv文件
                            top_1000_rows.to_csv('selected_pdbfile_by_pocket_similarity_nik.csv', index=None)

                            
                            # 基于recap规则打碎分子
                            fragments = []
                            for i, row in top_1000_rows.iterrows():
                                mol = Chem.MolFromSmiles(row[6])
                                if mol:
                                    fragment_list = Recap.RecapDecompose(mol).children
                                    #fragment_list = [[Chem.MolToSmiles(fragment)] for fragment in fragment_list]  # 将每个碎片单独成一个列表
                                    fragment_list = [fragment for fragment in fragment_list]
                                    for fragment in fragment_list:
                                        #fragments.append([row[6], Chem.MolToSmiles(Chem.MolFromSmiles(fragment[0])), fragment[0]])
                                        fragments.append([row[6], Chem.MolToSmiles(Chem.MolFromSmiles(fragment)), fragment])
                            
                            # 将碎片和原smiles另存为csv文件
                            fragments_df = pd.DataFrame(fragments, columns=['original_smiles', 'fragment_smiles', 'fragment'])
                            fragments_df.to_csv('fragments_file_by_pocket_similarity_nik.csv', index=None)
                            
                            # 将碎片库单独另存为csv文件
                            fragment_library = list(set(fragments_df['fragment'].tolist()))
                            fragment_library_df = pd.DataFrame(fragment_library, columns=['fragment'])
                            fragment_library_df.to_csv('fragment_library_file_by_pocket_similarity_nik.csv', index=None)
                            
                            # 随机重组碎片生成新的分子
                            new_molecules = []
                            new_molecules_smiles=[]
                            n=0
                            while n < 10000:  # 生成10000个新分子
                                selected_fragments = random.sample(fragment_library, 2)
                                #random.shuffle(fragment_library)
                                fragments_mols = [Chem.MolFromSmiles(smiles) for smiles in selected_fragments]
                                m1=fragments_mols[0]
                                m2=fragments_mols[1]
                                m = Chem.CombineMols(fragments_mols[0], fragments_mols[1])
                                mw = Chem.RWMol(m)
                                
                                #conn = []
                                for i in range(0,m1.GetNumAtoms()):
                                    if m1.GetAtomWithIdx(i).GetSymbol() == "*":  #判断是否是取代位点
                                        neighbor1_idx = m1.GetAtomWithIdx(i).GetNeighbors()[0].GetIdx() #取代位点邻位原子在分子中的序号
                                        for j in range(0,m2.GetNumAtoms()):
                                            if m2.GetAtomWithIdx(j).GetSymbol() == "*":
                                                neighbor2_idx = m2.GetAtomWithIdx(j).GetNeighbors()[0].GetIdx()
                                                mw.AddBond(neighbor1_idx,neighbor2_idx+m1.GetNumAtoms(),Chem.BondType.SINGLE)
                                                mw.RemoveAtom(j+m1.GetNumAtoms())
                                                mw.RemoveAtom(i)
                                                new_smiles=Chem.MolToSmiles(mw,canonical=True)
                                                mw=Chem.MolFromSmiles(new_smiles)
                                                
                                                # 计算分子质量
                                                mol_weight = Descriptors.MolWt(mw)
                                                logP = Crippen.MolLogP(mw)
                                                num_h_acceptors = Descriptors.NumHAcceptors(mw)
                                                num_h_donors = Descriptors.NumHDonors(mw)
                                                psa = Descriptors.TPSA(mw)
                                                
                                                #将新生成的分子进行五规则筛选并去重
                                                if mol_weight <= 500 and logP <= 5 and num_h_acceptors <= 10 and num_h_donors <= 5 and psa <= 140 and (new_smiles not in new_molecules_smiles):
                                                    #conn.append(Chem.MolToSmiles(mw)) #一个位点取代结束，进行存储
                                                    new_molecules.append([n+1,selected_fragments[0],selected_fragments[1], new_smiles])
                                                    new_molecules_smiles.append(new_smiles)
                                                    n += 1
                                                else:
                                                    n = n
                            
                            # 将新分子和smiles另存为csv文件
                            new_molecules_df = pd.DataFrame(new_molecules, columns=['id', 'frag1','frag2','new_smiles'])
                            new_molecules_df.to_csv('new_molecules_file_by_pocket_similarity_filter_NIK_10000.csv', index=None)
