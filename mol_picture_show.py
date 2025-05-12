# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 12:04:49 2022

@author: ll
"""


from rdkit.Chem import AllChem as ch

from rdkit.Chem import Draw as Draw

from rdkit import DataStructs
import pandas as pd
import numpy as np


#载入分子库

#suppl = ch.SDMolSupplier(r'E:\llinfile\dataset\SOS1_dataset\scaffold dataset\scaffold_set_100_remove_dup.sdf')
#data = pd.read_csv(r'SOS1_pdb_molecule.csv')
data = pd.read_csv(r'scaffold_set_smiles_processed.csv')
suppl = data["smiles"]

molss = [x for x in suppl if x is not None]
mols=[]
for m in molss:
    #s = ch.MolToSmiles(m)
    mol=ch.MolFromSmiles(m)
    mols.append(mol)

print(len(mols))

img=Draw.MolsToGridImage([m for m in mols], molsPerRow=5,subImgSize=(400,400))
print(img)
img.save(r"scaffold_set_smiles_nik.png")
img.show()

#print(result)

'''
#img=Draw.MolsToGridImage([m for m, sim in result], subImgSize=(500,500),legends=[mol.GetProp('Name') +': '+ str(sim)for mol, sim in result])
img=Draw.MolsToGridImage([m for m, sim in result], molsPerRow=5,subImgSize=(400,400),legends=[str(sim)for mol, sim in result])
print(img)

# kekulize=False, 

#img.save(r"D:\fragment_lib_sos1_3.png")
img.show()

'''




