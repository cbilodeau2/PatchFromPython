# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 11:06:24 2019
Takes a PDB files as an input, parses into patches and saves this as pickle
@author: camil
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib.colors as col
from scipy.optimize import curve_fit
import pandas as pd
from math import factorial
import random
import time
import Bio
import pickle

from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=1)


# Class defined for patches:
class Patch:
    def __init__(self,center,res_list,identification,pairlist):
        self.center = center
        self.res_list = res_list
        self.num = identification
        self.pairlist = pairlist

start = time.time()


# Load in PDB file:
structure_id = "fab"
filename = "fabA.pdb"
structure = parser.get_structure(structure_id,filename)

# Create list of all residues in protein:
residue_list = structure[0][" "].get_list()

# Define number of resiudes:
n_res = np.shape(residue_list)[0]

# Define patch (a class of objects) that surrounds each residue
patch_size = 10 # Number of residues per patch
all_patches = []

# For every residue, we define a patch
for centers in range(0,10): # n_res):
    print('Analyzing Patch', centers)
    patch_list = []
    res_cent = residue_list[centers]
    res_cent_cog = res_cog(res_cent)
    
    # Determine which residues should be included in the patch
    for non_centers in range(0,n_res):
        if np.shape(patch_list)[0]< patch_size:
            patch_list.append(residue_list[non_centers])
            
        elif np.shape(patch_list)[0]>= patch_size:
            
            # Calculate most distant residue in patch:
            patch_dist_list = [[np.linalg.norm(res_cog(patch_list[i])-res_cent_cog),i] for i in range(0,np.shape(patch_list)[0])]
            max_dist = max([patch_dist_list[i][0] for i in range(0,np.shape(patch_dist_list)[0])])

            # Calculate distance between center and non-center residue:
            dist  = np.linalg.norm(res_cent_cog-np.array(res_cog(residue_list[non_centers])))

            # If non-center residue is closer than farthest patch residue, 
            # Replace the patch residue with non-center residue:
            if dist<max_dist:
                for i in range(0,np.shape(patch_dist_list)[0]):
                    if patch_dist_list[i][0]==max_dist:
                        patch_list[i] = residue_list[non_centers]
    
    # Define final patch and add to patch list:
    new_patch = Patch(res_cog(residue_list[centers]),patch_list,centers,[]) 
    all_patches.append(new_patch)
    
    

end = time.time()
 
print('Time Elapsed:',end-start)

pickle_out = open("patches1.pickle","wb")
pickle.dump(all_patches, pickle_out)
pickle_out.close()


