# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 15:12:08 2019
Takes patches1.pickle as input (list of patches) and adds pairlist
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

# Calculate center of geometry for a given residue:
def res_cog(residue):
    coord = [residue.get_list()[i].get_coord() for i in range(0,np.shape(residue.get_list())[0])]
    cog = np.mean(coord,axis=0)
    return cog

# Class defined for patches:
class Patch:
    def __init__(self,center,res_list,identification,pairlist):
        self.center = center
        self.res_list = res_list
        self.num = identification
        self.pairlist = pairlist

class Pair:
    def __init__(self,residue1,residue2,atom1,atom2,distance):
        self.res1 = residue1
        self.res2 = residue2
        self.atom1 = atom1
        self.atom2 = atom2
        self.dist = distance

start = time.time()

pickle_in = open("patches1.pickle","rb")
all_patches = pickle.load(pickle_in)
i = 0
for patch in all_patches:
    i+=1
    print('Analyzing Patch',i)
    pairlist = []
    for residue1 in patch.res_list:
        for atom1 in residue1.get_list():
            if not (atom1.name[0]=='H'):
                for residue2 in patch.res_list:
                    for atom2 in residue2.get_list():
                        if not (atom2.name[0]=='H'):
                            #print(atom1)
                            dist = np.linalg.norm(atom1.get_coord()-atom2.get_coord())
                            pairlist.append(Pair(residue1,residue2,atom1,atom2,dist))
    patch.pairlist = pairlist



end = time.time()

pickle_out = open("patches2.pickle","wb")
pickle.dump(all_patches, pickle_out)
pickle_out.close()

print('Time Elapsed:',end-start)


