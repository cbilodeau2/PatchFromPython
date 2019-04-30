# -*- coding: utf-8 -*-
"""
Created on Thu Mar  7 15:42:50 2019

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
import keras

from keras.layers import Input, Dense
from keras.models import Model

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
        
def takeDistance(Pair):
    return Pair.dist

# This is the size of our encoded representations
encoding_dim = 32

# Load in patch representation
pickle_in = open("patches2.pickle","rb")
all_patches = pickle.load(pickle_in)

# Prune distances (only use closest n distances):
prune = 3000
i=0
for patch in all_patches:
    i+=1
    #print('Analyzing Patch',i)
    patch.pairlist.sort(key=takeDistance)
    patch.pairlist = patch.pairlist[0:prune]
 
# BEGIN AUTOENCODER: (From SampleAutoEncoder):    

input_img = Input(shape=(3000,))
# "encoded" is the encoded representation of the input
encoded = Dense(encoding_dim, activation='relu')(input_img)
# "decoded" is the lossy reconstruction of the input
decoded = Dense(3000, activation='sigmoid')(encoded)

# this model maps an input to its reconstruction
autoencoder = Model(input_img, decoded)

# create a placeholder for an encoded (x-dimensional) input
encoded_input = Input(shape=(encoding_dim,))
# retrieve the last layer of the autoencoder model
decoder_layer = autoencoder.layers[-1]
# create the decoder model
decoder = Model(encoded_input, decoder_layer(encoded_input))

autoencoder.compile(optimizer='adadelta', loss='binary_crossentropy')

#x_train = x_train.astype() / 10.
#x_test = x_test.astype('float32') / 10.
#x_train = x_train.reshape((len(x_train), np.prod(x_train.shape[1:])))
#x_test = x_test.reshape((len(x_test), np.prod(x_test.shape[1:])))
#print(x_train.shape)
#print(x_test.shape)

