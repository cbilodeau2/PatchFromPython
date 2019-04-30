# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 18:40:36 2019

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


#ax = plt.axes()
z_slice=70
plt.scatter(gridx,gridy,s=carbon_density[:,:,z_slice].T*100) #,carbon_density)
#plt.xlim([-15,-10])
#plt.ylim([-40,-30])
plt.show()

plt.scatter(gridx,gridy,s=nitrogen_density[:,:,z_slice].T*100) #,carbon_density)
plt.show()

plt.scatter(gridx,gridy,s=oxygen_density[:,:,z_slice].T*100) #,carbon_density)
plt.show()

plt.scatter(gridx,gridy,s=sulfur_density[:,:,z_slice].T*100) #,carbon_density)
plt.show()

#ax.plot(carbon_density[:,61,45])
#ax.plot(carbon_density[:,60,45])
#ax.plot(carbon_density[:,59,45])
#ax.plot(carbon_density[:,58,45])