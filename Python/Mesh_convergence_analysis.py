# -*- coding: utf-8 -*-
"""
Created on Thu Oct 08 12:08:15 2015

@author: Thomas
"""
#import xlrd
#import sys
import numpy as np
import matplotlib.pyplot as plt
import os
from thesis_functions import *
code_dir    = os.getcwd()
comsol_meshdir = code_dir[:-6] + "Data\\Comsol_mesh_convergence\\"
sample_dt   = 0.001
tnew        = np.arange( 0 , 0.5+sample_dt , sample_dt) 

#%%  
file_dict_c = {}
file_dict_c[0] = ".5s_20_8_4_HT.CSV"
file_dict_c[1] = ".5s_30_8_6_HT.CSV"
file_dict_c[2] = ".5s_50_10_3_HT.CSV"
file_dict_c[3] = ".5s_60_15_8_HT.CSV"

plt.figure(1)
for i in range(0,len(file_dict_c)):
    name    = file_dict_c[i]
    noflap  = ".5s_50_10_4_noflap.CSV"
    norot   = "flap10rot000.CSV"
    t,de    = read_comsolfiles(comsol_meshdir,norot,noflap,name,tnew)
    plt.plot(t*10,de)    
uplim = 1.5e-5
y_axis = [0,9,uplim,-uplim]
plt.xlabel('$\\frac{f}{f_{flap}}$')
plt.ylabel('$\Delta \epsilon$')
plt.axis(y_axis)
plt.legend(["20x8x4","30x8x6","50x10x3","60x15x8","80x12x4"],loc = 1)
plt.yticks(np.linspace(y_axis[2],y_axis[3], 3))
plt.ticklabel_format(axis='y',scilimits=(0,0))