# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 10:59:17 2015
@author: TMohren
Analyze data files for Msc. thesis 
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from thesis_functions import *
code_dir        = os.getcwd()
comsol_dir      = code_dir[:-6] + "Data\\Comsol_simulations\\"
matlab_dir      = code_dir[:-6] + "Data\\Matlab_simulations\\"
comsol_meshdir  = code_dir[:-6] + "Data\\Comsol_mesh_convergence\\"

#%% Load MATLAB results----------------------------------------------------------------------------------------
mat_dict    = {}
mat_dict[0] = 'MATLAB_300per0_flapf10.mat'
mat_dict[1] = 'MATLAB_030per0_flapf10.mat'
mat_dict[2] = 'MATLAB_300per1_flapf10.mat'
mat_dict[3] = 'MATLAB_030per1_flapf10.mat'
com_dict    = {}
com_dict[0] = 'flap10rot000.csv'
com_dict[1] = 'flap0rot300.csv'
com_dict[2] = 'flap10rot300.csv'
com_dict[3] = 'flap0rot030.csv'
com_dict[4] = 'flap10rot030.csv'
com_dict[5] = 'flap0rot300_per.csv'
com_dict[6] = 'flap10rot300_per.csv'
com_dict[7] = 'flap0rot030_per.csv'
com_dict[8] = 'flap10rot030_per.csv'

uplim       = 2e-5
freq_axis   = [0,2.5,0,uplim/2]
flap_name   = "flap10rot000.csv"
rot_name    = "flap0rot300.csv"
flaprot_name = "flap10rot300.csv"
fft_range   = [0,1000]
sample_dt   = 0.001
tnew        = np.arange( 0 , 1+sample_dt , sample_dt) 
t_nondim    = tnew*10
y_axis      = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]   
titleplot   = '-'
for i in range(4):
    #MATLAB part
    t,y         = read_matlabfiles( matlab_dir +mat_dict[i])
    freq,FFT    = do_fft(t,y,fft_range)
    flap_freq   = 10
    t_nondim    = t*flap_freq #+ 2
    f_nondim    = freq/flap_freq 
    add_plot(i+1,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    add_plot(i+1,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics =2)
    print "matlab"
    print  find_peak(tnew,y,fft_range)   
        
    #COMSOL part 
    t,y         = read_comsolfiles(comsol_dir ,com_dict[0],com_dict[2*i+1],com_dict[2*i+2],tnew)
    flapvec     = t*10 
    freq,FFT    = do_fft(t,y,fft_range)
    flap_freq   = 10
    t_nondim    = t*flap_freq
    f_nondim    = freq/flap_freq 
    add_plot(i+6,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    add_plot(i+6,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2)
    print "comsol"
    print  find_peak(tnew,y,fft_range)   
    
    plt.figure(i+1)
    plt.legend(["Euler-Lagrange","FEM"])

# fine tuning of axes 
uplim       = uplim/2
y_axis      = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]  
freq_axis   = [0,2.5,0,uplim/2]

plt.figure(3)
plt.subplot(211)
plt.axis(y_axis)
plt.yticks(np.linspace(y_axis[2],y_axis[3], 3))
plt.subplot(212)
plt.axis(freq_axis)
plt.yticks(np.linspace(freq_axis[2],freq_axis[3], 2))

plt.figure(4)
plt.subplot(211)
plt.axis(y_axis)
plt.yticks(np.linspace(y_axis[2],y_axis[3], 3))
plt.subplot(212)
plt.axis(freq_axis)
plt.yticks(np.linspace(freq_axis[2],freq_axis[3], 2))
plt.xticks(np.linspace(freq_axis[0],freq_axis[1], 6))

plt.figure(8)
plt.subplot(211)
plt.axis(y_axis)
plt.yticks(np.linspace(y_axis[2],y_axis[3], 3))
plt.subplot(212)
plt.axis(freq_axis)
plt.yticks(np.linspace(freq_axis[2],freq_axis[3], 2))

plt.figure(9)
plt.subplot(211)
plt.axis(y_axis)
plt.yticks(np.linspace(y_axis[2],y_axis[3], 3))
plt.subplot(212)
plt.axis(freq_axis)
plt.yticks(np.linspace(freq_axis[2],freq_axis[3], 2))