# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 08:30:32 2015
Analyze data for thesis, ook figures van maken 
@author: Thomas
"""

#%pyplot qt 
import sys,os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from thesis_functions import *
code_dir = os.getcwd()
data_dir = code_dir[:-6] + "data\\flapping_data_10_04\\"
matplotlib.rc('xtick', labelsize=30) 
matplotlib.rc('ytick', labelsize=30) 
matplotlib.rc('font', **font)
#%% --------------------------------------------------------------------------------------------
file_dict = {}
data_dir = "D:/Mijn_documenten/Dropbox/ThomasMohren/data/flapping_data_10_04/"
file_dict[0] = "fd_flapping_halfbridge"
file_dict[1] = "fd_flapping_hb_xrot3_flap10"
file_dict[2] = "fd_flapping_hb_xrot5_flap10"
file_dict[3] = "fd_flapping_hb_yrot3_flap10"
file_dict[4] = "fd_flapping_quarterbridge"
file_dict[5] = "fd_flapping_qb_xrot3_flap10"
file_dict[6] = "fd_flapping_qb_xrot5_flap10"
file_dict[7] = "fd_flapping_qb_yrot3_flap10"

file_dict[8] = "fd_flapping_hb_xrot3p_flap10"
file_dict[9] = "fd_flapping_hb_yrot3p_flap10"
file_dict[10] = "fd_flapping_qb_xrot3p_flap10"
file_dict[11] = "fd_flapping_qb_yrot3p_flap10"
lincols = ['r','b','g']
add_f   = [0,0.09,0.18]
    
#%% constant rotation X half bridge 

tnew        = np.arange(1,69,0.001)
uplim       = round(NV_2_strain(50+511,"half"),5)
freq_axis   = [0,2.5,0,uplim/2]
fft_range   = [40000,45000]
titleplot   = ""
peak_mat    = np.zeros([3,10,2])
        
for i in np.arange(0,3,1):
    file_name   = file_dict[i]
    SD_data     = pd.read_csv(data_dir + file_name + ".CSV")
    y           = datafile_analyze("A2",SD_data,tnew)
    y           = NV_2_strain(y+511,"half") 
    fft_range = [40000,45000]
    freq,FFT    = do_fft(tnew,y,fft_range)
    flap_freq   = find_f(tnew,y,fft_range)
    t_nondim    = tnew*flap_freq
    f_nondim    = freq/flap_freq
    y_axis      = [510,520,-uplim,uplim]    
    line = add_plot(1,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    plt.setp(line,color = lincols[i])
#
    for l in np.arange(0,10,1):
        fft_range = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        f_nondim    = freq/flap_freq + add_f[i]
        peak_mat[i,l,:] = find_peaks(tnew,y,fft_range)[0]
        line = add_plot(1,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])
        
plt.figure(1)
plt.subplot(212)
for i in range(3):
    peak_mean = [np.mean(peak_mat[i,:,0]),np.mean(peak_mat[i,:,1])]
    peak_std = [np.std(peak_mat[i,:,0]),np.std(peak_mat[i,:,1])]
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
    
#%% Periodic rotation X half bridge 
titleplot   = "Periodic rotation $\Omega_X$, Half Bridge"
tnew        = np.arange(1,69,0.001)
uplim       = round(NV_2_strain(50+511,"half"),5)
freq_axis   = [0,2.5,0,uplim/2]
fft_range   = [40000,45000]
titleplot   = ""
peak_mat    = np.zeros([2,10,2])
offpeak_mat    = np.zeros([4,10])
i_it        = [0,0,0,0,0,0,0,0,1]
ff_mat = []
for k in [0,8]:
    i           = i_it[k]
    file_name   = file_dict[k]
    SD_data     = pd.read_csv(data_dir + file_name + ".CSV")
    y            = datafile_analyze("A2",SD_data,tnew)
    y           = NV_2_strain(y+511,"half") 
    freq,FFT    = do_fft(tnew,y,fft_range)
    flap_freq   = find_f(tnew,y,fft_range)
    t_nondim    = tnew*flap_freq
    f_nondim    = freq/flap_freq
    y_axis      = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]    
    y_axis      = [725,735,-uplim,uplim]  
    line        = add_plot(2,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    plt.setp(line,color = lincols[i])
#
    for l in np.arange(0,10,1):
        fft_range = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        ff_mat.append(flap_freq)
        f_nondim    = freq/flap_freq + add_f[i]
        peak_mat[i,l,:] = find_peaks(tnew,y,fft_range)[0]
        line = add_plot(2,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])
        if(i==1):
            ff_mat.append(flap_freq)
            offpeak_mat[:,l] = find_offpeaks(tnew,y,fft_range)
            offpeak_mat[:2,l] = offpeak_mat[:2,l]
plt.figure(2)
plt.legend(["0 rps","3 rps"])

plt.subplot(212)
for i in range(2):
    peak_mean = [np.mean(peak_mat[i,:,0]),np.mean(peak_mat[i,:,1])]
    peak_std = [np.std(peak_mat[i,:,0]),np.std(peak_mat[i,:,1])]
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
    
f_f = np.mean(ff_mat)
plt.axvline((f_f+3)/f_f+add_f[1],color = 'k')
plt.axvline((f_f-3)/f_f+add_f[1],color = 'k')
plt.axvline((f_f*2+3)/f_f+add_f[1],color = 'k')
plt.axvline((f_f*2-3)/f_f+add_f[1],color = 'k')

#print offpeak_mat
p_mat       = offpeak_mat
p_top       = p_mat[:,p_mat[2,:].argsort()][:,5:]
print p_top
peak_mean   = np.mean(p_top,1)
peak_std    = np.std(p_top,1)
print peak_mean
print f_f
plt.figure(2)
plt.subplot(212)
plt.errorbar(peak_mean[0:2]/f_f+add_f[1],peak_mean[2:4],peak_std[2:4], fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    

#%% constant rotation Y half bridge 

tnew = np.arange(1,69,0.001)
uplim = round(NV_2_strain(50+511,"half"),5)
freq_axis = [0,2.5,0,uplim/2]
fft_range = [40000,45000]
titleplot= "Constant rotation $\Omega_Y$, Half Bridge"
titleplot= ""
peak_mat = np.zeros([2,10,2])
i_it = [0,0,0,1]

for k in [0,3]:
    i = i_it[k]
    file_name = file_dict[k]
    SD_data = pd.read_csv(data_dir + file_name + ".CSV")
    
    y = datafile_analyze("A2",SD_data,tnew)
    y = NV_2_strain(y+511,"half") 
  
    freq,FFT= do_fft(tnew,y,fft_range)
    flap_freq = find_f(tnew,y,fft_range)
    t_nondim = tnew*flap_freq
    f_nondim = freq/flap_freq
    y_axis = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]    #    
    y_axis      = [730,740,-uplim,uplim]  
    line = add_plot(3,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    plt.setp(line,color = lincols[i])
#
    for l in np.arange(0,10,1):
        fft_range = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        f_nondim    = freq/flap_freq + add_f[i]
        peak_mat[i,l,:] = find_peaks(tnew,y,fft_range)[0]
        line = add_plot(3,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])    
    
plt.subplot(212)
for i in range(2):
    peak_mean = [np.mean(peak_mat[i,:,0]),np.mean(peak_mat[i,:,1])]
    peak_std = [np.std(peak_mat[i,:,0]),np.std(peak_mat[i,:,1])]
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    

#%% Periodic rotation Y half bridge 

titleplot= "Periodic rotation $\Omega_Y$, Half Bridge"
tnew = np.arange(1,69,0.001)
uplim = round(NV_2_strain(50+511,"half"),5)
freq_axis = [0,2.5,0,uplim/2]
fft_range = [40000,45000]
titleplot= ""
peak_mat = np.zeros([2,10,2])
offpeak_mat    = np.zeros([4,10])
i_it = [0,0,0,0,0,0,0,0,0,1]
ff_mat = []

for k in [0,9]:
    i = i_it[k]
    file_name = file_dict[k]
    SD_data = pd.read_csv(data_dir + file_name + ".CSV")
    
    y = datafile_analyze("A2",SD_data,tnew)
    y = NV_2_strain(y+511,"half") 
  
    freq,FFT= do_fft(tnew,y,fft_range)
    flap_freq = find_f(tnew,y,fft_range)
    t_nondim = tnew*flap_freq
    f_nondim = freq/flap_freq
    y_axis      = [225,235,-uplim,uplim]  
    ff_mat = []

    line = add_plot(4,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    plt.setp(line,color = lincols[i])
#
    for l in np.arange(0,10,1):
        fft_range = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        f_nondim    = freq/flap_freq + add_f[i]
        peak_mat[i,l,:] = find_peaks(tnew,y,fft_range)[0]
        line = add_plot(4,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])    
        ff_mat.append(flap_freq)
        if(i==1):
            ff_mat.append(flap_freq)
            offpeak_mat[:,l] = find_double_offpeaks(tnew,y,fft_range)
            offpeak_mat[:2,l] = offpeak_mat[:2,l]#/flap_freq11a
    
plt.figure(4)
plt.subplot(212)
for i in range(2):
    peak_mean = [np.mean(peak_mat[i,:,0]),np.mean(peak_mat[i,:,1])]
    peak_std = [np.std(peak_mat[i,:,0]),np.std(peak_mat[i,:,1])]
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
    
f_f = np.mean(ff_mat)
plt.axvline((f_f+3)/f_f+add_f[1],color = 'k')
plt.axvline((f_f-3)/f_f+add_f[1],color = 'k')
plt.axvline((f_f*2+3)/f_f+add_f[1],color = 'k')
plt.axvline((f_f*2-3)/f_f+add_f[1],color = 'k')

p_mat       = offpeak_mat
p_top       = p_mat[:,p_mat[2,:].argsort()]
f_f         = np.mean(ff_mat)
peak_mean   = np.mean(p_top,1)
peak_std    = np.std(p_top,1)
print peak_mean
plt.subplot(212)
plt.errorbar(peak_mean[0:2]/f_f+add_f[1],peak_mean[2:4],peak_std[2:4], fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
