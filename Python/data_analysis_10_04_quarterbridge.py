# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 08:30:32 2015
@author: Thomas
Analysis of quarterbridge experiments
"""
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
lincols     = ['r','b','g']
add_f       = [0,0.09,0.18]

#%% ----------------------------------X_const_qb_ratio find 
# recurring parameters
tnew        = np.arange(1,69,0.001)
fft_range   = [40000,45000]
titleplot   = ""
# refreshed parameters 
uplim       = round(NV_2_strain(500+511,"quarter"),4)
freq_axis   = [0,2.5,0,uplim]
peak_mat    = np.zeros([3,10,2])
i           = 0
k_range     = [4,5,6]
for k in k_range:
    file_name   = file_dict[k]
    SD_data     = pd.read_csv(data_dir + file_name + ".CSV")
    y           = datafile_analyze("A0",SD_data,tnew)
    y           = NV_2_strain(y+511,"quarter") 
    freq,FFT    = do_fft(tnew,y,fft_range)
    flap_freq   = find_f(tnew,y,fft_range)
    t_nondim    = tnew*flap_freq
    f_nondim    = freq/flap_freq
    y_axis      = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]    
    line        = add_plot(1,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])    
    plt.setp(line,color = lincols[i])
#
    for l in np.arange(0,10,1):
        fft_range   = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        f_nondim    = freq/flap_freq + add_f[i]
        peak_mat[i,l,:] = find_peaks(tnew,y,fft_range)[0]
        line = add_plot(1,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])
    i = i +1 
ratio1 = np.mean(peak_mat[1,[peak_mat[1,:,0].argsort()[5:]],0])/np.mean(peak_mat[0,[peak_mat[0,:,0].argsort()[5:]],0])
ratio2 = np.mean(peak_mat[2,[peak_mat[2,:,0].argsort()[5:]],0])/np.mean(peak_mat[1,[peak_mat[1,:,0].argsort()[5:]],0])
plt.figure(1)
plt.subplot(212)
for i in range(3):
    p_mat       = peak_mat[i,:,:]
    p_top       = p_mat[p_mat[:,0].argsort()][5:]
    peak_mean   = np.mean(p_top,0)
    peak_std    = np.std(p_top,0)
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
#%% ----------------------X_const_compensated------------------------------------------------------------------------------
# refreshed parameters 
uplim       = round(NV_2_strain(50+511,"quarter"),5)
freq_axis   = [0,2.5,0,uplim]
peak_mat    = np.zeros([3,10,2])
i           = 0 
for k in k_range:
    file_name = file_dict[k]
    SD_data = pd.read_csv(data_dir + file_name + ".CSV")
    y = (datafile_analyze("A0",SD_data,tnew) - datafile_analyze("A1",SD_data,tnew))
    if (i ==1):
        y = y/ratio1
    elif ( i ==2):
        y = y/ratio2
    y           = NV_2_strain(y+511,"quarter") 
    freq,FFT    = do_fft(tnew,y,fft_range)
    flap_freq   = find_f(tnew,y,fft_range)
    t_nondim    = tnew*flap_freq
    f_nondim    = freq/flap_freq
    y_axis      = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]    
    y_axis      = [490,500,-uplim,uplim]  
    line        = add_plot(2,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    plt.setp(line,color = lincols[i])

    for l in np.arange(0,10,1):
        fft_range   = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        f_nondim    = freq/flap_freq + add_f[i]
        peak_mat[i,l,:] = find_peaks(tnew,y,fft_range)[0]
        line        = add_plot(2,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])
    i = i +1
plt.figure(2)
plt.subplot(212)
for i in range(3):
    p_mat       = peak_mat[i,:,:]
    p_top       = p_mat[p_mat[:,0].argsort()][5:]
    peak_mean   = np.mean(p_top,0)
    peak_std    = np.std(p_top,0)
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
#%% ----------------------------------Y_const_qb_ratio find 
uplim       = round(NV_2_strain(500+511,"quarter"),4)
freq_axis   = [0,2.5,0,uplim]
peak_mat    = np.zeros([2,10,2])
i           = 0
k_range     = [4,7]
for k in k_range:
    file_name = file_dict[k]
    SD_data = pd.read_csv(data_dir + file_name + ".CSV")
    y           = datafile_analyze("A0",SD_data,tnew)
    y           = NV_2_strain(y+511,"quarter") 
    freq,FFT    = do_fft(tnew,y,fft_range)
    flap_freq   = find_f(tnew,y,fft_range)
    t_nondim    = tnew*flap_freq
    f_nondim    = freq/flap_freq
    y_axis      = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]    
    y_axis      = [210,220,-uplim,uplim]  
    line        = add_plot(4,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    plt.setp(line,color = lincols[i])
    for l in np.arange(0,10,1):
        fft_range = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        f_nondim    = freq/flap_freq + add_f[i]
        peak_mat[i,l,:] = find_peaks(tnew,y,fft_range)[0]
        line = add_plot(4,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])
    i = i +1 
ratio = np.mean(peak_mat[1,[peak_mat[1,:,0].argsort()[5:]],0])/np.mean(peak_mat[0,[peak_mat[0,:,0].argsort()[5:]],0])

plt.figure(4)
plt.subplot(212)
for i in range(2):
    p_mat       = peak_mat[i,:,:]
    p_top       = p_mat[p_mat[:,0].argsort()][5:]
    peak_mean   = np.mean(p_top,0)
    peak_std    = np.std(p_top,0)
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
#%% ----------------------Y_const_compensated------------------------------------------------------------------------------
# refreshed parameters 
uplim       = round(NV_2_strain(50+511,"quarter"),5)
freq_axis   = [0,2.5,0,uplim]
peak_mat    = np.zeros([2,10,2])
i           = 0 
for k in k_range:
    file_name = file_dict[k]
    SD_data = pd.read_csv(data_dir + file_name + ".CSV")
    y = (datafile_analyze("A0",SD_data,tnew) - datafile_analyze("A1",SD_data,tnew))
    if (i ==1):
        y = y/ratio
    y           = NV_2_strain(y+511,"quarter") 
    freq,FFT    = do_fft(tnew,y,fft_range)
    flap_freq   = find_f(tnew,y,fft_range)
    t_nondim    = tnew*flap_freq
    f_nondim    = freq/flap_freq
    y_axis      = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]    
    y_axis      = [490,500,-uplim,uplim]  
    line        = add_plot(5,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    plt.setp(line,color = lincols[i])
    for l in np.arange(0,10,1):
        fft_range = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        f_nondim    = freq/flap_freq + add_f[i]
        peak_mat[i,l,:] = find_peaks(tnew,y,fft_range)[0]
        line        = add_plot(5,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])
    i = i +1
plt.figure(5)
plt.subplot(212)
for i in range(2):
    p_mat       = peak_mat[i,:,:]
    p_top       = p_mat[p_mat[:,0].argsort()][5:]
    peak_mean   = np.mean(p_top,0)
    peak_std    = np.std(p_top,0)
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
#%% ----------------------------------X_periodic_qb_ratio find 
# refreshed parameters 
uplim       = round(NV_2_strain(500+511,"quarter"),4)
freq_axis   = [0,2.5,0,uplim]
peak_mat    = np.zeros([2,10,2])
i           = 0
k_range     = [4,10]
for k in k_range:
    file_name = file_dict[k]
    SD_data = pd.read_csv(data_dir + file_name + ".CSV")
    y = datafile_analyze("A0",SD_data,tnew)
    y = NV_2_strain(y+511,"quarter") 
    freq,FFT    = do_fft(tnew,y,fft_range)
    flap_freq   = find_f(tnew,y,fft_range)
    t_nondim    = tnew*flap_freq
    f_nondim    = freq/flap_freq
    y_axis      = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]    
    y_axis      = [550,560,-uplim,uplim]  
    line        = add_plot(7,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    plt.setp(line,color = lincols[i])
    for l in np.arange(0,10,1):
        fft_range = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        f_nondim    = freq/flap_freq + add_f[i]
        peak_mat[i,l,:] = find_peaks(tnew,y,fft_range)[0]
        line = add_plot(7,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])
    i = i +1 
ratio = np.mean(peak_mat[1,[peak_mat[1,:,0].argsort()[5:]],0])/np.mean(peak_mat[0,[peak_mat[0,:,0].argsort()[5:]],0])
plt.figure(7)
plt.subplot(212)
for i in range(2):
    p_mat       = peak_mat[i,:,:]
    p_top       = p_mat[p_mat[:,0].argsort()][5:]
    peak_mean   = np.mean(p_top,0)
    peak_std    = np.std(p_top,0)
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
#%% ----------------------X_const_compensated------------------------------------------------------------------------------
# refreshed parameters 
uplim       = round(NV_2_strain(50+511,"quarter"),5)
freq_axis   = [0,2.5,0,uplim]
peak_mat    = np.zeros([2,10,2])
offpeak_mat = np.zeros([4,10])
i           = 0
for k in k_range:
    file_name = file_dict[k]
    SD_data = pd.read_csv(data_dir + file_name + ".CSV")
    y = (datafile_analyze("A0",SD_data,tnew) - datafile_analyze("A1",SD_data,tnew))
    if (i ==1):
        y = y/ratio
    y           = NV_2_strain(y+511,"quarter") 
    freq,FFT    = do_fft(tnew,y,fft_range)
    flap_freq   = find_f(tnew,y,fft_range)
    t_nondim    = tnew*flap_freq
    f_nondim    = freq/flap_freq
    y_axis      = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]    
    y_axis      = [490,500,-uplim,uplim]  
    line = add_plot(8,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    plt.setp(line,color = lincols[i])
    for l in np.arange(0,10,1):
        fft_range = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        f_nondim    = freq/flap_freq + add_f[i]
        if(i==1):
            offpeak_mat[:,l] = find_offpeaks(tnew,y,fft_range)
            offpeak_mat[:2,l] = offpeak_mat[:2,l]/flap_freq
        peak_mat[i,l,:] = find_peaks(tnew,y,fft_range)[0]
        line = add_plot(8,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])
    i = i +1
    
plt.figure(8)
plt.subplot(212)
for i in range(2):
    p_mat       = peak_mat[i,:,:]
    p_top       = p_mat[p_mat[:,0].argsort()]
    peak_mean   = np.mean(p_top,0)
    peak_std    = np.std(p_top,0)
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
p_mat       = offpeak_mat
p_top       = p_mat[:,p_mat[2,:].argsort()][:,5:]

peak_mean   = np.mean(p_top,1)
peak_std    = np.std(p_top,1)
print peak_mean
plt.figure(8)
plt.subplot(212)
plt.errorbar(peak_mean[0:2]+add_f[1],peak_mean[2:4],peak_std[2:4], fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    

#%% ----------------------------------Y_per_qb_ratio find 
# refreshed parameters 
uplim       = round(NV_2_strain(500+511,"quarter"),4)
freq_axis   = [0,2.5,0,uplim]
peak_mat    = np.zeros([2,10,2])
i           = 0
offpeak_mat = np.zeros([4,10])
k_range     = [4,11]

for k in k_range:
    file_name = file_dict[k]
    SD_data = pd.read_csv(data_dir + file_name + ".CSV")
    y = datafile_analyze("A0",SD_data,tnew)
    y = NV_2_strain(y+511,"quarter") 
    freq,FFT    = do_fft(tnew,y,fft_range)
    flap_freq   = find_f(tnew,y,fft_range)
    t_nondim    = tnew*flap_freq
    f_nondim    = freq/flap_freq
    y_axis      = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]    
    line        = add_plot(11,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    plt.setp(line,color = lincols[i])
    for l in np.arange(0,10,1):
        fft_range = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        f_nondim    = freq/flap_freq + add_f[i]
        peak_mat[i,l,:] = find_peaks(tnew,y,fft_range)[0]
        line = add_plot(11,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])
    i = i +1 
ratio = np.mean(peak_mat[1,[peak_mat[1,:,0].argsort()[5:]],0])/np.mean(peak_mat[0,[peak_mat[0,:,0].argsort()[5:]],0])

plt.figure(11)
plt.subplot(212)
for i in range(2):
    p_mat       = peak_mat[i,:,:]
    p_top       = p_mat[p_mat[:,0].argsort()][5:]
    peak_mean   = np.mean(p_top,0)
    peak_std    = np.std(p_top,0)
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
#%% ----------------------Y_per_compensated------------------------------------------------------------------------------
# refreshed parameters 
uplim       = round(NV_2_strain(50+511,"quarter"),5)
freq_axis   = [0,2.5,0,uplim]
offpeak_mat = np.zeros([4,10])
i           = 0 

for k in k_range:
    file_name   = file_dict[k]
    SD_data     = pd.read_csv(data_dir + file_name + ".CSV")
    y           = (datafile_analyze("A0",SD_data,tnew) - datafile_analyze("A1",SD_data,tnew))
    if (i ==1):
        y = y/ratio
    y           = NV_2_strain(y+511,"quarter") 
    freq,FFT    = do_fft(tnew,y,fft_range)
    flap_freq   = find_f(tnew,y,fft_range)
    t_nondim    = tnew*flap_freq
    f_nondim    = freq/flap_freq
    y_axis      = [t_nondim[fft_range[0]],t_nondim[fft_range[1]],-uplim,uplim]    
    y_axis      = [700,710,-uplim,uplim]  
    line        = add_plot(12,211,t_nondim,y,titleplot,y_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
    plt.setp(line,color = lincols[i])
    find_double_offpeaks(tnew,y,fft_range)
    for l in np.arange(0,10,1):
        fft_range = [17999+5000*l,22999+5000*l]
        freq,FFT    = do_fft(tnew,y,fft_range)
        flap_freq   = find_f(tnew,y,fft_range)
        f_nondim    = freq/flap_freq + add_f[i]
        if(i==1):
            offpeak_mat[:,l] = find_double_offpeaks(tnew,y,fft_range)
            offpeak_mat[:2,l] = offpeak_mat[:2,l]/flap_freq
        peak_mat[i,l,:]     = find_peaks(tnew,y,fft_range)[0]
        line                = add_plot(12,212,f_nondim,abs(FFT),titleplot,freq_axis,xlab = '$\\frac{f}{f_{flap}}$ [-]',ylab = '$\Delta \epsilon$',n_ytics=2,tickstyle='sci',figsize = [10,10])
        plt.setp(line,color = lincols[i])
    i = i +1
    
plt.figure(12)
plt.subplot(212)
for i in range(2):
    p_mat       = peak_mat[i,:,:]
    p_top       = p_mat[p_mat[:,0].argsort()][5:]
    peak_mean   = np.mean(p_top,0)
    peak_std    = np.std(p_top,0)
    plt.errorbar([1+add_f[i],2+add_f[i]],peak_mean,peak_std, fmt='o',linewidth=2,color='k', capsize=10,capthick=2, zorder=3)    
