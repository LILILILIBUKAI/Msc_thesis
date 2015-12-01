# -*- coding: utf-8 -*-
"""
Created on Thu Oct 08 17:01:05 2015

@author: Thomas
"""
import sys,os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy.fftpack import fft, fftfreq, fftshift
from scipy.interpolate import interp1d
import scipy as sp
from numpy import genfromtxt

plt.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.rc('text', usetex=True)
font = {'family' : 'normal',
        'size'   : 18}
plt.rc('font', **font)

def datafile_analyze(name,SD_data,tnew):
    y_raw = SD_data[name]
    t = SD_data["micros"]/1e6
    y_fun        = interp1d(t,y_raw,'linear')
    y             = y_fun(tnew) - np.mean(y_fun(tnew))
    return y

def read_comsolfiles(file_dir,flap_name,rot_name,flaprot_name,tnew):
    """Read 3 .CSV files, do some subtractions and return one vector out
    Parameters
    ----------
    file_dir : string
    """
    rot    = genfromtxt(file_dir+rot_name, delimiter=',',skip_header=5)
    flap     = genfromtxt(file_dir+flap_name,delimiter=',',skip_header =5)
    flaprot = genfromtxt(file_dir+flaprot_name,delimiter=',',skip_header =5)
    
    del_flap    = flap[:,1]- flap[:,2]
    del_rot     = rot[:,1] - rot[:,2]
    del_flaprot = flaprot[:,1] - flaprot[:,2]
    idf_func        = interp1d(flap[:,0],del_flap,'linear')
    idflap          = idf_func(tnew)
    idr_func        = interp1d(rot[:,0],del_rot,'linear')
    idrot           = idr_func(tnew)
    idfr_func       =interp1d(flaprot[:,0],del_flaprot,'linear')
    idflaprot       = idfr_func(tnew)
    flaprot_cor = idflaprot-idflap -idrot
    return tnew, flaprot_cor
    
def read_forcefiles(file_dir,filename):
    data =  genfromtxt(file_dir+filename, delimiter=',',skip_header=3)
    t = data[:,0]
    flapvec = t*10
    vec_toplot = data[:,[8,11,14]]
    return flapvec,vec_toplot
    

def read_matlabfiles(fullname):
    matdict = sp.io.loadmat(fullname) 
    t_vec = matdict['T']
    dstrain = matdict['strainy'][:,1,1]-matdict['strainy'][:,9,1]    
    return t_vec,dstrain

def do_fft(t,vec,fft_range,sample_dt=0.001):
    N = fft_range[1] - fft_range[0]    
    signal          = vec[fft_range[0]:fft_range[1]] - np.mean(vec[fft_range[0]:fft_range[1]])
    freq            = scipy.fftpack.fftfreq(N,sample_dt)
    FFT_full        = 2.0/N  *scipy.fft( signal )
    FFT             = abs(FFT_full)          
    phase           = np.angle(FFT_full)
    return freq,FFT
    
def find_f(t,vec,fft_range,sample_dt=0.001):
    freq,FFT = do_fft(t,vec,fft_range,sample_dt)
    maxind      = FFT[0:len(FFT)/2].argsort()[-1:]#[::-1]
    dominant_freq   = freq[maxind]
    return dominant_freq

def NV_2_strain(NV,wheatstone):
    K = 2.0
    amp_gain = 1000.0
    VG = (NV*5.0/1023.0-2.5)/amp_gain
    if wheatstone == "quarter":
        epsilon = (2.0/5.0*VG)/(K*(0.5-VG/5.0))
    elif wheatstone == "half":
        epsilon = 4.0/5.0*VG/K
    return epsilon 
    
def find_peaks(t,vec,fft_range,sample_dt=0.001):
    freq,FFT        = do_fft(t,vec,fft_range,sample_dt)
    maxind          = FFT[0:len(FFT)/2].argsort()[-10:][::-1]
    max2            = [i for i in maxind[1:] if i >= maxind[0]+10]
    maxvals         = [FFT[maxind[0]] , FFT[max2[0] ] ]
    maxfreqs        = [freq[maxind[0]] , freq[max2[0]]  ] 
    return maxvals,maxfreqs

def find_peak(t,vec,fft_range,sample_dt=0.001):
    freq,FFT        = do_fft(t,vec,fft_range,sample_dt)
    maxind          = FFT[0:len(FFT)/2].argsort()[-10:][::-1]
    maxvals         = [FFT[maxind[0]]  ] 
    maxfreqs        = [freq[maxind[0]]  ] 
    return maxvals,maxfreqs
    
def find_offpeaks(t,vec,fft_range,sample_dt=0.001):
    freq,FFT    = do_fft(t,vec,fft_range,sample_dt)
    maxind      = FFT[0:len(FFT)/2].argsort()[-100:][::-1]
    maxhigh     = [i for i in maxind[1:] if (   (i >= maxind[0]+8)& (i <= maxind[0]+25)   )]
    maxlow      = [i for i in maxind[1:] if (  (i >= maxind[0]-25) & (i <= maxind[0]-8)    )]
    return freq[maxlow[0]],freq[maxhigh[0]],FFT[maxlow[0]],FFT[maxhigh[0]]
     
def find_offpeaks(t,vec,fft_range,sample_dt=0.001):
    freq,FFT = do_fft(t,vec,fft_range,sample_dt)
    maxind      = FFT[0:len(FFT)/2].argsort()[-100:][::-1]
    maxhigh = [i for i in maxind[1:] if (   (i >= maxind[0]+8)& (i <= maxind[0]+25)   )]
    maxlow = [i for i in maxind[1:] if (  (i >= maxind[0]-25) & (i <= maxind[0]-8)    )]
    return freq[maxlow[0]],freq[maxhigh[0]],FFT[maxlow[0]],FFT[maxhigh[0]]
    
def find_double_offpeaks(t,vec,fft_range,sample_dt=0.001):
    freq,FFT        = do_fft(t,vec,fft_range,sample_dt)
    maxind          = FFT[0:len(FFT)/2].argsort()[-100:]
    maxhigh         = [i for i in maxind[1:] if (   (i >= maxind[-1]*2+10)& (i <= maxind[-1]*2+20)   )]
    maxlow          = [i for i in maxind[1:] if (  (i >= maxind[-1]*2-22) & (i <= maxind[-1]*2-10)    )]
    return freq[maxlow[0]],freq[maxhigh[0]],FFT[maxlow[0]],FFT[maxhigh[0]]
    
def add_plot(fig_nr,subplt,t,vec,tit='none',ax_vec='none',xlab = 'Time [s]',ylab = '$\Delta \epsilon_{yy}$',n_ytics=3,tickstyle='sci',figsize = [15,5]):
    """Create a subplot on a specified figure
    Parameters
    ----------
    fig_nr : int
        The figure on which to plot
    subplt : int
        The subplot code, e.g. 131
    t : vector
        The x axis vector    
    vec: vector
        The y axis vector        
    tit : string
        Title of the plot        
    ax_vec : vector of length 4 (default='none')
        Tries to modify the plot axis, only works if vector is correctly specified    
    xlab : string (default='Time [s]')
        Name of xlabel    
    ylab : string (default '$\Delta \epsilon_{yy}$')
        Name of ylabel        
    n_ytics : int (default=3)
        Sets number of ticks on ylabel, takes min and max from ax_vec
        Only works if ax_vec is correctly specified        
    tickstyle : string (default='sci')
        Sets the tickstyle for the y-axis 
    """
    plt.figure(fig_nr,figsize)
    plt.subplot(subplt)
    line = plt.plot(t,vec)
    plt.title(tit)
    try:
        plt.axis(ax_vec)
        plt.yticks(np.linspace(ax_vec[2],ax_vec[3], n_ytics))
        plt.xticks(np.linspace(ax_vec[0],ax_vec[1], 3))
    except:
        print 'Axis not specified'
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.ylabel(ylab)
    plt.ticklabel_format(style =tickstyle,axis='y',scilimits=(0,0))
    plt.tight_layout()
    return line
