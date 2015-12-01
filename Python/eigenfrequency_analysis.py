# -*- coding: utf-8 -*-
"""
Created on Mon Oct 05 12:32:44 2015
@author: Thomas
Analyze fling test, find eigenfrequencies
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from thesis_functions import *
code_dir    = os.getcwd()
file_dir    = code_dir[:-6] + "Data\\flapping_data_10_04\\"
font        = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 30}
tnew        = np.arange(1,69,0.001)
uplim       = round(NV_2_strain(200+511,"half"),4)
freq_axis   = [0,150,0,uplim/5]
fft_range   = [41400,42200]
y_axis      = [tnew[fft_range[0]],tnew[fft_range[1]],-uplim, uplim]
file_name   = "fd_wing_flingtest"
SD_data     = pd.read_csv(file_dir + file_name + ".CSV")
y_raw       = datafile_analyze("A1",SD_data,tnew)
y           = NV_2_strain(y_raw+511,"half") 
t_plot      = tnew-42.40
y_axis      = [0,0.8,-uplim, uplim]
add_plot(1,121,t_plot,y,"",y_axis,xlab = 't [s]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci')
freq,FFT    = do_fft(tnew,y,fft_range)
maxind      = FFT[0:len(FFT)/2].argsort()[-10:][::-1]
print freq[maxind]

add_plot(1,122,freq,FFT,"",freq_axis,xlab = 'Frequency [Hz]',ylab = '$\Delta \epsilon$',n_ytics=3,tickstyle='sci',figsize = [9,9])
plt.plot([22.5,22.5],[0,1])
plt.plot([122.5,122.5],[0,1])