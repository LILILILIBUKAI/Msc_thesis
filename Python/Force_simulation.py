# -*- coding: utf-8 -*-
"""
Created on Thu Aug 06 11:30:55 2015
@author: Thomas
Simulation of forces on a point mass during flapping and rotating 
"""
from thesis_functions import *
import numpy as np 
import matplotlib.pyplot as plt

def force_comp(t,f_f,v_r,per):
    """Compute forces on a point mass flapping and rotating
    Parameters
    ----------
    file_dir : string
    """
    r = 0.025
    amp_f= np.deg2rad(15)
    amp_rot = np.deg2rad(30)
    n = len(t)    
    
    phi = amp_f * np.cos(2*np.pi*f_f* t)
    phi_d1 = amp_f * -np.sin(2*np.pi*f_f*t)* 2*np.pi*f_f
    phi_d2 = amp_f * -np.cos(2*np.pi*f_f*t) * (2*np.pi*f_f)**2 
    f_r = max(v_r)
    if per == 1:
        theta = amp_rot * np.outer(v_r , np.sin(2*np.pi*f_r*t))
        omega = amp_rot * np.outer(v_r , np.cos(2*np.pi*f_r*t) * (2*np.pi*f_r))
        omega_d1 =  amp_rot * np.outer(v_r , np.cos(2*np.pi*f_r*t) * (2*np.pi*f_r)**2)
    elif per == 0:
        theta = f_r*2*np.pi * np.outer(v_r,t)
        omega = np.outer(2*np.pi*v_r,np.ones([n,1]))
        omega_d1 = np.zeros([3,n])

    x = r*np.sin(phi)
    y = r*np.cos(phi)
    z = np.zeros([1,n])
    p = np.vstack([x,y,z])

    vx = r*np.cos(phi)*phi_d1
    vy = r*-np.sin(phi)*phi_d1
    vz = np.zeros([1,n])
    v = np.vstack([vx,vy,vz])      
        
    ax = r* ( -np.sin(phi)*phi_d1**2 + np.cos(phi)*phi_d2)        
    ay = r* ( -np.cos(phi)*phi_d1**2 - np.sin(phi)*phi_d2)
    az = np.zeros([1,n])
    a_in = np.vstack([ax,ay,az])

    a_cor = np.zeros([3,n])
    a_cent = np.zeros([3,n])
    a_tan = np.zeros([3,n])
    
    for i in range(n):
        a_cor[:,i] = -2*np.cross(omega[:,i],v[:,i])
        a_cent[:,i] = -np.cross(omega[:,i],np.cross(omega[:,i],p[:,i]))
        a_tan[:,i] = -np.cross(omega_d1[:,i],p[:,i])
    acc = {'in':a_in}
    acc['cor'] = a_cor
    acc['cent'] = a_cent
    acc['tan'] = a_tan
    return phi,acc


rho_p   = 1180
m       = 0.02*0.05*1.27e-4*rho_p
f_f     = 10
dt      = 1.0/1000
t       = np.arange(0,0.2+dt,dt)+0.025
f       = (t-0.025)*f_f
v_r     = np.array([3,0,0])

combinations    = np.hstack( [ np.vstack([ 3*np.eye(3),3*np.eye(3)]), np.array([0,0,0,1,1,1])[:,None] ] )
namelist        = ['x0','y0','z0','x1','y1','z1']
acc_dict        = {}

for i in range(len(combinations)):
    phi,acc_dict[namelist[i]]= force_comp(t,10,combinations[i,0:3],combinations[i,3])

for per in range(2):
    plt.figure(per+1, [10,10])
    for c in range(3):
        print c
        plt.subplot(4,3,c+1)
        plt.plot(f,np.rad2deg(phi))
        ax_vec = [0,2,-15,15]
        plt.axis(ax_vec)         
        plt.xticks(np.linspace(ax_vec[0],ax_vec[1], 3))
        plt.yticks(np.linspace(ax_vec[2],ax_vec[3], 3))
        plt.tight_layout()
        plt.legend(['Stroke angle'])
        for r in np.arange(1,4):
            plt.subplot(4,3,r*3+c+1)
            plt.plot( f, m*acc_dict[namelist[c+3*per]]['in'][r-1],'g')
            plt.plot( f, m*acc_dict[namelist[c+3*per]]['tan'][r-1],'r')
            plt.plot( f, m*acc_dict[namelist[c+3*per]]['cor'][r-1],'y')
            plt.plot( f, m*acc_dict[namelist[c+3*per]]['cent'][r-1],'b')
            ax_vec = [0,2,-0.004,0.004]
            plt.axis(ax_vec)         
            plt.xticks(np.linspace(ax_vec[0],ax_vec[1], 3))
            plt.yticks(np.linspace(ax_vec[2],ax_vec[3], 3))
            plt.ticklabel_format(style ='sci',axis='y', scilimits=(0,0) )
            plt.tight_layout()
            total = m* (acc_dict[namelist[c+3*per]]['cent'][r-1] + \
            acc_dict[namelist[c+3*per]]['tan'][r-1]+\
            acc_dict[namelist[c+3*per]]['cor'][r-1]+\
            acc_dict[namelist[c+3*per]]['in'][r-1] )
            plt.plot( f, total ,'k')
    plt.legend(['Linear acceleration','Angular acceleration','Coriolis force','Centrifugal force','Total'])