# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 12:43:06 2016

@author: c0906755

read in data file
optionally remove binaries
plot cluster
"""

import numpy as np
import matplotlib.pyplot as plt
import glob as gb
import cubefuncs as cf
plt.rcParams.update({'font.size': 14, 'font.family':'serif','text.usetex':False})



"""
fname=raw_input("Filename:")

has_bins=raw_input("Remove binaries? (y/n)")

if has_bins.lower()=='y':
    has_bins=True
elif has_bins.lower()=='n':
    has_bins=False
else:
    print "ERROR"
    exit
if has_bins:
    bin_scale=float(raw_input("Scale to remove binaries (pc)?"))
else bin_scale=999
    
cols=raw_input("Columns for x and y positions:").split(',')
for i in range(len(cols)):
    cols[i]=int(cols[i])

skip=int(raw_input("Header rows to skip:"))
"""
"""
fname="./Data/Lupus3_pc.txt"
has_bins=True
cols=[0,1]
skip=0
"""

#fname="../Collaborations/160816_Koepferl/SMA_ALMA/set4.txt"
#has_bins=True
#bin_scale=0.03
#cols=[0,1]
#skip=0

"""fname="../Collaborations/160816_Koepferl/SMA_ALMA/Extra/hyp_181_sinks.txt"
has_bins=True
bin_scale=0.03
cols=[0,1]
skip=0"""

"""fname="./Data/HiGAL/all_tiles_xyz_pc.txt"
has_bins=True
bin_scale=0.03
cols=[0,1,2]fname="./Data/HiGAL/all_tiles_xyz_pc.txt"
has_bins=True
bin_scale=0.03
cols=[0,1,2]
skip=0
dims=3
skip=0
dims=3"""

#fname="./Data/HiGAL/all_tiles_xyz_pc.txt"
has_bins=True
bin_scale=0.03
cols=[0,1]
skip=0
dims=2
n_lim=10 #minimum number of stars for statistical significance

for fname in gb.glob('../161007CloudsToPoints/Experimenting/*.dat'):
    print fname
    if dims==2:
        x,y=np.loadtxt(fname,usecols=cols,skiprows=skip,unpack=True)
    elif dims==3:
        x,y,z=np.loadtxt(fname,usecols=cols,skiprows=skip,unpack=True)
    n=len(x)
    
    if dims==2:
    #find area
        max_dist=0
        for i in range(n):
            for j in range(i+1,n):
                d=np.sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2)
                if d>max_dist:
                    max_dist=d
        print max_dist
        area=np.pi*(max_dist*0.5)**2
        
        expected_bin=10**(1.86*np.log10(n/area)-0.56)
        print expected_bin
    
    
    fig,(ax1,ax2)=plt.subplots(2,1,sharex=True,sharey=True)
    
    ax1.plot(x,y,'bo',mfc='None',ls='None')
    ax1.set_title(fname)
    
    if has_bins:
        if dims==2:
            x,y,_,bin_stars=cf.remove_binaries(list(x),list(y),scale=bin_scale,expected=expected_bin)
            if len(x)>n_lim:
                np.savetxt('./Data/Cloud_play/'+fname+"_nobins.txt",np.vstack((x,y)).T,
                           header=fname+"%i systems merged, %i remain"%((len(bin_stars)/2.),len(x)))    
        elif dims==3:
            #expecting no chance bineries, remove all real binaries
            x,y,z,bin_stars=cf.remove_binaries(list(x),list(y),z=list(z),scale=bin_scale,expected=0)
            if len(x)>n_lim:
                np.savetxt('./Data/Cloud_play'+fname+"_nobins.txt",np.vstack((x,y,z)).T,
                           header=fname+"%i systems merged, %i remain"%((len(bin_stars)/2.),len(x)))
        ax2.plot(x,y,'bo',mfc='None',ls='None')   
        ax2.set_title("%i systems merged, %i remain"%((len(bin_stars)/2.),len(x)))
        plt.savefig('./Data/Cloud_play/'+fname+"_10_nobins.png")
    plt.close()
