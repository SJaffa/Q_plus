#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  5 19:26:37 2017

@author: sjaffa

Test new structure:
    Preprocess (with units)
    Build CG
    Build MST
    Derive stats
    Calculate PCs
    
    PROBLEMS WITH C
    
"""

def cluster_radius(x,y):
    '''
    Find radius of cluster
    '''
    R_cluster=np.amax(np.sqrt((x-np.mean(x))**2+(y-np.mean(y))**2))
    return R_cluster
    
def cluster_area(x,y):
    '''
    Calculate area of cluster
    '''
    R_cluster=cluster_radius(x,y)
    A_cluster=np.pi*pow(R_cluster,2)
    return A_cluster

def expected_vis_bin(x,y):
    '''
    Find how many visual binaries (not gravitationally bound systems) 
    we expect in projection of a cluster with this density.
    '''
    n=len(x)
    area=cluster_area(x,y)
    
    expected_bin=10**(1.86*np.log10(n/area)-0.56)
    return expected_bin
    
def measure_logN(n):
    return np.log10(n)
    
def measure_logR(all_edges):
    all_edges.sort()
    return np.log(all_edges[-1]/all_edges[4])
    
def measure_mbar(n,mst,area):
    return (np.mean(mst[mst!=0])*(n-1.))/pow((n*area),0.5)
    
def measure_sbar(all_edges,rad):
    return np.mean(all_edges)/rad
    
def measure_muMST(mst):
    return np.mean(mst[mst!=0])
    
def measure_stdMST(mst):
    return np.std(mst[mst!=0])
    
    
def seven_measures(n, mst,all_edges,area,radius):
    logN=measure_logN(n)
    logR=measure_logR(all_edges)
    mbar=measure_mbar(n,mst,area)
    sbar=measure_sbar(all_edges, radius)
    muMST=measure_muMST(mst)
    stdMST=measure_stdMST(mst)
    A=measure_A(all_edges)
    return [logN, logR, mbar, sbar, muMST, stdMST, A]
    
def measure_A(edges):
    edges.sort()
    smax=np.amax(edges)
    m=np.array(edges)/smax
    n=len(m)
    A=m[0]*2 + m[n-1]*((n-1) - (n-2))
    for i in range(1,n-2):
        A+=m[i]*((i+1) - (i-1))
    return A/(2.*n)
    
def find_square(grid,p1,p2):
    n=grid.shape[0]
    for i in range(n):
        if grid[i,0]['pc1']>p1:
            break
    for j in range(n):
        if grid[0,j]['pc2']>p2:
            break
    return i-1,j-1
                
def find_c(testD,testG,A_val):
    c_file="./AllStars/C_D%3.2fG%i.res"%(testD,testG)
    arr=np.loadtxt(c_file,ndmin=2)
    cs=[]
    ps=[]
    for i in range(len(arr)):
        cs.append(1./arr[i,0])
        m=arr[i,1]
        s=arr[i,2]
        ps.append(qf.normgaus(A_val,m,s))
    expect=sum(np.array(ps)*np.array(cs))/sum(ps)
    var=sum(np.array(ps)*pow((np.array(cs)-expect),2))/sum(ps)
    return expect,var
    
#==============================================================================
import numpy as np
import matplotlib.pyplot as plt
import qfuncs as qf
import glob as gb

plt.rcParams.update({'font.size': 14, 'font.family':'serif','text.usetex':False})
plot_std=False
thing='L'

namelist=gb.glob('/home/sjaffa/Dropbox/PhDwork/160101DCGstars/Data/Real/*_nobins.txt')
ucols=[0,1]
skip=0

if plot_std:
    fig,(ax2,ax3)=plt.subplots(1,2,sharex=True,sharey=True,figsize=(18,6))
else:
    fig=plt.figure()
    fig.set_size_inches(8.5,6)
    ax2=plt.gca()

ax2.set_title(r'Mean $\cal{'+thing+'}$')
if plot_std:
    ax3.set_title(r'Standard deviation of $\cal{'+thing+'}$')
    ax3.set_xlabel("PC 1")
    ax3.set_ylabel("PC 2")
ax2.set_xlabel("PC 1")
ax2.set_ylabel("PC 2")

grid=np.load('grid_DL')
means=grid['mean'+thing]
stds=grid['std'+thing]

min1=np.min(grid['pc1'])-0.01
max1=np.max(grid['pc1'])+0.01
min2=np.min(grid['pc2'])-0.01
max2=np.max(grid['pc2'])+0.01

im1=ax2.imshow(means.T,aspect='auto',interpolation='none',
               origin='lower',extent=(min1,max1,min2,max2),vmin=0)#,vmax=meanmax)
if plot_std:
    im2=ax3.imshow(stds.T,aspect='auto',interpolation='none',
                   origin='lower',extent=(min1,max1,min2,max2),vmin=0)#

fig.colorbar(im1,ax=ax2,pad=0)
if plot_std:
    fig.colorbar(im2,ax=ax3,pad=0)


for fname in namelist:
    print fname

    x,y=np.loadtxt(fname,usecols=ucols, skiprows=skip, unpack=True)
    
    #exp_bins=expected_vis_bin(x,y)
    
    #x,y,_,bin_stars=qf.remove_binaries(list(x),list(y),
    #                                   expected=exp_bins)
    
    xnorm,ynorm=qf.normstar(x,y)
    mst,all_edges=qf.justMST(xnorm,ynorm)

    all_edges.sort()
    res=seven_measures(len(xnorm), 
                       mst,
                       all_edges,
                       cluster_area(xnorm,ynorm),
                       cluster_radius(xnorm,ynorm))
    
    data=np.array([(res[0],
                   res[1], 
                   res[2], 
                   res[3], 
                   res[4], 
                   res[5])],
                   dtype=[('logN', '<f4'), 
                         ('logR', '<f4'), 
                         ('mbar', '<f4'), 
                         ('sbar', '<f4'), 
                         ('muMST', '<f4'), 
                         ('stdMST', '<f4')])
    #print data['mbar'],data['sbar'],data['mbar']/data['sbar']
    
    Q=res[2]/res[3]
    if Q>0.7:
        print "Q = %3.2f: this cluster may be centrally concentrated"%Q
    if Q>0.9:
        print "Q = %3.2f: this cluster may be centrally concentrated"%Q
        #continue
                         
    p1,p2=qf.transform_to_pc(data)
    
    l1,=ax2.plot(p1,p2,'.',markersize=10,label=fname[63:-11])
    print p1,p2
    plt.legend()
    
    

    #================================================================
    import sys
    fout=sys.stdout

    ii,jj=find_square(grid,p1,p2)
    if grid[ii,jj]['n']==0:
        fout.write("Outside parameter space:")
        fout.write("PC1=%3.2f, PC2=%3.2f\n"%(p1,p2))
        continue
        #cf.pooledmeanvar(thisstuff,newstuff)
    else:
        d_found=grid[ii,jj]['meanD']
        d_std=grid[ii,jj]['stdD']
        l_found=grid[ii,jj]['meanL']
        l_std=grid[ii,jj]['stdL']
    
    fout.write("[%3.2f,%3.2f]],\n"%(d_found,d_std))
    fout.write("[%3.2f,%3.2f]],\n"%(l_found,l_std))
        
    #look up C using A measure
    #find nearest D,G
    f_found=2**d_found
    l_found=round(l_found)
    f_found=round(f_found)
    d_found=np.log2(f_found)
    
    print d_found,l_found
    
    print res[-1]
    c_file="./AllStars/C_D%3.2fG%i.res"%(d_found,l_found)
    arr=np.loadtxt(c_file,ndmin=2)
    cs=[]
    ps=[]
    for i in range(len(arr)):
        cs.append(1./arr[i,0])
        m=arr[i,1]
        s=arr[i,2]
        ps.append(qf.normgaus(res[-1],m,s))
    expect=sum(np.array(ps)*np.array(cs))/sum(ps)
    var=sum(np.array(ps)*pow((np.array(cs)-expect),2))/sum(ps)
    
    print cs
    print ps
    
    print expect,var
    
    c_mean,c_std=find_c(d_found,l_found,res[-1])
    fout.write("[%3.2f,%3.2f]],\n"%(c_mean,c_std))
    fout.write("c=[%3.2f,%3.2f]\n"%(1./c_mean,(c_std/(c_mean**2))))

