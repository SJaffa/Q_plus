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
    
    
    WORKS IF I USE MAKEMST BUT NOT JUSTMST! WHY?
    SAVE AND IMPORT GRIDSEARCH
    ADD FINDING C
    RETURN D,L,C AND ERRORS
    
"""

def cluster_radius(x,y):
    '''
    Find radius of cluster
    '''
    n=len(x)
    max_dist=0
    for i in range(n):
        for j in range(i+1,n):
            d=np.sqrt((x[i]-x[j])**2 + (y[i]-y[j])**2)
            if d>max_dist:
                max_dist=d
    return max_dist

def cluster_area(x,y):
    '''
    Calculate area of cluster
    '''
    max_dist=cluster_radius(x,y)
        
    area=np.pi*(max_dist*0.5)**2
    return area

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
    
def measure_A(edges,smax):
    edges.sort()
    m=np.array(edges)/smax
    n=len(m)
    A=m[0]*2 + m[n-1]*((n-1) - (n-2))
    for i in range(1,n-2):
        A+=m[i]*((i+1) - (i-1))
    return A/(2.*n)
    
def makeMST(x,y):

    x,y=qf.normstar(x,y)
    n_entries=len(x)
    R_cluster=np.amax(np.sqrt((x-np.mean(x))**2+(y-np.mean(y))**2))
    A_cluster=np.pi*pow(R_cluster,2)
    #make complete graph
    mst,all_edges=qf.justMST(x,y)
    #make resultslist
    reslist=[]
    reslist.append(np.log10(n_entries))# log(n)
    all_edges.sort()
    reslist.append(np.log(all_edges[-1]/all_edges[4]))#Rmeasure
    reslist.append((np.mean(mst[mst!=0])*(n_entries-1.))/pow((n_entries*A_cluster),0.5))#mbar
    reslist.append(np.mean(all_edges)/R_cluster)#sbar
    reslist.append(np.mean(mst[mst!=0]))#mean mst
    reslist.append(np.std(mst[mst!=0]))#std mst
    reslist.append(measure_A(mst[mst!=0],all_edges[-1]))#Ameasure
    return reslist,x,y

#==============================================================================
import numpy as np
import matplotlib.pyplot as plt
import qfuncs as qf

fname='/home/sjaffa/Dropbox/PhDwork/160101DCGstars/Data/Real/Taurus_nobins.txt'
ucols=[0,1]
skip=0


x,y=np.loadtxt(fname,usecols=ucols, skiprows=skip, unpack=True)

#exp_bins=expected_vis_bin(x,y)

#x,y,_,bin_stars=qf.remove_binaries(list(x),list(y),
#                                   expected=exp_bins)

#xnorm,ynorm=qf.normstar(x,y)
#mst,all_edges=qf.justMST(xnorm,ynorm)

#res=seven_measures(len(x), 
#                   mst,
#                   all_edges,
#                   cluster_area(xnorm,ynorm),
#                   cluster_radius(xnorm,ynorm))

res,x,y=makeMST(x,y)
print res

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

Q=data['mbar']/data['sbar']
if Q>0.7:
    print "Q = %3.2f: this cluster may be centrally concentrated"%Q
    
                     
pcs=qf.transform_to_pc(data)

l1,=ax2.plot(pcs[0,0],pcs[1,0],'*',markersize=20,label=fname)
print fname
print pcs


