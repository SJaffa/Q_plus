# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 15:40:17 2016

@author: c0906755
"""
import numpy as np
import qfuncs as cf
import matplotlib.pyplot as plt
#from scipy.optimize import curve_fit
import glob as gb
import sys
#import annotateplot as ap
#import os

import warnings
warnings.filterwarnings("ignore")

plt.rcParams.update({'font.size': 10, 'font.family':'serif','text.usetex':False})

execfile('gridsearch.py')

D=np.array([ 1.       ,  1.5849625,  2.       ])
C=np.array([  22.,   3.,  2.])
G=np.array([ 4.,  5.,  6.])

#%%
llist=[]

datafile=sys.argv[1]



fout=sys.stdout

ucols=sys.argv[2]
ucols=map(int, ucols.strip('[]').split(','))

#fout.write("\n"+datafile+"\n")
print "%s\n"%datafile,
#print ucols
#try:
res,x,y=cf.makeMST(datafile,ucols=ucols)

if len(x)<10:
    print "Only %i stars"%len(x)

data=np.array([(res[4],
               res[5], 
               res[7], 
               res[8], 
               res[9], 
               res[10])],
              dtype=[('logN', '<f4'), 
                     ('logR', '<f4'), 
                     ('mbar', '<f4'), 
                     ('sbar', '<f4'), 
                     ('muMST', '<f4'), 
                     ('stdMST', '<f4')])

Q=data['mbar']/data['sbar']
if Q>0.9:
    print "Q = %3.2f: this cluster may be centrally concentrated"%Q
    #continue
                     
pcs=cf.transform_to_pc(data)

A_measure=res[6]
ax=plt.gca()
l1,=ax.plot(pcs[0,0],pcs[1,0],'*',markersize=20,label=datafile)#,mfc=datafile[1])
#af = ap.AnnoteFinder([pcs[0,0]],[pcs[1,0]], [dataname], ax=ax2)
#fig.canvas.mpl_connect('button_press_event', af)

#grid_info=np.load('pcgrid.npy')

llist.append(l1)

p1=pcs[0,0]
p2=pcs[1,0]

tmp=True
while tmp:

    print "PC1 = %3.2f, PC2 = %3.2f"%(p1,p2)
    if p2<-0.3:
        print "This cluster is on the edge of parameter space \nand may not be fractally substructured: PC2<-0.3"
    
    try:
        thissquare,_=cf.find_square(grid_info,p1,p2)
    except TypeError:
        fout.write("Outside parameter space:")
        fout.write("PC1=%3.2f, PC2=%3.2f\n"%(p1,p2))
    
    est=grid_info[thissquare]
    
    if est.N>2: 
        fout.write(est.printD())
        fout.write(est.printG())
        g_found=est.G.mean
        d_found=est.D.mean
    else:
        neighbs,newN=cf.find_neighbs(grid_info,thissquare)
        if newN==0:
            fout.write("Outside parameter space (empty super-square):")
            fout.write("PC1=%3.2f, PC2=%3.2f\n"%(p1,p2))
            exit
        #merge D
        thisN=grid_info[neighbs[0]].N
        thisDmean=grid_info[neighbs[0]].D.mean
        thisDvar=grid_info[neighbs[0]].D.std
        for i in range(1,len(neighbs)):
            thisstuff=[thisDmean,thisDvar,thisN]
            newstuff=[grid_info[neighbs[i]].D.mean,
                      grid_info[neighbs[i]].D.std,
                      grid_info[neighbs[i]].N]
            thisDmean,thisDvar,thisN=cf.pooledmeanvar(thisstuff,newstuff)
        fout.write("[%3.2f,%3.2f],"%(thisDmean,thisDvar))
        
        #merge G
        thisN=grid_info[neighbs[0]].N
        thisGmean=grid_info[neighbs[0]].G.mean
        thisGvar=grid_info[neighbs[0]].G.std
        for i in range(1,len(neighbs)):
            thisstuff=[thisGmean,thisGvar,thisN]
            newstuff=[grid_info[neighbs[i]].G.mean,
                      grid_info[neighbs[i]].G.std,
                      grid_info[neighbs[i]].N]
            thisGmean,thisGvar,thisN=cf.pooledmeanvar(thisstuff,newstuff)
        fout.write("[%3.2f,%3.2f],"%(thisGmean,thisGvar))
        g_found=thisGmean
        d_found=thisDmean
    #look up C using A measure
    #find nearest D,G
    f_found=2**d_found
    g_found=round(g_found)
    f_found=round(f_found)
    d_found=np.log2(f_found)
    #print A_measure
    c_mean,c_std=cf.find_c(d_found,g_found,A_measure)
    #fout.write("[%3.2f,%3.2f]],\n"%(c_mean,c_std))
    fout.write("[%3.2f,%3.2f]\n"%(1./c_mean,(c_std/(c_mean**2))))
    
    tmp=False

#except MemoryError:
#    continue

plt.legend(loc='upper left')
plt.show()