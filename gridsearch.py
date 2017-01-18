# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 17:46:20 2016

@author: c0906755
"""

import numpy as np
import qfuncs as cf
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14, 'font.family':'serif','text.usetex':True})

plot_std=False

def plot_rect(xmin,xmax,ymin,ymax,col='b',alpha=1):
    rect=plt.Rectangle((xmin, ymin), xmax-xmin, ymax-ymin,facecolor=col,alpha=alpha)
    plt.gca().add_patch(rect)
    
class liststat(object):
    def __init__(self):
        self.list=[]
    def populate(self,l):
        self.list=l
        self.mean=np.mean(l)
        self.std=np.std(l)
    
class grid_square(object):
    def __init__(self,pc1min,pc2min,pc1max,pc2max,npts,i,j):
        self.pc1min=pc1min
        self.pc1max=pc1max
        self.pc2min=pc2min
        self.pc2max=pc2max
        self.N=npts
        self.i=i
        self.j=j
        if self.N>0:
            self.D=liststat()
            self.C=liststat()
            self.G=liststat()
            self.pc1=liststat()
            self.pc2=liststat()
        else:
            self.D=-99
            self.C=-99
            self.G=-99
            self.pc1=-99
            self.pc2=-99
        
    def printD(self):
        return "[%3.2f,%3.2f]"%(self.D.mean,self.D.std)
        
    def printC(self):
        return "[%3.2f,%3.2f]"%(self.C.mean,self.C.std)
        
    def printG(self):
        return "[%3.2f,%3.2f]"%(self.G.mean,self.G.std)
    
    
D,C,G,_=cf.get_parameter_space()
        
fullres=np.load('fullres.npy')

full=fullres[fullres['flag']=='True'] #takes only present clusters

full['C'][full['C']>5]=4 #relabels c=infty

pcs=cf.transform_to_pc(full[['logN','logR','mbar','sbar','muMST','stdMST']])

alldat=np.zeros(full.shape,dtype=[('D', '<f8'),
                                  ('C', '<f8'), #changed C to floats so can inverse
                                  ('G', '<i4'),
                                  ('R', '<i4'), 
                                  ('logN', '<f4'), 
                                  ('logR', '<f4'), 
                                  ('mbar', '<f4'), 
                                  ('sbar', '<f4'), 
                                  ('muMST', '<f4'), 
                                  ('stdMST', '<f4'), 
                                  ('flag', '|S50'), 
                                  ('pc1','f4'), #added 
                                  ('pc2','f4')]) #added
for name in full.dtype.names:
    alldat[name]=full[name]
alldat['pc1']=pcs[0,:]
alldat['pc2']=pcs[1,:]

min1=np.min(pcs[0,:])-0.01
max1=np.max(pcs[0,:])+0.01
min2=np.min(pcs[1,:])-0.01
max2=np.max(pcs[1,:])+0.01

N=50
thing='D'
n_bad=0

n_pts=np.zeros((N-1,N-1))
means=np.ones((N-1,N-1))*-10
stds=np.ones((N-1,N-1))*-10

grid_info=[]

x=np.linspace(min1,max1,N)
y=np.linspace(min2,max2,N)


dx=x[1]-x[0]
dy=y[1]-y[0]

for i,j in cf.cartesian((range(N-1),range(N-1))):
    here=cf.select_data(alldat,d=np.mean([x[i],x[i+1]]),dstring='pc1',sep=dx/2.)
    here=cf.select_data(here,d=np.mean([y[j],y[j+1]]),dstring='pc2',sep=dy/2.)
    means[i,j]=np.mean(here[thing])  
    stds[i,j]=np.std(here[thing])
    n_pts[i,j]=len(here)
    grid_info.append(grid_square(x[i],y[j],x[i+1],y[j+1],len(here),i,j))
    
    if len(here)==0:
        continue

    else:
        grid_info[-1].D.populate(here['D'])
        grid_info[-1].C.populate(here['C'])
        grid_info[-1].G.populate(here['G'])
        grid_info[-1].pc1.populate(here['pc1'])
        grid_info[-1].pc2.populate(here['pc2'])

    
meanmax=np.max(means)
meanmin=np.min(means)
stdmax=np.max(stds)
stdmin=np.min(stds)

if plot_std:
    fig,(ax2,ax3)=plt.subplots(1,2,sharex=True,sharey=True,figsize=(18,6))
else:
    fig=plt.figure()
    fig.set_size_inches(8.5,6)
    ax2=plt.gca()

if thing=='G':
    thing='L'
ax2.set_title(r'Mean $\cal{'+thing+'}$')
if plot_std:
    ax3.set_title(r'Standard deviation of $\cal{'+thing+'}$')
    ax3.set_xlabel("PC 1")
    ax3.set_ylabel("PC 2")
ax2.set_xlabel("PC 1")
ax2.set_ylabel("PC 2")


im1=ax2.imshow(means.T,aspect='auto',interpolation='none',
               origin='lower',extent=(min1,max1,min2,max2),vmin=0)#,vmax=meanmax)
if plot_std:
    im2=ax3.imshow(stds.T,aspect='auto',interpolation='none',
                   origin='lower',extent=(min1,max1,min2,max2),vmin=0)#

fig.colorbar(im1,ax=ax2,pad=0)
if plot_std:
    fig.colorbar(im2,ax=ax3,pad=0)

#plt.savefig("grid_"+thing+".pdf")

#np.save('pcgrid',grid_info,allow_pickle=True)