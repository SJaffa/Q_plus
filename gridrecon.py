# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 15:40:17 2016

@author: c0906755
"""
import numpy as np
import qfuncs as cf
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import glob as gb
import sys
#import annotateplot as ap

plt.rcParams.update({'font.size': 10, 'font.family':'serif','text.usetex':False})

    
def find_square(grid,p1,p2):
    for i in range(len(grid)):
        gs=grid[i]
        if (gs.pc1min<p1 and gs.pc1max>p1):
            if (gs.pc2min<p2 and gs.pc2max>p2):
                return i,gs.N
    
def find_neighbs(grid,i):
    this=grid[i]
    dx=this.pc1max-this.pc1min
    dy=this.pc2max-this.pc2min
    xmid=this.pc1min+0.5*dx
    ymid=this.pc2min+0.5*dy
    try:
        n1,nn1=find_square(grid,xmid-dx,ymid+dy)
    except TypeError:
        n1=0
        nn1=0
        pass
    try:
        n2,nn2=find_square(grid,xmid,ymid+dy)
    except TypeError:
        n2=0
        nn2=0
        pass
    try:
        n3,nn3=find_square(grid,xmid+dx,ymid+dy)
    except TypeError:
        n3=0
        nn3=0
        pass
    try:
        n4,nn4=find_square(grid,xmid-dx,ymid)
    except TypeError:
        n4=0
        nn4=0
        pass
    try:
        n5,nn5=find_square(grid,xmid+dx,ymid)
    except TypeError:
        n5=0
        nn5=0
        pass
    try:
        n6,nn6=find_square(grid,xmid-dx,ymid-dy)
    except TypeError:
        n6=0
        nn6=0
        pass
    try:
        n7,nn7=find_square(grid,xmid,ymid-dy)
    except TypeError:
        n7=0
        nn7=0
        pass
    try:
        n8,nn8=find_square(grid,xmid+dx,ymid-dy)
    except TypeError:
        n8=0
        nn8=0
        pass
    mask_empty_squares=np.array([nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8])
    n_tot=sum(mask_empty_squares)
    return np.array([n1,n2,n3,n4,n5,n6,n7,n8])[mask_empty_squares>0],n_tot

def find_c(testD,testG,A_val):
    c_file="./AllStars/C_D%3.2fG%i.res"%(testD,testG)
    arr=np.loadtxt(c_file,ndmin=2)
    cs=[]
    ps=[]
    for i in range(len(arr)):
        cs.append(1./arr[i,0])
        m=arr[i,1]
        s=arr[i,2]
        ps.append(cf.normgaus(A_val,m,s))
    expect=sum(np.array(ps)*np.array(cs))/sum(ps)
    var=sum(np.array(ps)*pow((np.array(cs)-expect),2))/sum(ps)
    return expect,var


D=np.array([ 1.       ,  1.5849625,  2.       ])
C=np.array([  22.,   3.,  2.])
G=np.array([ 4.,  5.,  6.])

#testfiles=["AllStars/D%3.2fC%iG%iR%i.star"%(d,c,g,r) for (d,c,g,r) in cf.cartesian((D,C,G,R))]

#%%
llist=[]
#testfiles=[["./Data/Real/Chamaeleon.dat",'None'],
#           ["./Data/Real/ChamaeleonI_nobins.txt",'k']]

#testfiles=[["./Data/HiGAL/all_tiles_latlon.txt",'y'],
#           ["./Data/HiGAL/long_218.dat",'r'],
#           ["./Data/HiGAL/long_221.dat",'k'],
#           ["./Data/HiGAL/long_224.dat",'b']]

#testfiles=[["./Data/HiGAL/unbound_latlondis.dat",'y'],
#           ["./Data/HiGAL/prestellar_latlondis.dat",'r'],
#            ["./Data/HiGAL/protostellar_latlondis.dat",'b'],
#            ["./Data/HiGAL/all_tiles_latlon.txt",'k'],]

#testfiles=[["./Data/Real/Taurus_nobins.txt",'k'],
#           ["./Data/Real/ChamaeleonI_nobins.txt",'g'],
#           ["./Data/Real/Lupus3_nobins.txt",'b'],
#           ["./Data/Real/IC348_nobins.txt",'y']]
           
#testfiles=[ "./Data/Scott/ST1.2_nobins.txt",
#           "./Data/Scott/ST1.4_nobins.txt",
#           "./Data/Scott/ST1.6_nobins.txt",
#           "./Data/Scott/ST1.8_nobins.txt",
#           "./Data/Scott/ST2.0_nobins.txt",]

#testfiles=gb.glob('./Data/Koepferl/Nobins/Qsim*_nobins.txt')
#testfiles=gb.glob('/home/sjaffa/Dropbox/PhDwork/160101DCGstars/Data/Frac_clouds/D232C22G6R0_points.txt')

"""testfiles=['../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_10_100.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_10_500.dat',
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_10_1000.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_10_2000.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_10_4000.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_10_6000.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_original.dat',
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_5_100.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_5_500.dat',
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_5_1000.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_5_2000.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_5_4000.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_30_std_5_6000.dat', ]"""

#testfiles=['../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_100.dat', 
#'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_500.dat',
#'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_1000.dat', 
#'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_2000.dat', 
#'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_4000.dat', 
#'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_6000.dat', ]

#testfiles=["../Collaborations/160816_Koepferl/SMA_ALMA/set1_nobins.txt",
#           "../Collaborations/160816_Koepferl/SMA_ALMA/set2_nobins.txt",
#           "../Collaborations/160816_Koepferl/SMA_ALMA/set3_nobins.txt",
#           "../Collaborations/160816_Koepferl/SMA_ALMA/set4_nobins.txt",
#           '../Collaborations/160816_Koepferl/SMA_ALMA/Extra/hyp_181_sinks_nobins.txt']

#testfiles=gb.glob('./Data/Clouds/OrionA_*.dat')       
#testfiles=gb.glob('./Data/Clouds/Taurus160_*.dat')
#testfiles=['../Collaborations/160816_Koepferl/SMA_ALMA/Extra/hyp_181_sinks_nobins.txt']

#testfiles=gb.glob('../Collaborations/160818_Smilgys/cluster*.dat')
#testfiles=gb.glob('../Collaborations/160816_Koepferl/Q_analysis/Q*.txt')
#testfiles=gb.glob('/home/sjaffa/Dropbox/PhDwork/Collaborations/161127_Boily/ForSJ/r15_Ser*_nobins.txt')

#testfiles=gb.glob('/home/sjaffa/Dropbox/PhDwork/170103QplusClouds/HiGal/slice_*n2000.pp')

testfiles=gb.glob('/home/sjaffa/Dropbox/PhDwork/Collaborations/170109_Urquhart/G*')

fout=sys.stdout
names=[]
ii=0
for datafile in testfiles:
    dataname=datafile[60:]
    ucols=[1,2]
    #fout.write("\n"+datafile+"\n")
    print "['%s',"%dataname,
    #print ucols
    try:
        res,x,y=cf.makeMST(datafile,ucols=ucols)
        if len(x)<10:
                    print "Only %i stars"%len(x)
                    continue
        names.append(datafile[40:-8])
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
        if Q>0.7:
            print "Q = %3.2f: this cluster may be centrally concentrated"%Q
            #continue
                             
        pcs=cf.transform_to_pc(data)
        
        A_measure=res[6]
        l1,=ax2.plot(pcs[0,0],pcs[1,0],'*',markersize=20,label=dataname)#,mfc=datafile[1])
        #af = ap.AnnoteFinder([pcs[0,0]],[pcs[1,0]], [dataname], ax=ax2)
        #fig.canvas.mpl_connect('button_press_event', af)
        
        #grid_info=np.load('pcgrid.npy')
        
        llist.append(l1)
        
        p1=pcs[0,0]
        p2=pcs[1,0]
        
        try:
            thissquare,_=find_square(grid_info,p1,p2)
        except TypeError:
            fout.write("Outside parameter space:")
            fout.write("PC1=%3.2f, PC2=%3.2f\n"%(p1,p2))
            continue
        
        est=grid_info[thissquare]

        if est.N>2: 
            fout.write(est.printD())
            fout.write(est.printG())
            g_found=est.G.mean
            d_found=est.D.mean
        else:
            neighbs,newN=find_neighbs(grid_info,thissquare)
            if newN==0:
                fout.write("Outside parameter space (empty super-square):")
                fout.write("PC1=%3.2f, PC2=%3.2f\n"%(p1,p2))
                continue
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
        c_mean,c_std=find_c(d_found,g_found,A_measure)
        fout.write("[%3.2f,%3.2f]],\n"%(c_mean,c_std))
        #fout.write("c=[%3.2f,%3.2f]\n"%(1./c_mean,(c_std/(c_mean**2))))
        
        ii+=1
        
    except MemoryError:
        continue


#leg=plt.legend(llist,['Taurus','Chamaeleon I','Lupus 3','IC 348'])
#leg=plt.legend(llist,['v1.2','v1.4','v1.6','v1.8','v2.0'])
#leg=plt.legend(llist,['Set 1','Set 2','Set 3','Set 4', 'hyp_181_sinks'])
#leg=plt.legend(llist,[r'$l$ 217-224',r'$l=218 \pm 1$',r'$l=221 \pm 1$',r'$l=224 \pm 1$'])
#leg=plt.legend(llist,['Unbound','Prestellar','Protostellar','All'])
#leg=plt.legend(llist,names,loc='upper left')
#leg=plt.legend(llist,['With binaries','Without binaries'],loc='upper right')
#leg=plt.legend(llist,['142','132'])
#leg=plt.legend(llist,['100','500','1000','2000','4000','6000','original (338)'])
#leg=plt.legend(llist,['Series 2','Series 3'])
#ax2.add_artist(leg)
plt.legend(loc='upper left')

#plt.savefig("../Collaborations/160816_Koepferl/SMA_ALMA/G_grid.pdf")
