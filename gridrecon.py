# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 15:40:17 2016

@author: c0906755
"""
import numpy as np
import cubefuncs as cf
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from matplotlib import use
#use('pdf')
import glob as gb
import sys

plt.rcParams.update({'font.size': 14, 'font.family':'serif','text.usetex':True})

def pltlabels(ax):
    #ax.set_xlabel("PC 1")
    ax.set_ylabel("PC 2")
    #ax.legend()
    
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
    #plt.figure()
    for i in range(len(arr)):
        cs.append(1./arr[i,0])
        m=arr[i,1]
        s=arr[i,2]
        ps.append(cf.normgaus(A_val,m,s))
        #plt.errorbar(m,1./arr[i,0],xerr=s)
    #plt.vlines(A_val,0,1)
    #print cs
    #print ps
    #plt.figure()
    #plt.plot(cs,ps/sum(ps),'ko')
    #plt.xlim((-0.1,0.6))
    popt,pcov=curve_fit(cf.normgaus,cs,ps/sum(ps),[.3,0.2],bounds=(0,2))
    #print popt,
    #fit=cf.normgaus(np.linspace(0,2,50),*popt)
    #print fit
   
    #plt.plot(np.linspace(0,2,50),fit,'k:')
    #plt.title("D%3.2f G%i A=%3.2f"%(testD,testG,A_val))
    #plt.show()
    return popt#cs[np.argmax(ps)],max(ps)/sum(ps)#popt


dtomks={'1.00': 's',#o
        '1.58': '^',#f
        '2.00': 'o',#o
        '2.32': 's',#f
        '2.58': '^'}#o
    
dtomfc={'1.00': 'None',#o
        '1.58': 'b',#f
        '2.00': 'None',#o
        '2.32': 'b',#f
        '2.58': 'None'}#o
        
ctomks={'1.00': 's',#o
        '2.00': '^',#f
        '3.00': 'o',#o
        '22.00': 's'}#f
    
ctomfc={'1.00': 'None',#o
        '2.00': 'b',#f
        '3.00': 'None',#o
        '22.00': 'b'}#f
        
gtomks={'3.00': 's',#o
        '4.00': '^',#f
        '5.00': 'o',#o
        '6.00': 's',#f
        '7.00': '^',#o
        '8.00': 'o'}#f
    
gtomfc={'3.00': 'None',#o
        '4.00': 'b',#f
        '5.00': 'None',#o
        '6.00': 'b',#f
        '7.00': 'None',#o
        '8.00': 'b'}#f


"""meanparams=np.load('pcmeans.npy')
mp=meanparams[meanparams['flag']=='True']


fig,(ax1,ax2,ax3)=plt.subplots(3,1,sharex=True,sharey=True)
fig.subplots_adjust(hspace=0)
fig.set_size_inches(8.27,11.69)

ax3.set_xlabel("PC 1")
for d in D:
    data=mp[abs(mp['D']-d)<0.1]
    ax1.scatter(data['pc1'][:,0],data['pc2'][:,0],
            facecolor=[dtomfc["%3.2f"%t] for t in data['D']],
            marker=dtomks["%3.2f"%d],s=50,
            label=r'$\cal{D}=$%3.2f'%d,
            edgecolor='b')
pltlabels(ax1)

#plt.figure(2)
for c in C:
    data=mp[abs(mp['C']-c)<0.1]
    ax2.scatter(data['pc1'][:,0],data['pc2'][:,0],
            facecolor=[ctomfc["%3.2f"%t] for t in data['C']],
            marker=ctomks["%3.2f"%c],s=50,
            label=r'$\cal{C}=$%3.2f'%c,
            edgecolor='b')
pltlabels(ax2)
 

for g in G:
    data=mp[abs(mp['G']-g)<0.1]
    ax3.scatter(data['pc1'][:,0],data['pc2'][:,0],
            facecolor=[gtomfc["%3.2f"%t] for t in data['G']],
            marker=gtomks["%3.2f"%g],s=50,
            label=r'$\cal{G}=$%3.2f'%g,
            edgecolor='b')
pltlabels(ax3)

ax1.set_yticks(ax1.get_yticks()[:-2])

#plt.savefig('paramspace.pdf')"""

#D,C,G,_=cf.get_parameter_space()
#R=range(10)

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

testfiles=['../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_100.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_500.dat',
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_1000.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_2000.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_4000.dat', 
'../161007CloudsToPoints/Experimenting/D158C22G6R10_npix_20_std_5_6000.dat', ]

#testfiles=["../Collaborations/160816_Koepferl/SMA_ALMA/set1_nobins.txt",
#           "../Collaborations/160816_Koepferl/SMA_ALMA/set2_nobins.txt",
#           "../Collaborations/160816_Koepferl/SMA_ALMA/set3_nobins.txt",
#           "../Collaborations/160816_Koepferl/SMA_ALMA/set4_nobins.txt",
#           '../Collaborations/160816_Koepferl/SMA_ALMA/Extra/hyp_181_sinks_nobins.txt']

#testfiles=gb.glob('./Data/Clouds/OrionA_*.dat')       
#testfiles=gb.glob('./Data/Clouds/Taurus160_*.dat')
#testfiles=['../Collaborations/160816_Koepferl/SMA_ALMA/Extra/hyp_181_sinks_nobins.txt']
#for (d,c,g) in [(1.00,22,6)]:#cf.cartesian((D,C,G)):
#    print d,g,c
#    testfiles=["AllStars/D%3.2fC%iG%iR%i.star"%(d,c,g,11)]# for r in range(10)]
#    outf="./in-out-grid/cor_D%3.2fC%iG%i_2.txt"%(d,c,g)
#fout=open(outf,'w')
fout=sys.stdout
#testfiles=gb.glob('../Collaborations/160818_Smilgys/cluster*.dat')
#testfiles=gb.glob('../Collaborations/160816_Koepferl/Q_analysis/Q*.txt')
names=[]
#ax2=plt.gca()
#cols=['k','g','b','r','y']
ii=0
for datafile in testfiles:
    for ucols in [[0,1]]:#,[0,2],[1,2]]:
        #fout.write("\n"+datafile+"\n")
        print "['%s',"%datafile,
        #print ucols
        try:
            res,x,y=cf.makeMST(datafile,ucols=ucols)
            if len(x)<10:
                        print "Only %i stars"%len(x)
                        continue
            #plt.figure()
            #plt.plot(x,y,mfc='None',mec='b',marker='o',ls='None')
            #plt.xlabel('sx_source')
            #plt.ylabel('sy_source')
            #plt.axis('equal')
            #plt.xlim((28,37))
            #plt.ylim((27,34))
            #plt.title(datafile[0])
            #print res
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
                                 
            pcs=cf.transform_to_pc(data)
            
            A_measure=res[6]
            #plt.figure(1)
            #l1,=ax1.plot(pcs[0,0],pcs[1,0],'*',markersize=20,label=datafile)#,mfc='None')
            #plt.figure(2)
            l1,=ax2.plot(pcs[0,0],pcs[1,0],'*',markersize=20,label=datafile)#,mfc=datafile[1])
            #plt.figure(3)
            #ax3.plot(pcs[0,0],pcs[1,0],'*',markersize=20)#,mfc='None')
            
            grid_info=np.load('pcgrid.npy')
            
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
            #fout.write(str(est.N))
            if est.N>2: 
                fout.write(est.printD())#+"\n")
                #print est.printC()
                fout.write(est.printG())#+"\n")
                g_found=est.G.mean
                d_found=est.D.mean
            else:
                #fout.write("Zero cluster square \n")
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
                #fout.write(str(thisN)+"\n")
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
                #print thisN
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
#fout.close()
#leg=plt.legend(llist,['Taurus','Chamaeleon I','Lupus 3','IC 348'])
#leg=plt.legend(llist,['v1.2','v1.4','v1.6','v1.8','v2.0'])
#leg=plt.legend(llist,['Set 1','Set 2','Set 3','Set 4', 'hyp_181_sinks'])
#leg=plt.legend(llist,[r'$l$ 217-224',r'$l=218 \pm 1$',r'$l=221 \pm 1$',r'$l=224 \pm 1$'])
#leg=plt.legend(llist,['Unbound','Prestellar','Protostellar','All'])
#leg=plt.legend(llist,names,loc='upper left')
#leg=plt.legend(llist,['With binaries','Without binaries'],loc='upper right')
#leg=plt.legend(llist,['142','132'])
leg=plt.legend(llist,['100','500','1000','2000','4000','6000','original (338)'])
ax2.add_artist(leg)
#plt.legend()

#plt.savefig("../Collaborations/160816_Koepferl/SMA_ALMA/G_grid.pdf")
