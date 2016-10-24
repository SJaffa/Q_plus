import numpy as np
from scipy.sparse.csgraph import minimum_spanning_tree
import re
import multiprocessing as mp
#import glob
import matplotlib.pyplot as plt
from scipy.spatial import KDTree

Rmax = 10

def stars(d, c, gmax):
    L0 = 2.
    S = 2.
    r = 0
    R = Rmax
    #d = np.log2(dp)
    #c = 2.+(1./cp)

    while r < R:
        #print r,
        np.random.seed(r)
        all_limits = makelims(d, np.ceil(gmax))
        star_poss = []
        Nbrem = 0
        for lims in all_limits:
            for l in lims:
                Lg = l[1]-l[0]
                g = (np.log(L0)-np.log(Lg))/np.log(S)
                Nb = S**((3-c)*(np.ceil(gmax)-g))
                if g == np.ceil(gmax) and gmax%1!=0:
                    Nb = Nb*(gmax%1)
                Nbrange = int(Nb)
                Nbrem += Nb%1.
                if Nbrem > 1.:
                    Nbrem -= 1.
                    Nbrange += 1

                for i in np.arange(Nbrange):
                    x = np.random.uniform(l[0],l[1])
                    y = np.random.uniform(l[2],l[3])
                    z = np.random.uniform(l[4],l[5])
                    star_poss.append([x, y, z])

        #print len(star_poss)
        star_poss = cull2sphere(star_poss)
        star_poss = rotate_view(star_poss)
        descname = "D%3.2fC%3.2fG%3.2f"%(d,c,gmax)
        #check for zero mass clusters
        if len(star_poss) == 0:
            r += 1
            R += 1
        else:
            #SAVE POSITIONS OF STARS TO FILE
            starname = "PStars/"+descname+"R%i.star"%r
            np.savetxt(starname, star_poss, fmt="%8.5f",
                delimiter='\t', header=str(len(star_poss))+"\n")
            #save filename to filenamelist
            filesname = "PStars/"+descname+".starnames"
            with open(filesname, "a") as namefile:
                namefile.write(starname+"\n")
            r += 1
    return filesname

def conditions3D():
    xmin = -1.
    xmax = 1.
    ymin = -1.
    ymax = 1.
    zmin = -1
    zmax = 1
    limits = [[xmin, xmax, ymin, ymax, zmin, zmax]]
    return limits
    
def conditions2D():
    xmin = -1.
    xmax = 1.
    ymin = -1.
    ymax = 1.
    limits = [[xmin, xmax, ymin, ymax]]
    return limits

def cube3D(Xmin, Xmax, Ymin, Ymax, Zmin, Zmax, rem=0, dim=3.0, nDiv=2):
    prob = pow(nDiv, dim-3)
    Xmid = np.mean((Xmin, Xmax))
    Ymid = np.mean((Ymin, Ymax))
    Zmid = np.mean((Zmin, Zmax))
    lims = []
    posslims = np.array([[Xmin, Xmid, Ymin, Ymid, Zmin, Zmid],
        [Xmin, Xmid, Ymid, Ymax, Zmin, Zmid],
        [Xmid, Xmax, Ymin, Ymid, Zmin, Zmid],
        [Xmid, Xmax, Ymid, Ymax, Zmin, Zmid],
        [Xmin, Xmid, Ymin, Ymid, Zmid, Zmax],
        [Xmin, Xmid, Ymid, Ymax, Zmid, Zmax],
        [Xmid, Xmax, Ymin, Ymid, Zmid, Zmax],
        [Xmid, Xmax, Ymid, Ymax, Zmid, Zmax]])
    np.random.shuffle(posslims)
    prange = int(prob*8)
    rem += (prob*8.)%1.
    if rem > 1.:
        rem -= 1.
        prange += 1
    for j in np.arange(prange):
        lims.append(posslims[j])
    return lims, rem
    
def cube2D(Xmin, Xmax, Ymin, Ymax, rem=0, dim=2.0, nDiv=2):
    prob = pow(nDiv, dim-2)
    Xmid = np.mean((Xmin, Xmax))
    Ymid = np.mean((Ymin, Ymax))
    lims = []
    posslims = np.array([[Xmin, Xmid, Ymin, Ymid],
        [Xmin, Xmid, Ymid, Ymax],
        [Xmid, Xmax, Ymin, Ymid],
        [Xmid, Xmax, Ymid, Ymax]])
    np.random.shuffle(posslims)
    prange = int(prob*4)
    rem += (prob*4.)%1.
    if rem > 1.:
        rem -= 1.
        prange += 1
    for j in np.arange(prange):
        lims.append(posslims[j])
    return lims, rem

def makelims(F,G):
    generation=0
    current_limits=conditions3D()
    rem=np.random.rand()
    all_lims=[current_limits]
    while generation<G:
        l_loop=range(len(current_limits))

        new_limits=[]
        for i in l_loop:
            new_lims,rem=cube3D(*current_limits[i],rem=rem,dim=F)
            new_limits.extend(new_lims)
        current_limits=new_limits
        all_lims.append(current_limits)
        generation+=1
    return all_lims
    
def makelims_2d(F,G):
    generation=0
    current_limits=conditions2D()
    rem=np.random.rand()
    all_lims=[current_limits]
    while generation<G:
        l_loop=range(len(current_limits))

        new_limits=[]
        for i in l_loop:
            new_lims,rem=cube2D(*current_limits[i],rem=rem,dim=F)
            new_limits.extend(new_lims)
        current_limits=new_limits
        all_lims.append(current_limits)
        generation+=1
    return all_lims

def rotate_view(stars):
    relev=np.random.uniform(-90,90)
    razim=np.random.uniform(0,360)
    phi=np.deg2rad(relev)
    theta=np.deg2rad(razim)
    T1=np.array([[-np.sin(phi), 0, -np.cos(phi)],
                   [0, 1, 0],
                   [np.cos(phi), 0, -np.sin(phi)]])
    T2=np.array([[np.cos(theta), np.sin(theta), 0],
                   [-np.sin(theta), np.cos(theta), 0],
                   [0, 0, 1]])
    T=np.dot(T1,T2)
    poss=[]
    for (X,Y,Z) in stars:
        p=np.array([X,Y,Z])
        p_rot=np.dot(T,p)
        poss.append(p_rot)
    return poss

def cull2sphere(poss):
    star_poss_2=[]
    for s in poss:
        if s[0]**2+s[1]**2+s[2]**2<1.:
            star_poss_2.append(s)
    return star_poss_2

def writelist(out,l):
    out.writelines(["%s\t" % item  for item in l])
    out.write("\n")

#def A_measure(edges):
#    edges.sort()
#    cum_edges=np.cumsum(edges)
#    return ((len(cum_edges)-1)*cum_edges[-1] 
#        - 2*np.sum(cum_edges[:-1]))/(len(cum_edges)*cum_edges[-1])
        
def A_measure(edges,smax):
    eta=1
    edges.sort()
    m=np.array(edges)/smax
    n=len(m)
    A=m[0]*2**eta + m[n-1]*(n-1**eta - (n-2)**eta)
    for i in range(1,n-2):
        A+=m[i]*((i+1)**eta - (i-1)**eta)
    return A/(2.*n**eta)

def justMST(x,y,z=[8]):
    if len(z)==1: #2D data
        grid=np.zeros((len(x),len(y)))
        all_edges=[]
        for i in range(len(x)):
          for j in range(len(x)):
            if i!=j:
              if grid[j][i]==0:
                grid[i][j]=np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)
                all_edges.append(grid[i][j])   
        Tcsr = minimum_spanning_tree(grid,overwrite=True)
        mst=Tcsr.toarray().astype(float)
    
    elif len(z)==len(x):
        length=len(x)
        grid=np.zeros((length,length))
        all_edges=[]
        for i in range(length):
          for j in range(length):
            if i!=j:
              if grid[j][i]==0:
                grid[i][j]=np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2+(z[i]-z[j])**2)
                all_edges.append(grid[i][j])   
        #make mst
        Tcsr = minimum_spanning_tree(grid,overwrite=True)
        mst=Tcsr.toarray().astype(float)
    else:
        print "Error: x,y, and z must have the same length"
        quit()
    return mst, all_edges
    
    
def makeMST(datafile,dim=2,ucols=[0,1]):
    #d,c,g,r=getname(datafile)
    d,c,g,r=0,0,0,0
    #read in data
    with open(datafile, 'r') as infile:
        x,y=np.loadtxt(infile,skiprows=0,unpack=True,usecols=ucols)
    #basis statistics
    x,y=normstar(x,y)
    n_entries=len(x)
    R_cluster=np.amax(np.sqrt((x-np.mean(x))**2+(y-np.mean(y))**2))
    A_cluster=np.pi*pow(R_cluster,2)
    #make complete graph
    if dim==2:
        mst,all_edges=justMST(x,y)
    elif dim==3:
        mst,all_edges=justMST(x,y,z)
    #make resultslist
    reslist=[]
    reslist.append(d)#D
    reslist.append(c)#C
    reslist.append(g)#G
    reslist.append(r)#R
    reslist.append(np.log10(n_entries))# log(n)
    all_edges.sort()
    reslist.append(np.log(all_edges[-1]/all_edges[4]))#Rmeasure
    reslist.append(A_measure(mst[mst!=0],all_edges[-1]))#Ameasure
    reslist.append((np.mean(mst[mst!=0])*(n_entries-1.))/pow((n_entries*A_cluster),0.5))#mbar
    reslist.append(np.mean(all_edges)/R_cluster)#sbar
    reslist.append(np.mean(mst[mst!=0]))#mean mst
    reslist.append(np.std(mst[mst!=0]))#std mst
    reslist.append(np.mean(all_edges))#mean all
    reslist.append(np.std(all_edges))#std all
    reslist.append(reslist[7]/reslist[8]) #Q parameter
    return reslist,x,y

def mst(d, c, g):
    fname = "PStars/D%3.2fC%3.2fG%3.2f.starnames"%(d, c, g)
    outlist = []
    pool = mp.Pool()
    namefile = open(fname, "r")
    names = namefile.read().split()
    namefile.close()
    for datafile in names:
        #outlist.append(makeMST(datafile))
        outlist.append(pool.apply_async( makeMST, args=(datafile, )))
    pool.close()
    pool.join()
    outlist = [output.get() for output in outlist]
    outfile = "PStars/D%3.2fC%3.2fG%3.2f.res"%(d, c, g)
    with open(outfile, 'a') as fileout:
        for i in outlist:
            writelist(fileout, i) 
    return outfile

def n_stars(d, c, g):
    #d = np.log2(dp)
    #c = 2.+(1./cp)
    alpha = pow(2., (c+d-3.))
    n0 = pow(2., (3.-c)*g)
    if (alpha-1. < 0.01):
        return n0*(g+1.)
    else:   
        return n0*(pow(alpha, (g+1.))-1.)/(alpha-1.)
    #return Ntot

def analyse(d, c, g):
    resfile = "PStars/D%3.2fC%3.2fG%3.2f.res"%(d, c, g)
    alldata = np.loadtxt(resfile, ndmin=2)

    names = ["Fractal Dimension", "Density Scaling Exponent",
             "Maximum Generation", "Random seed", "Log10(N_stars)",
             "r measure", "A measure", "Normalised mean edge length of MST",
             "Normalised correlation length", "Mean of MST edges",
             "Standard deviation of MST edges", "Mean of all edges",
             "Standard deviation of all edges"]

    statfile = "PStars/D%3.2fC%3.2fG%3.2f.stat"%(d, c, g)
    outfile = open(statfile, 'w')
    outfile.write(statfile+"\n")
    for i in range(len(names)):
        m = np.mean(alldata[:, i])
        s = np.std(alldata[:, i])
        outfile.write("#%s\t%f\t%f\n"%(names[i], m, s))
    outfile.close()

def readstats(fname):
    #read stats from file
    with open(fname,'r') as sfile:
        sfile.readline()
        slist = []
        p=re.compile('\D+')
        for i in range(13):
            l=sfile.readline()
            #name=l[:-18]
            parts=p.split(l)
            if i==4:
                m=int(parts[2])+float(parts[3])*10**-(len(parts[3]))
                s=int(parts[4])+float(parts[5])*10**-(len(parts[5]))
            else:
                m=int(parts[1])+float(parts[2])*10**-(len(parts[2]))
                s=int(parts[3])+float(parts[4])*10**-(len(parts[4]))
            slist.append([m,s,i])
        
    return slist
    
def getname(fname):
    #Find parameters in filename
    p=re.compile('\D+')
    parts=p.split(fname)
    #print fname
    d=float(parts[1])+(float(parts[2])*10**-(len(parts[2])))
    c=int(parts[3])#+(float(parts[4])*10**-(len(parts[4])))
    g=int(parts[4])#+(float(parts[6])*10**-(len(parts[6])))
    #print len(parts)
    if len(parts)>8:
        r=int(parts[7])
    else:
        r=9999
    return d,c,g,r

def normgaus(val, mu, sig):
    #find prob of val being drawn from gaussian mu, sig
    if sig<0.001:
        sig=0.001
    normfac=1./(pow(2.*np.pi, 0.5)*sig)
    ex=np.exp(-((val-mu)**2)/(2.*(sig**2)))
    return normfac*ex
    
def gaus(val, mu, sig, A):
    #find prob of val being drawn from gaussian mu, sig
    ex=A*np.exp(-((val-mu)**2)/(2.*(sig**2)))
    return ex
    
def covar(X,Y):
    return sum((X - np.mean(X))*(Y - np.mean(Y)))/(len(X)-1.)
    
def potential(x,y):
    n=len(x)
    GPE=0

    for i in range(n):
        for j in range(i+1,n):
            r=np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)
            GPE = GPE + (1./r)
    GPE = 5*GPE/(3*(n**2))
    return GPE

def get_data(params,fullstats,xname,yname):
    meanx=[]
    meany=[]
    stdx=[]
    stdy=[]
    #thing=[]
    for [dp,c,g] in params:
        d=np.log2(dp)
        print d,c,g  ,  
        data=fullstats[(abs(fullstats['D']-d)<0.1) & (fullstats['C']==c) & (fullstats['G']==g)]
        #data=fullstats[(abs(fullstats['D']-d)<0.1) & (fullstats['G']==g)]
        #data=fullstats[(abs(fullstats['D']-d)<0.1)]
        print len(data)
        if len(data)>0:
            #thing.append(c)
            if (xname=='D' or xname =='C' or xname == 'G'):
                meanx.append(float(np.mean(data[xname])))
                stdx.append(0)
            else:
                meanx.append(float(np.mean(data[xname]['mu'])))
                stdx.append(float(np.mean(data[xname]['sig'])))
            if (yname=='D' or yname =='C' or yname == 'G'):
                meany.append(float(np.mean(data[yname])))
                stdy.append(0)
            else:
                meany.append(float(np.mean(data[yname]['mu'])))
                stdy.append(float(np.mean(data[yname]['sig'])))
        
    return [meanx,stdx],[meany,stdy]#,thing
    
def fit(val,name,stats):
    #print val,name,stats
    return normgaus(val, stats[name]['mu'], stats[name]['sig'])

def pooledmeanvar(x,y):
    '''calculates mean and variance of hte combination of two data sets, x and y'''
    xmean = x[0]
    xvar  = x[1]
    ymean = y[0]
    yvar  = y[1]
    xn    = x[2] #number of points
    yn    = y[2]
    
    xyn    = xn + yn
    xymean = (xn * xmean + yn * ymean) / xyn
    xyvar  = np.sqrt(((xn * xvar**2 + yn * yvar**2)/xyn) + ((xn * yn)/(xyn**2))*((xmean - ymean)**2))
    
    return xymean, xyvar, xyn

def normstar(x,y,z=np.array([0])):
    ''' reads in 2d or 3d star cluster and normalises
    positions to radius of 1 centred on mean position'''
    
    xmean=np.mean(x)
    ymean=np.mean(y)
    
    x=x-xmean #centre on mean position
    y=y-ymean
    posmax=max(max(abs(x)),max(abs(y))) #find maximum distance from mean
    
    if len(z)>1:
        #3D cluster
        zmean=np.mean(z)
        z=z-zmean
        posmax=max(posmax,max(abs(z)))
    
    #Scale all positions to maximum distance from centre.
    invmax=1./posmax
    
    x=x*invmax
    y=y*invmax
    z=z*invmax
    
    if len(z)>1:
        return x,y,z
    else:
        return x,y

def get_dicts():
    dtocol={'1.00': 'k',
        '1.58': 'r',
        '2.00': 'g',
        '2.32': 'b',
        '2.58': 'y'}
        
    dtomks={'1.00': 's',
            '1.58': '^',
            '2.00': 'o',
            '2.32': 'v',
            '2.58': 'h'}
        
    gtocol={'3.00': 'k',
            '4.00': 'r',
            '5.00': 'g',
            '6.00': 'b',
            '7.00': 'y',
            '8.00': 'c'}
            
    gtomks={3: 's',
            4: '^',
            5: 'o',
            6: 'h',
            7: 'v',
            8: 'p'}
            
    ctocol={'0.00': 'k',
            '1.00': 'y',
            '2.00': 'g',
            '3.00': 'b',
            '22.00': 'r'}
            
    ctomks={0: 'h',
            1: '^',
            2: 'o',
            3: 's',
            22: 'v'}
    return dtocol,dtomks,ctocol,ctomks,gtocol,gtomks
    
def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n / arrays[0].size
    out[:,0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m,1:])
        for j in xrange(1, arrays[0].size):
            out[j*m:(j+1)*m,1:] = out[0:m,1:]
    return out
    
def select_data(full,d=99,c=99,g=99,dstring='D',cstring='C',gstring='G',sep=0.1):
    if d!=99:
        full=full[(abs(full[dstring]-d)<sep)]
    if c!=99:
        full=full[(abs(full[cstring]-c)<sep)]
    if g!=99:
        full=full[(abs(full[gstring]-g)<sep)]
    return full
    
def read_all_data(save_array=False):
    #Define parameter space
    D,C,G,_=get_parameter_space()
    names=['D','C','G','R','logN','logR','mbar','sbar','muMST','stdMST']
    #structured array
    stats=np.zeros((len(D)*len(C)*len(G)*100),
                   dtype=[('D','f8'),
                          ('C','i4'),
                          ('G','i4'),
                          ('R','i4'),
                          ('logN','f4'),
                          ('logR','f4'),
                          ('mbar','f4'),
                          ('sbar','f4'),
                          ('muMST','f4'),
                          ('stdMST','f4'),
                          ('flag','S50')])
    paramsets=cartesian((D,C,G))
    for n in range(len(paramsets)):
        d,c,g=paramsets[n]
        fname="./AllStars/D%3.2fC%iG%i.res"%(d,c,g)
        try:
            tmp=np.loadtxt(fname,usecols=[0,1,2,3,4,5,7,8,9,10])#,unpack=True)
            for t in range(len(names)):
                stats[n*100:(n*100)+100][names[t]]=tmp[:,t]
            stats[n*100:(n*100)+100]['flag']='True'
        except IOError:
            stats[n*100:(n*100)+100]['D']=d
            stats[n*100:(n*100)+100]['C']=c
            stats[n*100:(n*100)+100]['G']=g
            if n_stars(d,c,g)<20:
                stats[n*100:(n*100)+100]['flag']='Small N'
            elif n_stars(d,c,g)>8000:
                stats[n*100:(n*100)+100]['flag']='Large N'

    #log all R
    stats[:]['logR']=np.log(stats[:]['logR'])
    if save_array:
        np.save('fullres',stats)#
    return stats

def read_all_means():
    D,C,G,F=get_parameter_space()
    names=['D','C','G','R','logN','logR','mbar','sbar','muMST','stdMST']
    paramsets=cartesian((D,C,G))
    
    stats=read_all_data()
    
    meanstats=np.zeros(len(paramsets),
               dtype=[('D','f8',(2)),
                      ('C','i4',(2)),
                      ('G','i4',(2)),
                      ('R','i4',(2)),
                      ('logN','f4',(2)),
                      ('logR','f4',(2)),
                      ('mbar','f4',(2)),
                      ('sbar','f4',(2)),
                      ('muMST','f4',(2)),
                      ('stdMST','f4',(2)),
                      ('flag','S50')])

    for n in range(len(paramsets)):
        d,c,g=paramsets[n]
        thisset=select_data(stats,d=d,c=c,g=g)
        for p in names[:3]:
            meanstats[n][p]=(np.mean(thisset[p]),np.std(thisset[p]))
        meanstats[n]['flag']=thisset[0]['flag']
        if thisset[0]['flag']=='True':
            for p in names[3:]:
                meanstats[n][p]=(np.mean(thisset[p]),np.std(thisset[p]))
    #np.save('meanres',meanstats)            
    return meanstats
    

    
def get_fiducial_values():
    d=1.58
    c=3
    g=5
    return d,c,g
    
def get_col_names():
    return ['D','C','G','R','logN','logR','mbar','sbar','muMST','stdMST']
    
def calculate_parameter_means():
    D,C,G,_=get_parameter_space()

    full=read_all_data()
    data=full[full['flag']=='True']
    names=get_col_names()
    
    means=np.zeros(1,
                   dtype=[('logN','f8'),
                          ('logR','f8'),
                          ('mbar','f8'),
                          ('sbar','f8'),
                          ('muMST','f8'),
                          ('stdMST','f8')])

    stds=np.zeros_like(means)

    for n in names[4:]:
        print n
        means[n]=np.mean(data[n])
        stds[n]=np.std(data[n])
    #np.savez('fullmeans',means=means,stds=stds)
    return means

def mean_adjust_data(save_array=False):   
    alldata=read_all_data()
    
    data=alldata[alldata['flag']=='True']
    names=data.dtype.names
    
    numnames=list(names[4:-1])
    nums=data[numnames]
    numarray=nums.view('f4').reshape(nums.shape + (-1,))
    
    tmp=np.load('fullmeans.npz')
    means=tmp['means']
    stds=tmp['stds']
    meaned=np.zeros_like(data)
    
    for i in range(len(means[0])):
        meaned[numnames[i]]=(numarray[:,i]-means[0][i])/stds[0][i]

    xmeaned=meaned[list(names[4:-1])]
    if save_array:
        np.save('full_adjusted',xmeaned)
    return xmeaned
    
def transform_to_pc(tdata,method='cov'):
                     
    #subtract means
    if method=='cov':
        tmp=np.load('fullmeans.npy')
        meanarr=struct_to_array(tmp)
        
        numarray=struct_to_array(tdata)
        meaned=np.zeros_like(numarray,dtype=float)
        
        for i in range(len(meanarr)):
            meaned[i,:]=(numarray[i,:]-meanarr[i])
            
        #multiply by eigenvecors
        eig=np.load('eigenfull.npz')
        evects=eig['evects'][:,:2]
        newData=np.dot(evects.T,meaned)
    elif method=='corr':
        tmp=np.load('new_fullmeans.npz')
        meanarr=struct_to_array(tmp['means'])
        stdarr=struct_to_array(tmp['stds'])
        
        numarray=struct_to_array(tdata)
        meaned=np.zeros_like(numarray,dtype=float)
        
        for i in range(len(meanarr)):
            meaned[i,:]=(numarray[i,:]-meanarr[i])/stdarr[i]
            
        #multiply by eigenvecors
        eig=np.load('new_eigen.npz')
        evects=eig['evects'][:,:2]
        newData=np.dot(evects.T,meaned)
    return newData
    
def struct_to_array(struct):
    tmp=[]
    [ tmp.append(struct[col]) for col in struct.dtype.names]
    arr=np.array(tmp,dtype=float)
    return arr
    
def get_parameter_space():
    F = np.array([2.0, 3.0, 4.0, 5.0]) #ax0
    C = np.array([2.0, 3.0, 22.0]) #ax1
    G = np.array([3.0, 4.0, 5.0, 6.0, 7.0, 8.0]) #ax2
    D=np.log2(F)
    
    return D,C,G,F

def plot_cluster(d,c,g,r):
    data=np.loadtxt("./AllStars/D%3.2fC%iG%iR%i.star"%(d,c,g,r),usecols=[0,1])
    plt.figure()
    plt.plot(data[:,0],data[:,1],'k.')
    return data

def find_binaries(x,y,z=[],v=0.03):
    if len(z)<1:
        data=np.vstack((x,y)).T
    else:
        data=np.vstack((x,y,z)).T
        #print data.T
        KDT=KDTree(data)
        #print "Tree made"
        binaries=list(KDT.query_pairs(v))
        #print binaries
        try:
            return binaries#[0]
        except IndexError:
            return []

def remove_binaries(x,y,z=[],scale=0.03,expected=0):
    if len(z)<1: #2D
        bins=find_binaries(x,y,v=scale)
        deleted_stars=[]
        while len(bins)>expected:
            #print bins
            i,j = bins[0]
            xa=x.pop(i)
            xb=x.pop(j-1)
            ya=y.pop(i)
            yb=y.pop(j-1)
            xn=(xa+xb)*0.5
            yn=(ya+yb)*0.5
            x.append(xn)
            y.append(yn)
            deleted_stars.extend([[xa,ya],[xb,yb]])
            bins=find_binaries(x,y,v=scale)
    else: #3D
        bins=find_binaries(x,y,z=z,v=scale)
        deleted_stars=[]
        while len(bins)>expected:
            #print bins
            i,j = bins[0] #index of binary pair
            xa=x.pop(i) #remove both from x,y,z
            xb=x.pop(j-1)
            ya=y.pop(i)
            yb=y.pop(j-1)
            za=z.pop(i)
            zb=z.pop(j-1)
            xn=(xa+xb)*0.5 #position of new "system"
            yn=(ya+yb)*0.5
            zn=(za+zb)*0.5
            x.append(xn) #add new system
            y.append(yn)
            z.append(zn)
            deleted_stars.extend([[xa,ya,za],[xb,yb,zb]])
            bins=find_binaries(x,y,z=z,v=scale)
    return x,y,z,deleted_stars