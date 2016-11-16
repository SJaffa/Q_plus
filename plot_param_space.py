# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 14:09:35 2016

@author: c0906755
"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

plt.rcParams.update({'font.size': 6, 'font.family':'serif','text.usetex':False})

def plot_3d_error(ax,fx,fy,fz,xerror,yerror,zerror,label=None,mark='o'):
    a=0.3
    ax.plot([fx+xerror, fx-xerror], [fy, fy], [fz, fz], marker="_",c='k',alpha=a)
    ax.plot([fx, fx], [fy+yerror, fy-yerror], [fz, fz], marker="_",c='k',alpha=a)
    ax.plot([fx, fx], [fy, fy], [fz+zerror, fz-zerror], marker="_",c='k',alpha=a)
    ax.plot([fx], [fy], [fz], linestyle="None", marker=mark,label=label)

def plot_3d_line(ax,datalist,label):
    ds=[]
    cs=[]
    gs=[]        
    for (name,d,g,cinv) in datalist:
        ds.append(d[0])
        cs.append(cinv[0])
        gs.append(g[0])
    ax.plot(ds,gs,cs,'-',label=label) 

fig = plt.figure(dpi=100)
ax = fig.add_subplot(111, projection='3d')

if 0:
    real_clusters=[["Taurus",[1.00,0.00],[8.00,0.00],[0.84,0.25]],
                   ["Chamaeleon I",[1.46,0.23],[5.20,0.40],[0.84,0.29]],
                   ["Lupus 3",[1.00,0.00],[7.00,0.34],[0.00,1.02]],
                   ["IC 348",[2.16,0.10],[3.38,0.00],[0.00,1.03]]]
                  
    for (name,d,g,cinv) in real_clusters:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                         label=name,mark='^')

if 0:
    koepferl_data=[["Set 1",[1.03,0.12],[4.95,0.55],[0.00,1.02]],
                   ["Set 2",[1.04,0.14],[4.93,0.65],[0.00,1.07]],
                   ["Set 3",[1.00,0.00],[6.29,0.35],[0.84,1.02]],
                   ["Set 4",[1.00,0.00],[6.38,0.25],[0.00,1.10]],
                   ["hyp_181_sinks",[1.00,0.00],[5.60,0.49],[0.00,0.92]]]
    
    for (name,d,g,cinv) in koepferl_data:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='o')
                      
if 0:
    koepferl_points=[['142',[1.00,0.00],[4.33,0.47],[0.99,0.57]],
                    ['132',[1.19,0.27],[3.33,0.47],[0.83,0.28]]]

    for (name,d,g,cinv) in koepferl_points:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='o')
                      
if 0:
    koepferl_sim=[['52_F',[1.00,0.00],[6.33,0.00],[0.99,0.55]],
                ['52_F',[1.00,0.00],[5.83,0.42],[0.00,1.02]],
                ['52_F',[1.00,0.00],[6.40,0.00],[0.97,0.87]],
                ['42_B',[1.00,0.00],[7.33,0.00],[0.00,0.96]],
                ['42_B',[1.00,0.00],[7.33,0.00],[0.00,1.04]],
                ['42_B',[1.00,0.00],[7.00,0.00],[0.00,1.02]],
                ['162_F',[1.00,0.00],[7.73,0.39],[0.00,1.04]],
                ['162_F',[1.00,0.00],[7.10,0.25],[0.00,0.91]],
                ['162_F',[1.00,0.00],[7.57,0.31],[0.00,1.03]],
                ['82_T',[1.00,0.00],[6.50,0.00],[0.00,0.91]],
                ['82_T',[1.00,0.00],[6.33,0.00],[0.00,0.91]],
                ['82_T',[1.00,0.00],[6.67,0.00],[0.00,0.91]],
                ['122_F',[1.00,0.00],[6.20,0.40],[0.00,0.97]],
                ['122_F',[1.00,0.00],[6.25,0.43],[0.00,0.99]],
                ['122_F',[1.00,0.00],[6.56,0.39],[0.00,0.91]],
                ['152_T',[1.00,0.00],[8.00,0.00],[0.00,1.06]],
                ['52_B',[1.00,0.00],[6.18,0.56],[1.00,0.60]],
                ['52_B',[1.00,0.00],[7.00,0.00],[0.00,1.01]],
                ['52_B',[1.00,0.00],[8.00,0.00],[0.00,1.03]],
                ['142_B',[1.00,0.00],[7.75,0.43],[0.00,1.03]],
                ['142_B',[1.00,0.00],[7.64,0.39],[0.00,1.02]],
                ['142_B',[1.00,0.00],[7.10,0.25],[0.00,0.96]],
                ['112_B',[1.00,0.00],[7.10,0.25],[0.00,0.91]],
                ['112_B',[1.00,0.00],[6.67,0.47],[0.00,0.93]],
                ['112_B',[1.00,0.00],[6.67,0.47],[0.00,0.92]],
                ['132_B',[1.00,0.00],[7.00,0.00],[0.00,0.91]],
                ['132_B',[1.00,0.00],[7.00,0.00],[0.00,0.92]],
                ['132_B',[1.00,0.00],[6.96,0.43],[0.00,0.93]],
                ['92_F',[1.00,0.00],[6.38,0.25],[0.00,1.02]],
                ['92_F',[1.00,0.00],[6.64,0.26],[0.00,0.91]],
                ['92_F',[1.00,0.00],[6.11,0.24],[0.06,1.15]],
                ['102_B',[1.00,0.00],[6.33,0.47],[0.00,1.05]],
                ['102_B',[1.00,0.00],[7.00,0.32],[0.00,0.91]],
                ['102_B',[1.00,0.00],[7.00,0.00],[0.00,0.91]],
                ['152_F',[1.00,0.00],[7.10,0.25],[0.00,0.91]],
                ['152_F',[1.00,0.00],[7.67,0.47],[0.00,1.03]],
                ['152_F',[1.00,0.00],[6.67,0.47],[0.00,0.92]],
                ['152_B',[1.00,0.00],[7.75,0.43],[0.00,1.03]],
                ['152_B',[1.00,0.00],[7.86,0.27],[0.00,1.02]],
                ['152_B',[1.00,0.00],[7.41,0.35],[0.00,1.04]],
                ['82_F',[1.00,0.00],[6.11,0.24],[0.00,1.05]],
                ['82_F',[1.00,0.00],[6.64,0.26],[0.00,0.91]],
                ['30_B',[1.00,0.00],[6.67,0.00],[0.00,0.91]],
                ['30_B',[1.00,0.00],[6.00,0.00],[0.00,1.06]],
                ['30_B',[1.00,0.00],[6.33,0.00],[0.00,1.04]],
                ['92_B',[1.00,0.00],[6.33,0.47],[0.00,1.09]],
                ['92_B',[1.00,0.00],[7.33,0.46],[0.00,0.93]],
                ['92_B',[1.00,0.00],[7.00,0.00],[0.00,0.91]],
                ['180_B',[1.00,0.00],[7.86,0.27],[0.00,1.02]],
                ['180_B',[1.00,0.00],[7.40,0.39],[0.00,1.10]],
                ['180_B',[1.00,0.00],[7.75,0.43],[0.00,1.02]],
                ['132_F',[1.00,0.00],[6.33,0.47],[0.00,0.99]],
                ['132_F',[1.00,0.00],[7.00,0.00],[0.00,0.91]],
                ['132_F',[1.00,0.00],[7.00,0.82],[0.00,0.91]],
                ['122_T',[1.00,0.00],[6.33,0.00],[0.00,0.92]],
                ['122_T',[1.00,0.00],[7.33,0.00],[0.00,0.91]],
                ['122_T',[1.00,0.00],[6.33,0.00],[0.00,0.91]],
                ['32_T',[1.00,0.00],[7.00,0.00],[0.00,0.91]],
                ['32_T',[1.00,0.00],[6.50,0.00],[0.00,0.94]],
                ['32_T',[1.00,0.00],[6.00,0.00],[0.97,0.50]],
                ['102_T',[1.00,0.00],[6.00,0.00],[0.00,0.91]],
                ['102_T',[1.00,0.00],[6.67,0.00],[0.00,0.91]],
                ['102_T',[1.00,0.00],[5.87,0.41],[0.00,0.91]],
                ['180_F',[1.00,0.00],[7.86,0.27],[0.00,1.02]],
                ['180_F',[1.00,0.00],[7.40,0.39],[0.00,1.10]],
                ['180_F',[1.00,0.00],[7.75,0.43],[0.00,1.02]],
                ['102_F',[1.00,0.00],[6.86,0.44],[0.00,0.91]],
                ['102_F',[1.00,0.00],[6.86,0.44],[0.00,0.91]],
                ['102_F',[1.00,0.00],[7.11,0.37],[0.00,0.91]],
                ['172_F',[1.00,0.00],[7.75,0.43],[0.00,1.03]],
                ['172_F',[1.00,0.00],[6.67,0.47],[0.00,0.91]],
                ['172_F',[1.00,0.00],[6.96,0.43],[0.00,0.92]],
                ['132_T',[1.00,0.00],[5.93,0.19],[0.00,0.91]],
                ['132_T',[1.00,0.00],[6.67,0.00],[0.00,0.91]],
                ['92_T',[1.00,0.00],[7.33,0.00],[0.00,0.91]],
                ['92_T',[1.00,0.00],[6.33,0.00],[0.00,0.91]],
                ['172_B',[1.00,0.00],[7.32,0.46],[0.00,0.97]],
                ['172_B',[1.00,0.00],[7.75,0.43],[0.00,1.02]],
                ['172_B',[1.00,0.00],[7.75,0.43],[0.00,1.02]],
                ['31_T',[1.00,0.00],[6.00,0.00],[0.00,1.07]],
                ['31_T',[1.00,0.00],[6.20,0.00],[0.00,1.08]],
                ['31_T',[1.00,0.00],[6.20,0.00],[0.00,0.91]],
                ['28_B',[1.00,0.00],[4.88,0.56],[0.00,1.09]],
                ['32_B',[1.00,0.00],[7.00,0.00],[0.00,0.91]],
                ['32_B',[1.00,0.00],[6.33,0.00],[1.00,0.61]],
                ['32_B',[1.00,0.00],[6.00,0.00],[0.96,0.47]],
                ['112_F',[1.00,0.00],[7.40,0.37],[0.00,0.91]],
                ['112_F',[1.00,0.00],[7.00,0.00],[0.00,0.91]],
                ['112_F',[1.00,0.00],[7.25,0.43],[0.00,0.91]],
                ['72_F',[1.00,0.00],[6.11,0.24],[0.00,1.06]],
                ['72_F',[1.00,0.00],[6.31,0.38],[0.00,1.09]],
                ['72_F',[1.00,0.00],[6.50,0.00],[0.00,0.91]],
                ['142_F',[1.00,0.00],[6.67,0.47],[0.00,0.91]],
                ['142_F',[1.00,0.00],[7.58,0.37],[0.00,1.04]],
                ['142_F',[1.00,0.00],[7.00,0.00],[0.00,0.91]],
                ['42_F',[1.00,0.00],[4.86,0.99],[0.81,1.03]],
                ['42_F',[1.00,0.00],[6.50,0.00],[0.00,0.92]],
                ['112_T',[1.00,0.00],[6.33,0.00],[0.00,0.93]],
                ['112_T',[1.00,0.00],[7.33,0.00],[0.00,0.91]],
                ['112_T',[1.00,0.00],[6.00,0.29],[0.00,0.91]],
                ['72_B',[1.00,0.00],[6.64,0.26],[0.00,0.91]],
                ['72_B',[1.00,0.00],[6.29,0.35],[0.00,1.10]],
                ['72_B',[1.00,0.00],[6.75,0.00],[0.00,0.91]],
                ['122_B',[1.00,0.00],[7.10,0.25],[0.00,0.91]],
                ['122_B',[1.00,0.00],[6.67,0.47],[0.00,0.91]],
                ['122_B',[1.00,0.00],[6.33,0.47],[0.98,0.84]],
                ['162_B',[1.00,0.00],[7.86,0.27],[0.00,1.03]],
                ['162_B',[1.00,0.00],[7.75,0.43],[0.00,1.02]],
                ['162_B',[1.00,0.00],[7.40,0.39],[0.00,1.02]],
                ['31_B',[1.00,0.00],[6.00,0.00],[0.00,1.11]],
                ['31_B',[1.00,0.00],[6.67,0.00],[0.00,0.91]],
                ['29_T',[1.00,0.00],[4.83,0.00],[0.00,1.02]],
                ['82_B',[1.00,0.00],[7.00,0.00],[0.00,0.91]],
                ['82_B',[1.00,0.00],[7.20,0.32],[0.00,0.92]],
                ['82_B',[1.00,0.00],[7.25,0.00],[0.00,0.91]],
                ['30_T',[1.00,0.00],[6.00,0.00],[0.00,1.07]],
                ['30_T',[1.00,0.00],[6.00,0.00],[0.00,1.04]],
                ['30_T',[1.00,0.00],[5.25,0.25],[0.00,1.13]],
                ['62_B',[1.00,0.00],[7.29,0.27],[0.00,0.91]],
                ['62_B',[1.00,0.00],[7.25,0.35],[0.00,0.93]],
                ['62_B',[1.00,0.00],[7.33,0.00],[0.00,1.07]],
                ['142_T',[1.00,0.00],[7.33,0.00],[0.00,0.91]],
                ['142_T',[1.00,0.00],[7.00,0.00],[0.00,0.91]],
                ['29_B',[1.00,0.00],[6.00,0.00],[0.00,0.94]],
                ['29_B',[1.00,0.00],[4.88,0.00],[0.00,1.07]],
                ['62_F',[1.00,0.00],[6.33,0.41],[1.00,0.79]],
                ['62_F',[1.00,0.00],[6.40,0.00],[0.96,0.89]],
                ['62_F',[1.00,0.00],[7.00,0.00],[0.00,0.91]],]
                
    for (name,d,g,cinv) in koepferl_sim:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='s')
if 0:
    higal_longsets=[["All tiles",[1.58,0.00],[6.00,0.00],[0.75,0.18]],
                    ['long_218',[1.44,0.25],[5.25,0.43],[1.00,0.75]],
                    ['long_221',[1.12,0.23],[5.80,1.17],[0.94,0.44]],
                    ['long_224',[2.27,0.12],[4.17,0.37],[0.63,0.07]]]
                    
    for (name,d,g,cinv) in higal_longsets:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='s')
if 0:
    higal_types=[['Unbound',[1.00,0.00],[7.41,0.29],[0.67,1.11]],
                 ['Prestellar',[1.41,0.26],[6.43,0.73],[0.79,0.23]],
                 ['Protostellar',[1.00,0.00],[7.71,0.00],[0.83,0.24]],
                 ['All',[1.58,0.00],[6.00,0.00],[0.75,0.18]]]
                      
    for (name,d,g,cinv) in higal_types:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='o')

if 0:
    smilgys=[['cluster12',[1.00,0.00],[5.00,0.00],[0.00,0.98]],
             ['cluster7',[1.00,0.00],[6.50,0.50],[0.00,0.91]],
             ['cluster16',[1.00,0.00],[4.38,0.29],[0.00,1.10]],
             ['cluster36',[1.06,0.16],[3.83,0.43],[0.95,0.89]],
             ['cluster3',[1.23,0.28],[6.20,1.17],[0.97,0.49]],
             ['cluster8',[1.00,0.00],[4.33,0.47],[0.00,1.02]],
             ['cluster20',[1.08,0.20],[3.57,0.49],[0.00,1.10]],
             ['cluster18',[1.36,0.28],[3.12,0.33],[0.00,0.91]],
             ['cluster5',[1.00,0.00],[6.42,0.34],[0.00,1.02]],
             ['cluster19',[1.00,0.00],[4.33,0.47],[0.00,1.06]],
             ['cluster4',[1.00,0.00],[8.00,0.00],[0.00,1.02]],
             ['cluster13',[1.00,0.00],[5.67,0.94],[0.00,0.97]],
             ['cluster22',[1.00,0.00],[4.83,0.00],[0.00,0.94]],
             ['cluster9',[1.19,0.27],[5.67,0.47],[0.00,0.93]],
             ['cluster10',[1.25,0.29],[3.86,0.83],[0.00,1.06]],]
            
    for (name,d,g,cinv) in smilgys:
        #if g[0]>4:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='o')
        #else:
        #    plt.plot([d[0]],[g[0]],[cinv[0]],'ko',alpha=0.3)
if 0:
    sills=[['DR 21',[1.58,0.00],[6.86,0.35],[0.00,0.60]]]
    
    for (name,d,g,cinv) in sills:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='o')
                      
                      
if 0:
    taurus070=[['100',[2.21,0.23],[3.20,0.40],[0.75,0.18]],#
               ['500',[2.25,0.22],[4.10,0.30],[0.69,0.12]],#
               ['1000',[1.72,0.20],[5.83,0.69],[0.00,0.97]],#
               ['2000',[1.87,0.22],[6.31,0.58],[0.00,0.63]],#
               ['3000',[2.13,0.16],[5.58,0.49],[0.00,0.71]],#
               ['4000',[2.10,0.13],[5.70,0.40],[0.65,0.57]],#
               ['5000',[2.00,0.00],[6.00,0.00],[0.63,0.30]],]#
       
    for (name,d,g,cinv) in taurus070:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='o')
    plot_3d_line(ax,taurus070,'070')
if 0:
    taurus160=[['100',[2.02,0.33],[3.67,0.75],[0.00,0.91]],
               ['500',[1.85,0.26],[5.33,0.67],[0.00,0.93]],
               ['1000',[1.77,0.21],[5.82,0.83],[1.00,0.64]],
               ['2000',[1.95,0.14],[6.12,0.33],[0.00,0.70]],
               ['4000',[1.79,0.21],[6.50,0.50],[0.67,0.53]],]
    plot_3d_line(ax,taurus160,'160')  
    for (name,d,g,cinv) in taurus160:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='o')
                 
if 0:
    taurus250=[['100',[1.00,0.00],[6.67,0.47],[0.98,0.87]],
               ['500',[1.05,0.15],[7.58,0.00],[0.82,0.23]],
               ['1000',[1.00,0.00],[7.75,0.43],[1.02,0.75]],
               ['2000',[1.29,0.29],[7.50,0.50],[0.00,1.12]],
               ['3000',[1.44,0.00],[8.00,0.00],[0.24,0.33]],
               ['4000',[1.58,0.00],[7.80,0.32],[0.24,0.33]],
               ['5000',[1.58,0.00],[7.67,0.47],[0.24,0.33]],
               ['10000',[1.58,0.00],[7.75,0.35],[0.24,0.33]],]
    plot_3d_line(ax,taurus250,'250')
    for (name,d,g,cinv) in taurus250:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='s')
              
if 0:
    taurus350=[['100',[1.00,0.00],[6.17,0.90],[0.94,0.43]],
               ['500',[1.17,0.26],[6.86,0.35],[1.01,0.77]],
               ['1000',[1.13,0.21],[7.69,0.39],[1.02,0.78]],
               ['2000',[1.44,0.25],[7.25,0.43],[0.58,0.21]],
               ['4000',[1.58,0.00],[7.67,0.33],[0.24,0.33]],
               ['6000',[1.63,0.00],[7.38,0.46],[0.54,0.14]],]
    plot_3d_line(ax,taurus350,'350')
    for (name,d,g,cinv) in taurus350:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='^')

if 0:
    taurus500=[['100',[1.00,0.00],[6.17,0.90],[0.95,0.45]],
               ['500',[1.21,0.27],[6.69,0.60],[1.01,0.72]],
               ['1000',[1.00,0.00],[8.00,0.00],[1.02,0.75]],
               ['2000',[1.29,0.29],[7.50,0.50],[0.00,1.12]],
               ['4000',[1.58,0.00],[7.67,0.33],[0.24,0.33]],
               ['6000',[1.60,0.09],[7.63,0.56],[0.24,0.33]]]
    plot_3d_line(ax,taurus500,'500')
    for (name,d,g,cinv) in taurus500:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='*')

if 0:
    orionA=[['100',[1.00,0.00],[5.00,0.00],[0.85,0.28]],
            ['500',[1.55,0.29],[6.00,0.58],[0.00,0.93]],
            ['1000',[1.58,0.00],[6.78,0.42],[0.00,0.61]],
            ['2000',[1.78,0.25],[6.50,0.65],[0.00,0.68]],
            ['3000',[1.58,0.00],[7.00,0.00],[0.65,0.58]],
            ['4000',[1.79,0.21],[6.50,0.50],[0.63,0.28]],
            ['5000',[1.76,0.15],[6.75,0.52],[0.59,0.23]],]

    for (name,d,g,cinv) in orionA:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='o')
   

if 1:
    star_original=[['Cluster',[1.58,0],[6,0],[0,0]]]
    for (name,d,g,cinv) in star_original:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='*')     
if 1:
    cloud_10=[['100',[1.25,0.29],[4.86,0.35],[0.84,0.27]],
                ['500',[1.55,0.29],[6.00,0.58],[0.00,0.98]],
                ['1000',[1.69,0.18],[5.75,0.43],[0.85,0.31]],
                ['2000',[1.46,0.23],[7.20,0.40],[0.00,0.74]],
                ['4000',[1.72,0.20],[7.00,0.82],[0.60,0.24]],
                ['6000',[1.96,0.11],[6.10,0.27],[0.55,0.16]],
                ['original (336)',[1.99,0.18],[4.22,0.63],[0.97,0.78]],]

    for (name,d,g,cinv) in cloud_10:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='o')
                      
if 1:
    cloud_5=[['std_5_100',[1.00,0.00],[5.00,0.00],[0.84,0.28]],
             ['std_5_500',[1.58,0.00],[5.82,0.39],[0.00,1.11]],
             ['std_5_1000',[1.77,0.21],[5.82,0.83],[0.82,0.27]],
             ['std_5_2000',[1.78,0.25],[6.50,0.65],[0.44,0.74]],
             ['std_5_4000.dat',[1.58,0.00],[7.33,0.47],[0.59,0.21]],]

    for (name,d,g,cinv) in cloud_5:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='s')


if 1:
    cloud_pix_20=[['100',[1.25,0.29],[4.86,0.35],[0.84,0.27]],
                  ['500',[1.66,0.17],[5.40,0.49],[0.74,0.17]],
                  ['1000',[1.58,0.00],[6.55,0.50],[0.00,0.61]],
                  ['2000',[1.87,0.22],[6.31,0.58],[0.62,0.61]],
                  ['4000',[1.64,0.14],[7.13,0.62],[0.61,0.26]],
                  ['6000',[2.00,0.00],[6.00,0.00],[0.55,0.16]],]
    plot_3d_line(ax,cloud_pix_20,'20x20 pixels')
    for (name,d,g,cinv) in cloud_pix_20:
        plot_3d_error(ax,d[0],g[0],cinv[0],d[1],g[1],cinv[1],
                      label=name,mark='^')

                      
plt.legend()

#ax.set_xlim((0.5,3))
#ax.set_ylim((2,9))
#ax.set_zlim((-0.1,0.6))

ax.set_xlabel(r'$\cal{D}$')
ax.set_ylabel(r'$\cal{G}$')
ax.set_zlabel(r'$\cal{C}^{-1}$')

plt.show()