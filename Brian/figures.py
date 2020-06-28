#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 14:16:08 2020

@author: ronaldo
"""
from brian2 import *
from generateNetwork import *
import matplotlib.pyplot as plt
import sys
import os

# Raster Plot
def rasterPlot(xValues,yValues,duration,figure,N,saveFigure):
    
    ticks=[0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,
       8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,
       16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,
       24.5,25.5,26.5,27.5,28.5]   

    areasName=['V1','V2','V4','DP','MT',
           '8m','5','8l','TEO','2','F1',
           'STPc','7A','46d','10','9/46v',
           '9/46d','F5','TEpd','PBr','7m','7B',
           'F2','STPi','PROm','F7','8B','STPr','24c']   

    plt.figure()
    plt.plot(xValues, 1.0*yValues/(4*400), '.',markersize=1)
    plt.plot([0, duration/ms], np.arange(N+1).repeat(2).reshape(-1, 2).T, 'k-')
    plt.ylabel('Area')
    plt.yticks(np.arange(N))
    plt.xlabel('time (ms)')
    ylim(0,N)
    yticks(ticks[:N],areasName[:N])
    xlim(0,duration/ms)

    # Save figure
    if saveFigure== 'yes':
        plt.savefig(path+'/figures/figure'+figure+'_'+str(N)+'areas.png')     

    plt.show()

    return 0

# Plot for Mnaximum firing rate
def firingRatePlot(maxrategood,maxratebad,figure,N,saveFigure):
    
    areasName=['V1','V2','V4','DP','MT',
           '8m','5','8l','TEO','2','F1',
           'STPc','7A','46d','10','9/46v',
           '9/46d','F5','TEpd','PBr','7m','7B',
           'F2','STPi','PROm','F7','8B','STPr','24c']   
    
    
    plt.semilogy(range(N), maxratebad, '.-',color = (0.4660, 0.6740, 0.1880),linewidth = 1,label='weak GBA')
    plt.semilogy(range(N), maxrategood,'.-', color = (0.4940,0.1840,0.5560),linewidth = 1,label='strong GBA')
    plt.xticks(range(N), areasName)
    plt.xticks(rotation=90) 
    plt.xlabel('Area')
    plt.ylabel('Maximum Firing Rate (Hz)')
    plt.legend(('weak GBA', 'strong GBA'))
    
    # Save figure
    if saveFigure== 'yes':
        plt.savefig(path+'/figures/figure'+figure+'_'+str(N)+'areas.png')    
    
    plt.show()    
    
    return 0


# Arguments
N=int(sys.argv[1])
figure=str(sys.argv[2])
saveData=str(sys.argv[3])
saveFigure=str(sys.argv[4])
useSaveData=str(sys.argv[4])

# Path
path = os.path.abspath(os.getcwd())
 
# Simulations

# Figure 5B
if figure== '5B':
    
    file= path+'/files/spikes_figure'+figure+'_'+str(N)+'areas.txt'
    duration=800*ms
    
    if useSaveData!='yes':
        monitors = run_network(N,'asynchronous','weak',duration/ms)
        xValues=monitors.t/ms
        yValues=monitors.i
        # Save data
        if saveData== 'yes':
            np.savetxt(file,np.column_stack([monitors.t/ms,monitors.i]))     
    else:
        temp=np.loadtxt(file)
        xValues=temp[:,0]
        yValues=temp[:,1]

    rasterPlot(xValues,yValues,duration,figure,N,saveFigure)

# Figure 5C
elif figure== '5C':
    
    file= path+'/files/spikes_figure'+figure+'_'+str(N)+'areas.txt'
    duration=800*ms
    
    if useSaveData!='yes':
        monitors = run_network(N,'asynchronous','strong',duration/ms)
        xValues=monitors.t/ms
        yValues=monitors.i
        # Save data
        if saveData== 'yes':
            np.savetxt(file,np.column_stack([monitors.t/ms,monitors.i]))     
    else:
        temp=np.loadtxt(file)
        xValues=temp[:,0]
        yValues=temp[:,1]

    rasterPlot(xValues,yValues,duration,figure,N,saveFigure)

# Figure 5E
elif figure== '5E':
    
    file1= path+'/files/spikes_figure5B_'+str(N)+'areas.txt'
    file2= path+'/files/spikes_figure5C_'+str(N)+'areas.txt'
    
    duration=800*ms
    
    if useSaveData!='yes':
        monitorsBad = run_network(N,'asynchronous','weak',duration/ms)
        bad=np.column_stack([monitorsBad.t/ms,monitorsBad.i])
        if saveData== 'yes':
            np.savetxt(file1,bad) 
    else:
        bad=np.loadtxt(file1)
        
    if useSaveData!='yes':
        monitorsGood = run_network(N,'asynchronous','strong',duration/ms)
        good=np.column_stack([monitorsGood.t/ms,monitorsGood.i])
        if saveData== 'yes':
            np.savetxt(file2,good) 
    else:
        good=np.loadtxt(file2)

    maxratebad, maxrategood= firingRate(N,good,bad,duration)
    firingRatePlot(maxrategood,maxratebad,figure,N,saveFigure)
    
# Figure 6A
elif figure== '6A':
    
    file= path+'/files/spikes_figure'+figure+'_'+str(N)+'areas.txt'
    duration=440*ms
    
    if useSaveData!='yes':
        monitors = run_network(N,'synchronous','weak',duration/ms)
        xValues=monitors.t/ms
        yValues=monitors.i
        # Save data
        if saveData== 'yes':
            np.savetxt(file,np.column_stack([monitors.t/ms,monitors.i]))     
    else:
        temp=np.loadtxt(file)
        xValues=temp[:,0]
        yValues=temp[:,1]

    rasterPlot(xValues,yValues,duration,figure,N,saveFigure)


# Figure 6B
elif figure== '6B':    
    
    file= path+'/files/spikes_figure'+figure+'_'+str(N)+'areas.txt'
    duration=440*ms
    
    if useSaveData!='yes':
        monitors = run_network(N,'synchronous','strong',duration/ms)
        xValues=monitors.t/ms
        yValues=monitors.i
        # Save data
        if saveData== 'yes':
            np.savetxt(file,np.column_stack([monitors.t/ms,monitors.i]))     
    else:
        temp=np.loadtxt(file)
        xValues=temp[:,0]
        yValues=temp[:,1]

    rasterPlot(xValues,yValues,duration,figure,N,saveFigure)

# Figure 6D
elif figure== '6D':
    
    file1= path+'/files/spikes_figure6A_'+str(N)+'areas.txt'
    file2= path+'/files/spikes_figure6B_'+str(N)+'areas.txt'
    
    duration=440*ms
    
    if useSaveData!='yes':
        monitorsBad = run_network(N,'synchronous','weak',duration/ms)
        bad=np.column_stack([monitorsBad.t/ms,monitorsBad.i])
        if saveData== 'yes':
            np.savetxt(file1,bad) 
    else:
        bad=np.loadtxt(file1)
        
    if useSaveData!='yes':
        monitorsGood = run_network(N,'synchronous','strong',duration/ms)
        good=np.column_stack([monitorsGood.t/ms,monitorsGood.i])
        if saveData== 'yes':
            np.savetxt(file2,good) 
    else:
        good=np.loadtxt(file2)
    

    maxratebad, maxrategood= firingRate(N,good,bad,duration)
    firingRatePlot(maxrategood,maxratebad,figure,N,saveFigure)
    


