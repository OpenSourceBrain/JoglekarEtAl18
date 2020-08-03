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
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Utils/') 
import utils


# Arguments
N=int(sys.argv[1])
figure=str(sys.argv[2])
saveData=str(sys.argv[3])
saveFigure=str(sys.argv[4])
useSaveData=str(sys.argv[5])

# Path
path = os.path.abspath(os.getcwd())
 
# Simulations

# Figure 5B
if figure== '5B':
    
    file= path+'/files/spikes_figure'+figure+'_'+str(N)+'areas.txt'
    duration=800
    
    if useSaveData!='yes':
        monitors = run_network(N,'asynchronous','weak',duration)
        xValues=monitors.t/ms
        yValues=monitors.i
        # Save data
        if saveData== 'yes':
            np.savetxt(file,np.column_stack([monitors.t/ms,monitors.i]))     
    else:
        temp=np.loadtxt(file)
        xValues=temp[:,0]
        yValues=temp[:,1]

    utils.rasterPlot(xValues,yValues,duration,figure,N,saveFigure,path)

# Figure 5C
elif figure== '5C':
    
    file= path+'/files/spikes_figure'+figure+'_'+str(N)+'areas.txt'
    duration=800
    
    if useSaveData!='yes':
        monitors = run_network(N,'asynchronous','strong',duration)
        xValues=monitors.t/ms
        yValues=monitors.i
        # Save data
        if saveData== 'yes':
            np.savetxt(file,np.column_stack([monitors.t/ms,monitors.i]))     
    else:
        temp=np.loadtxt(file)
        xValues=temp[:,0]
        yValues=temp[:,1]

    utils.rasterPlot(xValues,yValues,duration,figure,N,saveFigure,path)

# Figure 5E
elif figure== '5E':
    
    file1= path+'/files/spikes_figure5B_'+str(N)+'areas.txt'
    file2= path+'/files/spikes_figure5C_'+str(N)+'areas.txt'
    
    duration=800
    
    if useSaveData!='yes':
        monitorsBad = run_network(N,'asynchronous','weak',duration)
        bad=np.column_stack([monitorsBad.t/ms,monitorsBad.i])
        if saveData== 'yes':
            np.savetxt(file1,bad) 
    else:
        bad=np.loadtxt(file1)
        
    if useSaveData!='yes':
        monitorsGood = run_network(N,'asynchronous','strong',duration)
        good=np.column_stack([monitorsGood.t/ms,monitorsGood.i])
        if saveData== 'yes':
            np.savetxt(file2,good) 
    else:
        good=np.loadtxt(file2)

    maxratebad,_= utils.firingRate(N,bad,duration)
    maxrategood,_= utils.firingRate(N,good,duration)
    utils.firingRatePlot(maxrategood,maxratebad,figure,N,saveFigure,path)
    
# Figure 6A
elif figure== '6A':
    
    file= path+'/files/spikes_figure'+figure+'_'+str(N)+'areas.txt'
    duration=440
    
    if useSaveData!='yes':
        monitors = run_network(N,'synchronous','weak',duration)
        xValues=monitors.t/ms
        yValues=monitors.i
        # Save data
        if saveData== 'yes':
            np.savetxt(file,np.column_stack([monitors.t/ms,monitors.i]))     
    else:
        temp=np.loadtxt(file)
        xValues=temp[:,0]
        yValues=temp[:,1]

    utils.rasterPlot(xValues,yValues,duration,figure,N,saveFigure,path)


# Figure 6B
elif figure== '6B':    
    
    file= path+'/files/spikes_figure'+figure+'_'+str(N)+'areas.txt'
    duration=440
    
    if useSaveData!='yes':
        monitors = run_network(N,'synchronous','strong',duration)
        xValues=monitors.t/ms
        yValues=monitors.i
        # Save data
        if saveData== 'yes':
            np.savetxt(file,np.column_stack([monitors.t/ms,monitors.i]))     
    else:
        temp=np.loadtxt(file)
        xValues=temp[:,0]
        yValues=temp[:,1]

    utils.rasterPlot(xValues,yValues,duration,figure,N,saveFigure,path)

# Figure 6D
elif figure== '6D':
    
    file1= path+'/files/spikes_figure6A_'+str(N)+'areas.txt'
    file2= path+'/files/spikes_figure6B_'+str(N)+'areas.txt'
    
    duration=440
    
    if useSaveData!='yes':
        monitorsBad = run_network(N,'synchronous','weak',duration)
        bad=np.column_stack([monitorsBad.t/ms,monitorsBad.i])
        if saveData== 'yes':
            np.savetxt(file1,bad) 
    else:
        bad=np.loadtxt(file1)
        
    if useSaveData!='yes':
        monitorsGood = run_network(N,'synchronous','strong',duration)
        good=np.column_stack([monitorsGood.t/ms,monitorsGood.i])
        if saveData== 'yes':
            np.savetxt(file2,good) 
    else:
        good=np.loadtxt(file2)
    

    maxratebad,_= utils.firingRate(N,bad,duration)
    maxrategood,_= utils.firingRate(N,good,duration)
    utils.firingRatePlot(maxrategood,maxratebad,figure,N,saveFigure,path)
    


