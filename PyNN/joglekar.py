#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 21:08:32 2020

@author: ronaldo
"""
#
#
# Work in progress
#
#

from pyNN.utility import get_script_args, Timer, ProgressBar
from pyNN.random import NumpyRNG, RandomDistribution
from importlib import import_module
from scipy.io import loadmat
import numpy as np
import os

# Testing with Brian
from pyNN.neuroml import *

''' To Do
1) Compare increment in current with increment in voltage (create two programs in Brian and compare)
2) Check how to add input 
3) Original code uses a delta function for synapse model, maybe it is 
    possible to use an alpha function (with a very small tau) 
    (see implementation of Brunel in PyNN)
4) Check if it is running in NEST
5) Convert anatomical data to pickle or other format (not .mat)
6) Create function to choose between synchronous or asynchronous and weak or strong GBA                                 
7) Noise in integrate and fire model (stochastic model)                                 
'''

# Just for test
setup()



def setConnections():
    path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Matlab/' #path to .mat files    
    #==================== Fln ==========================#
    dataTemp = loadmat(path +'subgraphData.mat')
    conn=dataTemp['flnMat']
    # Remove area LIP (30). It is not used in Joglekar's paper
    conn=conn[0:29,0:29]
    connBin=(conn>0)*1  
    #==================== distances ====================#
    dataTemp = loadmat(path +'subgraphWiring29.mat')
    # wiring is a 30x30 distance matrix with values given in mm
    dist=dataTemp['wiring']
    # Remove area LIP (30). It is not used in Joglekar's paper
    dist=dist[0:29,0:29]
    #==================== hierarchy ====================#
    hierVals = loadmat(path +'hierValspython.mat')
    hierValsnew = hierVals['hierVals'][:]
    #hierarchy normalized. 
    hierNorm=hierValsnew/max(hierValsnew)

    return conn,connBin,dist,hierNorm


#============================= Parameters ====================================# 

R            = 50              # membrane resitance [MOhm]
tauMemE      = 20.0            # excitatory neuron membrane time constant [ms]
tauMemI      = 10.0            # inhibirtory neuron membrane time constant [ms]
CmE          = tauMemE/R       # Capacitance excitatory neurons [nF]
CmI          = tauMemI/R       # Capacitance inhibitory neurons [nF]
tauRef       = 2.0             # refractory time [ms]
Vl           = -70.0           # resting potential [mV]
Vt           = -50.0           # threshold [mV]
Vr           = -60.0           # reset [mV]  
Nareas       = 29              # number of cortical areas
NE           = 40             # number of excitatory neurons per area
NI           = 10             # number of inhibitory neurons per area
prob         = 0.1             # connection probability (local and interareal)
velocity     = 3.5             # axonal conduction velocity [m/s] or [mm/ms]
conn, connBin, dist, hier =setConnections()    # Define fln and distances

paramsE = {'cm': CmE, 
           'tau_m': tauMemE, 
           'v_rest': Vl, 
           'v_thresh': Vt, 
           'tau_refrac': tauRef,
           'v_reset': Vr}


paramsI = {'cm': CmI, 
           'tau_m': tauMemI, 
           'v_rest': Vl, 
           'v_thresh': Vt, 
           'tau_refrac': tauRef,
           'v_reset': Vr}

#===================== Parameters from Brunel  ===============================#
tauSyn      = 0.1     # synaptic time constant [ms]
# synaptic weights, scaled for alpha functions, such that
# for constant membrane potential, charge J would be deposited
fudge = 0.00041363506632638  # ensures dV = J at V=0

# excitatory weight: JE = J_eff / tauSyn * fudge
JE = (J_eff/tauSyn)*fudge

#===================== Parameters for  connections ===========================#

# random connections
connector = FixedProbabilityConnector(p_connect=prob) 

# synaptic weights for asynchronous regime and strong GBA  
wEE = 0.01       # excitatory -> excitatory (intra-areal)
wEI = 0.075      # inhibitory -> excitatory (intra-areal)
wIE = 0.05       # excitatory -> inhibitory (intra-areal)
wII = 0.075      # inhibitory -> inhibitory (intra-areal)
muEE = 0.05      # excitatory -> excitatory (inter-areal)
muIE = 19/4      # excitatory -> inhibitory (inter-areal)

# Alpha parameters (related to hierarchy)
alpha=.68 # it changes for synchronous 

# delays for intra-areal connections 
delay_local = 2  # [ms]

# Lists to store populations ans projectios
pE=[]        # List of excitatory populations
pI=[]        # List of inhibitory populations
intraEE=[]   # List of EE intra-areal projections
intraEI=[]   # List of EI intra-areal projections
intraIE=[]   # List of IE intra-areal projections
intraII=[]   # List of II intra-areal projections
interEE=[]   # List of EE inter-areal projections
interIE=[]   # List of IE inter-areal projections

#=========== Create populations and define intra-areal connections ===========#
for i in range(Nareas):
#========================= Create Populations ================================#
    pE.append(Population(NE, IF_cond_exp,paramsE,label="Exc"+str(i+1)))
    pI.append(Population(NI, IF_cond_exp,paramsI,label="Inh"+str(i+1)))
#======================= Synapses for intra-areal ============================#
    synEE=StaticSynapse(weight=(1+alpha*hier[i][0])*wEE, delay=delay_local)
    synEI=StaticSynapse(weight=-wEI, delay=delay_local)
    synIE=StaticSynapse(weight=(1+alpha*hier[i][0])*wIE, delay=delay_local)
    synII=StaticSynapse(weight=-wII, delay=delay_local)  
#======================== Intra-Areal Connections ============================#
    intraEE.append(Projection(pE[i], pE[i], connector, synEE))
    intraEI.append(Projection(pI[i], pE[i], connector, synEI))
    intraIE.append(Projection(pE[i], pI[i], connector, synIE))
    intraII.append(Projection(pI[i], pI[i], connector, synII))

#========================== Inter-areal connections ==========================#
for i in range(Nareas): # source
    
    targets=np.where(connBin[i,:]==1)[0] # targets for source i
    
    for j in targets:
        
        # delay for inter-areal connections 
        delayMean=dist[i,j]/velocity  # mean of a gaussian distribution
        delayVar=0.1*delayMean        # std for a gaussian distribution (10% of the mean)
        delayInter = RandomDistribution('normal', [delayMean, delayVar], rng=NumpyRNG(seed=4242))
        
        # Synapse types for inter-areal projections
        synInterEE=StaticSynapse(weight=(1+alpha*hier[j][0])*conn[i,j]*muEE, delay=delayInter)
        synInterIE=StaticSynapse(weight=(1+alpha*hier[j][0])*conn[i,j]*muIE, delay=delayInter)
        
        # Inter-areal projections from i to j
        interEE.append(Projection(pE[i], pE[j], connector, synInterEE))
        interIE.append(Projection(pE[i], pI[j], connector, synInterIE))
        
    
run(100.0)

