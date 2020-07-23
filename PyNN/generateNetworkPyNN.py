#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 29 09:45:18 2020

@author: ronaldo
"""

import scipy.io
import numpy as np
import random as pyrand
import time
import os

from pyNN.utility import get_script_args, Timer, ProgressBar
from pyNN.random import NumpyRNG, RandomDistribution
from importlib import import_module

import matplotlib.pyplot as plt

    
def spikeData(pop):
    
    indexSpike=[]
    timeSpike=[]
    
    spikes =  pop.get_data('spikes', gather=False)
    spiketrains = spikes.segments[0].spiketrains
    
    for spiketrain in spiketrains:
        source_id = spiketrain.annotations['source_id']
        source_index = spiketrain.annotations['source_index']

        for t in spiketrain:
            indexSpike.append(source_index)
            timeSpike.append(t.magnitude.tolist())
    
    # Index and time
    return np.column_stack((timeSpike,indexSpike))            


def computeFiringRate(spikeMonitor,group,time_end,dt):       

    # array to store firing rate and time
    firing_rate=np.zeros((2,round((time_end+dt)/dt)))
    
    # vetor com tempos
    firing_rate[0,:]=np.arange(0,int((time_end+dt)/dt),int(dt/dt))
           
    tempList=[]
    
    for idx in range(0,len(group)):
        
        tempList.extend(np.where(spikeMonitor[:,1]==group[idx])[0].tolist()[:])
    
    tempArray=np.copy(np.asarray(tempList))

    tempArray=np.copy(spikeMonitor[tempArray])
    times=np.sort(tempArray[:,0])  

    # tempos dos disparos e quantidade de neuronios por tempo
    frTemp=np.unique(times,return_counts=True)

    # indices dos tempos dos disparos
    frTemp2=(np.round(frTemp[0]/dt)).astype(int)

    # Firing Rate
    firing_rate[1,frTemp2]=frTemp[1]* (1.0/(dt/(10**3)))/len(group)

    return firing_rate

def sliding_window(serie,width,dt):
    # width in ms   
    width_dt = int(width / 2 / dt)*2 + 1
    used_width = width_dt * dt
    window = np.ones(width_dt)
    
    return np.convolve(serie, window * 1. / sum(window), mode='same')


def gen_params(n_areas,regime, gba, duration):    
    
    
    para= {'N'         : 2000,    # Number of neurons 
           'NAreas'    : n_areas, # Number of areas
           'Ne'        : 0.8,     # fraction of excitatory neurons
           'Vr'        : -70.,    # Membrane potential rest
           'Vreset'    : -60.,    # Membrane potential reset
           'Vt'        : -50.,    # Membrane potential threshold
           'taumE'     : 20.,     # Membrane time constant (excitatory neurons) 
           'taumI'     : 10.,     # Membrane time constant (inhibitory neurons)  
           'tref'      : 2.,      # refractory time
           'probIntra' : .1,      # connection probability intra area
           'probInter' : .1,      # connection probability inter areas   
           'sigma'     : 3.,      # Noise 
           'alpha'     : 4.,      # gradient 
           'dlocal'    : 2.,      # delays local 
           'speed'     : 3.5,     # axonal conduction velocity
           'lrvar'     : 0.1,     # standard deviation delay long range
           'path'      : os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Matlab/', #path to .mat files,
           'tau_syn_e' :  0.1,
           'dt'        :  0.1 
           }
    
    #hierarchy values file 
    hierVals = scipy.io.loadmat(para["path"]+'hierValspython.mat')
    hierValsnew = hierVals['hierVals'][:]
    hier=hierValsnew/max(hierValsnew)#hierarchy normalized. 
    para['hier']=hier[:para["NAreas"]]

    #fln values file 
    flnMatp = scipy.io.loadmat(para["path"]+'efelenMatpython.mat')
    conn=flnMatp['flnMatpython'][:][:] #fln values..Cij is strength from j to i 
    para['conn']=conn[:para["NAreas"],:para["NAreas"]]

    distMatp = scipy.io.loadmat(para["path"]+'subgraphWiring29.mat')
    distMat=distMatp['wiring'][:][:] #distances between areas values..
    delayMat = distMat/para['speed']
    para['delayMat']=delayMat[:para["NAreas"],:para["NAreas"]]
    
    # Membrane resitance
    R=50
    para['R']=R
    para['duration']    = duration;   
    if regime=='asynchronous':  
    
        # general for assynchronous regime
        para['VextE']    = 14.2
        para['VextI']    = 14.7  
        para['muIE']     = .19/4 
        para['wII']      = .075
        para['wEE']      = .01
        para['wIE']      = .075
        
        if gba=='weak':
            para['wEI']      = .0375
            para['muEE']     = .0375
            para['currval']  = 15       
        
        elif gba=='strong':    
            para['wEI']      = .05
            para['muEE']     = .05    
            para['currval']  = 6.3 
      
        para['currdur']       = 150
        
    elif regime=='synchronous':
        
        # general for synchronous regime
        para['muIE']     = .19 
        para['wII']      = .3
        para['wEE']      = .04
        para['wIE']      = .3
        
        if gba=='weak':
            para['VextI']    = 14.0  
            para['VextE']    = 15.4
            para['wEI']      = .56
            para['muEE']     = .16
        elif gba=='strong':    
            para['VextI']    = 14.0  
            para['VextE']    = 16.0
            para['wEI']      = .98
            para['muEE']     = .25    
        
        para['currdur']      = 8
        para['currval']      = 10.1


    # Parameter to make synapses equivalent to delta function
    # Check the notebooks
    para['coeffE']=0.67957046
    para['coeffI']=2*para['coeffE']
    
    return para


def network(para,simulator):
    
    
    if simulator == "-nest":
        sim = import_module("pyNN.nest")    
    elif simulator == "-neuron":   
        sim = import_module("pyNN.neuron")
    elif simulator == "-brian":
        sim = import_module("pyNN.brian")

    # initialize pyNN simulation
    sim.setup(timestep=para["dt"])
    
    
    # Parameters for excitatory neurons
    E_parameters={
                "tau_m":para["taumE"], 
                "cm":para["taumE"]/para["R"], 
                "v_rest":para["Vr"], 
                "v_reset":para["Vreset"], 
                "v_thresh":para["Vt"], 
                "tau_refrac":para["tref"],
                "tau_syn_E":para["tau_syn_e"],
                "tau_syn_I":para["tau_syn_e"],
                "i_offset":(para["VextE"]/para["R"])}
    
    # Paramters for inhibitory neurons
    I_parameters={
                "tau_m":para["taumI"], 
                "cm":para["taumI"]/para["R"], 
                "v_rest":para["Vr"], 
                "v_reset":para["Vreset"], 
                "v_thresh":para["Vt"], 
                "tau_refrac":para["tref"],
                "tau_syn_E":para["tau_syn_e"],
                "tau_syn_I":para["tau_syn_e"],
                "i_offset":(para["VextI"]/para["R"])}
    
    ############ All Excitatory and inhibitory neurons ########################
    
    # number of excitatory neurons in one network
    NE_net=int(para['N']*para['Ne'])
    # Total number of excitatory neurons
    NE=int(para['NAreas']*NE_net)
    
    # number of excitatory neurons in one network
    NI_net=int(para['N']*(1-para['Ne']))
    # Total number of inhibitory neurons
    NI=int((para['NAreas']*NI_net))
    
    popE =sim.Population(NE, sim.IF_curr_alpha,E_parameters,label="popE")
    popI =sim.Population(NI, sim.IF_curr_alpha,I_parameters,label="popI")
    
    ################################ Noise ####################################
    
    # Noise on excitatory neurons
    stdNoiseE = (para["sigma"]/para["R"])*(para["taumE"]**0.5)/(para["dt"]**0.5)
    for i in range(NE):
        noise = sim.NoisyCurrentSource(mean=0, stdev=stdNoiseE, start=0.0, stop=para["duration"])
        popE[i].inject(noise)
        popE[i].inject(noise)
    
    # Noise on inhibitory neurons    
    stdNoiseI = (para["sigma"]/para["R"])*(para["taumI"]**0.5)/(para["dt"]**0.5)
    for i in range(NI):
        noise = sim.NoisyCurrentSource(mean=0, stdev=stdNoiseI, start=0.0, stop=para["duration"])
        popI[i].inject(noise)
        popI[i].inject(noise)
    
    # paramters for initial conditions
    kernelseed  = 5456532
    rng = NumpyRNG(kernelseed, parallel_safe=True)
    uniformDistr = RandomDistribution('uniform', low=para["Vr"], high=para["Vt"], rng=rng)
    sim.initialize(popE, v=uniformDistr)
    sim.initialize(popI, v=uniformDistr)
    
    # Separate population in population views
    popEList=[]
    popIList=[]
    
    # Store projections
    EE=[]
    IE=[]
    EI=[]
    II=[]
    EElongRange=[]
    IElongRange=[]
    
    for i in range(para['NAreas']):
        
        # store sub populations in lists 
        popEList.append(popE[(i*NE_net):((i+1)*NE_net)])
        popIList.append(popI[(i*NI_net):((i+1)*NI_net)])
     
        #### Synapses
        
        # Weights for recurrent connections
        wEE_alpha = (((1+para['alpha']*para["hier"][i])*para["wEE"])/para['coeffE'])[0]      #[nA]
        wIE_alpha = (((1+para['alpha']*para["hier"][i])*para["wIE"])/para['coeffI'])[0]      #[nA]  
        wEI_alpha = (para["wEI"]/para['coeffE'])*-1                          #[nA]
        wII_alpha = (para["wII"]/para['coeffI'])*-1                          #[nA]
    
        # Connections    
        EE_connections = sim.Projection(popEList[i], popEList[i], sim.FixedProbabilityConnector(p_connect=para["probIntra"]),
                                         sim.StaticSynapse(weight=wEE_alpha, delay=para["dlocal"]))
        IE_connections = sim.Projection(popEList[i], popIList[i], sim.FixedProbabilityConnector(p_connect=para["probIntra"]),
                                         sim.StaticSynapse(weight=wIE_alpha, delay=para["dlocal"]))
        EI_connections = sim.Projection(popIList[i], popEList[i], sim.FixedProbabilityConnector(p_connect=para["probIntra"]),
                                         sim.StaticSynapse(weight=wEI_alpha, delay=para["dlocal"]))
        II_connections = sim.Projection(popIList[i], popIList[i], sim.FixedProbabilityConnector(p_connect=para["probIntra"]),
                                        sim.StaticSynapse(weight=wII_alpha, delay=para["dlocal"]))
        
        # Store projections in lists
        EE.append(EE_connections)
        IE.append(IE_connections)
        EI.append(EI_connections)
        II.append(II_connections)

    
    # Long Range connections
    for i in range(para['NAreas']):
        for j in range(para['NAreas']):
            
            if i!=j:
                
                # Weights
                wEE_alphaLR =  ((1 + para['alpha']*para["hier"][j]) * para['muEE'] * para["conn"][j,i])[0]/para['coeffE']
                wIE_alphaLR =  ((1 + para['alpha']*para["hier"][j]) * para['muIE'] * para["conn"][j,i])[0]/para['coeffI']
                
                # Delay
                # Mean for delay distribution 
                meanlr = para["delayMat"][j,i]
                # Standard deviation for delay distribution
                varlr =  para['lrvar']*meanlr
                dLR=RandomDistribution('normal', [meanlr, varlr], rng=NumpyRNG(seed=4242))
                
                # Connections
                EE_connectionsLR = sim.Projection(popEList[i], popEList[j], sim.FixedProbabilityConnector(p_connect=para["probInter"]),
                                         sim.StaticSynapse(weight=wEE_alphaLR, delay=dLR))
                IE_connectionsLR = sim.Projection(popEList[i], popIList[j], sim.FixedProbabilityConnector(p_connect=para["probInter"]),
                                         sim.StaticSynapse(weight=wIE_alphaLR, delay=dLR))
       
                # Store projections in list
                EElongRange.append(EE_connectionsLR)
                IElongRange.append(IE_connectionsLR)
    
    
    # Stimulus
    amplitude =para['currval']/para['R'] #[mV]
    pulse = sim.DCSource(amplitude=amplitude, start=300.0, stop=300.0+para['currdur'])
    pulse.inject_into(popEList[0])
    
    
    # Record data 
    popE.record('spikes') 
    popI.record('spikes') 
    
    # Run
    sim.run(para['duration'])
        
    # Store spikes
    spikesE_in = popE.get_data()
    spikesI_in = popI.get_data()
    
    
    # Generate array with spike data for popE
    spkpopE=spikeData(popE)
    spkpopI=spikeData(popI)
    
    return spkpopE, spkpopI

def run_network(n_areas,regime,gba,duration,simulator):
    
    # Parameters
    para = gen_params(n_areas,regime,gba,duration) 
    # Run Network 
    E_spike,_= network(para,simulator)
    
    return E_spike

def firingRate(N,goodprop,badprop,duration):

    binsize = 10
    stepsize =  1  
    
    # Store maximum firing rate for each area
    maxratebad  = np.empty([N,1])
    maxrategood = np.empty([N,1])
    
    # sort net spikes
    netspikebad = len(badprop)
    netspikegood= len(goodprop)
    
    badpropsorted = badprop[badprop[:,1].argsort(),]
    goodpropsorted = goodprop[goodprop[:,1].argsort(),]
            
    netbinno = int( 1+(duration)-(binsize))
    popratebad = np.empty([N,netbinno ])
    poprategood = np.empty([N,netbinno ])
            
     
    countbad = 0
    countgood = 0#for each spike. 
            
    monareaktimeallbad = []
    monareaktimeallgood = []
                    
    for u in range(N):
        monareaktimebad = []
        monareaktimegood = []
        
        while((countbad < netspikebad) and (badpropsorted[countbad,1]<1600*(u+1)) ):
          monareaktimebad.append(badpropsorted[countbad,0])#append spike times for each area.
          countbad = countbad + 1
          
        while((countgood < netspikegood) and (goodpropsorted[countgood,1]<1600*(u+1)) ):
          monareaktimegood.append(goodpropsorted[countgood,0])#append spike times for each area.
          countgood = countgood + 1
          
        valsbad = np.histogram(monareaktimebad, bins=int(duration/stepsize))
        valsgood = np.histogram(monareaktimegood, bins=int(duration/stepsize))
        
        valszerobad = valsbad[0]
        valszerogood = valsgood[0]
    
        astep = binsize/1
        
        valsnewbad = np.zeros(netbinno)
        valsnewgood = np.zeros(netbinno)
        
        acount = 0
        while acount < netbinno:        
            valsnewbad[acount] = sum(valszerobad[int(acount):int(acount+astep)])
            valsnewgood[acount] = sum(valszerogood[int(acount):int(acount+astep)])
            acount=acount+1
    
        valsratebad = valsnewbad*((1000/binsize) /(1600) ) # divide by no of neurons per E pop. 
        valsrategood = valsnewgood*((1000/binsize) /(1600) )    
        popratebad[u,:], poprategood[u,:] = valsratebad, valsrategood    
        
        #compute population firing rates. 
        
        maxratebad[u,0] = max(valsratebad[int(len(valsratebad)/3):])
        maxrategood[u,0] = max(valsrategood[int(len(valsrategood)/3):])
        
            
    return maxratebad, maxrategood

    