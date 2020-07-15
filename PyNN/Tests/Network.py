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

import matplotlib.pyplot as plt

from pyNN.utility import get_script_args, Timer, ProgressBar
from pyNN.random import NumpyRNG, RandomDistribution
import sys


if sys.argv[1] == "-nest":
    from pyNN.nest import *    
elif sys.argv[1] == "-neuron":   
    from pyNN.neuron import *
elif sys.argv[1] == "-brian":
    from pyNN.brian import *
   
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
    
###############################################################################
############################## Parameters #####################################
###############################################################################
NE            = 1600
NI            = 400 
R             = 50.0         # membrane resitance [MOhm]
tauE_m        = 20.0         # excitatory neuron membrane time constant [ms]
tauI_m        = 10.0         # inhibitory neuron membrane time constant [ms]
cE_m          = tauE_m/R     # Capacitance excitatory neurons [nF]
cI_m          = tauI_m/R     # Capacitance excitatory neurons [nF]
tauRef        = 2.0          # refractory time [ms]
Vrest         = -70.0        # resting potential [mV]
Vt            = -50.0        # threshold [mV]
Vreset        = -60.0        # reset [mV]  
tau_syn_e     =  0.1         # time constant for synapse [ms] 

# External input for assynchronous behavior
VextE         = 14.2         # External input in excitatory neurons [mV]
VextI         = 14.7         # External input in excitatory neurons [mV]

# Noise
sigmaV=3.0   # [mV]

# Coefficient (delta to alpha synapse)
coeffE=0.67957046
coeffI=2*coeffE

# Synaptic weights [mV]
wEE = .01 
wIE = .075 
wEI = -.0375 
wII = -.075 

# Connection probability
plocal = 0.1

# delay local
d= 2.0    #[ms]
###############################################################################
############################## PyNN  ##########################################
###############################################################################
# timestep
dt=0.1
# initialize pyNN simulation
setup(timestep=dt)
    
# duration
duration=1000

# Excitatory popuation
popE =Population(NE, IF_curr_alpha,IF_curr_alpha.default_parameters,label="popE")

# Inhibitory popuation
popI =Population(NI, IF_curr_alpha,IF_curr_alpha.default_parameters,label="popI")

# Set parameter for excitatory population
popE.set(tau_m=tauE_m) 
popE.set(cm=cE_m) 
popE.set(v_rest=Vrest) 
popE.set(v_reset=Vreset) 
popE.set(v_thresh=Vt) 
popE.set(tau_refrac=tauRef)
popE.set(tau_syn_E=tau_syn_e)
popE.set(tau_syn_I=tau_syn_e)
popE.set(i_offset=(VextE/R))

# Set parameter for inhibitory population
popI.set(tau_m=tauI_m) 
popI.set(cm=cI_m) 
popI.set(v_rest=Vrest) 
popI.set(v_reset=Vreset) 
popI.set(v_thresh=Vt) 
popI.set(tau_refrac=tauRef)
popI.set(tau_syn_E=tau_syn_e)
popI.set(tau_syn_I=tau_syn_e)
popI.set(i_offset=(VextI/R))

# reset weigths for alpha synapse
wEE_alpha = (wEE/coeffE)      #[nA]
wIE_alpha = (wIE/coeffI)      #[nA]  
wEI_alpha = (wEI/coeffE)      #[nA]
wII_alpha = (wII/coeffI)      #[nA]

# Noise on excitatory neurons
stdNoiseE = (sigmaV/R)*(tauE_m**0.5)/(dt**0.5)
for i in range(NE):
    noise = NoisyCurrentSource(mean=0, stdev=stdNoiseE, start=0.0, stop=duration)
    popE[i].inject(noise)

# Noise on inhibitory neurons    
stdNoiseI = (sigmaV/R)*(tauI_m**0.5)/(dt**0.5)
for i in range(NI):
    noise = NoisyCurrentSource(mean=0, stdev=stdNoiseI, start=0.0, stop=duration)
    popI[i].inject(noise)

# Synapses
EE_connections = Projection(popE, popE, FixedProbabilityConnector(p_connect=plocal),
                                    StaticSynapse(weight=wEE_alpha, delay=d))
IE_connections = Projection(popE, popI, FixedProbabilityConnector(p_connect=plocal),
                                    StaticSynapse(weight=wIE_alpha, delay=d))
EI_connections = Projection(popI, popE, FixedProbabilityConnector(p_connect=plocal),
                                    StaticSynapse(weight=wEI_alpha, delay=d))
II_connections = Projection(popI, popI, FixedProbabilityConnector(p_connect=plocal),
                                    StaticSynapse(weight=wII_alpha, delay=d))

# initial conditions
kernelseed  = 5456532
rng = NumpyRNG(kernelseed, parallel_safe=True)
uniformDistr = RandomDistribution('uniform', low=Vrest, high=Vt, rng=rng)
initialize(popE, v=uniformDistr)
initialize(popI, v=uniformDistr)

# Record
popE.record('spikes') 
popI.record('spikes') 
popE.record(['v'])  # record other variables from first two neurons
popI.record(['v'])  # record other variables from first two neurons


# Run
run(duration)
    
# Store spikes
spikesE_in = popE.get_data()
spikesI_in = popI.get_data()

############### Generate array with spike data for popE
spkpopE=spikeData(popE)
spkpopI=spikeData(popI)


################# Firing Rate
frE= computeFiringRate(spkpopE,range(NE),duration,dt)
frWindowE=sliding_window(frE[1,:],10,dt) # window of 10 ms

frI= computeFiringRate(spkpopI,range(NI),duration,dt)
frWindowI=sliding_window(frI[1,:],10,dt) # window of 10 ms


# Plot
plt.subplot(211)
plt.plot(spkpopE[:,0],spkpopE[:,1],'k.',markersize=2)
plt.plot(spkpopI[:,0],NE+spkpopI[:,1],'r.',markersize=2)
plt.ylabel('# Neuron',fontsize=15)
plt.subplot(212)
plt.plot(frE[0,:],frWindowE,'k')
plt.plot(frI[0,:],frWindowI,'r')
plt.ylabel('Firing Rate (Hz)',fontsize=15)
plt.xlabel('Time (ms)',fontsize=15)
plt.show()


