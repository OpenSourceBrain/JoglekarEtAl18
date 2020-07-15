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
import sys

import matplotlib.pyplot as plt

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

# timestep
dt=0.1
# initialize pyNN simulation
setup(timestep=dt)
    
# duration
duration=500

# number of areas
nAreas=2

# Axonal velocity
speed=3.5

# Alpha
alpha=4

# mu
muIE     = .19/4 
muEE     = .0375

# path to .mat files
path='/home/ronaldo/github/JoglekarEtAl18/Matlab/'
#hierarchy values file 
hierVals = scipy.io.loadmat(path+'hierValspython.mat')
hierValsnew = hierVals['hierVals'][:]
hier=hierValsnew/max(hierValsnew)#hierarchy normalized. 
hier=np.squeeze(hier[:nAreas])

#fln values file 
flnMatp = scipy.io.loadmat(path+'efelenMatpython.mat')
conn=flnMatp['flnMatpython'][:][:] #fln values..Cij is strength from j to i 
conn=conn[:nAreas,:nAreas]

distMatp = scipy.io.loadmat(path+'subgraphWiring29.mat')
distMat=distMatp['wiring'][:][:] #distances between areas values..
delayMat = distMat/speed
delay=delayMat[:nAreas,:nAreas]



###############################################################################
############################## PyNN  ##########################################
###############################################################################


# Parameters for excitatory neurons
E_parameters={
            "tau_m":tauE_m, 
            "cm":cE_m, 
            "v_rest":Vrest, 
            "v_reset":Vreset, 
            "v_thresh":Vt, 
            "tau_refrac":tauRef,
            "tau_syn_E":tau_syn_e,
            "tau_syn_I":tau_syn_e,
            "i_offset":(VextE/R)}
# Paramters for inhibitory neurons
I_parameters={
            "tau_m":tauI_m, 
            "cm":cI_m, 
            "v_rest":Vrest, 
            "v_reset":Vreset, 
            "v_thresh":Vt, 
            "tau_refrac":tauRef,
            "tau_syn_E":tau_syn_e,
            "tau_syn_I":tau_syn_e,
            "i_offset":(VextI/R)}

############ All Excitatory and inhibitory neurons ############################
popE =Population(nAreas*NE, IF_curr_alpha,E_parameters,label="popE")
popI =Population(nAreas*NI, IF_curr_alpha,I_parameters,label="popI")

################################ Noise ########################################

# Noise on excitatory neurons
stdNoiseE = (sigmaV/R)*(tauE_m**0.5)/(dt**0.5)
for i in range(nAreas*NE):
    noise = NoisyCurrentSource(mean=0, stdev=stdNoiseE, start=0.0, stop=duration)
    popE[i].inject(noise)
    popE[i].inject(noise)

# Noise on inhibitory neurons    
stdNoiseI = (sigmaV/R)*(tauI_m**0.5)/(dt**0.5)
for i in range(nAreas*NI):
    noise = NoisyCurrentSource(mean=0, stdev=stdNoiseI, start=0.0, stop=duration)
    popI[i].inject(noise)
    popI[i].inject(noise)

# paramters for initial conditions
kernelseed  = 5456532
rng = NumpyRNG(kernelseed, parallel_safe=True)
uniformDistr = RandomDistribution('uniform', low=Vrest, high=Vt, rng=rng)
initialize(popE, v=uniformDistr)
initialize(popI, v=uniformDistr)

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
for i in range(nAreas):
    
    # store sub populations in lists 
    popEList.append(popE[(i*NE):((i+1)*NE)])
    popIList.append(popI[(i*NI):((i+1)*NI)])


    #### Synapseos
    
    # Weights for recurrent connections
    wEE_alpha = (((1+alpha*hier[i])*wEE)/coeffE)      #[nA]
    wIE_alpha = (((1+alpha*hier[i])*wIE)/coeffI)      #[nA]  
    wEI_alpha = (wEI/coeffE)                          #[nA]
    wII_alpha = (wII/coeffI)                          #[nA]

    # Connections    
    EE_connections = Projection(popEList[i], popEList[i], FixedProbabilityConnector(p_connect=plocal),
                                     StaticSynapse(weight=wEE_alpha, delay=d))
    IE_connections = Projection(popEList[i], popIList[i], FixedProbabilityConnector(p_connect=plocal),
                                     StaticSynapse(weight=wIE_alpha, delay=d))
    EI_connections = Projection(popIList[i], popEList[i], FixedProbabilityConnector(p_connect=plocal),
                                     StaticSynapse(weight=wEI_alpha, delay=d))
    II_connections = Projection(popIList[i], popIList[i], FixedProbabilityConnector(p_connect=plocal),
                                    StaticSynapse(weight=wII_alpha, delay=d))
    
    # Store projections in lists
    EE.append(EE_connections)
    IE.append(IE_connections)
    EI.append(EI_connections)
    II.append(II_connections)
                                

    # Record data 
    popEList[i].record('spikes') 
    popIList[i].record('spikes') 


# Long Range connections
for i in range(nAreas):
    for j in range(nAreas):
        
        if i!=j:
            
            # Weights
            wEE_alphaLR =  ((1 + alpha * hier[j]) * muEE * conn[j,i])/coeffE
            wIE_alphaLR =  ((1 + alpha * hier[j]) * muIE * conn[j,i])/coeffI
            
            # Delay
            dLR=RandomDistribution('normal', [delay[j,i], 0.1*delay[j,i]], rng=NumpyRNG(seed=4242))
            
            # Connections
            EE_connectionsLR = Projection(popEList[i], popEList[j], FixedProbabilityConnector(p_connect=plocal),
                                     StaticSynapse(weight=wEE_alphaLR, delay=dLR))
            IE_connectionsLR = Projection(popEList[i], popIList[j], FixedProbabilityConnector(p_connect=plocal),
                                     StaticSynapse(weight=wIE_alphaLR, delay=dLR))
   
            # Store projections in list
            EElongRange.append(EE_connectionsLR)
            IElongRange.append(IE_connectionsLR)


# Stimulus
amplitudeV =10.1 #[mV]
pulse = DCSource(amplitude=amplitudeV/R, start=300.0, stop=320.0)
pulse.inject_into(popEList[0])

# Run
run(duration)
    
# Store spikes
spikesE_V1_in = popEList[0].get_data()
spikesI_V1_in = popIList[0].get_data()
spikesE_V2_in = popEList[1].get_data()
spikesI_V2_in = popIList[1].get_data()


############### Generate array with spike data for popE
spkpopE_V1=spikeData(popEList[0])
spkpopI_V1=spikeData(popIList[0])

spkpopE_V2=spikeData(popEList[1])
spkpopI_V2=spikeData(popIList[1])


################# Firing Rate
frE_V1= computeFiringRate(spkpopE_V1,range(NE),duration,dt)
frWindowE_V1=sliding_window(frE_V1[1,:],10,dt) # window of 10 ms

frI_V1= computeFiringRate(spkpopI_V1,range(NI),duration,dt)
frWindowI_V1=sliding_window(frI_V1[1,:],10,dt) # window of 10 ms

frE_V2= computeFiringRate(spkpopE_V2,range(NE,2*NE),duration,dt)
frWindowE_V2=sliding_window(frE_V2[1,:],10,dt) # window of 10 ms

frI_V2= computeFiringRate(spkpopI_V2,range(NI,2*NI),duration,dt)
frWindowI_V2=sliding_window(frI_V2[1,:],10,dt) # window of 10 ms

# Plot
plt.subplot(221)
plt.title('V1',fontsize=15)
plt.plot(spkpopE_V1[:,0],spkpopE_V1[:,1],'k.',markersize=2)
plt.plot(spkpopI_V1[:,0],NE+spkpopI_V1[:,1],'r.',markersize=2)
plt.ylabel('# Neuron',fontsize=15)
plt.subplot(223)
plt.plot(frE_V1[0,:],frWindowE_V1,'k')
plt.plot(frI_V1[0,:],frWindowI_V1,'r')
plt.ylabel('Firing Rate (Hz)',fontsize=15)
plt.xlabel('Time (ms)',fontsize=15)

plt.subplot(222)
plt.title('V2',fontsize=15)
plt.plot(spkpopE_V2[:,0],spkpopE_V2[:,1]-NE,'k.',markersize=2)
plt.plot(spkpopI_V2[:,0],NE+(spkpopI_V2[:,1]-NI),'r.',markersize=2)
plt.ylabel('# Neuron',fontsize=15)
plt.subplot(224)
plt.plot(frE_V2[0,:],frWindowE_V2,'k')
plt.plot(frI_V2[0,:],frWindowI_V2,'r')
plt.ylabel('Firing Rate (Hz)',fontsize=15)
plt.xlabel('Time (ms)',fontsize=15)



plt.show()



