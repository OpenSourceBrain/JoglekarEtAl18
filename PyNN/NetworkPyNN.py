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

# Testing with neuroml
from pyNN.neuron import *


def plot_spiketrains(segment):
    for spiketrain in segment.spiketrains:
        y = np.ones_like(spiketrain) * spiketrain.annotations['source_id']
        plt.plot(spiketrain, y, 'b.')
        plt.ylabel(segment.name)
        plt.setp(plt.gca().get_xticklabels(), visible=False)

def plot_signal(signal, index, colour='b'):
    label = "Neuron %d" % signal.annotations['source_ids'][index]
    plt.plot(signal.times, signal[:, index], colour, label=label)
    plt.ylabel("%s (%s)" % (signal.name, signal.units._dimensionality.string))
    plt.setp(plt.gca().get_xticklabels(), visible=False)
    plt.legend()
    
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
tau_syn_e     = 0.1          # time constant for synapse [ms] 

# External input for assynchronous behavior
VextE         = 14.2         # External input in excitatory neurons [mV]
VextI         = 14.7         # External input in excitatory neurons [mV]

# Noise
sigmaV=3.0   # [mV]

# Coefficient
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

# initialize pyNN simulation
setup(timestep=0.1)
    
# duration
duration=500

# Excitatory popuation
popE =Population(NE, IF_curr_alpha,IF_curr_alpha.default_parameters,label="popE")

# Excitatory popuation
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
stdNoiseE = (sigmaV/R)*(tauE_m**0.5)   
for i in range(NE):
    noise = NoisyCurrentSource(mean=0, stdev=stdNoiseE, start=0.0, stop=duration, dt=0.1)
    popE[i].inject(noise)

# Noise on inhibitory neurons    
stdNoiseI = (sigmaV/R)*(tauI_m**0.5)
for i in range(NI):
    noise = NoisyCurrentSource(mean=0, stdev=stdNoiseI, start=0.0, stop=duration, dt=0.1)
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
kernelseed  = 1    
rng = NumpyRNG(kernelseed, parallel_safe=True)
uniformDistr = RandomDistribution('uniform', low=Vrest, high=Vt, rng=rng)
initialize(popE, v=uniformDistr)
initialize(popI, v=uniformDistr)

# Record
popE.record('spikes') 
popI.record('spikes') 

# Run
run(duration)
    
# Store spikes
spikesE_in = popE.get_data()
spikesI_in = popI.get_data()

import matplotlib.pyplot as plt
fig_settings = {
    'lines.linewidth': 0.5,
    'axes.linewidth': 0.5,
    'axes.labelsize': 'small',
    'legend.fontsize': 'small',
    'font.size': 8
}

plt.rcParams.update(fig_settings)
plt.figure(1, figsize=(6, 8))

plot_spiketrains(spikesE_in.segments[0])
plot_spiketrains(spikesI_in.segments[0])


# for array in spikesE_in.segments[0].analogsignals:
#     for i in range(array.shape[1]):
#         plt.subplot(10,1,i+1)
#         plot_signal(array, i)

# plt.xlabel("time (%s)" % array.times.units._dimensionality.string)
# plt.setp(plt.gca().get_xticklabels(), visible=True)


plt.show()




##### Criar uma funãço para isso
# spike_data = {}
# spike_data['senders'] = []
# spike_data['times'] = []
# index_offset = 1

# for pop in [popE , popI]:
#     spikes =  pop.get_data('spikes', gather=False)
#     #print(spikes.segments[0].all_data)
#     num_rec = len(spikes.segments[0].spiketrains)
#     print("Extracting spike info (%i) for %i cells in %s"%(num_rec,pop.size,pop.label))
#     #assert(num_rec==len(spikes.segments[0].spiketrains))
#     for i in range(num_rec):
#         ss = spikes.segments[0].spiketrains[i]
#         for s in ss:
#             index = i+index_offset
#             #print("Adding spike at %s in %s[%i] (cell %i)"%(s,pop.label,i,index))
#             spike_data['senders'].append(index)
#             spike_data['times'].append(s)
#     index_offset+=pop.size
