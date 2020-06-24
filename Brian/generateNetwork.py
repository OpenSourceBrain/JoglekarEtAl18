# -*- coding: utf-8 -*-
""" to create figures for spiking network models in joglekar et al neuron 2018
"""
from __future__ import division
from brian2 import *
prefs.codegen.target = 'auto'


import scipy.io
import numpy as np
import random as pyrand
import time
import os
 
rnd_seed = 1
pyrand.seed(324823+rnd_seed)
numpy.random.seed(324823+rnd_seed)


def gen_params(n_areas,regime, gba, duration):    
    
    
    para= {'N'         : 2000,    # Number of neurons 
           'NAreas'    : n_areas, # Number of areas
           'Ne'        : 0.8,     # fraction of excitatory neurons
           'Vr'        : -70.*mV, # Membrane potential rest
           'Vreset'    : -60.*mV, # Membrane potential reset
           'Vt'        : -50.*mV, # Membrane potential threshold
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
           'path'      : os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Matlab/' #path to .mat files
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
    
    
    
    para['k']=400
    
    # Membrane resitance
    R=50*Mohm
    
    para['duration']    = duration*ms;   
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
            para['currval']  = (300*pA*R)/mV 
        
        elif gba=='strong':    
            para['wEI']      = .05
            para['muEE']     = .05    
            para['currval']  = (126*pA*R)/mV 
      
        para['currdur']       = 1500
        
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
        
        para['currdur']      = 80
        para['currval']      = (202*pA*R)/mV

    return para


def equations():
    
    # Equations
    eqsE = Equations('''
    dV/dt=(-(V-Vr) + stimulus(t,i) + Vext )*(1./tau) + 
    (sigma*(1./tau)**0.5)*xi : volt (unless refractory)
    
    Vext : volt    
    tau: second
    sigma : volt
    Vr:volt
    
    ''' )

    eqsI = Equations('''
      dV/dt=(-(V-Vr) + Vext )*(1./tau) + 
      (sigma*(1./tau)**0.5)*xi : volt (unless refractory)
     
    Vext : volt    
    tau: second
    sigma : volt
    Vr:volt  
      
      ''')
    
    return eqsE,eqsI

def setStimulus(para):
    
    # Stimulus
    netsteps = round(para['duration']/defaultclock.dt)
    
    a1 = np.zeros([3000,1]) #input given to v1 for fixed duration. 
    a2 = para['currval']*np.ones([para['currdur'],1])
    a3 = np.zeros([  int(netsteps - 3000 - para['currdur']) , 1])
    aareaone = np.vstack((a1,a2,a3)) 


    timelen = len(aareaone)
    excotherareas = para['k']*4*(para['NAreas']-1)
    aareaonenet = np.tile(aareaone,(1,para['k']*4))
    arest = np.zeros([timelen, excotherareas])
    netarr = np.hstack((aareaonenet,arest))
    
    stimulus = TimedArray(netarr*mV, dt=defaultclock.dt)
    
    return stimulus



def network(para):
    
    # Equations
    eqsE,eqsI=equations()
    
    # Stimulus
    stimulus=setStimulus(para)
    
    
    # Total number of excitatory neurons
    NE=int(para['NAreas']*para['N']*para['Ne'])
    
    # Total number of inhibitory neurons
    NI=int((para['NAreas']*para['N'])-NE)
    
    # Parameters
    paraVt     = para['Vt']
    paraVreset = para['Vreset']
    
    # Neuron groups
    E = NeuronGroup(N=NE, method='euler',model=eqsE, threshold='V > paraVt', reset='V=paraVreset', refractory=para['tref']*ms)
    I = NeuronGroup(N=NI, method='euler',model=eqsI, threshold='V > paraVt', reset='V=paraVreset', refractory=para['tref']*ms)

    #E I across areas
    Exc, Inh = [], []
    Exc = [ E[y*(para['k']*4):(y+1)*(para['k']*4)] for y in range(para['NAreas'])]
    Inh = [ I[z*(para['k']):(z+1)*(para['k'])] for z in range(para['NAreas'])] 
    
    # List to store connections
    Exc_C_loc  = [None]*para['NAreas']
    Inh_C_loc  = [None]*para['NAreas']
    EtoI_C_loc = [None]*para['NAreas']
    ItoE_C_loc = [None]*para['NAreas']

    Exc_C_lr_fromi =[]
    EtoI_C_lr_fromi =[]

    #set up synaptic connections 
    h = 0
    while h < para['NAreas']:
      #print(h)  #local. 
      Exc_C_loc[h] = Synapses(Exc[h], Exc[h], 'w:volt', delay = para["dlocal"]*ms, on_pre='V+=w')  
      Inh_C_loc[h] = Synapses(Inh[h], Inh[h], 'w:volt', delay = para["dlocal"]*ms, on_pre='V+= w ')  
      EtoI_C_loc[h] = Synapses(Exc[h], Inh[h],'w:volt', delay = para["dlocal"]*ms, on_pre='V+= w ')    
      ItoE_C_loc[h] = Synapses(Inh[h], Exc[h],'w:volt', delay = para["dlocal"]*ms, on_pre='V+= w ') 
          
      Exc_C_loc[h].connect(p =  para["probIntra"]) 
      Inh_C_loc[h].connect(p =  para["probIntra"]) 
      EtoI_C_loc[h].connect(p = para["probIntra"]) 
      ItoE_C_loc[h].connect(p = para["probIntra"]) 
      
      Exc_C_loc[h].w = (1+para["alpha"]*para["hier"][h])*para['wEE']*mV
      Inh_C_loc[h].w = -para['wII']*mV
      EtoI_C_loc[h].w = (1+para["alpha"]*para["hier"][h])*para['wIE']*mV
      ItoE_C_loc[h].w = -para['wEI']*mV
      
      j = 0 #long range to j. 
      while j < para['NAreas']:
        if j!= h:  
            #print j
            exc_lr_itoj, etoi_lr_itoj = None, None
    
            exc_lr_itoj = Synapses(Exc[h], Exc[j], 'w:volt', on_pre='V+= w ') 
            etoi_lr_itoj = Synapses(Exc[h], Inh[j], 'w:volt', on_pre='V+= w ')
                    
            exc_lr_itoj.connect(p = para["probInter"])    
            etoi_lr_itoj.connect(p = para["probInter"])  
            
            exc_lr_itoj.w =  (1 + para["alpha"] * para["hier"][j]) * para['muEE'] * para["conn"][j,h]*mV
            etoi_lr_itoj.w = (1 + para["alpha"] * para["hier"][j]) * para['muIE'] * para["conn"][j,h]*mV
            
            # Mean for delay distribution 
            meanlr = para["delayMat"][j,h]
            # Standard deviation for delay distribution
            varlr =  para['lrvar']*meanlr
            
            exc_lr_itoj.delay = np.random.normal(meanlr,varlr,len(exc_lr_itoj.w))*ms
            etoi_lr_itoj.delay = np.random.normal(meanlr,varlr,len(etoi_lr_itoj.w))*ms
            
            Exc_C_lr_fromi.append(exc_lr_itoj)
            EtoI_C_lr_fromi.append(etoi_lr_itoj)
            
        j = j + 1       
      h = h + 1

    
    # Initial conditions
    E.V = para['Vr'] + rand(len(E)) * (para['Vt'] - para['Vr'])
    E.tau   = para['taumE']*ms
    E.Vext  = para['VextE']*mV
    E.sigma = para['sigma']*mV
    E.Vr = para['Vr']
    
    I.V = para['Vr'] + rand(len(I)) * (para['Vt'] - para['Vr'])
    I.tau   = para['taumI']*ms
    I.Vext  = para['VextI']*mV
    I.sigma = para['sigma']*mV
    I.Vr = para['Vr']

    
    # Monitors 
    monitorsE = SpikeMonitor(E)
    monitorsI = SpikeMonitor(I)

    # Network    
    net = Network(E,I,Exc_C_loc,EtoI_C_loc,ItoE_C_loc,Inh_C_loc,Exc_C_lr_fromi,EtoI_C_lr_fromi,monitorsE,monitorsI)

    net.store()
    print("net stored")
    net.run(para['duration'], report='text')
    
    
    return monitorsE, monitorsI


def run_network(n_areas,regime,gba,duration):
    
    # Parameters
    para = gen_params(n_areas,regime,gba,duration) 
    # Run Network - Brian
    monitor_spike, monitor_v = network(para)
    
    return monitor_spike


def firingRate(N,goodprop,badprop,duration):

    binsize = 10*ms
    stepsize =  1*ms  
    
    # Store maximum firing rate for each area
    maxratebad  = np.empty([N,1])
    maxrategood = np.empty([N,1])
    
    # sort net spikes
    netspikebad = len(badprop)
    netspikegood= len(goodprop)
    
    badpropsorted = badprop[badprop[:,1].argsort(),]
    goodpropsorted = goodprop[goodprop[:,1].argsort(),]
            
    netbinno = int( 1+(duration/ms)-(binsize/ms))
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
    
        astep = binsize/(1*ms)
        
        valsnewbad = np.zeros(netbinno)
        valsnewgood = np.zeros(netbinno)
        
        acount = 0
        while acount < netbinno:        
            valsnewbad[acount] = sum(valszerobad[int(acount):int(acount+astep)])
            valsnewgood[acount] = sum(valszerogood[int(acount):int(acount+astep)])
            acount=acount+1
    
        valsratebad = valsnewbad*((1000*ms/binsize) /(1600) ) # divide by no of neurons per E pop. 
        valsrategood = valsnewgood*((1000*ms/binsize) /(1600) )    
        popratebad[u,:], poprategood[u,:] = valsratebad, valsrategood    
        
        #compute population firing rates. 
        
        maxratebad[u,0] = max(valsratebad[int(len(valsratebad)/3):])
        maxrategood[u,0] = max(valsrategood[int(len(valsrategood)/3):])
        
            
    return maxratebad, maxrategood

    
