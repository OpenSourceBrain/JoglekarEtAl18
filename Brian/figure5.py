# -*- coding: utf-8 -*-
""" to create figures for spiking network models in joglekar et al neuron 2018
"""
from __future__ import division
from brian2 import *
prefs.codegen.target = 'auto'

import matplotlib.pyplot as plt

import scipy.io
import numpy as np
import random as pyrand
import time
 

def gen_params(regime, gba, duration):    
    
    
    para= {'N'         : 2000,    # Number of neurons 
           'NAreas'    : 3,       # Number of areas
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
           'alpha'     : 4,       # gradient 
           'dlocal'    : 2 ,      # delays local 
           'speed'     : 3.5,     # axonal conduction velocity
           'lrvar'     : 0.1      # standard deviation delay long range
           }

    R=50*Mohm
    para['duration']    = duration*ms;   
    if regime=='assynchronous':  
    
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
            para['currval']       = (300*pA*R)/mV
        
        elif gba=='strong':    
            para['wEI']      = .05
            para['muEE']     = .05    
            para['currval']       = (126*pA*R)/mV
      
        para['currdur']       = 1500
        
    elif regime=='synchronous':
        
        # general for assynchronous regime
        para['muIE']     = .19 
        para['wII']      = .3
        para['wEE']      = .04
        para['wIE']      = .3
        
        if gba=='weak':
            para['VextI']    = 14.0  
            para['VextE']    = 15.4
            para['wEI']      = .056
            para['muEE']     = .016
        elif gba=='strong':    
            para['VextI']    = 14.0  
            para['VextE']    = 16.0
            para['wEI']      = .98
            para['muEE']     = .25    
        
        para['currdur']       = 80
        para['currval']       = (200*pA*R)/mV

    return para


def equations():

    eqs =Equations('''
    dV/dt=(-(V-Vr) + Vext )*(1./tau) + stimulus(t,i)*(1./tau)+ 
        (sigma*(1./tau)**0.5)*xi : volt (unless refractory)

    Vext : volt    
    tau: second
    sigma : volt
    Vr:volt
    
    ''' )

    return eqs


def setConnections(para):
    
    sources=[] # to store the sources
    targets=[] # to store the targets
    wExc=[]    # to store the excitatory synaptic weights
    wInh=[]    # to store the inhibitory synaptic weights 
    delays=[]  # to store the delays
    
    # all neurons 
    index=np.arange(para["N"]*para["NAreas"])
    # neurons separated by networks
    nets=np.split(index,para["NAreas"])
    
    # Array with all excitatory neurons
    # Array with all inhibitory neurons
    tempN=np.split(nets,[int(para["Ne"]*para["N"]),para["N"]],axis=1)[0:2]
    neuronsE=tempN[0].flatten()
    neuronsI=tempN[1].flatten()
    
    #hierarchy values file 
    hierVals = scipy.io.loadmat('hierValspython.mat')
    hierValsnew = hierVals['hierVals'][:]
    hier=hierValsnew/max(hierValsnew)#hierarchy normalized. 
    hier=hier[:para["NAreas"]]

    #fln values file 
    flnMatp = scipy.io.loadmat('efelenMatpython.mat')
    conn=flnMatp['flnMatpython'][:][:] #fln values..Cij is strength from j to i 
    conn=conn[:para["NAreas"],:para["NAreas"]]

    distMatp = scipy.io.loadmat('subgraphWiring29.mat')
    distMat=distMatp['wiring'][:][:] #distances between areas values..
    delayMat = distMat/para['speed']
    delayMat=delayMat[:para["NAreas"],:para["NAreas"]]
    
    for areaSource in range(para["NAreas"]):
         
        for neuron in nets[areaSource]:
            
            # Temp for weights
            tempWeights=[] 
            ######################### Recurrent connections ###################           
            netTemp=nets[areaSource][:]
            # to avoid autapses
            netTemp=np.delete(netTemp,np.where(nets[areaSource]==neuron))
            # compute connections
            tempConn=np.random.rand(len(netTemp))<para["probIntra"]
            # save sources
            sources.extend([neuron]*len(netTemp[tempConn]))
            # save targets
            targets.extend(netTemp[tempConn])
            
            ##################### Set Weights #################################
            if neuron<int(para["Ne"]*(nets[0][-1]+1)+nets[areaSource][0]): # Excitatory
                
                # array to store excitatory weights for local connections
                wExcTemp=np.ones(len(netTemp[tempConn]))
                # E->E
                _,idxE,_=np.intersect1d(netTemp[tempConn],neuronsE,return_indices=True)
                wExcTemp[idxE]=np.ones(len(wExcTemp[idxE]))*(1+para["alpha"]*hier[areaSource])*para["wEE"]
                # E->I 
                _,idxI,_=np.intersect1d(netTemp[tempConn],neuronsI,return_indices=True)
                wExcTemp[idxI]=np.ones(len(wExcTemp[idxI]))*(1+para["alpha"]*hier[areaSource])*para["wIE"]
                # local excitatory weights
                wExc.extend(wExcTemp)
                # set at zero the inhibitory connections between these neurons     
                wInh.extend([0]*len(netTemp[tempConn]))
                # delays for excitatory local connections
                delays.extend(np.ones(len(netTemp[tempConn]))*para["dlocal"])
                ############################### Long Range ###################
                # delete area source from array of areas
                netsTarget=np.delete(nets,areaSource,0) 
                # connections 
                temp=np.random.rand(np.shape(netsTarget)[0],np.shape(netsTarget)[1])<para["probInter"]
                # target neurons
                targets.extend(netsTarget[temp])
                # source neuron
                sources.extend([neuron]*len(netsTarget[temp]))
                # number of target neurons
                nNeuronsTarget=np.sum(temp*1,1)
                # Weight for long range connections
                Weights=conn[np.delete(np.arange(para["NAreas"]),areaSource,0),areaSource] #### REVER ESSA ESTRUTURA talvez [,] seja melhor
                # Store weights for long range connections
                [tempWeights.extend((1+para["alpha"]*hier[i])*Weights[i].repeat(nNeuronsTarget[i])) for i in range(nNeuronsTarget.size)][0]
                # Convert list to array
                tempWeights=np.asarray(tempWeights)
                ############## Set Weights according to target ################
                # Return idx for connection where the target is excitatory
                _,idxE,_=np.intersect1d(netsTarget[temp],neuronsE,return_indices=True)
                # Return idx for connection where the target is inhibitory
                _,idxI,_=np.intersect1d(netsTarget[temp],neuronsI,return_indices=True)
                # Adjust long range connections EE
                tempWeights[idxE]=tempWeights[idxE]*para["muEE"]
                # Adjust long range connections EI
                tempWeights[idxI]=tempWeights[idxI]*para["muIE"]
                # Store weights
                wExc.extend(tempWeights)
                # Set at zero weights for inhibitory connections
                wInh.extend([0]*np.sum(nNeuronsTarget))
                # delays for long range connections
                [delays.extend(np.random.normal(delayMat[i,areaSource],para['lrvar'] *delayMat[i,areaSource],nNeuronsTarget[i])) for i in range(nNeuronsTarget.size)][0]
                ################################################################
                    
            else:
                 # array to store inhibitory weights for local connections
                wInhTemp=np.ones(len(netTemp[tempConn]))
                # I->E
                _,idxE,_=np.intersect1d(netTemp[tempConn],neuronsE,return_indices=True)
                wInhTemp[idxE]=np.ones(len(wInhTemp[idxE]))*para["wEI"]*-1
                # I->I 
                _,idxI,_=np.intersect1d(netTemp[tempConn],neuronsI,return_indices=True)
                wInhTemp[idxI]=np.ones(len(wInhTemp[idxI]))*para["wII"]*-1
                # set at zero the excitatory connections between these neurons
                wExc.extend([0]*len(netTemp[tempConn]))
                # local inhibitory weights
                wInh.extend(wInhTemp) 
                # delays for inhibitory local connections
                #delays.extend(np.random.randn(len(netTemp[tempConn]))*par["delayStdRR"]+par["delayMeanRR"])
                delays.extend(np.ones(len(netTemp[tempConn]))*para["dlocal"])
                ###################################################################        

    return sources,targets,wExc,wInh,delays


def setVextTau(para):
    
    tau=[]
    Vext=[]
    
    for i in range(para['NAreas']):
        tau.extend(np.ones(int(para['N']*para['Ne']))*para['taumE']) 
        tau.extend(np.ones(int(para['N']-para['N']*para['Ne']))*para['taumI'])
        
        Vext.extend(np.ones(int(para['N']*para['Ne']))*para['VextE']) 
        Vext.extend(np.ones(int(para['N']-para['N']*para['Ne']))*para['VextI'])
 
    return tau,Vext


def createStimulus(para):
    
    netsteps = round(para['duration']/defaultclock.dt)
    
    a1 = np.zeros([3000,1]) #input given to v1 for fixed duration. 
    a2 = para['currval']*np.ones([para['currdur'],1])
    a3 = np.zeros([  int(netsteps - 3000 - para['currdur']) , 1])
    aareaone = np.vstack((a1,a2,a3)) 

    timelen = len(aareaone)
    aareaonenet = np.tile(aareaone,(1,int(para['N']*para['Ne'])))
    
    
    arest = np.zeros([timelen, int((para['NAreas']*para['N'])-para['N']*para['Ne'])])
    netarr = np.hstack((aareaonenet,arest))
    
    stim = TimedArray(netarr*mV, dt=defaultclock.dt)

    return stim

def setMonitors(monitors,para):

    # all neurons 
    index=np.arange(para["N"]*para["NAreas"])
    # neurons separated by networks
    nets=np.split(index,para["NAreas"])
    # Array with all excitatory neurons
    # Array with all inhibitory neurons
    tempN=np.split(nets,[int(para["Ne"]*para["N"]),para["N"]],axis=1)[0:2]
    neuronsE=tempN[0].flatten()
    neuronsI=tempN[1].flatten()
    
    _,idxE,_=np.intersect1d(monitors.i,neuronsE,return_indices=True)
    _,idxI,_=np.intersect1d(monitors.i,neuronsI,return_indices=True)
    
    # Arrays to store neuron index and spiking time
    monE=np.zeros((2,len(idxE)))
    monI=np.zeros((2,len(idxI)))
    
    # neuron index              
    monE[0,:] = monitors.i[idxE]           
    monI[0,:] = monitors.i[idxI]              
    # spinking time              
    monE[1,:] = monitors.t[idxE]           
    monI[1,:] = monitors.t[idxI]              
    
    return monE,monI



def plotUtils(monE, para):
    # Set the index of neurons for plotting
    
     # all neurons 
    index=np.arange(para["N"]*para["NAreas"])
    # neurons separated by networks
    nets=np.split(index,para["NAreas"])
    # number of inhibitory neurons
    Ni=para["N"]-(para["N"]*para["Ne"])
    
    for i in range(1,para["NAreas"]):
        _,idxE,_=np.intersect1d(monE[0,:],nets[i],return_indices=True)
        monE[0,idxE]=monE[0,idxE]-((i*Ni)-1)
    
    return monE


def network(regime,gba,duration):
    
    para = gen_params(regime,gba,duration) 
    eqs = equations()
    sources,targets,wExc,wInh,delays=setConnections(para)
    taum,VextInput =setVextTau(para)
    stimulus=createStimulus(para)
    
    Ntotal=int(para['N']*para['NAreas'])
    paraVt     = para['Vt']
    paraVreset = para['Vreset']
    P = NeuronGroup(Ntotal, method='euler', model=eqs, threshold='V > paraVt', reset='V=paraVreset', refractory=para['tref']*ms)
    
    # Excitatory Synapses (RR and LR)
    CE = Synapses(P, P, model='w : volt',on_pre='V+=w')
    CE.connect(i=sources, j=targets)
    CE.w=wExc*mV
    CE.delay=delays*ms
    
    # Inhibitory Synapses (RR)
    CI = Synapses(P, P, model='w : volt',on_pre='V+=w')
    CI.connect(i=sources, j=targets)
    CI.w=wInh*mV
    CI.delay=delays*ms
    
    # Initial conditions
    P.V     = para['Vr'] + rand(len(P)) * (para['Vt'] - para['Vr'])
    P.tau   = taum*ms
    P.Vext  = VextInput*mV
    P.sigma = para['sigma']*mV
    P.Vr = para['Vr']
    
    #monitor system behavior -- spikes and state variables. 
    monitors = SpikeMonitor(P)
    monitorstatev = StateMonitor(P,'V',record=True)
    
    run(para['duration'], report='text')

    mE,mI=setMonitors(monitors,para)
    mE=plotUtils(mE, para)
    
    return mE,mI,stimulus

mE, _,stimulus = network('assynchronous','weak',1000)
mE2, _,_ = network('assynchronous','strong',1000)


################################### Plot #####################################
plt.figure()
subplot(131)
plt.plot(mE[1,:]/ms, mE[0,:], '.',markersize=1)
plt.plot([0, max(mE[1,:]/ms)], np.arange(3+1).repeat(2).reshape(-1, 2).T*1600, 'k-')

subplot(132)
plt.plot(mE2[1,:]/ms, mE2[0,:], '.',markersize=1)
plt.plot([0, max(mE[1,:]/ms)], np.arange(3+1).repeat(2).reshape(-1, 2).T*1600, 'k-')
