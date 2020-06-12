# -*- coding: utf-8 -*-
""" to create figures for spiking network models in joglekar et al neuron 2018
"""
from __future__ import division
from brian2 import *
prefs.codegen.target = 'auto'

import matplotlib.pyplot as plt

import scipy.io
import numpy as np
import numpy.random
import random as pyrand

#hierarchy values file 
hierVals = scipy.io.loadmat('/Users/maddy/anaconda/hierValspython.mat')
hierValsnew = hierVals['hierVals'][:]
netwParams_hier=hierValsnew/max(hierValsnew)#hierarchy normalized. 

#fln values file 
flnMatp = scipy.io.loadmat('/Users/maddy/anaconda/efelenMatpython.mat')
flnMat=flnMatp['flnMatpython'][:][:] #fln values..Cij is strength from j to i 

distMatp = scipy.io.loadmat('/Users/maddy/anaconda/subgraphWiring29.mat')
distMat=distMatp['wiring'][:][:] #distances between areas values..


from brian2 import defaultclock
 
rnd_seed = 1
pyrand.seed(324823+rnd_seed)
numpy.random.seed(324823+rnd_seed)

#generate parameters for networks
def gen_params(extra_params=dict()):    
    # SINGLE NEURON PARAMETERS
    para= {'Vr' : -70.*mV, 'Vreset' : -60.*mV, 'Vt' : -50.*mV, 
    'taum' : 20. * ms, 'tref' : 2.*ms, 
    'taumI': 10. * ms,
    'k' : 400, 
    'p' : .1, 'pintarea': .1, #connection density local and long range
    'N_area' : 29, 'isFB' : True,
    'pE': -.0 , 'pI': .0, #interareal E to E is pintarea + pE but E to I is pintarea + pI
    'sigmaval' : 3.*mV, #noise 
    'muI'    :  14.1*mV, 'muE' : 14.1*mV , 
    'omegaEEsp' :  0.04*mV,#local strengths. 
    'omegaIEsp' : .3*mV, 
    'omegaEIsp' : 0.24*mV,  
    'omegaIIsp' : .3*mV,
    'muEEsp'    : 0.06*mV, 'muIEsp'    : 0.1*mV,#long range strengths 
    'alpha'     : .68,
    'dee' : 0.*ms, 'die': 0.*ms, 'dei' : 0.*ms, 'dii': 0.*ms,'dlocal':0.,#delays local modified later
    'speed': 3.5, 
     }

    # PARAMETERS INTRODUCED FROM EXTERNAL CODE
    for key in extra_params:
        para[key] = extra_params[key] # overwrite the old value
        
    return para

para = gen_params() 

binsize,stepsize = 5*ms,1*ms #this is not used - but take one ms steps.   

#local delay, conduction speed.
para['dlocal'],para['speed'] = 2.,3.5

#current strengths and durations. for asynchronous and synchronous propagation cases

#for synchronous regime. 
#synchronous regime weak gba -- set currval = 0 for background
para['muE'], para['muI'],para['alpha'],duration,currdur,currval = 15.4*mV, 14.*mV,(4./.68)*.68, 440*ms,80,10.1 
para['muEEsp'],para['omegaEIsp'],para['muIEsp'] = 0.16*mV, 0.56*mV, .19*mV 

#synchronous regime strong gba -- set currval = 0 for background
#para['muE'], para['muI'],para['alpha'],duration,currdur,currval = 16.*mV, 14.*mV,(4./.68)*.68, 440*ms,80,10.1
#para['muEEsp'],para['omegaEIsp'],para['muIEsp'] = 0.25*mV, 0.98*mV, .19*mV 


#for asynchronous regime. 
#para['muIEsp'],para['omegaIIsp'],para['omegaEEsp'],para['omegaIEsp'] =.19/4*mV, .075*mV,.01*mV, .075*mV
#currdur,currval = 1500, 6.3

#asynchronous regime weak gba -- set currval = 0 for background
#para['omegaEIsp'], para['muEEsp'] = .0375*mV, .0375*mV 
#para['muI'],para['muE'] = 14.7*mV, 14.2*mV 

#asynchronous regime strong gba -- set currval = 0 for background
#para['omegaEIsp'], para['muEEsp'] = .05*mV, .05*mV 
#para['muI'],para['muE'] = 14.7*mV, 14.2*mV 

##if no feedback 
##flnMat = np.tril(flnMat)

delaylrgaus, lrvar = True, .1

netsteps = round(duration/defaultclock.dt)

arealen = para['N_area']#no of areas. 

#"""
a1 = np.zeros([3000,1]) #input given to v1 for fixed duration. 
a2 = currval*np.ones([currdur,1])
a3 = np.zeros([  int(netsteps - 3000 - currdur) , 1])
aareaone = np.vstack((a1,a2,a3)) 

""" #try input to area 2
currvalarea2 = 15.8
a2area2 = currvalarea2*np.ones([currdur,1])
aareaarea2 = np.vstack((a1, a2area2, a3))
timelen = len(aareaone)
excotherareasone, excotherareastwo = para['k']*4*8, para['k']*4*19
aareaonenet, aareaarea2net = np.tile(aareaone,(1,para['k']*4)), np.tile(aareaarea2,(1,para['k']*4))
arestone, aresttwo = np.zeros([timelen, excotherareasone]), np.zeros([timelen, excotherareastwo])
netarr = np.hstack((aareaonenet,arestone,aareaarea2net,aresttwo)) #new for fig. 
"""

timelen = len(aareaone)
excotherareas = para['k']*4*(arealen-1)
aareaonenet = np.tile(aareaone,(1,para['k']*4))
arest = np.zeros([timelen, excotherareas])
netarr = np.hstack((aareaonenet,arest))

inputtoE1 = TimedArray(netarr*mV, dt=defaultclock.dt)
Inpcur = inputtoE1

#put in parameters
paraVr, paraVt, paraVreset, paramuE, paramuI, parataum, parataumI, parasigmaval = para['Vr'], para['Vt'], para['Vreset'], para['muE'], para['muI'], para['taum'], para['taumI'], para['sigmaval']
paraalpha, paraomegaEEsp, paraomegaEIsp, paraomegaIEsp, paraomegaIIsp = para['alpha'], para['omegaEEsp'], para['omegaEIsp'], para['omegaIEsp'], para['omegaIIsp']
plocal, plongr = para['p'], para['pintarea']
paramuEEsp, paramuIEsp = para['muEEsp'], para['muIEsp']
dlocal = para['dlocal']

#system eqn
eqs = Equations('''
dV/dt=(-(V-paraVr) + inputtoE1(t,i) + paramuE )*(1./parataum) + (parasigmaval*(1./parataum)**0.5)*xi : volt (unless refractory)
''' )

eqsI = Equations('''
dV/dt=(-(V-paraVr) + paramuI )*(1./parataumI) + (parasigmaval*(1./parataumI)**0.5)*xi : volt (unless refractory)
''')

#E I populations
E = NeuronGroup(N=para['k']*4*arealen, method='euler', model=eqs, threshold='V > paraVt', reset='V=paraVreset', refractory=para['tref'])
I = NeuronGroup(N=para['k']*arealen, method='euler',model=eqsI, threshold='V > paraVt', reset='V=paraVreset', refractory=para['tref'])

#E I across areas
Exc, Inh = [], []
Exc = [ E[y*(para['k']*4):(y+1)*(para['k']*4)] for y in range(arealen)]
Inh = [ I[z*(para['k']):(z+1)*(para['k'])] for z in range(arealen)] 

delayMat = distMat/para['speed']

Exc_C_loc, Inh_C_loc, EtoI_C_loc, ItoE_C_loc = [None]*arealen, [None]*arealen, [None]*arealen, [None]*arealen 

Exc_C_lr_fromi, EtoI_C_lr_fromi =[], []

#set up synaptic connections 
h = 0
while h < arealen:
  print h  #local. 
  Exc_C_loc[h] = Synapses(Exc[h], Exc[h], 'w:volt', delay = dlocal*ms, on_pre='V+=w')  
  Inh_C_loc[h] = Synapses(Inh[h], Inh[h], 'w:volt', delay = dlocal*ms, on_pre='V+= w ')  
  EtoI_C_loc[h] = Synapses(Exc[h], Inh[h], 'w:volt', delay = dlocal*ms, on_pre='V+= w ')    
  ItoE_C_loc[h] = Synapses(Inh[h], Exc[h], 'w:volt', delay = dlocal*ms, on_pre='V+= w ') 
      
  Exc_C_loc[h].connect(p = plocal) 
  Inh_C_loc[h].connect(p = plocal) 
  EtoI_C_loc[h].connect(p = plocal) 
  ItoE_C_loc[h].connect(p = plocal) 
  
  Exc_C_loc[h].w = (1+paraalpha*netwParams_hier[h])*paraomegaEEsp
  Inh_C_loc[h].w = -paraomegaIIsp
  EtoI_C_loc[h].w = (1+paraalpha*netwParams_hier[h])*paraomegaIEsp
  ItoE_C_loc[h].w = -paraomegaEIsp
  
  j = 0 #long range to j. 
  while j < arealen:
    if j!= h:  
        print j
        exc_lr_itoj, etoi_lr_itoj = None, None

        exc_lr_itoj = Synapses(Exc[h], Exc[j], 'w:volt', on_pre='V+= w ') 
        etoi_lr_itoj = Synapses(Exc[h], Inh[j], 'w:volt', on_pre='V+= w ')
                
        exc_lr_itoj.connect(p = plongr)    
        etoi_lr_itoj.connect(p = plongr)  
        
        exc_lr_itoj.w =  (1 + paraalpha * netwParams_hier[j]) * paramuEEsp * flnMat[j,h]
        etoi_lr_itoj.w = (1 + paraalpha * netwParams_hier[j]) * paramuIEsp * flnMat[j,h]
        
        meanlr, varlr = delayMat[j,h], lrvar*delayMat[j,h]
        exc_lr_itoj.delay = np.random.normal(meanlr,varlr,len(exc_lr_itoj.w))*ms
        etoi_lr_itoj.delay = np.random.normal(meanlr,varlr,len(etoi_lr_itoj.w))*ms
        
        Exc_C_lr_fromi.append(exc_lr_itoj)
        EtoI_C_lr_fromi.append(etoi_lr_itoj)
        
    j = j + 1       
  h = h + 1

#monitor system behavior -- spikes and state variables. 
monitors = SpikeMonitor(E)
monitorsI = SpikeMonitor(I)
monitorstatev = [StateMonitor(pp,'V',record=True) for pp in Exc]
monitorstatevI = [StateMonitor(ppp,'V',record=True) for ppp in Inh]

E.V = para['Vr'] + rand(len(E)) * (para['Vt'] - para['Vr'])
I.V = para['Vr'] + rand(len(I)) * (para['Vt'] - para['Vr'])

print "before net created"

net = Network(E,I,Exc_C_loc,EtoI_C_loc,ItoE_C_loc,Inh_C_loc,Exc_C_lr_fromi,EtoI_C_lr_fromi,monitors,monitorsI)#,monitorstatev,monitorstatevI

print "net created"
net.store()
print "net stored"
net.run(duration, report='text')


#monitor population firing rates. 

maxrate = np.empty([arealen,1])
meanrate = np.empty([arealen,1])   

netspike = len(monitors.i)
allspike = np.empty([netspike,2])
#monitors save spikes -- neuron numbers and spike times. 
allspike[:,0]=monitors.t/ms
allspike[:,1]=monitors.i

#sort spikes
allspikesorted = allspike[allspike[:,1].argsort(),]

netbinno = int( 1+(duration/ms)-(binsize/ms))
poprate = np.empty([arealen,netbinno ])

u = 0 #for areas. 
count = 0#for each spike. 

monareaktimeall = []
while u<arealen:
    monareaktime = []
    while((count < netspike) and (allspikesorted[count,1]<para['k']*4*(u+1)) ):
      monareaktime.append(allspikesorted[count,0])#append spike times. for each area.
      count = count + 1
    vals = []
    vals = numpy.histogram(monareaktime, bins=int(duration/stepsize))
    valszero = vals[0]

    astep = binsize/(1*ms)
    valsnew = np.zeros(netbinno )
    acount = 0
    while acount < netbinno:        
        valsnew[acount] = sum(valszero[acount:acount+astep]) 
        acount=acount+1

    valsrate = valsnew*((1000*ms/binsize) /(4.*para['k']) ) # divide by no of neurons per E pop. 
    poprate[u,:] = valsrate    
   
    #mean and peak firing rates
    maxrate[u,0] = max(valsrate[int(len(valsrate)/3):])
    meanrate[u,0] = np.mean(valsrate[int(len(valsrate)/3):int(3*len(valsrate)/5)]) 

#    print u
    monareaktimeall.append(monareaktime)
    u = u+1
    
    
#for I population firing rates (used only when deciding background rate )
meanrateI = np.empty([arealen,1])   
poprate_I = np.empty([arealen,netbinno ])
netspikeI = len(monitorsI.i)
allspikeI = np.empty([netspikeI,2])
allspikeI[:,0]=monitorsI.t/ms
allspikeI[:,1]=monitorsI.i
allspikesortedI = allspikeI[allspikeI[:,1].argsort(),]
u = 0 #for areas. 
count = 0#for each spike. 
monareaktimeall_I = []
while u<arealen:
    monareaktime_I = []
    while((count < netspikeI) and (allspikesortedI[count,1]<para['k']*(u+1)) ):
      monareaktime_I.append(allspikesortedI[count,0])#append spike times. for each area.
      count = count + 1
    vals = []
    vals = numpy.histogram(monareaktime_I, bins=int(duration/stepsize))
    
    valszero = vals[0]
    valsnew = np.zeros(netbinno )
    acount = 0
    while acount < netbinno:        
        valsnew[acount] = sum(valszero[acount:acount+astep]) 
        acount=acount+1
        
    valsrate = valsnew*((1000*ms/binsize) /(1.*para['k']) ) 
    poprate_I[u,:] = valsrate
    meanrateI[u,0] = np.mean(valsrate[int(len(valsnew)/3):int(3*len(valsnew)/5)])
#    meanrateI[u,0] = np.mean(valsrate[int(len(valsnew)/3):])
#    meanrateI[u,0] = np.mean(vals[0][int(len(vals[0])/3):])
    print u
    monareaktimeall_I.append(monareaktime_I)
    u = u+1
    
"""
plt.figure()
plt.plot(range(arealen), maxrate)
plt.xlabel('Area')
plt.ylabel('max rate for area')
plt.show()
"""
#print maxrate
#print "mean E:", meanrate, mean(meanrate)
#print "mean I:", meanrateI, mean(meanrateI)
#print mean(meanrate), mean(meanrateI)

plt.figure()
plt.plot(monitors.t/ms, 1.0*monitors.i/(4*para['k']), '.',markersize=1)
plt.plot([0, duration/ms], np.arange(arealen+1).repeat(2).reshape(-1, 2).T, 'k-')
plt.ylabel('Area')
plt.yticks(np.arange(arealen))
plt.xlabel('time (ms)')
ylim(0,arealen)
yticks([0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5,10.5,11.5,12.5,13.5,14.5,15.5,16.5,17.5,18.5,19.5,20.5,21.5,22.5,23.5,24.5,25.5,26.5,27.5,28.5],['V1','V2','V4','DP','MT','8m','5','8l','TEO','2','F1','STPc','7A','46d','10','9/46v','9/46d','F5','TEpd','PBr','7m','7B','F2','STPi','PROm','F7','8B','STPr','24c'])
xlim(300,360)
plt.title('muI/E/alph, wIE, wII, wEE, muIE, wEI, muEE, cur, val = %2.1f,%2.1f,%1.f, %1.3f, %1.3f,%1.2f, %1.3f, %1.2f, %1.2f, %2.f,%2.1f' %(1000*para['muI'],1000*para['muE'],para['alpha'], 1000*para['omegaIEsp'],1000*para['omegaIIsp'],1000*para['omegaEEsp'], 1000*para['muIEsp'], 1000*para['omegaEIsp'], 1000*para['muEEsp'], currdur/10,currval ) )
plt.show()

"""
allspikebad = allspike #weak gba 
#np.save('/Users/maddy/Dropbox/rishicode_maddy070115/allspikebad',allspikebad)
#np.save('/Users/maddy/Dropbox/rishicode_maddy070115/meanratebad',meanrate)
#np.save('/Users/maddy/Dropbox/rishicode_maddy070115/meanrateIbad',meanrateI)
#np.save('/Users/maddy/Dropbox/rishicode_maddy070115/maxratebad',maxrate)
badprop=load('/Users/maddy/Dropbox/rishicode_maddy070115/allspikebad.npy')

allspikegood = allspike #strong gba 
np.save('/Users/maddy/Dropbox/rishicode_maddy070115/allspikegood',allspikegood)
#np.save('/Users/maddy/Dropbox/rishicode_maddy070115/meanrategood',meanrate)
#np.save('/Users/maddy/Dropbox/rishicode_maddy070115/meanrateIgood',meanrateI)
#np.save('/Users/maddy/Dropbox/rishicode_maddy070115/maxrategood',maxrate)
goodprop=load('/Users/maddy/Dropbox/rishicode_maddy070115/allspikegood.npy')
"""
