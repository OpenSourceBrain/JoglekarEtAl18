# -*- coding: utf-8 -*-
"""
consciousness figure - save data
"""
from __future__ import division
from brian2 import *
##from brian2.only import *
prefs.codegen.target = 'auto'

import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import numpy.random
import random as pyrand

from brian2 import defaultclock


def testcons(seedval, curr):

  netwParams_hier = np.load('./netwParams_hiervals.npy')
  distMat = np.load('./distMatval.npy')
  flnMat = np.load('./flnMatval.npy')

  def gen_params(extra_params=dict()):    

      para= {'Vr' : -70.*mV, 'Vreset' : -60.*mV, 'Vt' : -50.*mV, 
      'taum' : 20. * ms, 'tref' : 2.*ms, 
      'taumI': 10. * ms,
      'k' : 400, 
      'p' : .1, 'pintarea': .1,
      'N_area' : 29, 'isFB' : True, 
      'sigmaval' : 3.*mV, 
       }

      # PARAMETERS INTRODUCED/MODIFYED FROM EXTERNAL CODE
      for key in extra_params:
          para[key] = extra_params[key] # overwrite the old value

      return para

  para = gen_params() #adjust below since orig values were calculated with k 100

  binsize = 10*ms 
  para['dlocal'],para['speed'] = 2.,3.5

  duration = 550*ms 

  para['alpha'] = 4.  

  para['muIEsp'],para['omegaIIsp'],para['omegaEEsp'],para['omegaIEsp'] =.19/4*mV, .075*mV,.01*mV, .075*mV

  para['omegaEIsp'], para['muEEsp'] = .05*mV, .05*mV #good. async
  para['muI'],para['muE'] = 14.7*mV, 14.2*mV

  #used
  #rnd =1 : case1 rnd=2: case2 ... rnd=10:case3 ... rnd=5:case4 ..rnd=18:case5
  #currlist = [0.0] + np.arange(3.5,6.1,1./6).tolist()
  
  rnd_seed = seedval
  pyrand.seed(324823+rnd_seed)
  seed(324823 + rnd_seed)
  numpy.random.seed(324823+rnd_seed)

  currval,currdur = curr, 1500 
  #flnMat = np.load('./flnMatshuf2.npy')
  #flnMat = np.tril(flnMat)

  netsteps = round(duration/defaultclock.dt)
  arealen = para['N_area']

  a1 = np.zeros([3000,1]) #try changing to see if signal propagates to higher area later
  a2 = currval*np.ones([currdur,1])
  a3 = np.zeros([  int(netsteps - 3000 - currdur) , 1])
  aareaone = np.vstack((a1,a2,a3)) #this changed. 

  #"""
  #give input to v1
  timelen = len(aareaone)
  excotherareas = para['k']*4*(arealen-1)
  aareaonenet = np.tile(aareaone,(1,para['k']*4))
  arest = np.zeros([timelen, excotherareas])
  netarr = np.hstack((aareaonenet,arest))
  #"""
  inputtoE1 = TimedArray(netarr*mV, dt=defaultclock.dt)
  Inpcur = inputtoE1

  paraVr, paraVt, paraVreset, paramuE, paramuI, parataum, parataumI, parasigmaval = para['Vr'], para['Vt'], para['Vreset'], para['muE'], para['muI'], para['taum'], para['taumI'], para['sigmaval']
  paraalpha, paraomegaEEsp, paraomegaEIsp, paraomegaIEsp, paraomegaIIsp = para['alpha'], para['omegaEEsp'], para['omegaEIsp'], para['omegaIEsp'], para['omegaIIsp']
  plocal, plongr = para['p'], para['pintarea']
  paramuEEsp, paramuIEsp = para['muEEsp'], para['muIEsp']
  dlocal = para['dlocal']


  eqs = Equations('''
  dV/dt=(-(V-paraVr) + inputtoE1(t,i) + paramuE )*(1./parataum) + (parasigmaval*(1./parataum)**0.5)*xi : volt (unless refractory)
  ''' )

  eqsI = Equations('''
  dV/dt=(-(V-paraVr) + paramuI )*(1./parataumI) + (parasigmaval*(1./parataumI)**0.5)*xi : volt (unless refractory)
  ''')

  E = NeuronGroup(N=para['k']*4*arealen, method='euler', model=eqs, threshold='V > paraVt', reset='V=paraVreset', refractory=para['tref'])
  I = NeuronGroup(N=para['k']*arealen, method='euler',model=eqsI, threshold='V > paraVt', reset='V=paraVreset', refractory=para['tref'])

  Exc, Inh = [], []
  Exc = [ E[y*(para['k']*4):(y+1)*(para['k']*4)] for y in range(arealen)]
  Inh = [ I[z*(para['k']):(z+1)*(para['k'])] for z in range(arealen)] 

  delayMat = distMat/para['speed']

  Exc_C_loc, Inh_C_loc, EtoI_C_loc, ItoE_C_loc = [None]*arealen, [None]*arealen, [None]*arealen, [None]*arealen 

  Exc_C_lr_fromi, EtoI_C_lr_fromi =[], []


  h = 0
  while h < arealen:
    #print h   #local. 

    Exc_C_loc[h] = Synapses(Exc[h], Exc[h], 'w:volt', delay = dlocal*ms, on_pre='V+=w')  
    Inh_C_loc[h] = Synapses(Inh[h], Inh[h], 'w:volt', delay = dlocal*ms, on_pre='V+= w ')  
    EtoI_C_loc[h] = Synapses(Exc[h], Inh[h], 'w:volt', delay = dlocal*ms, on_pre='V+= w ')    
    ItoE_C_loc[h] = Synapses(Inh[h], Exc[h], 'w:volt', delay = dlocal*ms, on_pre='V+= w ') 

    Exc_C_loc[h].connect(p = plocal) #this step is taking longest time. rate determining. 
    Inh_C_loc[h].connect(p = plocal) 
    EtoI_C_loc[h].connect(p = plocal) 
    ItoE_C_loc[h].connect(p = plocal) 

    Exc_C_loc[h].w = (1+paraalpha*netwParams_hier[h])*paraomegaEEsp
    Inh_C_loc[h].w = -paraomegaIIsp
    EtoI_C_loc[h].w = (1+paraalpha*netwParams_hier[h])*paraomegaIEsp
    ItoE_C_loc[h].w = -paraomegaEIsp

    m = 0 #long range to m. 
    while m < arealen:
      if m!= h:  
          exc_lr_itoj, etoi_lr_itoj = None, None
  #        print m
          exc_lr_itoj = Synapses(Exc[h], Exc[m], 'w:volt', on_pre='V+= w ') 
          etoi_lr_itoj = Synapses(Exc[h], Inh[m], 'w:volt', on_pre='V+= w ')

          exc_lr_itoj.connect(p = plongr)  #long time.   
          etoi_lr_itoj.connect(p = plongr)  

          exc_lr_itoj.w =  (1 + paraalpha * netwParams_hier[m]) * paramuEEsp * flnMat[m,h]
          etoi_lr_itoj.w = (1 + paraalpha * netwParams_hier[m]) * paramuIEsp * flnMat[m,h]

          meanlr, varlr = delayMat[m,h], .1*delayMat[m,h]
          exc_lr_itoj.delay = np.random.normal(meanlr,varlr,len(exc_lr_itoj.w))*ms
          etoi_lr_itoj.delay = np.random.normal(meanlr,varlr,len(etoi_lr_itoj.w))*ms

          Exc_C_lr_fromi.append(exc_lr_itoj)
          EtoI_C_lr_fromi.append(etoi_lr_itoj)

      m = m + 1       
    h = h + 1

  monitors = SpikeMonitor(E)

  # Setup the network, and run it
  E.V = para['Vr'] + rand(len(E)) * (para['Vt'] - para['Vr'])
  I.V = para['Vr'] + rand(len(I)) * (para['Vt'] - para['Vr'])

  print "before net created"

  net = Network(E,I,Exc_C_loc,EtoI_C_loc,ItoE_C_loc,Inh_C_loc,Exc_C_lr_fromi,EtoI_C_lr_fromi,monitors)

  print "net created"

  net.run(duration, report='text')

  print"hey"

  maxrate = np.empty([arealen,1])
  meanrate = np.empty([arealen,1])   
  netspike = len(monitors.i)
  allspike = np.empty([netspike,2])
  allspike[:,0]=monitors.t/ms
  allspike[:,1]=monitors.i
  allspikesorted = allspike[allspike[:,1].argsort(),]

  netbinno = int( 1+(duration/ms)-(binsize/ms))
  poprate = np.empty([arealen,netbinno ])

  u = 0 #for areas. 
  count = 0#for each spike. 
  stepsize = 1*ms
  monareaktimeall = []
  while u<arealen:
      monareaktime = []

      while((count < netspike) and (allspikesorted[count,1]<1600*(u+1)) ):
        monareaktime.append(allspikesorted[count,0])#append spike times. for each area.
        count = count + 1

      vals= []
      vals = numpy.histogram(monareaktime, bins=int(duration/stepsize))    
      valszero = vals[0]  #now valsbad[0] is a big vector of 500 points. binsize/stepsize = aplus say.
      astep = binsize/(1*ms)
      valsnew = np.zeros(netbinno)    
      acount = 0
      while acount < netbinno:        
          valsnew[acount] = sum(valszero[acount:acount+astep])
          acount=acount+1

      valsrate = valsnew*((1000*ms/binsize) /(1600) ) #new divide by no of neurons per E pop. 
      poprate[u,:] = valsrate   
      maxrate[u,0] = max(valsrate[int(len(valsrate)/3):])
      monareaktimeall.append(monareaktime)
      u = u+1
      #np.save('./poprate_curr_'+str(round(curr,2))+'_seed_'+str(seedval),poprate)                                             
      #print poprate                                                                                                                       
      #print "done"                                                                                                                        

