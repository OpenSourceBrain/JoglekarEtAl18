"""
Created on Sun Dec  4 13:34:31 2016
@author: maddy
"""
#analyze consciousness data for various random seeds. 
import matplotlib.pyplot as plt
import scipy.io
import numpy as np
import numpy.random
import random as pyrand

from brian2 import *

netbinno = 1 + 550 - 10
popratenet = np.zeros([29, netbinno])

currlist = [0.0] + np.arange(3.5,6.1,1./6).tolist()                            
s_list = [1,2,10,5,18]
diffratecur = np.zeros([29,17])

for currind in range(len(currlist)):

    popratenet = np.zeros([29, netbinno])

    for seedval in s_list:
       curr = currlist[currind]
       #pop = np.load('./poprate_orig/poprate_curr_'+str(round(curr,2))+'_seed_'+str(seedval)+'.npy')
       #pop = np.load('./poprate_curr_'+str(round(curr,2))+'_seed_'+str(seedval)+'.npy')
       #pop = np.load('./poprate_flnshuf2/poprate_curr_'+str(round(curr,2))+'_seed_'+str(seedval)+'.npy')

       for a in range(29):
           diffratecur[a,currind] = diffratecur[a,currind] + ( max(pop[a,300:]) - mean(pop[a,100:285]) )*(1./len(s_list))

diffratecur = ( diffratecur.transpose() - diffratecur[:,0]).transpose()       

list1rate = ['6','9','10','11','15','17','19','22','23','24','25','27','28']

#np.save('./nonnormalizeddiffrate_nofb',diffratecur)
#np.save('./nonnormalizeddiffrate_flnshuf2',diffratecur)

for a in range(29):
#  if str(a) not in list1rate:
    diffratecur[a,:] = diffratecur[a,:]/max(diffratecur[a,:]) 

#np.save('./diffratecur_clust_orig',diffratecur)
#np.save('./diffratecur_clust_orig_nofb',diffratecur)
#np.save('./diffratecur_clust_orig_flnshuf2',diffratecur)

