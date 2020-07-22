#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 10:25:17 2020

@author: ronaldo
"""


import numpy as np
import os
import scipy.io
import matplotlib.pyplot as plt

# Load data
path=os.path.abspath(os.path.join(os.getcwd(), os.pardir))+'/Matlab/subgraphData.mat'
flnMat = scipy.io.loadmat(path)
conn=flnMat['flnMat'][:][:]

# Ventral areas
ventarea = np.array([0,1,2,8,18])
tempConn =conn[ventarea,:]
conn=tempConn[:,ventarea]

# Adjust bins
t=np.logspace(-6,0,13)
bins = t[:-1] + np.diff(t)/2
valores=[]
valores.extend(bins.tolist())

valores.append(bins[-1]+1.42302495)

# Plot data
plt.hist(conn[conn>0],bins=np.array(valores),ec='black',color='violet') 
plt.gca().spines['right'].set_color('none')
plt.gca().spines['top'].set_color('none')
plt.ylim([0, 5])
plt.xlim([10**-6, 3])
plt.ylabel('FLN: ventral stream')
plt.xscale('log')

plt.show()