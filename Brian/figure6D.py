#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 14:18:42 2020

@author: ronaldo
"""
from generateNetwork import *
import matplotlib.pyplot as plt

_,f1,_,_ = run_network(29,'synchronous','weak',600)
_,f2,_,_ = run_network(29,'9ynchronous','strong',600)


plt.semilogy(np.max(f1,axis=1),'-o',color='green')
plt.semilogy(np.max(f2,axis=1),'-o',color='purple')
plt.show()