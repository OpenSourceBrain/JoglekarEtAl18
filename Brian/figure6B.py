#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 14:16:08 2020

@author: ronaldo
"""
from generateNetwork import *
import matplotlib.pyplot as plt

mE1,_,_,_ = run_network(29,'synchronous','strong',600)

################################### Plot #####################################
plt.figure()
# raster plot
plt.plot(mE1[1,:], mE1[0,:], '.',markersize=1)
#  black line
plt.plot([0, max(mE1[1,:])], np.arange(3+1).repeat(2).reshape(-1, 2).T*1600, 'k-')
xlim(250,500)
plt.show()
