#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 18:04:51 2019


Code to simulate predator-prey model for Fussmann study (Science 2000)
Stationary simulations (fixed dilution rate)

@author: ThomasMBury
"""



# import python libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# import EWS function
import sys
sys.path.append('../../../early_warnings')
from ews_compute import ews_compute


#---------------------
# Directory for data output
#–----------------------

# Name of directory within data_export
dir_name = 'ews_stat'

if not os.path.exists('data_export/'+dir_name):
    os.makedirs('data_export/'+dir_name)




#--------------------------------
# Global parameters
#–-----------------------------


# Simulation parameters
dt = 0.01
t0 = 0
tmax = 200 # make large (to get idealised statistics from stationary distribution)
tburn = 100 # burn-in period
seed = 1 # random number generation seed
dbif1 = 0.150983 # first Hopf bifurcation (from Mma bif file)
dbif2 = 0.969538 # second Hopf bifurcation (from Mma bif file)



# EWS parameters
dt2 = 1 # spacing between time-series for EWS computation
rw = 0.25 # rolling window
bw = 0.1 # bandwidth
lags = [1,2,4] # autocorrelation lag times
ews = ['var','ac','sd','cv','skew','kurt','smax','aic','cf'] # EWS to compute
ham_length = 40 # number of data points in Hamming window
ham_offset = 0.5 # proportion of Hamming window to offset by upon each iteration
pspec_roll_offset = 10 # offset for rolling window when doing spectrum metrics






















