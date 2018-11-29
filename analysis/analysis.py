#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 12:19:21 2018

@author: Thomas Bury

Analysis of chemostat data from Fussmann et al.

"""

# import python libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# import EWS function
import sys
sys.path.append('../../early_warnings')
from ews_compute import ews_compute
from ews_spec import pspec_welch, pspec_metrics


# import data
raw = pd.read_excel('../data/raw_fussmann_2000.xls',header=[1])

# round delta column to 2d.p
raw['meandelta'] = raw['meandelta'].apply(lambda x: round(x,2))


## shift day# to start at 0
# function to take list and subtract minimum element
def zero_shift(array):
    return array-min(array)

# unique delta values
deltaVals = raw['meandelta'].unique()

# loop through delta values
for d in deltaVals: 
    # shift time values to start at 0
    raw.loc[ raw['meandelta']==d,'day#'] = zero_shift(
            raw.loc[ raw['meandelta']==d,'day#']) 


## index dataframe by meandelta and day#
raw.set_index(['meandelta','day#'], inplace=True)            

# compute number of data points for each value of delta
series_lengths=pd.Series(index=deltaVals)
series_lengths.index.name= 'meandelta'
for d in deltaVals:
    series_lengths.loc[d] = len(raw.loc[d])
    
# only keep delta values with over 40 data points (so power spec can be computed)
deltaValsFilt = series_lengths[ series_lengths>40 ].index

               
               
#–--------------------          
## Some plots of trajectories
#-----------------------

# Plot all Chlorella trajectories
raw['Chlorella'].unstack(level=0).plot(
        title = 'Chlorella trajectories')

# Plot all Brachionus trajectories
raw['Brachionus'].unstack(level=0).plot(
        title = 'Brachionus trajectories')


#--------------------------
# Compute EWS
#------------------------

# Make dataframe indexed by delta with columns for each EWS
df_ews_chlor = pd.DataFrame([])
df_ews_brach = pd.DataFrame([])

# Chlorella EWS compute
for d in deltaValsFilt:
    series = raw.loc[d,'Chlorella']
    # plug series into ews_compute - no rolling window (rw=1)
    df_ews = ews_compute(series,
                         roll_window = 1,
                         smooth = True,
                         band_width = 0.5,
                         ews = ['var','ac','smax','aic','cf','cv'],
                         lag_times = [1,2],
                         w_cutoff = 0.7
                         )
    # Final entry of dataframe gives overall EWS for time-series (no rollwindow)
    series_ews = df_ews.iloc[-1]
    series_ews.name = d
    # Add to dataframe 
    df_ews_chlor = df_ews_chlor.append(series_ews)
    
# Brachionus EWS compute
for d in deltaValsFilt:
    series = raw.loc[d,'Brachionus']
    # plug series into ews_compute - no rolling window (rw=1)
    df_ews = ews_compute(series,
                         roll_window = 1,
                         smooth = True,
                         band_width = 0.2,
                         ews = ['var','ac','smax','aic','cf','cv'],
                         lag_times = [1,2],
                         w_cutoff = 0.7
                         )
    # Final entry of dataframe gives overall EWS for time-series (no rollwindow)
    series_ews = df_ews.iloc[-1]
    series_ews.name = d
    # Add to dataframe 
    df_ews_brach = df_ews_brach.append(series_ews)   
    

#--------------------------
## Compute power spectrum functions
#---------------------------
    









#----------------
## EWS plots
#----------------
    
    
# Make plot of smoothing
df_ews[['State variable','Smoothing']].plot(title='Early warning signals')

# Plot of EWS metrics
fig1, axes = plt.subplots(nrows=5, ncols=1, sharex=True, figsize=(6,6))
df_ews_chlor[['Variance']].plot(ax=axes[0],title='Early warning signals')
df_ews_brach[['Variance']].plot(ax=axes[0],secondary_y=True)
df_ews_chlor[['Coefficient of variation']].plot(ax=axes[1])
df_ews_brach[['Coefficient of variation']].plot(ax=axes[1],secondary_y=True)
df_ews_chlor[['Lag-1 AC']].plot(ax=axes[2])
df_ews_brach[['Lag-1 AC']].plot(ax=axes[2],secondary_y=True)
df_ews_chlor[['Smax']].plot(ax=axes[3])
df_ews_brach[['Smax']].plot(ax=axes[3],secondary_y=True)
df_ews_chlor[['AIC hopf']].plot(ax=axes[4])
df_ews_brach[['AIC hopf']].plot(ax=axes[4],secondary_y=True)




