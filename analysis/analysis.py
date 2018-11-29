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


#-----------------------
# Parameters
#–-----------------------

# EWS computation parmaeters
band_width = 0.5 # band width of Gaussian smoothing
ham_length = 40 # length of Hamming window
w_cutoff = 0.5 # cutoff of higher frequencies
ews = ['var','ac','smax','aic','cf','cv'] # EWS to compute
lag_times = [1,2] # lag times for autocorrelation computation





#------------------
## Data curation
#-------------------

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
deltaValsFilt = series_lengths[ series_lengths > 40 ].index

               
               
#–--------------------          
## Plot of trajctories
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

## Chlorella EWS

# Set up a list to store output dataframes of EWS for each delta
appended_ews = []

# Loop through delta values
for d in deltaValsFilt:
    series = raw.loc[d,'Chlorella']
    # plug series into ews_compute - no rolling window (rw=1)
    df_temp = ews_compute(series,
                         roll_window = 1,
                         smooth = True,
                         band_width = band_width,
                         ews = ews,
                         lag_times = lag_times,
                         w_cutoff = w_cutoff
                         )
    # include a column in the dataframe for delta value
    df_temp['meandelta'] = d*np.ones([len(series)])
    # add DataFrame to list
    appended_ews.append(df_temp)
    
# concatenate EWS DataFrames - use delta value and time as indices
df_ews_chlor = pd.concat(appended_ews).set_index('meandelta',append=True).reorder_levels([1,0])


## Brachionus EWS

# Set up a list to store output dataframes of EWS for each delta
appended_ews = []

# Loop through delta values
for d in deltaValsFilt:
    series = raw.loc[d,'Brachionus']
    # plug series into ews_compute - no rolling window (rw=1)
    df_temp = ews_compute(series,
                         roll_window = 1,
                         smooth = True,
                         band_width = band_width,
                         ews = ews,
                         lag_times = lag_times,
                         w_cutoff = w_cutoff
                         )
    # include a column in the dataframe for delta value
    df_temp['meandelta'] = d*np.ones([len(series)])
    # add DataFrame to list
    appended_ews.append(df_temp)
    
# concatenate EWS DataFrames - use delta value and time as indices
df_ews_brach = pd.concat(appended_ews).set_index('meandelta',append=True).reorder_levels([1,0])





#--------------------------
## Compute power spectrum functions
#---------------------------

## Chlorella pspec

# Append power spectrum for each delta value into appeded_pspec
appended_pspec = []
for d in deltaValsFilt:
    pspec_temp = pspec_welch(df_ews_chlor.loc[d,'Residuals'].values,
                             1, ham_length=ham_length, w_cutoff=w_cutoff)
    # turn into dataframe and reset index
    pspec_temp = pspec_temp.to_frame().reset_index()
    # add a column with delta value
    pspec_temp['meandelta'] = d*np.ones(len(pspec_temp))
    # add to appeded ews
    appended_pspec.append(pspec_temp)
    
#  Concatenate dataframes of power spectra and set index as delta value and w
df_pspec_chlor = pd.concat(appended_pspec).set_index(['meandelta','Frequency'])


## Brachionus pspec

# Append power spectrum for each delta value into appeded_pspec
appended_pspec = []
for d in deltaValsFilt:
    pspec_temp = pspec_welch(df_ews_brach.loc[d,'Residuals'].values,
                             1, ham_length=ham_length, w_cutoff=w_cutoff)
    # turn into dataframe and reset index
    pspec_temp = pspec_temp.to_frame().reset_index()
    # add a column with delta value
    pspec_temp['meandelta'] = d*np.ones(len(pspec_temp))
    # add to appeded ews
    appended_pspec.append(pspec_temp)
  
#  Concatenate dataframes of power spectra and set index as delta value and w
df_pspec_brach = pd.concat(appended_pspec).set_index(['meandelta','Frequency'])



## Plot of single power spectrum along with its nonlinear fits
#
## Execute the function pspec_metrics to compute the AIC weights and fitting parameters
#spec_ews = pspec_metrics(pspec, ews=['smax', 'cf', 'aic', 'aic_params'])
## Define the power spectrum models
#def fit_fold(w,sigma,lam):
#    return (sigma**2 / (2*np.pi))*(1/(w**2+lam**2))
#        
#def fit_hopf(w,sigma,mu,w0):
#    return (sigma**2/(4*np.pi))*(1/((w+w0)**2+mu**2)+1/((w-w0)**2 +mu**2))
#        
#def fit_null(w,sigma):
#    return sigma**2/(2*np.pi)* w**0
#
## Make plot
#w_vals = np.linspace(-max(pspec.index), max(pspec.index), 100)
#fig2 = plt.figure(2)
#pspec.plot(label='Measured')
#plt.plot(w_vals, fit_fold(w_vals, spec_ews['Params fold']['sigma'], spec_ews['Params fold']['lam']),label='Fold (AIC='+str(round(spec_ews['AIC fold'],2))+')')
#plt.plot(w_vals, fit_hopf(w_vals, spec_ews['Params hopf']['sigma'], spec_ews['Params hopf']['mu'], spec_ews['Params hopf']['w0']),label='Hopf (AIC='+str(round(spec_ews['AIC hopf'],2))+')')
#plt.plot(w_vals, fit_null(w_vals, spec_ews['Params null']['sigma']),label='Null (AIC='+str(round(spec_ews['AIC null'],2))+')')
#plt.ylabel('Power')
#plt.legend()
#plt.title('Power spectrum and fits at time t='+str(t_pspec))
#








#----------------
## Plots of EWS against delta value
#----------------
        
# Make plot of smoothing for some delta value
df_ews_chlor.loc[1.37,['State variable','Smoothing']].plot(title='Early warning signals')

# Create dataframe with summary EWS (final entry in time of df_ews)
ews_summary_chlor = pd.DataFrame([])
for d in deltaValsFilt:
    series_temp = df_ews_chlor.loc[d].iloc[-1]
    series_temp.name = d
    ews_summary_chlor = ews_summary_chlor.append(series_temp)
    
ews_summary_brach = pd.DataFrame([])
for d in deltaValsFilt:
    series_temp = df_ews_brach.loc[d].iloc[-1]
    series_temp.name = d
    ews_summary_brach = ews_summary_brach.append(series_temp)


# Plot of EWS metrics
fig1, axes = plt.subplots(nrows=5, ncols=1, sharex=True, figsize=(6,6))
ews_summary_chlor[['Variance']].plot(ax=axes[0],title='Early warning signals')
ews_summary_brach[['Variance']].plot(ax=axes[0],secondary_y=True)
ews_summary_chlor[['Coefficient of variation']].plot(ax=axes[1])
ews_summary_brach[['Coefficient of variation']].plot(ax=axes[1],secondary_y=True)
ews_summary_chlor[['Lag-1 AC']].plot(ax=axes[2])
ews_summary_brach[['Lag-1 AC']].plot(ax=axes[2],secondary_y=True)
ews_summary_chlor[['Smax']].plot(ax=axes[3])
ews_summary_brach[['Smax']].plot(ax=axes[3],secondary_y=True)
ews_summary_chlor[['AIC hopf']].plot(ax=axes[4])
ews_summary_brach[['AIC hopf']].plot(ax=axes[4],secondary_y=True)




