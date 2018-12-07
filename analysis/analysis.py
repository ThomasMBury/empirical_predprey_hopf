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
import seaborn as sns

# import EWS function
import sys
sys.path.append('../../early_warnings')
from ews_compute import ews_compute
from ews_spec import pspec_welch, pspec_metrics


#-----------------------
# Parameters
#–-----------------------

# EWS computation parmaeters
band_width = 0.2 # band width of Gaussian smoothing
ham_length = 40 # length of Hamming window
w_cutoff = 0.5 # cutoff of higher frequencies
ews = ['var','ac','smax','aic','aic_params','cf','cv'] # EWS to compute
lag_times = [1,2] # lag times for autocorrelation computation





#------------------
## Data curation
#-------------------

# import data
raw = pd.read_excel('../data/raw_fussmann_2000.xls',header=[1])

# round delta column to 2d.p
raw['meandelta'] = raw['meandelta'].apply(lambda x: round(x,2))

# rename delta column
raw.rename(columns={'meandelta':'Delta'}, inplace=True)

## shift day# to start at 0
# function to take list and subtract minimum element
def zero_shift(array):
    return array-min(array)

# unique delta values
deltaVals = raw['Delta'].unique()

# loop through delta values
for d in deltaVals: 
    # shift time values to start at 0
    raw.loc[ raw['Delta']==d,'day#'] = zero_shift(
            raw.loc[ raw['Delta']==d,'day#']) 


## index dataframe by Delta and day#
raw.set_index(['Delta','day#'], inplace=True)            

# compute number of data points for each value of delta
series_lengths=pd.Series(index=deltaVals)
series_lengths.index.name= 'Delta'
for d in deltaVals:
    series_lengths.loc[d] = len(raw.loc[d])
    
# only keep delta values with over 40 data points (so power spec can be computed)
deltaValsFilt = series_lengths[ series_lengths > 40 ].index

               
               
#–--------------------          
## Plot of all trajctories
#-----------------------

## Plot all Chlorella trajectories
#raw['Chlorella'].unstack(level=0).plot(
#        title = 'Chlorella trajectories')
#
## Plot all Brachionus trajectories
#raw['Brachionus'].unstack(level=0).plot(
#        title = 'Brachionus trajectories')


#--------------------------
# Compute EWS
#------------------------

## Chlorella EWS

# Set up a list to store output dataframes of EWS for each delta
appended_ews = []

# Loop through delta values
print('Chlorella')
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
    df_temp['Delta'] = d*np.ones([len(series)])
    # add DataFrame to list
    appended_ews.append(df_temp)
    
    # Print complete    
    print('Delta = '+str(d)+' complete.')
    
# concatenate EWS DataFrames - use delta value and time as indices
df_ews_chlor = pd.concat(appended_ews).set_index('Delta',append=True).reorder_levels([1,0])


## Brachionus EWS
print('\nBrachionus')

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
    df_temp['Delta'] = d*np.ones([len(series)])
    # add DataFrame to list
    appended_ews.append(df_temp)
    
    # Print complete    
    print('Delta = '+str(d)+' complete.')
    
# concatenate EWS DataFrames - use delta value and time as indices
df_ews_brach = pd.concat(appended_ews).set_index('Delta',append=True).reorder_levels([1,0])



## Create dataframes with summary EWS (final entry in time of df_ews)
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
    pspec_temp['Delta'] = d*np.ones(len(pspec_temp))
    # add to appeded ews
    appended_pspec.append(pspec_temp)
    
#  Concatenate dataframes of power spectra and set index as delta value and w
df_pspec_chlor = pd.concat(appended_pspec).set_index(['Delta','Frequency'])


## Brachionus pspec

# Append power spectrum for each delta value into appeded_pspec
appended_pspec = []
for d in deltaValsFilt:
    pspec_temp = pspec_welch(df_ews_brach.loc[d,'Residuals'].values,
                             1, ham_length=ham_length, w_cutoff=w_cutoff)
    # turn into dataframe and reset index
    pspec_temp = pspec_temp.to_frame().reset_index()
    # add a column with delta value
    pspec_temp['Delta'] = d*np.ones(len(pspec_temp))
    # add to appeded ews
    appended_pspec.append(pspec_temp)
  
#  Concatenate dataframes of power spectra and set index as delta value and w
df_pspec_brach = pd.concat(appended_pspec).set_index(['Delta','Frequency'])



##--------------------------------
### Grid plots of power spectra
##--------------------------------
#
#
## Remove indexing for plotting
#plotdf_pspec_chlor = df_pspec_chlor.reset_index()
#plotdf_pspec_brach = df_pspec_brach.reset_index()
#
## Plot properties
#small_size = 8
#med_size = 10
#big_size = 12
#plt.rc('font', size=small_size)          # controls default text sizes
#plt.rc('axes', titlesize=med_size)     # fontsize of the axes title
#plt.rc('axes', labelsize=med_size)    # fontsize of the x and y labels
#plt.rc('xtick', labelsize=small_size)    # fontsize of the tick labels
#plt.rc('ytick', labelsize=small_size)    # fontsize of the tick labels
#plt.rc('legend', fontsize=small_size)    # legend fontsize
#plt.rc('figure', titlesize=med_size)  # fontsize of the figure title
#
## Chlorella plot grid
#plot_pspec_chlor = sns.FacetGrid(plotdf_pspec_chlor, 
#                  col='Delta',
#                  col_wrap=3,
#                  sharey=False,
#                  aspect=1.5,
#                  size=1.8
#                  )
#plot_pspec_chlor = plot_pspec_chlor.map(plt.plot, 'Frequency', 'Power spectrum')
## Change y labels
#axes = plot_pspec_chlor.axes
#for ax in axes[::3]:
#    ax.set_ylabel('Power')
## Change titles
#i=0
#for ax in axes:
#    ax.set_title('Delta = '+str(deltaValsFilt[i]))
#    i=i+1
#
#
### Brachionus plot grid
#plot_pspec_brach = sns.FacetGrid(plotdf_pspec_brach, 
#                  col='Delta',
#                  col_wrap=3,
#                  sharey=False,
#                  aspect=1.5,
#                  size=1.8
#                  )
#plot_pspec_brach = plot_pspec_brach.map(plt.plot, 'Frequency', 'Power spectrum')
## Change y labels
#axes = plot_pspec_brach.axes
#for ax in axes[::3]:
#    ax.set_ylabel('Power')
## Change titles
#i=0
#for ax in axes:
#    ax.set_title('Delta = '+str(deltaValsFilt[i]))
#    i=i+1


#-------------------------
## Compute nonlinear fits for each power spectrum and display as grid plot
#-----------------------------
    
# Initialise a DataFrame
df_spec_metrics = pd.DataFrame([])
# Loop over delta values
for d in deltaValsFilt:
    # compute spectral metrics using pspec_metrics
    series_temp1 = pspec_metrics(
        df_pspec_chlor.loc[d,'Power spectrum'],
        ews=['smax', 'aic', 'aic_params'])
    series_temp2 = pspec_metrics(
        df_pspec_brach.loc[d,'Power spectrum'],
        ews=['smax', 'aic', 'aic_params'])
    # add delta value and species as an entry
    series_temp1['Delta'] = d
    series_temp2['Delta'] = d
    series_temp1['Species'] = 'Chlor'
    series_temp2['Species'] = 'Brach'
    # add to a DataFrame
    df_spec_metrics = df_spec_metrics.append(series_temp1, ignore_index=True)
    df_spec_metrics = df_spec_metrics.append(series_temp2, ignore_index=True)
  
# Index by species and delta
df_spec_metrics.set_index(['Species','Delta'], inplace=True)    

    

# Define the power spectrum models
def fit_fold(w,sigma,lam):
    return (sigma**2 / (2*np.pi))*(1/(w**2+lam**2))
        
def fit_hopf(w,sigma,mu,w0):
    return (sigma**2/(4*np.pi))*(1/((w+w0)**2+mu**2)+1/((w-w0)**2 +mu**2))
        
def fit_null(w,sigma):
    return sigma**2/(2*np.pi)* w**0


# List of data frames to append
append_df = []

# Frequency values from pspec computation
wVals = df_pspec_chlor.index.levels[1]

# Frequency values to use when plotting nonlinear fits
wVals_dense = np.linspace(min(wVals), max(wVals), 10*len(wVals))
#wVals = df_pspec_chlor.index.levels[1]

# Loop over delta values
for d in deltaValsFilt:
    for species in ['Chlor', 'Brach']:
    
        # Fold fit values
        pspec_fold = fit_fold(wVals_dense, df_spec_metrics.loc[species ,d]['Params fold']['sigma'],
                 df_spec_metrics.loc[species, d]['Params fold']['lam'])
        # Hopf fit values
        pspec_hopf = fit_hopf(wVals_dense, df_spec_metrics.loc[species ,d]['Params hopf']['sigma'],
                 df_spec_metrics.loc[species ,d]['Params hopf']['mu'],
                 df_spec_metrics.loc[species ,d]['Params hopf']['w0'])
        # Null fit values
        pspec_null = fit_null(wVals_dense, df_spec_metrics.loc[species ,d]['Params null']['sigma'])
        
        # Create dictionary for dataframe
        dic = {'Species': [species for i in range(len(wVals_dense))],
                          'Delta': d*np.ones(len(wVals_dense)),
                          'Frequency': wVals_dense,
                          'Fold fit': pspec_fold,
                          'Hopf fit': pspec_hopf,
                          'Null fit': pspec_null
                          }
    
        df_temp = pd.DataFrame(data=dic)
        
        # Set the index to frequency and species
        df_temp.set_index(['Species', 'Delta','Frequency'], inplace=True)
        
        # Put empirical power spectra in its own Dataframe indexed by frequency and species
        spec_empirical = df_pspec_chlor.loc[d] if species=='Chlor' else df_pspec_brach.loc[d]
        spec_empirical['Species'] = [species for i in wVals]
        spec_empirical['Delta'] = d*np.ones(len(wVals))
        spec_empirical.set_index(['Species','Delta'], inplace=True, append=True)
        spec_empirical = spec_empirical.reorder_levels(['Species','Delta','Frequency'])
        
        
        # Concatenate with empirical spectrum
        df_temp = pd.concat(
                [df_temp, spec_empirical],
                axis=1)
        
               
        # Add to appended list of dataframes to append
        append_df.append(df_temp)
        

# Concatenate all dataframes
df_pspec_fits = pd.concat(append_df, axis=0)



# Chlorella grid plot

g = sns.FacetGrid(df_pspec_fits.loc['Chlor'].reset_index(level=['Delta','Frequency']), 
                  col='Delta',
                  col_wrap=3,
                  sharey=False,
                  aspect=1.5,
                  size=1.8
                  )
# Axes properties
axes = g.axes
for ax in axes[::3]:
    ax.set_ylabel('Power')
axes[6].set_ylim(top=1.1*max(df_pspec_fits.loc[('Chlor',0.69),'Power spectrum']))
axes[7].set_ylim(top=1.1*max(df_pspec_fits.loc[('Chlor',0.95),'Power spectrum']))

# Plots
g.map(plt.plot, 'Frequency', 'Fold fit', linewidth=1)
g.map(plt.plot, 'Frequency', 'Hopf fit', color='r',linewidth=1)
g.map(plt.plot, 'Frequency', 'Null fit', color='g',linewidth=1)
g.map(plt.plot, 'Frequency', 'Power spectrum', color='k', linewidth=1)


# Brachionus grid plot

g = sns.FacetGrid(df_pspec_fits.loc['Brach'].reset_index(level=['Delta','Frequency']), 
                  col='Delta',
                  col_wrap=3,
                  sharey=False,
                  aspect=1.5,
                  size=1.8
                  )
g.map(plt.plot, 'Frequency', 'Fold fit', linewidth=1)
g.map(plt.plot, 'Frequency', 'Hopf fit', color='r',linewidth=1)
g.map(plt.plot, 'Frequency', 'Null fit', color='g',linewidth=1)
g.map(plt.plot, 'Frequency', 'Power spectrum', color='k', linewidth=1)


# Axes properties
axes = g.axes
for ax in axes[::3]:
    ax.set_ylabel('Power')










#----------------
## Plots of EWS against delta value
#----------------
        
## Make plot of smoothing for some delta value
#df_ews_chlor.loc[1.37,['State variable','Smoothing']].plot(title='Early warning signals')

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



# print out hopf parameters
for d in deltaValsFilt:
    print(df_spec_metrics.loc['Chlor',:]['Params hopf'].loc[d])









