#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 11 2019

@author: Thomas Bury

Early warning signal anlaysis using bootstrapping with the Fussmann dataset

"""

# Import python libraries
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Import EWS function
import sys
sys.path.append('../../early_warnings')
from ews_compute import ews_compute

# Import bootstrap module
sys.path.append('../../ews_bootstrap')
from roll_bootstrap import roll_bootstrap, mean_ci



# Name of directory within data_export
dir_name = 'bootstrap_block20_span05'

if not os.path.exists('../data_export/'+dir_name):
    os.makedirs('../data_export/'+dir_name)


#-----------------------
# Parameters
#–-----------------------

# EWS computation parameters
span = 0.5 # span used for Loess filtering of time-series
ham_length = 40 # length of Hamming window
ham_offset = 0.5 # offset of Hamming windows
w_cutoff = 0.8 # cutoff of higher frequencies
ews = ['var','ac','smax','aic'] # EWS to compute
lags = [1,2,10] # lag times for autocorrelation computation (lag of 10 to show decreasing AC where tau=T/2)
sweep = False # whether to sweep over initialisation parameters during AIC fitting


# Bootstrapping parameters
block_size = 40 # size of blocks used to resample time-series
bs_type = 'Stationary' # type of bootstrapping
n_samples = 2 # number of bootstrapping samples to take
roll_offset = 20 # rolling window offset (note: don't need a rolling window for Fussmann data)







#------------------
## Data import and curation
#-------------------

# Import raw data
raw = pd.read_excel('../data/raw_fussmann_2000.xls',header=[1])

# Round delta column to 2d.p
raw['meandelta'] = raw['meandelta'].apply(lambda x: round(x,2))

# Rename column labels
raw.rename(columns={'meandelta':'Delta', 'day#':'Time'}, inplace=True)

## Shift day# to start at 0

# Function to take list and subtract minimum element
def zero_shift(array):
    return array-min(array)

# Unique delta values
deltaVals = raw['Delta'].unique()

# Loop through delta values
for d in deltaVals: 
    # Shift time values to start at 0
    raw.loc[ raw['Delta']==d,'Time'] = zero_shift(
            raw.loc[ raw['Delta']==d,'Time']) 

## Index dataframe by Delta and Time
raw.set_index(['Delta','Time'], inplace=True)            

# Export trajectories as a csv file
raw_traj = raw[['Chlorella','Brachionus']]
raw_traj.to_csv("../data_export/series_data.csv", index_label=['Delta','Time'])               
               
               
# Compute number of data points for each value of delta
series_lengths=pd.Series(index=deltaVals)
series_lengths.index.name= 'Delta'
for d in deltaVals:
    series_lengths.loc[d] = len(raw.loc[d])
    
# Only consider the delta values for which the corresponding trajecotories have over 25 data points
deltaValsFilt = series_lengths[ series_lengths > 25 ].index

       






#-------------------------------------
# Compute bootstrapped samples of time-series data
#–----------------------------------

# Compute bootstrapped residuals of each time-series

# Initialise list for DataFrames of bootstrapped time-series
list_df_samples = []

# Loop over dilution rates
for d in deltaValsFilt:
    # Loop over species
    for species in ['Chlorella','Brachionus']:
        
        # Time-series to work with
        series = raw_traj.loc[d][species]
        
        # Compute bootstrapped series   
        df_samples_temp = roll_bootstrap(series,
                           span = span,
                           roll_window = 1,
                           n_samples = n_samples,
                           bs_type = bs_type,
                           block_size = block_size
                           )

        # Add column for dilution rate and species
        df_samples_temp['Species'] = species
        df_samples_temp['Dilution rate'] = d
        
        # Remove indexing
        df_samples_temp.reset_index(inplace=True)
        
        # Add DF to list
        list_df_samples.append(df_samples_temp)
        
# Concatenate dataframes of samples
df_samples = pd.concat(list_df_samples)

# Drop the Time column (redundant as no rolling window used)
df_samples.drop('Time', axis=1, inplace=True)

# Set a sensible index
df_samples.set_index(['Dilution rate', 'Species', 'Sample', 'Wintime'], inplace=True)

# Sort the index
df_samples.sort_index(inplace=True)


        
        
        
#----------------------------------------------
# Compute EWS for each bootstrapped time-series
#-------------------------------------------------

# List to store EWS DataFrames
list_df_ews = []
# List to store power spectra DataFrames of one of the samples
list_pspec = []

# Sample values
sampleVals = np.array(df_samples.index.levels[2])


# Loop through dilution rate
for d in deltaValsFilt:
    
    # Loop through species
    for species in ['Chlorella', 'Brachionus']:
        
        # Loop through sample values
        for sample in sampleVals:
            
            # Compute EWS of sample series
            series_temp = df_samples.loc[(d, species, sample, ),'x']
            
            ews_dic = ews_compute(series_temp,
                              roll_window = 1, 
                              smooth = 'None',
                              ews = ews,
                              lag_times = lags,
                              upto='Full',
                              sweep=sweep,
                              ham_length=ham_length,
                              ham_offset=ham_offset)
            
            ## The DataFrame of EWS
            df_ews_temp = ews_dic['EWS metrics']
            
            # Include columns for dilution rate, species and sample number
            df_ews_temp['Dilution rate'] = d
            df_ews_temp['Species'] = species
            df_ews_temp['Sample'] = sample

            # Drop NaN values
            df_ews_temp = df_ews_temp.dropna()        
            
            # Append list_df_ews
            list_df_ews.append(df_ews_temp)
            
            
            ## The DataFrame of power spectra
            df_pspec_temp = ews_dic['Power spectrum'][['Empirical']].dropna()
            # Add columns for species and dilution rate
            df_pspec_temp['Dilution rate'] = d
            df_pspec_temp['Species'] = species
            df_pspec_temp['Sample'] = sample
            
            # Append list
            list_pspec.append(df_pspec_temp)
    
    # Print update
    print('EWS for d = %.2f complete' % d)
        
# Concatenate EWS DataFrames. Index [Dilution rate, species , Sample]
df_ews_boot = pd.concat(list_df_ews).reset_index(drop=True).set_index(['Dilution rate','Species', 'Sample'])
# Sort the index 
df_ews_boot.sort_index(inplace=True)

# Concatenate power spectrum DataFrames
df_pspec_boot = pd.concat(list_pspec).reset_index().set_index(['Dilution rate','Species','Sample','Frequency'])
# Drop the time column
df_pspec_boot.drop('Time', inplace=True, axis=1)
df_pspec_boot.sort_index(inplace=True)



#--------------------------
# Compute confidence intervals of EWS
#-----------------------------












#----------------
## Some plots of EWS
#----------------


# Smax
df_ews_boot.xs(('Chlorella'), level=('Species'))['Smax'].unstack(level='Sample').plot()

# AIC hopf
df_ews_boot.xs(('Chlorella'), level=('Species'))['AIC hopf'].unstack(level='Sample').plot()

# Power spectrum
df_pspec_boot.loc[0.04, 'Brachionus']['Empirical'].unstack(level='Sample').plot()

df_pspec_boot.loc[0.04, 'Chlorella']['Empirical'].unstack(level='Sample').plot()





#-----------------------
# Export data for plotting in MMA
#-----------------------
    

# EWS data
ews_export = ['Variance','Lag-1 AC','Lag-2 AC','Lag-10 AC','Smax','AIC fold','AIC hopf','AIC null']
df_ews_boot[ews_export].to_csv('../data_export/'+dir_name+'/ews.csv')


# Export empirical pspec data for plotting in MMA
df_pspec_boot.to_csv('../data_export/'+dir_name+'/pspec.csv')










