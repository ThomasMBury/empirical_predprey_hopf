#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 12:19:21 2018

@author: Thomas Bury

Analysis of chemostat data from Fussmann et al.

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




#-----------------------
# Parameters
#–-----------------------

# EWS computation parmaeters
band_width = 0.1 # band width of Gaussian smoothing
ham_length = 40 # length of Hamming window
ham_offset = 0.5 # offset of Hamming windows
w_cutoff = 0.8 # cutoff of higher frequencies
ews = ['var','ac','smax','aic','aic_params','cf','cv'] # EWS to compute
lag_times = [1,2] # lag times for autocorrelation computation



#------------------
## Data import and curation
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

# export trajectories as a csv file
raw_traj = raw[['Chlorella','Brachionus']]
raw_traj.to_csv("../data_export/series_data.csv", index_label=['Delta','day#'])               
               
               
# compute number of data points for each value of delta
series_lengths=pd.Series(index=deltaVals)
series_lengths.index.name= 'Delta'
for d in deltaVals:
    series_lengths.loc[d] = len(raw.loc[d])
    
# only consider the delta values for which the corresponding trajecotories have over 25 data points
deltaValsFilt = series_lengths[ series_lengths > 25 ].index

           

#--------------------------
# Compute EWS
#------------------------


# Set up a list to store output dataframes for each delta
appended_ews = []
appended_pspec = []

# Loop through species
for species in ['Chlorella','Brachionus']:
    # Update status
    print('\n'+species)
    # Loop through delta values
    for d in deltaValsFilt:
        series = raw.loc[d,species]
        # plug series into ews_compute - no rolling window (rw=1)
        ews_dic = ews_compute(series,
                             roll_window = 1,
                             smooth = True,
                             band_width = band_width,
                             ews = ews,
                             lag_times = lag_times,
                             ham_length = ham_length,
                             ham_offset = ham_offset,
                             w_cutoff = w_cutoff
                             )
        # DataFrame of EWS metrics
        df_ews_temp = ews_dic['EWS metrics']
        # DataFrame for power spectra and their fits
        df_pspec_temp = ews_dic['Power spectrum']
        
        # Include a column for delta value and species in each DataFrame
        df_ews_temp['Delta'] = d
        df_ews_temp['Species'] = species
        
        df_pspec_temp['Delta'] = d
        df_pspec_temp['Species'] = d
        
        # add DataFrames to list
        appended_ews.append(df_ews_temp)
        appended_pspec.append(df_pspec_temp)
        
        # Update status    
        print('Delta = '+str(d)+' complete.')
    
# concatenate EWS DataFrames - use species, delta and time as indices
df_ews_full = pd.concat(appended_ews).reset_index().set_index(['Species','Delta','Time'])
# concatenate pspec DataFrames - use species, delta, time and frequency as indices
df_pspec = pd.concat(appended_pspec).reset_index().set_Index(['Species','Delta','Time','Frequency'])

## Reduce dataframes by getting rid of NaN cells and dropping the Time index
df_ews = df_ews_full.dropna().reset_index(level=2, drop=True) 
df_pspec = df_pspec.reset_index(level=2, drop=True)



#--------------------------
## Grid plot of all trajectories and smoothing
#--------------------------



## Chlorella trajectories plot
g = sns.FacetGrid(df_ews_full.reset_index(level=['Delta','Time']), 
                  col='Delta',
                  col_wrap=3,
                  sharey=False,
                  aspect=1.5,
                  size=1.8
                  )
# Plots
plt.rc('axes', titlesize=10) 
g.map(plt.plot, 'Time', 'State variable', color='tab:blue', linewidth=1)
g.map(plt.plot, 'Time', 'Smoothing', color='tab:orange', linewidth=1)
# Axes properties
axes = g.axes
# Assign plot label
plot_traj_chlor=g



## Brachionus trajectories plot
g = sns.FacetGrid(df_ews_brach_full.reset_index(level=['Delta','Time']), 
                  col='Delta',
                  col_wrap=3,
                  sharey=False,
                  aspect=1.5,
                  size=1.8
                  )
# Plots
plt.rc('axes', titlesize=10) 
g.map(plt.plot, 'Time', 'State variable', color='tab:blue', linewidth=1)
g.map(plt.plot, 'Time', 'Smoothing', color='tab:orange', linewidth=1)
# Axes properties
axes = g.axes
for ax in axes[::3]:
    ax.set_ylabel('Brachionus')
# Assign plot label
plot_traj_chlor=g


#----------------
## Plots of EWS against delta value
#----------------
        


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





#---------------------------------
## Power spectra visualisation
#--------------------------------

# Limits for x-axis
xmin = -2.5
xmax = 2.5

## Chlorella

g = sns.FacetGrid(df_pspec_chlor.reset_index(level=['Delta','Frequency']), 
                  col='Delta',
                  col_wrap=3,
                  sharey=False,
                  aspect=1.5,
                  size=1.8
                  )
# Plots
plt.rc('axes', titlesize=10) 
g.map(plt.plot, 'Frequency', 'Empirical', color='k', linewidth=1)
g.map(plt.plot, 'Frequency', 'Fit fold', color='b', linestyle='dashed', linewidth=1)
g.map(plt.plot, 'Frequency', 'Fit hopf', color='r', linestyle='dashed', linewidth=1)
g.map(plt.plot, 'Frequency', 'Fit null', color='g', linestyle='dashed', linewidth=1)
# Axes properties
axes = g.axes
# Global axes properties
for i in range(len(axes)):
    ax=axes[i]
    d=deltaValsFilt[i]
    ax.set_ylim(bottom=0, top=1.1*max(df_pspec_chlor.loc[d]['Empirical'].loc[xmin:xmax].dropna()))
    ax.set_xlim(left=xmin, right=xmax)
    ax.set_xticks([-2,-1,0,1,2])
    ax.set_title('Delta = %.2f' % deltaValsFilt[i])
    # AIC weights
    xpos=0.7
    ypos=0.9
    ax.text(xpos,ypos,
            '$w_f$ = %.1f' % df_ews_chlor.loc[d]['AIC fold'],
            fontsize=9,
            color='b',
            transform=ax.transAxes)  
    ax.text(xpos,ypos-0.12,
            '$w_h$ = %.1f' % df_ews_chlor.loc[d]['AIC hopf'],
            fontsize=9,
            color='r',
            transform=ax.transAxes)
    ax.text(xpos,ypos-2*0.12,
            '$w_n$ = %.1f' % df_ews_chlor.loc[d]['AIC null'],
            fontsize=9,
            color='g',
            transform=ax.transAxes)
# Y labels
for ax in axes[::3]:
    ax.set_ylabel('Power')
    
# Specific Y limits
axes[1].set_ylim(top=0.004)
axes[8].set_ylim(top=0.04)
#axes[3].set_ylim(top=3)

# Assign to plot label
pspec_plot_chlor=g





## Brachionus grid plot
# x-axes limits
g = sns.FacetGrid(df_pspec_brach.reset_index(level=['Delta','Frequency']), 
                  col='Delta',
                  col_wrap=3,
                  sharey=False,
                  aspect=1.5,
                  size=1.8
                  )
# Plots
plt.rc('axes', titlesize=10) 
g.map(plt.plot, 'Frequency', 'Empirical', color='k', linewidth=1)
g.map(plt.plot, 'Frequency', 'Fit fold', color='b', linestyle='dashed', linewidth=1)
g.map(plt.plot, 'Frequency', 'Fit hopf', color='r', linestyle='dashed', linewidth=1)
g.map(plt.plot, 'Frequency', 'Fit null', color='g', linestyle='dashed', linewidth=1)
# Axes properties
axes = g.axes
# Global axes properties
for i in range(len(axes)):
    ax=axes[i]
    d=deltaValsFilt[i]
    ax.set_ylim(bottom=0, top=1.1*max(df_pspec_brach.loc[d]['Empirical'].loc[xmin:xmax].dropna()))
    ax.set_xlim(left=xmin, right=xmax)
    ax.set_xticks([-2,-1,0,1,2])
    ax.set_title('Delta = %.2f' % deltaValsFilt[i])
    # AIC weights
    xpos=0.7
    ypos=0.9
    ax.text(xpos,ypos,
            '$w_f$ = %.1f' % df_ews_brach.loc[d]['AIC fold'],
            fontsize=9,
            color='b',
            transform=ax.transAxes)  
    ax.text(xpos,ypos-0.12,
            '$w_h$ = %.1f' % df_ews_brach.loc[d]['AIC hopf'],
            fontsize=9,
            color='r',
            transform=ax.transAxes)
    ax.text(xpos,ypos-2*0.12,
            '$w_n$ = %.1f' % df_ews_brach.loc[d]['AIC null'],
            fontsize=9,
            color='g',
            transform=ax.transAxes)
# Y labels
for ax in axes[::3]:
    ax.set_ylabel('Power')
    
# Specific Y limits
axes[1].set_ylim(top=1)
axes[2].set_ylim(top=1)
axes[3].set_ylim(top=3)

# Assign to plot label
pspec_plot_brach=g




#-----------------------
# Export data and plots
#-----------------------
    

# Export pspec plots
pspec_plot_chlor.savefig("../figures/pspec_grid_chlor.png", dpi=200)
pspec_plot_brach.savefig("../figures/pspec_grid_brach.png", dpi=200)

# Export EWS data for plotting in MMA
cols=['Variance','Coefficient of variation','Lag-1 AC','Lag-2 AC','Smax','AIC fold','AIC hopf','AIC null']
df_ews_chlor[cols].to_csv("../data_export/ews_chlor.csv")
df_ews_brach[cols].to_csv("../data_export/ews_brach.csv")

# Export empirical pspec data for plotting in MMA
df_pspec_chlor['Empirical'].dropna().to_frame().to_csv('../data_export/pspec_chlor.csv',index_label=['Delta','Frequency'])
df_pspec_brach['Empirical'].dropna().to_frame().to_csv('../data_export/pspec_brach.csv',index_label=['Delta','Frequency'])










