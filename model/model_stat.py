#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 18:04:51 2019


Code to simulate predator-prey model for Fussmann study (Science 2000)
Stationary simulations (fixed dilution rate)

@author: ThomasMBury
"""



#------------------------
# Import modules
#–----------------------

# Default python modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# EWS module
import sys
sys.path.append('../../early_warnings')
from ews_compute import ews_compute


#----------------------
# Useful functions
#-----------------------

# Apply operation to column of DataFrame in place
def apply_inplace(df, field, fun):
    """ Apply function to a column of a DataFrame in place."""
    return pd.concat([df.drop(field, axis=1), df[field].apply(fun)], axis=1)



#---------------------
# Directory for data output
#–----------------------

# Name of directory within data_export
dir_name = 'ews_stat_sweep_long'

if not os.path.exists('data_export/'+dir_name):
    os.makedirs('data_export/'+dir_name)




#--------------------------------
# Global parameters
#–-----------------------------


# Simulation parameters
dt = 0.01
t0 = 0
tmax = 4000 # make large (to get idealised statistics from stationary distribution)
tburn = 500 # burn-in period
seed = 1 # random number generation seed
dbif1 = 0.150983 # first Hopf bifurcation (from Mma bif file)
dbif2 = 0.969538 # second Hopf bifurcation (from Mma bif file)
dl = 0.005 # low delta value
dh = 1.6 # high delta value



# EWS parameters
dt2 = 1 # spacing between time-series for EWS computation
rw = 1 # rolling window (compute EWS using full time-series)
bw = 1 # bandwidth (take the whole dataset as stationary)
lags = [1,2,3,6] # autocorrelation lag times
ews = ['var','ac','sd','cv','skew','kurt','smax','aic','cf'] # EWS to compute
ham_length = 40 # number of data points in Hamming window
ham_offset = 0.5 # proportion of Hamming window to offset by upon each iteration
pspec_roll_offset = 20 # offset for rolling window when doing spectrum metrics


#----------------------------------
# Model
#----------------------------------


# Function for model dynamics (variables [n,c,r,b])
def de_fun(state, control, params):
    '''
    Inputs:
        state: array of state variables [n,c,r,b]
        control: control parameter that is to be varied
        params: list of parameter values [ni,bc,kc,bb,kb,epsilon,m,lamda]
    Output:
        array of gradient vector (derivative)
    '''
    [n,c,r,b] = state
    [ni,bc,kc,bb,kb,epsilon,m,lamda] = params
    d = control
    
    # Gradients for each variable to increment by
    n_grad = d*(ni-n) - bc*n*c/(kc+n)
    c_grad = bc*n*c/(kc+n) - bb*c*b/((kb+c)*epsilon) - d*c
    r_grad = bb*c*r/(kb+c) - (d+m+lamda)*r
    b_grad = bb*c*r/(kb+c) - (d+m)*b
            
    return np.array([n_grad, c_grad, r_grad, b_grad])
    
    
   
# System parameters
    
ni=80 # nitrogen inflow concentration
bc=3.3 # maximum birth rate of Chlorella
kc=4.3  # half saturation constant of Chlorella
bb=2.25 # maximum birth rate of Brachionus
kb=15   # half-saturation constant of Brachionus
epsilon=0.25    # assimilation efficacy of Brachionus
m=0.055 # mortality of Brachionus
lamda=0.4   # decay of fecundity of Brachionus

# Parameter list
params = [ni,bc,kc,bb,kb,epsilon,m,lamda]

# Control parameter values
deltaVals = np.arange(dl, dh, 0.005)


# Noise parameters
sigma_n = 0 # amplitude for N
sigma_c = 0.01 # amplitude for Chlorella
sigma_r = 0 # amplitude for R
sigma_b = 0.02 # amplitude for Brachionus

# Initial conditions
n0 = 2
c0 = 5
r0 = 1
b0 = 2




#
#
##--------------------------------------------
## Simulate (stationary) realisations of model for each delta value
##-------------------------------------------
#
#
#
#
#
### Implement Euler Maryuyama for stocahstic simulation
#
#
## Set seed
#np.random.seed(seed)
#
## Initialise a list to collect trajectories
#list_traj_append = []
#
## loop over delta values
#print('\nBegin simulations \n')
#for d in deltaVals:
#    
#    # Initialise array to store time-series data
#    t = np.arange(t0,tmax,dt) # Time array
#    x = np.zeros([len(t), 4]) # State array
#
#    
#    # Create brownian increments (s.d. sqrt(dt))
#    dW_n_burn = np.random.normal(loc=0, scale=sigma_n*np.sqrt(dt), size = int(tburn/dt))
#    dW_n = np.random.normal(loc=0, scale=sigma_n*np.sqrt(dt), size = len(t)) 
#    
#    dW_c_burn = np.random.normal(loc=0, scale=sigma_c*np.sqrt(dt), size = int(tburn/dt))
#    dW_c = np.random.normal(loc=0, scale=sigma_c*np.sqrt(dt), size = len(t))
#  
#    dW_r_burn = np.random.normal(loc=0, scale=sigma_r*np.sqrt(dt), size = int(tburn/dt))
#    dW_r = np.random.normal(loc=0, scale=sigma_r*np.sqrt(dt), size = len(t))
#    
#    dW_b_burn = np.random.normal(loc=0, scale=sigma_b*np.sqrt(dt), size = int(tburn/dt))
#    dW_b = np.random.normal(loc=0, scale=sigma_b*np.sqrt(dt), size = len(t))
#    
#    # Noise vectors
#    dW_burn = np.array([dW_n_burn,
#                        dW_c_burn,
#                        dW_r_burn,
#                        dW_b_burn]).transpose()
#    
#    dW = np.array([dW_n, dW_c, dW_r, dW_b]).transpose()
# 
#    # IC as a state vector
#    x0 = np.array([n0,c0,r0,b0])
#    
#    # Run burn-in period on initial condition
#    for i in range(int(tburn/dt)):
#        # Update with Euler Maruyama
#        x0 = x0 + de_fun(x0, d, params)*dt + dW_burn[i]
#        # Make sure that state variable remains >= 0 
#        x0 = [np.max([k,0]) for k in x0]
#        
#        
#    # Initial condition post burn-in period
#    x[0]=x0
#    
#    # Run simulation
#    for i in range(len(t)-1):
#        x[i+1] = x[i] + de_fun(x[i], d, params)*dt + dW[i]
#        # make sure that state variable remains >= 0 
#        x[i+1] = [np.max([k,0]) for k in x[i+1]]
#            
#    # Store series data in a DataFrame
#    data = {'Delta': d,
#                'Time': t,
#                'Nitrogen': x[:,0],
#                'Chlorella': x[:,1],
#                'Reproducing Brachionus': x[:,2],
#                'Brachionus': x[:,3]}
#    df_temp = pd.DataFrame(data)
#    # Append to list
#    list_traj_append.append(df_temp)
#    
#    print('Simulation with d='+str(d)+' complete')
#
##  Concatenate DataFrame from each realisation
#df_traj = pd.concat(list_traj_append)
#df_traj.set_index(['Delta','Time'], inplace=True)
#
#
## Coarsen time-series to have spacing dt2 (for EWS computation)
#df_traj_filt = df_traj.loc[::int(dt2/dt)]
#
#
## Make units of Chlor and Brach consistent with Fussmann experiments
#
#nC = 50000/1000000 # conversion factor to 10^6 cells/ml of Chlorella
#nB = 5 # conversion factor to females/ml of Brachiouns
#
#df_traj_filt = apply_inplace(df_traj_filt, 'Chlorella',
#                             lambda x: nC*x)
#df_traj_filt = apply_inplace(df_traj_filt, 'Reproducing Brachionus',
#                             lambda x: nB*x)
#df_traj_filt = apply_inplace(df_traj_filt, 'Brachionus',
#                             lambda x: nB*x)
#
#
## Export simulation data to SD Memory
#df_traj_filt.to_csv('/Volumes/SDMemory/Datasets/fussmann/sim_data.csv')
#


# Import simulation data from SD Memory
df_traj_filt = pd.read_csv('/Volumes/SDMemory/Datasets/fussmann/sim_data.csv')
df_traj_filt.set_index(['Delta','Time'], inplace=True)

# Delta Values
deltaVals = np.array(df_traj_filt.index.levels[0])


#----------------------
## Execute ews_compute for each delta value and each variable
#---------------------


# Set up a list to store output dataframes from ews_compute
# We will concatenate them at the end
appended_ews = []
appended_pspec = []

# loop through realisation number
print('\nBegin EWS computation\n')
for d in deltaVals:
    # loop through sate variable
    for var in ['Chlorella', 'Brachionus']:
        
        ews_dic = ews_compute(df_traj_filt.loc[d][var], 
                          roll_window = rw, 
                          band_width = bw,
                          lag_times = lags, 
                          ews = ews,
                          ham_length = ham_length,
                          ham_offset = ham_offset,
                          pspec_roll_offset = pspec_roll_offset
                          )
        
        # The DataFrame of EWS
        df_ews_temp = ews_dic['EWS metrics']
        # The DataFrame of power spectra
        df_pspec_temp = ews_dic['Power spectrum']
        
        # Include a column in the DataFrames for delta value and variable
        df_ews_temp['Delta'] = d
        df_ews_temp['Variable'] = var
        
        df_pspec_temp['Delta'] = d
        df_pspec_temp['Variable'] = var
                
        # Add DataFrames to list
        appended_ews.append(df_ews_temp)
        appended_pspec.append(df_pspec_temp)
        
    # Print status every realisation
    print('EWS for delta =  '+str(d)+' complete')


# Concatenate EWS DataFrames. Index [Delta, Variable, Time]
df_ews_full = pd.concat(appended_ews).reset_index().set_index(['Delta','Variable','Time'])
# Concatenate power spectrum DataFrames. Index [Realisation number, Variable, Time, Frequency]
df_pspec = pd.concat(appended_pspec).reset_index().set_index(['Delta','Variable','Time','Frequency'])


# Refine DataFrame to just have EWS data (no time dependence)
df_ews = df_ews_full.dropna().reset_index(level=2, drop=True).reorder_levels(['Variable', 'Delta'])
df_pspec = df_pspec.reset_index(level=2, drop=True).reorder_levels(['Variable', 'Delta','Frequency'])



#--------------------------
## Grid plot of all trajectories
#--------------------------

# Set up frame and axes
g = sns.FacetGrid(df_ews_full.reset_index(), 
                  col='Delta',
                  hue='Variable',
                  palette='Set1',
                  col_wrap=3,
                  sharey=False,
                  aspect=1.5,
                  height=1.8
                  )
# Set plot title size
plt.rc('axes', titlesize=10)
# Plot state variable
g.map(plt.plot, 'Time', 'State variable', linewidth=1)
## Plot smoothing
#g.map(plt.plot, 'Time', 'Smoothing', color='tab:orange', linewidth=1)
# Axes properties
axes = g.axes
# Assign plot label
plot_traj = g

# Axes properties
axes = g.axes
for i in range(len(axes)):
    ax=axes[i]
    d=deltaVals[i]
    ax.set_ylim(bottom=0, top=60)

    
## Export plot
#plot_traj.savefig("../figures/empirical_series2.png", dpi=200)


#----------------
## Plots of EWS against delta value
#----------------

# Plot of EWS metrics
fig1, axes = plt.subplots(nrows=5, ncols=1, sharex=True, figsize=(6,6))
df_ews.loc['Chlorella'][['Variance']].plot(ax=axes[0],title='Early warning signals')
df_ews.loc['Brachionus'][['Variance']].plot(ax=axes[0],secondary_y=True)
df_ews.loc['Chlorella'][['Coefficient of variation']].plot(ax=axes[1])
df_ews.loc['Brachionus'][['Coefficient of variation']].plot(ax=axes[1],secondary_y=True)
df_ews.loc['Chlorella'][['Lag-1 AC']].plot(ax=axes[2])
df_ews.loc['Brachionus'][['Lag-1 AC']].plot(ax=axes[2],secondary_y=True)
df_ews.loc['Chlorella'][['Smax']].plot(ax=axes[3])
df_ews.loc['Brachionus'][['Smax']].plot(ax=axes[3],secondary_y=True)
df_ews.loc['Chlorella'][['AIC hopf']].plot(ax=axes[4], ylim=(0,1.1))
df_ews.loc['Brachionus'][['AIC hopf']].plot(ax=axes[4],
          secondary_y=True, ylim=(0,1.1))



#---------------------------------
## Power spectra visualisation
#--------------------------------

# Limits for x-axis
xmin = -np.pi
xmax = np.pi

## Chlorella
species='Chlorella'
g = sns.FacetGrid(df_pspec.loc[species].reset_index(level=['Delta','Frequency']), 
                  col='Delta',
                  col_wrap=3,
                  sharey=False,
                  aspect=1.5,
                  height=1.8
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
    d=deltaVals[i]
#    ax.set_ylim(bottom=0, top=1.1*max(df_pspec.loc[species,d]['Empirical'].loc[xmin:xmax].dropna()))
    ax.set_xlim(left=xmin, right=xmax)
    ax.set_xticks([-3,-2,-1,0,1,2,3])
    ax.set_title('Delta = %.2f' % deltaVals[i])
    # AIC weights
    xpos=0.7
    ypos=0.9
    ax.text(xpos,ypos,
            '$w_f$ = %.1f' % df_ews.loc[species,d]['AIC fold'],
            fontsize=9,
            color='b',
            transform=ax.transAxes)  
    ax.text(xpos,ypos-0.12,
            '$w_h$ = %.1f' % df_ews.loc[species,d]['AIC hopf'],
            fontsize=9,
            color='r',
            transform=ax.transAxes)
    ax.text(xpos,ypos-2*0.12,
            '$w_n$ = %.1f' % df_ews.loc[species,d]['AIC null'],
            fontsize=9,
            color='g',
            transform=ax.transAxes)
# Y labels
for ax in axes[::3]:
    ax.set_ylabel('Power')
    
## Specific Y limits
#for ax in axes[:4]:
#    ax.set_ylim(top=0.004)
#for ax in axes[6:9]:
#    ax.set_ylim(top=0.25)
# Assign to plot label
pspec_plot_chlor=g



## Brachionus
species='Brachionus'
g = sns.FacetGrid(df_pspec.loc[species].reset_index(level=['Delta','Frequency']), 
                  col='Delta',
                  col_wrap=3,
                  sharey=False,
                  aspect=1.5,
                  height=1.8
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
    d=deltaVals[i]
#    ax.set_ylim(bottom=0, top=1.1*max(df_pspec.loc[species,d]['Empirical'].loc[xmin:xmax].dropna()))
    ax.set_xlim(left=xmin, right=xmax)
    ax.set_xticks([-3,-2,-1,0,1,2,3])
    ax.set_title('Delta = %.2f' % deltaVals[i])
    # AIC weights
    xpos=0.7
    ypos=0.9
    ax.text(xpos,ypos,
            '$w_f$ = %.1f' % df_ews.loc[species,d]['AIC fold'],
            fontsize=9,
            color='b',
            transform=ax.transAxes)  
    ax.text(xpos,ypos-0.12,
            '$w_h$ = %.1f' % df_ews.loc[species,d]['AIC hopf'],
            fontsize=9,
            color='r',
            transform=ax.transAxes)
    ax.text(xpos,ypos-2*0.12,
            '$w_n$ = %.1f' % df_ews.loc[species,d]['AIC null'],
            fontsize=9,
            color='g',
            transform=ax.transAxes)
# Y labels
for ax in axes[::3]:
    ax.set_ylabel('Power')
    
## Specific Y limits
#for ax in axes[:4]:
#    ax.set_ylim(top=0.004)
#for ax in axes[6:9]:
#    ax.set_ylim(top=0.25)
# Assign to plot label
pspec_plot_brach=g




#------------------------------------
## Export data / figures
#-----------------------------------


## Export EWS data

# Chlorella ews
df_ews_chlor = df_ews.loc['Chlorella']
df_ews_chlor.to_csv('data_export/'+dir_name+'/ews_chlor.csv')

# Brachionus EWS data
df_ews_brach = df_ews.loc['Brachionus']
df_ews_brach.to_csv('data_export/'+dir_name+'/ews_brach.csv')


## Export power spectrum (empirical data)

# Chlorella pspecs
df_pspec_chlor = df_pspec.loc['Chlorella','Empirical'].dropna()
df_pspec_chlor.to_csv('data_export/'+dir_name+'/pspec_chlor.csv')


# Brachionus pspecs
df_pspec_brach = df_pspec.loc['Brachionus', 'Empirical'].dropna()
df_pspec_brach.to_csv('data_export/'+dir_name+'/pspec_brach.csv')




