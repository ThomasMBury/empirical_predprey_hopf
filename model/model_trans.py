#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16, 2019

@author: ThomasMBury

Code to simulate predator-prey model for Fussmann study (Science 2000)
Transient simulations (vary the dilution rate over time)

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
    return pd.concat([df.drop(field, axis=1), df[field].apply(fun)], axis=1)


#---------------------
# Directory for data output
#–----------------------

# Name of directory within data_export
dir_name = 'ews_1'

if not os.path.exists('data_export/'+dir_name):
    os.makedirs('data_export/'+dir_name)


#--------------------------------
# Global parameters
#–-----------------------------


# Simulation parameters
dt = 0.01
t0 = 0
tmax = 200
tburn = 100 # burn-in period
numSims = 2
seed = 1 # random number generation seed
dl = 0.08 # min value of dilution rate (control)
dh = 0.18 # max value of dilution rate (control)
dbif = 0.15 # bifurcation value


# EWS parameters
dt2 = 1 # spacing between time-series for EWS computation
rw = 0.25 # rolling window
bw = 0.1 # bandwidth
lags = [1,2,4] # autocorrelation lag times
ews = ['var','ac','sd','cv','skew','kurt','smax','aic','cf'] # EWS to compute
ham_length = 40 # number of data points in Hamming window
ham_offset = 0.5 # proportion of Hamming window to offset by upon each iteration
pspec_roll_offset = 10 # offset for rolling window when doing spectrum metrics


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
d=1.5 # dilution rate (control parameter)
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


# Noise parameters
sigma_n = 0 # amplitude for N
sigma_c = 0.02 # amplitude for Chlorella
sigma_r = 0 # amplitude for R
sigma_b = 0.02 # amplitude for Brachionus

# Initial conditions
n0 = 2
c0 = 5
r0 = 1
b0 = 2
# State vector
x0 = np.array([n0,c0,r0,b0])


#--------------------------------------------
# Simulate (transient) realisations of model
#-------------------------------------------



# Initialise array to store time-series data
t = np.arange(t0,tmax,dt) # Time array
x = np.zeros([len(t), 4]) # State array

# Set up control parameter d, that increases linearly in time from dl to dh
d = pd.Series(np.linspace(dl,dh,len(t)),index=t)
# Time at which bifurcation occurs
tbif = d[d > dbif].index[1]

## Implement Euler Maryuyama for stocahstic simulation


# Set seed
np.random.seed(seed)

# Initialise a list to collect trajectories
list_traj_append = []

# Loop over simulations
print('\nBegin simulations \n')
for j in range(numSims):
    
    
    # Create brownian increments (s.d. sqrt(dt))
    dW_n_burn = np.random.normal(loc=0, scale=sigma_n*np.sqrt(dt), size = int(tburn/dt))
    dW_n = np.random.normal(loc=0, scale=sigma_n*np.sqrt(dt), size = len(t)) 
    
    dW_c_burn = np.random.normal(loc=0, scale=sigma_c*np.sqrt(dt), size = int(tburn/dt))
    dW_c = np.random.normal(loc=0, scale=sigma_c*np.sqrt(dt), size = len(t))
  
    dW_r_burn = np.random.normal(loc=0, scale=sigma_r*np.sqrt(dt), size = int(tburn/dt))
    dW_r = np.random.normal(loc=0, scale=sigma_r*np.sqrt(dt), size = len(t))
    
    dW_b_burn = np.random.normal(loc=0, scale=sigma_b*np.sqrt(dt), size = int(tburn/dt))
    dW_b = np.random.normal(loc=0, scale=sigma_b*np.sqrt(dt), size = len(t))
    
    # Noise vectors
    dW_burn = np.array([dW_n_burn,
                        dW_c_burn,
                        dW_r_burn,
                        dW_b_burn]).transpose()
    
    dW = np.array([dW_n, dW_c, dW_r, dW_b]).transpose()
    
    # Run burn-in period on x0
    for i in range(int(tburn/dt)):        
        
        x0 = x0 + de_fun(x0, dl, params)*dt + dW_burn[i]
        
    # Initial condition post burn-in period
    x[0]=x0
    
    # Run simulation
    for i in range(len(t)-1):
        x[i+1] = x[i] + de_fun(x[i], d.iloc[i], params)*dt + dW[i]
        # make sure that state variable remains >= 0 
        x[i+1] = [np.max([k,0]) for k in x[i+1]]
            
    # Store series data in a DataFrame
    data = {'Realisation number': (j+1)*np.ones(len(t)),
                'Time': t,
                'Nitrogen': x[:,0],
                'Chlorella': x[:,1],
                'Reproducing Brachionus': x[:,2],
                'Brachionus': x[:,3]}
    df_temp = pd.DataFrame(data)
    # Append to list
    list_traj_append.append(df_temp)
    
    print('Simulation '+str(j+1)+' complete')

#  Concatenate DataFrame from each realisation
df_traj = pd.concat(list_traj_append)
df_traj.set_index(['Realisation number','Time'], inplace=True)





#----------------------
## Execute ews_compute for each realisation and each variable
#---------------------

# Filter time-series to have time-spacing dt2
df_traj_filt = df_traj.loc[::int(dt2/dt)]

# Set up a list to store output dataframes from ews_compute
# We will concatenate them at the end
appended_ews = []
appended_pspec = []

# loop through realisation number
print('\nBegin EWS computation\n')
for i in range(numSims):
    # loop through sate variable
    for var in ['Chlorella', 'Brachionus']:
        
        ews_dic = ews_compute(df_traj_filt.loc[i+1][var], 
                          roll_window = rw, 
                          band_width = bw,
                          lag_times = lags, 
                          ews = ews,
                          ham_length = ham_length,
                          ham_offset = ham_offset,
                          pspec_roll_offset = pspec_roll_offset,
                          upto=tbif)
        
        # The DataFrame of EWS
        df_ews_temp = ews_dic['EWS metrics']
        # The DataFrame of power spectra
        df_pspec_temp = ews_dic['Power spectrum']
        
        # Include a column in the DataFrames for realisation number and variable
        df_ews_temp['Realisation number'] = i+1
        df_ews_temp['Variable'] = var
        
        df_pspec_temp['Realisation number'] = i+1
        df_pspec_temp['Variable'] = var
                
        # Add DataFrames to list
        appended_ews.append(df_ews_temp)
        appended_pspec.append(df_pspec_temp)
        
    # Print status every realisation
    if np.remainder(i+1,1)==0:
        print('EWS for realisation '+str(i+1)+' complete')


# Concatenate EWS DataFrames. Index [Realisation number, Variable, Time]
df_ews = pd.concat(appended_ews).reset_index().set_index(['Realisation number','Variable','Time'])
# Concatenate power spectrum DataFrames. Index [Realisation number, Variable, Time, Frequency]
df_pspec = pd.concat(appended_pspec).reset_index().set_index(['Realisation number','Variable','Time','Frequency'])


## Compute ensemble statistics of EWS over all realisations (mean, pm1 s.d.)
#ews_names = ['Variance', 'Lag-1 AC', 'Lag-2 AC', 'Lag-4 AC', 'AIC fold', 'AIC hopf', 'AIC null', 'Coherence factor']
#
#df_ews_means = df_ews[ews_names].mean(level='Time')
#df_ews_deviations = df_ews[ews_names].std(level='Time')



#-------------------------
# Plots to visualise EWS
#-------------------------

# Realisation number to plot
plot_num = 1
var = 'Chlorella'
## Plot of trajectory, smoothing and EWS of var (x or y)
fig1, axes = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(6,6))
df_ews.loc[plot_num,var][['State variable','Smoothing']].plot(ax=axes[0],
          title='Early warning signals for a single realisation')
df_ews.loc[plot_num,var]['Variance'].plot(ax=axes[1],legend=True)
df_ews.loc[plot_num,var][['Lag-1 AC','Lag-2 AC','Lag-4 AC']].plot(ax=axes[1], secondary_y=True,legend=True)
df_ews.loc[plot_num,var]['Smax'].dropna().plot(ax=axes[2],legend=True)
df_ews.loc[plot_num,var]['Coherence factor'].dropna().plot(ax=axes[2], secondary_y=True, legend=True)
df_ews.loc[plot_num,var][['AIC fold','AIC hopf','AIC null']].dropna().plot(ax=axes[3],legend=True)


## Define function to make grid plot for evolution of the power spectrum in time
def plot_pspec_grid(tVals, plot_num, var):
    
    g = sns.FacetGrid(df_pspec.loc[plot_num,var].loc[t_display].reset_index(), 
                  col='Time',
                  col_wrap=3,
                  sharey=False,
                  aspect=1.5,
                  height=1.8
                  )

    g.map(plt.plot, 'Frequency', 'Empirical', color='k', linewidth=2)
    g.map(plt.plot, 'Frequency', 'Fit fold', color='b', linestyle='dashed', linewidth=1)
    g.map(plt.plot, 'Frequency', 'Fit hopf', color='r', linestyle='dashed', linewidth=1)
    g.map(plt.plot, 'Frequency', 'Fit null', color='g', linestyle='dashed', linewidth=1)
    # Axes properties 
    axes = g.axes
    # Set y labels
    for ax in axes[::3]:
        ax.set_ylabel('Power')
        # Set y limit as max power over all time
        for ax in axes:
            ax.set_ylim(top=1.05*max(df_pspec.loc[plot_num,var]['Empirical']), bottom=0)
#            ax.set_yscale('log')
    return g

#  Choose time values at which to display power spectrum
t_display = df_pspec.index.levels[2][::1].values

plot_pspec_x = plot_pspec_grid(t_display,1,'Chlorella')
#plot_pspec_y = plot_pspec_grid(t_display,1,'y')


##------------------------------------
### Export data / figures
##-----------------------------------
#
## Export power spectrum evolution (grid plot)
#plot_pspec_x.savefig('figures/pspec_evol_x.png', dpi=200)
#plot_pspec_y.savefig('figures/pspec_evol_y.png', dpi=200)
#
### Export the first 5 realisations to see individual behaviour
## EWS DataFrame (includes trajectories)
#df_ews.loc[:5].to_csv('data_export/'+dir_name+'/ews_singles.csv')
## Power spectrum DataFrame (only empirical values)
#df_pspec.loc[:5,'Empirical'].dropna().to_csv('data_export/'+dir_name+'/pspecs.csv',
#            header=True)
#
## Export ensemble statistics
#df_ews_means.to_csv('data_export/'+dir_name+'/ews_ensemble_mean.csv')
#df_ews_deviations.to_csv('data_export/'+dir_name+'/ews_ensemble_std.csv')








