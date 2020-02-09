#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 18:04:51 2019


Code to simulate predator-prey model for Fussmann study (Science 2000)
Tune the noise in the model to match that in the data

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
import ewstools


#----------------------
# Useful functions
#-----------------------

# Apply operation to column of DataFrame in place
def apply_inplace(df, field, fun):
    """ Apply function to a column of a DataFrame in place."""
    return pd.concat([df.drop(field, axis=1), df[field].apply(fun)], axis=1)





#--------------------------------
# Global parameters
#–-----------------------------


# Simulation parameters
dt = 0.01
t0 = 0
tmax = 500 # make large (to get idealised statistics from stationary distribution)
tburn = 200 # burn-in period
seed = 1 # random number generation seed
dbif1 = 0.150983 # first Hopf bifurcation (from Mma bif file)
dbif2 = 0.969538 # second Hopf bifurcation (from Mma bif file)



# EWS parameters
dt2 = 1 # spacing between time-series for EWS computation
rw = 1 # rolling window (compute EWS using full time-series)
bw = 1 # bandwidth (take the whole dataset as stationary)
lags = [1,2,3,6] # autocorrelation lag times
ews = ['var','ac'] # EWS to compute

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
d = 0.1 # dilution rate

# Parameter list
params = [ni,bc,kc,bb,kb,epsilon,m,lamda]

# Control parameter values
deltaVals = np.array([0.1,1.05])

# Noise parameter (demograhpic noise)
sigma = 0.05
sigma_n = sigma
sigma_c = sigma
sigma_r = sigma
sigma_b = sigma


# Initial conditions
n0 = 2
c0 = 5
r0 = 1
b0 = 2


#--------------------------------------------
# Simulate (stationary) realisations of model for each delta value
#-------------------------------------------


## Implement Euler Maryuyama for stocahstic simulation


# Set seed
np.random.seed(seed)


# Initialise a list to collect trajectories
list_traj_append = []

# loop over delta values
print('\nBegin simulations \n')
for d in deltaVals:

    
        
    # Initialise array to store time-series data
    t = np.arange(t0,tmax,dt) # Time array
    x = np.zeros([len(t), 4]) # State array
    
    
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
     
    # IC as a state vector
    x0 = np.array([n0,c0,r0,b0])
    
    # Run burn-in period on initial condition
    for i in range(int(tburn/dt)):
        # Update with Euler Maruyama
        x0 = x0 + de_fun(x0, d, params)*dt + np.dot(np.sqrt(x0),dW_burn[i])
        # Make sure that state variable remains >= 0 
        x0 = [np.max([k,0]) for k in x0]
        
        
    # Initial condition post burn-in period
    x[0]=x0
    
    # Run simulation
    for i in range(len(t)-1):
        x[i+1] = x[i] + de_fun(x[i], d, params)*dt + np.dot(np.sqrt(x[i]),dW[i])
        # make sure that state variable remains >= 0 
        x[i+1] = [np.max([k,0]) for k in x[i+1]]
            
    # Store series data in a DataFrame
    data = {'Delta': d,
                'Time': t,
                'Nitrogen': x[:,0],
                'Chlorella': x[:,1],
                'Reproducing Brachionus': x[:,2],
                'Brachionus': x[:,3]}
    df_temp = pd.DataFrame(data)
    # Append to list
    list_traj_append.append(df_temp)
    print('Simulation with d='+str(d)+' complete')
    

#  Concatenate DataFrame from each realisation
df_traj = pd.concat(list_traj_append)
df_traj.set_index(['Delta','Time'], inplace=True)


# Coarsen time-series to have spacing dt2 (for EWS computation)
df_traj_filt = df_traj.loc[::int(dt2/dt)]


# Make units of Chlor and Brach consistent with Fussmann experiments

nC = 50000/1000000 # conversion factor to 10^6 cells/ml of Chlorella
nB = 5 # conversion factor to females/ml of Brachiouns

df_traj_filt = apply_inplace(df_traj_filt, 'Chlorella',
                             lambda x: nC*x)
df_traj_filt = apply_inplace(df_traj_filt, 'Reproducing Brachionus',
                             lambda x: nB*x)
df_traj_filt = apply_inplace(df_traj_filt, 'Brachionus',
                             lambda x: nB*x)



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
        
        ews_dic = ewstools.core.ews_compute(df_traj_filt.loc[d][var], 
                          roll_window = rw, 
                          band_width = bw,
                          lag_times = lags, 
                          ews = ews
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

# Refine DataFrame to just have EWS data (no time dependence)
df_ews = df_ews_full.dropna().reset_index(level=2, drop=True).reorder_levels(['Variable', 'Delta'])

# Variance
df_ews['Variance']




  
##----------------
### Plots of EWS against delta value
##----------------
#
## Plot of EWS metrics
#fig1, axes = plt.subplots(nrows=1, ncols=1, sharex=True, figsize=(6,6))
#df_ews.loc['Chlorella'][['Variance']].plot(ax=axes[0],title='Early warning signals')
#df_ews.loc['Brachionus'][['Variance']].plot(ax=axes[0],secondary_y=True)
#

