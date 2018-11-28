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
               
               
## Some plots

# Plot all Chlorella trajectories
raw['Chlorella'].unstack(level=0).plot(
        title = 'Chlorella trajectories')

# Plot all Brachionus trajectories
raw['Brachionus'].unstack(level=0).plot(
        title = 'Brachionus trajectories')


                      
                      
                      








