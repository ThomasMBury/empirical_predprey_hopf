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

# index by meandelta (identifies experiment setting)
raw.set_index(['meandelta','day#'], inplace=True)

# shift day# to start at 0

# function to take list and subtract minimum element
def zero_shift(array):
    return array-min(array)

# loop through delta values
deltaVals = raw.index.levels[0]
for d in deltaVals:
    # extract data for specific delta value
    df_temp = raw.loc[d]
    # remove time as index
    df_temp.reset_index(inplace=True)
    # shift time column to start at zero
    df_temp['day#'].apply(lambda x: x-min(df_temp['day#']))
    # put back in dataframe
    raw.loc[d] = df_temp


           
         



# delta values in data



                      
                      
                      








