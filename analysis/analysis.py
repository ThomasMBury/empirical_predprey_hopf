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

# index by time (day number) and by meandelta (identifies experiment setting)
raw.set_index(['meandelta','day#'], inplace=True)
                     
# start each time-series from day 0
raw


# delta values in data
deltaVals = raw.index.levels[0]


                      
                      
                      








