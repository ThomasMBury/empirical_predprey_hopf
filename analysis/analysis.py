#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 27 12:19:21 2018

@author: Thomas Bury

Analysis of Fussmann chemostat data using pandas

"""

# import python libraries
import numpy as np
import pandas as pd


# import data
raw = pd.read_excel('../data/raw_fussmann_2000.xls',header=[1])

# index by time (day number) and by meandelta (identifies experiment setting)
raw2 = raw.set_index(['meandelta','day#'])
                     
# 








