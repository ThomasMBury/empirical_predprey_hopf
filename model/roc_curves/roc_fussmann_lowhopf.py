#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 11:58:49 2020

Compute roc curves for EWS in the Fussmann predator prey model (without bootstrapping)
EWS prior to low Hopf bifurcation
Export data in a df for plotting in MMA

@author: tbury
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sklearn.metrics as metrics



# Function to compute ROC data from truth and indicator vals
# and return a df.
def roc_compute(truth_vals, indicator_vals):
    
    # Compute ROC curve and threhsolds using sklearn
    fpr, tpr, thresholds = metrics.roc_curve(truth_vals,indicator_vals)
    
    # Compute AUC (area under curve)
    auc = metrics.auc(fpr, tpr)
    
    # Put into a DF
    dic_roc = {'fpr':fpr, 'tpr':tpr, 'thresholds':thresholds, 'auc':auc}
    df_roc = pd.DataFrame(dic_roc)

    return df_roc


# Import Kendall tau values for Fold, Flip and Null trajectories
stateVar = 'Chlorella'
columns=['Variable','Variance','Lag-1 AC',
         'Lag-3 AC','Lag-6 AC',
         'Kurtosis','Skewness',
         'Smax']
df_ktau_hopf = pd.read_csv('../data_export/ews_trans_lowhopf2/ktau.csv',usecols=columns)[columns]
df_ktau_null = pd.read_csv('../data_export/ews_null_lowhopf/ktau.csv',usecols=columns)[columns]

# Use Ktau for Chlorella
df_ktau_hopf = df_ktau_hopf[df_ktau_hopf['Variable']==stateVar]
df_ktau_null = df_ktau_null[df_ktau_null['Variable']==stateVar]


# Add column for truth value
df_ktau_hopf['Truth']=1
df_ktau_null['Truth']=0

# Combine dataframes to provide data to compute ROC
df_roc_data = pd.concat([df_ktau_hopf,df_ktau_null],axis=0,ignore_index=True)



## ROC curve data for EWS prior to Fold
df_roc_var = roc_compute(df_roc_data['Truth'],df_roc_data['Variance'])
df_roc_ac1 = roc_compute(df_roc_data['Truth'],df_roc_data['Lag-1 AC'])
df_roc_ac3 = roc_compute(df_roc_data['Truth'],df_roc_data['Lag-3 AC'])
df_roc_ac6 = roc_compute(df_roc_data['Truth'],df_roc_data['Lag-6 AC'])
df_roc_smax = roc_compute(df_roc_data['Truth'],df_roc_data['Smax'])


## Plot ROC curve for EWS prior to Hopf bifurcation

plt.title('EWS prior to H1')
plt.plot(df_roc_var['fpr'], df_roc_var['tpr'], 'b', label = 'Variance, AUC = %0.2f' % df_roc_var['auc'][0])
plt.plot(df_roc_ac6['fpr'], df_roc_ac6['tpr'], 'g', label = 'Lag-6 AC, AUC = %0.2f' % df_roc_ac6['auc'][0])
plt.plot(df_roc_smax['fpr'], df_roc_smax['tpr'], 'k', label = 'Smax, AUC = %0.2f' % df_roc_smax['auc'][0])
#plt.plot(df_roc_sd['fpr'], df_roc_sd['tpr'], 'c', label = 'Standard deviation, AUC = %0.2f' % df_roc_sd['auc'][0])

plt.legend(loc = 'lower right')
plt.plot([0, 1], [0, 1],'r--')
plt.xlim([-0.02, 1])
plt.ylim([0, 1.02])
plt.ylabel('True Positive Rate')
plt.xlabel('False Positive Rate')
plt.show()




### Concatenate ROC dataframes for export
#df_roc = pd.concat([df_roc_var, df_roc_ac1,df_roc_ac3, df_roc_ac6, df_roc_smax], axis=0)
#df_roc.to_csv('data/roc_plotdata.csv')

    






