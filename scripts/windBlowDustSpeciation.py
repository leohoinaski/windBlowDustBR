#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:08:54 2024

@author: leohoinaski
"""
import pandas as pd
import os
import numpy as np

def speciate(windBlowDustFolder,FdustD):
    #windBlowDustFolder = os.path.dirname(os.path.dirname(os.getcwd()))+'/windBlowDustBR'
    spc = pd.read_csv(windBlowDustFolder+'/inputs/tables/weigth_perc_PM_CMAQ.csv')
    spc = spc[~spc['SPECIES_NAME'].isnull()]
    FdustDNew = np.zeros([FdustD.shape[0],spc.shape[0],FdustD.shape[1],FdustD.shape[2]]) 
    for index, row in spc.iterrows():
        print(index)
        FdustDNew[:,index,:,:] = FdustD*row['WP_MEAN']/100
        
    return FdustDNew