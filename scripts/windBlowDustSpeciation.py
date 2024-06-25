#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:08:54 2024

-------------------------windBlowDustSpeciation.py-----------------------------

classe utilizada para fazer a especiação química das emissões do windblowdust. 
Utilizamos o speciate para elabora a planilha weigth_perc_PM_CMAQ.csv com a 
especiação do material particulado emitido pelo solo.  


@author: leohoinaski

"""


import pandas as pd
import numpy as np


def speciate(windBlowDustFolder,FdustD):
    """
    função para a especiação química das emissões do windblowdust

    Parameters
    ----------
    windBlowDustFolder : path
        caminho para a pasta do módulo windblowdust.
    FdustD : np.array
        matriz com as emissões de partículas

    Returns
    -------
    FdustDNew : np.array
        matriz com as emissões especiadas.

    """
    print('=====STARTING windBlowDustSpeciation.py=====' )
    # abrindo csv com os perfis de especiação 
    spc = pd.read_csv(windBlowDustFolder+'/inputs/tables/weigth_perc_PM_CMAQ.csv')
    
    # usa todas as linhas que não tiver null
    spc = spc[~spc['SPECIES_NAME'].isnull()]
    
    # inicializa a matriz com as emissões especiadas
    FdustDNew = np.zeros([FdustD.shape[0],spc.shape[0],FdustD.shape[1],FdustD.shape[2]]) 
    
    # loop para cada espécie
    for index, row in spc.iterrows():
        
        print(index)
        # preenchendo a matriz
        FdustDNew[:,index,:,:] = FdustD*row['WP_MEAN']/100
        
        
    return FdustDNew