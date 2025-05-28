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
import geopandas as gpd

def contribution_areas(grids,lat,lon):
    
    # abrindo o shp com regioes com ferro
    shp_iron = gpd.read_file('/home/lcqar/MMA/windBlowDustBR/mnt/sdb1/inputs/MiningBR/BRASIL.shp')
    shp_iron = shp_iron[shp_iron['SUBS'].isin(['FERRO', 'MINÉRIO DE FERRO'])]
    shp_iron = gpd.GeoDataFrame(geometry=[shp_iron.unary_union], crs=shp_iron.crs)
        
    # cria um gdf com os grids
    df = pd.DataFrame({'geometry':grids})
    gdf = gpd.GeoDataFrame(df, crs=shp_iron.crs)
    gdf.set_geometry('geometry', inplace=True)
    
    # calcula a porcentagem da área coberta
    gdf['prct_iron'] = gdf.geometry.intersection(shp_iron.geometry.iloc[0]).area / gdf.geometry.area
        
    array_iron = gdf['prct_iron'].to_numpy().reshape((lat.shape[1]-1,lon.shape[0]-1)).transpose()
    
    contribution = np.zeros((2,lat.shape[0]-1, lon.shape[1]-1))
    
    contribution[0,:,:] = 1-array_iron
    
    contribution[1,:,:] = array_iron
    
    return contribution

def speciate(windBlowDustFolder,FdustD,grids,lat,lon,contribution):
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
    spc = pd.read_csv(windBlowDustFolder+'/inputs/tables/weigth_perc_PM_All_CMAQ.csv')
    
    # usa todas as linhas que não tiver null
    spc = spc[~spc['SPECIES_NAME'].isnull()]
    
    # inicializa a matriz com as emissões especiadas
    FdustDNew = np.zeros([FdustD.shape[0],spc.shape[0],FdustD.shape[1],FdustD.shape[2]]) 
    
    if type(contribution) == list:
        contribution = contribution_areas(grids,lat,lon)
    
    # loop para cada espécie
    for index, row in spc.iterrows():
        
        lista_tipos = []
        
        for type_mine in range(contribution.shape[0]):
            
            lista_tipos.append(contribution[type_mine,:,:]*FdustD*(row[type_mine+7]/100))                                 
            
        FdustDNew[:,index,:,:] = np.nansum(np.stack(lista_tipos), axis=0)  
        
    return FdustDNew, contribution

