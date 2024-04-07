#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 16:04:37 2024

Referências:
    https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016MS000823
    https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2010JD014649

@author: leohoinaski
"""

import regridMAPBIOMAS as regMap
import soilPrep as sp
import metPrep as mp
import windBlowDustCalc as wbd
import netCDFcreator as ncCreate
import os
import numpy as np
#import geopandas as gpd
import netCDF4 as nc
import wrf
import pandas as pd
from datetime import timedelta
import ismember
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', default=0, action='count')
    parser.add_argument('windBlowDustFolder')
    parser.add_argument('mcipPath')
    parser.add_argument('wrfoutFolder')
    parser.add_argument('domain')
    parser.add_argument('GDNAM')
    parser.add_argument('YEAR', type=int)
    parser.add_argument('RESET_GRID', type=int) # 0 OU 1
    args = parser.parse_args()
    
    windBlowDustFolder = args.windBlowDustFolder
    mcipPath = args.mcipPath
    wrfoutFolder = args.wrfoutFolder
    domain = args.domain
    domain = 'd0'+domain
    GDNAM  = args.GDNAM
    YEAR = args.YEAR
    RESET_GRID = args.RESET_GRID
    mcipMETCRO3Dpath = mcipPath+'/METCRO3D_'+GDNAM+'.nc'
    if RESET_GRID == 0:
        RESET_GRID=False
    else:
        RESET_GRID=True
    
    PM25 = {
      "Unit": '$\g.S^{-1}$',
      "tag":'PMFINE',
      "range":[0,2.5] # micrometers
    }
    
    PMC = {
      "Unit": '$\g.S^{-1}$',
      "tag":'PMC',
      "range":[2.5,10] # micrometers
    }
    
    PM10 = {
      "Unit": '$\g.S^{-1}$',
      "tag":'PM10',
      "range":[0,10] # micrometers
    }
    
    PM1 = {
      "Pollutant": "$NO_{2}$",
      "Unit": '$\g.S^{-1}$',
      "tag":'PMULTRAFINE',
      "range":[0,1] # micrometers
    }
    
    ALL = {
      "Unit": '$\g.S^{-1}$',
      "tag":'AllFractions',
      "fractions":['PMFINE','PMC','PM10'] # micrometers
    }


    idSoils = [23,30,25] #4.1. Praia, Duna e Areal  4.3. Mineração 4.4. Outras Áreas não Vegetadas
    dx = 0.1 # dx para integração das emissões das frações.
    Fractions = [PM25,PMC] # Lista com tipo de emissão por diâmetro. 
                            #Não precisa incluir o PM10 se já tiver PM25 e PM10
       
    inputFolder = windBlowDustFolder+'/inputs'
    tablePath = windBlowDustFolder+'/inputs/tables'
    outfolder = windBlowDustFolder+'/Outputs/'+GDNAM
    
    if os.path.isdir(outfolder):
        print('You have the outputs folder')
    else:
        os.makedirs(outfolder, exist_ok=True)
    
    ds = nc.Dataset(mcipMETCRO3Dpath)
    datesTimeMCIP = ncCreate.datePrepCMAQ(ds)
    file = [i for i in os.listdir(wrfoutFolder) if os.path.isfile(os.path.join(wrfoutFolder,i)) and \
             'wrfout_'+domain+'_'+str((datesTimeMCIP.datetime[0]- timedelta(days=1)).year).zfill(4)+'-'+\
                 str((datesTimeMCIP.datetime[0]- timedelta(days=1)).month).zfill(2)+'-'+\
                     str((datesTimeMCIP.datetime[0]- timedelta(days=1)).day).zfill(2) in i]
    wrfoutPath = wrfoutFolder+'/'+file[0]
    ds = nc.Dataset(wrfoutPath)
    datesTime = ncCreate.datePrepWRF(pd.to_datetime(wrf.extract_times(ds,wrf.ALL_TIMES)))
    lia, loc = ismember.ismember(np.array(datesTime.datetime), np.array(datesTimeMCIP.datetime))
    av,al,alarea,lat,lon,domainShp = regMap.main(wrfoutPath,GDNAM,inputFolder,outfolder,YEAR,idSoils,RESET_GRID)
    print(datesTime.iloc[lia,:])
    for EmisD  in Fractions:
        Dmax = np.max(EmisD['range'])
        Dmin = np.min(EmisD['range'])
        diam = np.arange(np.min(EmisD['range']),np.max(EmisD['range'])+dx,dx)
        FdustTotal=[]
        diamSelect = diam[(Dmin<=diam) & (Dmax>=diam)]
        for jj,diameters in enumerate(diamSelect):
            print(diameters)
            clayRegrid,sRef = sp.main(inputFolder,outfolder,domainShp,GDNAM,lat,lon,diameters,RESET_GRID)
            ustar,ustarT,ustarTd,avWRF,ustarWRF = mp.main(wrfoutPath,tablePath,av,al,diameters,clayRegrid,lia)
            Fdust,Fhd,Fhtot,Fvtot = wbd.wbdFlux(avWRF,alarea,sRef,ustarWRF,ustarT,ustarTd)
            RESET_GRID = False
            FdustTotal.append(Fdust)
        FdustTotal = np.stack(FdustTotal)
        FdustD = np.trapz(FdustTotal,dx=dx, axis=0)*1000 # converte para g VERIFICAR!
        ncCreate.createNETCDFtemporal(outfolder,'windBlowDust_',FdustD,datesTime.iloc[lia,:],mcipMETCRO3Dpath,EmisD)
        if EmisD==PM25:
            FdustFINE = FdustD
        elif EmisD==PMC:
            FdustCOARSE = FdustD
        else:
            print('You have selected an awkward fraction')
        try:
            FdustPM10 = FdustFINE+FdustCOARSE
            ncCreate.createNETCDFtemporal(outfolder,'windBlowDust_',FdustPM10,datesTime[lia],mcipMETCRO3Dpath,PM10)
        except:
            print('You do not have the fractions required for PM10')   
            
    
    FdustALL = [FdustFINE,FdustCOARSE,FdustPM10]
    FdustALL = np.stack(FdustALL,axis=0)
    ncCreate.createNETCDFtemporal(outfolder,'windBlowDust_',FdustALL,datesTime[lia],mcipMETCRO3Dpath,ALL)
    
