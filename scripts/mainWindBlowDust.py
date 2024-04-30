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
import windBlowDustSpeciation as wbds

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
  "Unit": '$\g.S^{-1}$',
  "tag":'PMULTRAFINE',
  "range":[0,1] # micrometers
}

ALL = {
  "Unit": '$\g.S^{-1}$',
  "tag":'AllFractions',
  "fractions":['PMFINE','PMC','PM10'] # micrometers
}

domain = 'd01'
GDNAM = 'BR_2019'
RESET_GRID = False
year = 2020


rootFolder =  os.path.dirname(os.path.dirname(os.getcwd()))
wrfoutFolder = rootFolder+'/BR_2019'
#wrfoutFolder='/home/lcqar/CMAQ_REPO/data/WRFout/BR/WRFd01_BR_20x20'
#wrfoutFolder='/home/WRFout/share/Congonhas/2021/d01'
#mcipPath='/home/artaxo/CMAQ_REPO/data/mcip/'+GDNAM
#mcipMETCRO3Dpath = mcipPath+'/METCRO3D_'+GDNAM+'.nc'
mcipMETCRO3Dpath = wrfoutFolder+'/METCRO3D_BR_2019.nc'
windBlowDustFolder = os.path.dirname(os.getcwd())
print(windBlowDustFolder)
#wrfoutFolder='/home/lcqar/CMAQ_REPO/data/WRFout/BR/WRFd01_BR_20x20'
#mcipMETCRO3Dpath ='/home/lcqar/CMAQ_REPO/data/mcip/BR_2019/METCRO3D_BR_2019.nc'


idSoils = [23,30,25] #4.1. Praia, Duna e Areal  4.3. Mineração 4.4. Outras Áreas não Vegetadas
dx = 0.1
Fractions = [PM25,PMC] # Lista com tipo de emissão por diâmetro. 
                        #Não precisa incluir o PM10 se já tiver PM25 e PM10

inputFolder = os.path.dirname(os.getcwd())+'/inputs'
print('Inputs folder = ' + inputFolder)
tablePath = os.path.dirname(os.getcwd())+'/inputs/tables'
print('Tables folder = ' + tablePath)
outfolder = os.path.dirname(os.getcwd())+'/Outputs/'+GDNAM
print('Outputs folder = ' + outfolder)


if os.path.isdir(outfolder):
    print('You have the outputs folder')
else:
    os.makedirs(outfolder, exist_ok=True)


ds = nc.Dataset(mcipMETCRO3Dpath)
datesTimeMCIP = ncCreate.datePrepCMAQ(ds)
# file = [i for i in os.listdir(wrfoutFolder) if os.path.isfile(os.path.join(wrfoutFolder,i)) and \
#          'wrfout_'+domain+'_'+str(datesTimeMCIP.year[0]).zfill(4)+'-'+\
#              str(datesTimeMCIP.month[0]).zfill(2)+'-'+\
#                  str(datesTimeMCIP.day[0]).zfill(2) in i]
file = [i for i in os.listdir(wrfoutFolder) if os.path.isfile(os.path.join(wrfoutFolder,i)) and \
         'wrfout_'+domain+'_'+str((datesTimeMCIP.datetime[0]- timedelta(days=1)).year).zfill(4)+'-'+\
             str((datesTimeMCIP.datetime[0]- timedelta(days=1)).month).zfill(2)+'-'+\
                 str((datesTimeMCIP.datetime[0]- timedelta(days=1)).day).zfill(2) in i]
wrfoutPath = wrfoutFolder+'/'+file[0]
ds = nc.Dataset(wrfoutPath)
datesTime = ncCreate.datePrepWRF(pd.to_datetime(wrf.extract_times(ds,wrf.ALL_TIMES)))
lia, loc = ismember.ismember(np.array(datesTime.datetime), np.array(datesTimeMCIP.datetime))
av,al,alarea,lat,lon,domainShp = regMap.main(wrfoutPath,GDNAM,inputFolder,outfolder,year,idSoils,RESET_GRID)

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
    ncCreate.createNETCDFtemporal(outfolder,'windBlowDust_',FdustD,datesTime[lia],mcipMETCRO3Dpath,EmisD)
    if EmisD==PM25:
        FdustFINE = FdustD
        FdustFINESpec = wbds.speciate(windBlowDustFolder, FdustFINE)
    elif EmisD==PMC:
        FdustCOARSE = FdustD
        FdustCOARSEpec = wbds.speciate(windBlowDustFolder, FdustFINE)
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
FdustSpeciated = FdustFINESpec + FdustCOARSEpec
ncCreate.createNETCDFtemporalSpeciated(windBlowDustFolder,outfolder,'windBlowDust_',FdustSpeciated,datesTime[lia],mcipMETCRO3Dpath)

# #%%
# import matplotlib.pyplot as plt
# import matplotlib.colors as colors
# shape_path= rootFolder+'/shapefiles/BR_regions.shp'   
# borderShape = gpd.read_file(shape_path)

# fig, ax = plt.subplots()
# ax.pcolor(lon,lat, np.nanmean(ustar[:, :, :].data,axis=0))
# borderShape.boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
# ax.set_title('ustar')

# fig, ax = plt.subplots()
# ax.pcolor(lon,lat, np.nanmean(ustarWRF[:, :, :].data,axis=0))
# borderShape.boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
# ax.set_title('ustarWRF')

# fig, ax = plt.subplots()
# ax.pcolor(lon,lat,np.nansum(avWRF,axis=0))
# borderShape.boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
# ax.set_title('avWRF')

# fig, ax = plt.subplots()
# ax.pcolor(lon,lat,sRef[0,:, :])
# borderShape.boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
# ax.set_title('sRef')

# fig, ax = plt.subplots()
# ax.pcolor(lon,lat,np.log(alarea[0,:, :]))
# borderShape.boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
# ax.set_title('alarea')

# fig, ax = plt.subplots()
# ax.pcolor(lon,lat,np.log(clayRegrid[0,:,:]))
# borderShape.boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
# ax.set_title('clayRegrid')

# fig, ax = plt.subplots()
# ax.pcolor(lon,lat,np.nanmean(Fvtot,axis=0))
# borderShape.boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
# ax.set_title('Fvtot')

# fig, ax = plt.subplots()
# ax.pcolor(lon,lat,np.nanmean(Fhtot,axis=0))
# borderShape.boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
# ax.set_title('Fhtot')

# fig, ax = plt.subplots()
# pcm = ax.pcolor(lon,lat,np.nansum(FdustD[:, :, :], axis=0),norm=colors.LogNorm())
# borderShape.boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
# ax.set_xticks([])
# ax.set_yticks([])
# cbar = fig.colorbar(pcm, ax=ax,fraction=0.04, pad=0.02,
#                         #extend='both', 
#                         #ticks=bounds,
#                         #spacing='uniform',
#                         orientation='horizontal',)

# cbar.ax.set_xlabel(EmisD['tag']+'\nWind blow Dust emission\n (g/s)', rotation=0,fontsize=8)
# ax.set_frame_on(False)
# cbar.ax.tick_params(labelsize=6) 
# fig.tight_layout()
# #ax.set_title('FdustD')

# fig, ax = plt.subplots()
# ax.scatter(ustarWRF.flatten(),FdustD.flatten())
# ax.set_title('FdustD vs ustar')
# ax.set_yscale('log')

# fig, ax = plt.subplots()
# ax.scatter(ustarWRF.flatten(),Fvtot.flatten())
# ax.set_title('Fvtot vs ustar')
# ax.set_yscale('log')

# fig, ax = plt.subplots()
# ax.scatter(ustarWRF.flatten(),Fhtot.flatten())
# ax.set_title('Fhtot vs ustar')
# ax.set_yscale('log')
