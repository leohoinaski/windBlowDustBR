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
import geopandas as gpd
import netCDF4 as nc
import wrf
import pandas as pd
from datetime import timedelta
import ismember
rootFolder =  os.path.dirname(os.path.dirname(os.getcwd()))
wrfoutFolder=rootFolder+'/BR_2019'
#wrfoutFolder='/mnt/sdb1/BR_2019'
domain = 'd01'
mcipMETCRO3Dpath =rootFolder+'/BR_2019/METCRO3D_BR_2019.nc'
#mcipMETCRO3Dpath ='/mnt/sdb1/BR_2019/METCRO3D_BR_2019.nc'

GRDNAM = 'BR_2019'
RESET_GRID = False
inputFolder = os.path.dirname(os.getcwd())+'/inputs'
tablePath = os.path.dirname(os.getcwd())+'/inputs/tables'
outfolder = os.path.dirname(os.getcwd())+'/outputs'
year = 2021
idSoils = [23,30,25] #4.1. Praia, Duna e Areal  4.3. Mineração 4.4. Outras Áreas não Vegetadas
EmisD = [1,2.5,10]
#EmisD = [5]
dx = 0.5


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
av,al,alarea,lat,lon,domainShp = regMap.main(wrfoutPath,GRDNAM,inputFolder,outfolder,year,idSoils,RESET_GRID)

diam = np.arange(0,np.max(EmisD)+dx,dx)

for ii, D in enumerate(EmisD):
    FdustTotal=[]
    diamSelect = diam[D>=diam]
    for jj,diameters in enumerate(diamSelect):
        print(diameters)
        clayRegrid,sRef = sp.main(inputFolder,outfolder,domainShp,GRDNAM,lat,lon,diameters,RESET_GRID)
        ustar,ustarT,ustarTd,avWRF = mp.main(wrfoutPath,tablePath,av,al,diameters,clayRegrid,lia)
        Fdust = wbd.wbdFlux(avWRF,alarea,sRef,ustar,ustarT,ustarTd)
        RESET_GRID = False
        FdustTotal.append(Fdust)
    FdustTotal = np.stack(FdustTotal)
    FdustD = np.trapz(FdustTotal,dx=dx, axis=0)
    ncCreate.createNETCDFtemporal(outfolder,'windBlowDust_',FdustD,datesTime,mcipMETCRO3Dpath,D)

import matplotlib.pyplot as plt
fig, ax = plt.subplots(4, 2)
ax[0, 0].pcolor(lon,lat, np.nanmean(ustar[:, :, :],axis=0))
ax[0, 1].pcolor(lon,lat,av[:, :])
ax[1, 0].pcolor(lon,lat,sRef[:, :])
ax[1, 1].pcolor(lon,lat,alarea[0,:, :])
# ax[2, 0].pcolor(lon,lat,np.nanmean(Fhd,axis=0))
# ax[2, 1].pcolor(lon,lat,np.nanmean(Fvtot,axis=0))
# ax[3, 0].pcolor(lon,lat,np.nanmean(Fhtot,axis=0))
ax[3, 1].pcolor(lon,lat,np.nansum(FdustD[:, :, :], axis=0))
for ii in range(0,4):
    for jj in range(0,2):
        shape_path= rootFolder+'/shapefiles/BR_regions.shp'   
        borderShape = gpd.read_file(shape_path)
        borderShape[borderShape['NM_MUN']=='Sul'].boundary.plot(edgecolor='black',linewidth=0.5,ax=ax[ii,jj])
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.pcolor(lon,lat,np.log(np.nansum(FdustD[ :, :, :], axis=0)))
borderShape[borderShape['NM_MUN']=='Sul'].boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
