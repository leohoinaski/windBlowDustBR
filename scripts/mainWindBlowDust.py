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


rootFolder =  os.path.dirname(os.path.dirname(os.getcwd()))
wrfoutFolder='/media/leohoinaski/HDD/SC_2019'
#wrfoutFolder='/mnt/sdb1/SC_2019'
domain = 'd02'
mcipMETCRO3Dpath ='/media/leohoinaski/HDD/SC_2019/METCRO3D_SC_2019.nc'
#mcipMETCRO3Dpath ='/mnt/sdb1/SC_2019/METCRO3D_SC_2019.nc'

GRDNAM = 'SC_2019'
RESET_GRID = False
inputFolder = os.path.dirname(os.getcwd())+'/inputs'
tablePath = os.path.dirname(os.getcwd())+'/inputs/tables'
outfolder = os.path.dirname(os.getcwd())+'/outputs'
year = 2021
idSoils = [23,30,25] #4.1. Praia, Duna e Areal  4.3. Mineração 4.4. Outras Áreas não Vegetadas
EmisD = [1,2.5,10]

for ii, D in enumerate(EmisD):
    ds = nc.Dataset(mcipMETCRO3Dpath)
    datesTime = ncCreate.datePrepCMAQ(ds)
    file = [i for i in os.listdir(wrfoutFolder) if os.path.isfile(os.path.join(wrfoutFolder,i)) and \
             'wrfout_'+domain+'_'+str(datesTime.year[0]).zfill(4)+'-'+\
                 str(datesTime.month[0]).zfill(2)+'-'+\
                     str(datesTime.day[0]).zfill(2) in i][0]
    wrfoutPath = wrfoutFolder+'/'+file
    ds = nc.Dataset(wrfoutPath)
    datesTime = ncCreate.datePrepWRF(pd.to_datetime(wrf.extract_times(ds,wrf.ALL_TIMES)))
    av,al,alarea,lat,lon,domainShp = regMap.main(wrfoutPath,GRDNAM,inputFolder,outfolder,year,idSoils,RESET_GRID)
    clayRegrid,sRef = sp.main(inputFolder,outfolder,domainShp,GRDNAM,lat,lon,D,RESET_GRID)
    ustar,ustarT,ustarTd,avWRF = mp.main(wrfoutPath,tablePath,av,al,alarea,D,clayRegrid)
    Fdust = wbd.wbdFlux(avWRF,alarea,sRef,ustar,ustarT,ustarTd)
    ncCreate.createNETCDFtemporal(outfolder,'windBlowDust_',Fdust,datesTime,mcipMETCRO3Dpath,D)
    RESET_GRID = False
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots(4, 2)
# ax[0, 0].pcolor(lon,lat, np.nanmean(ustar[:, :, :],axis=0))
# ax[0, 1].pcolor(lon,lat,av[:, :])
# ax[1, 0].pcolor(lon,lat,sRef[:, :])
# ax[1, 1].pcolor(lon,lat,al[ :, :])
# ax[2, 0].pcolor(lon,lat,np.nanmean(Fhd,axis=0))
# ax[2, 1].pcolor(lon,lat,np.nanmean(Fvtot,axis=0))
# ax[3, 0].pcolor(lon,lat,np.nanmean(Fhtot,axis=0))
# ax[3, 1].pcolor(lon,lat,np.nansum(np.nansum(Fdust[:, :, :, :], axis=0), axis=0))
# for ii in range(0,4):
#     for jj in range(0,2):
#         shape_path= rootFolder+'/shapefiles/BR_regions.shp'   
#         borderShape = gpd.read_file(shape_path)
#         borderShape[borderShape['NM_MUN']=='Sul'].boundary.plot(edgecolor='black',linewidth=0.5,ax=ax[ii,jj])
# import matplotlib.pyplot as plt
# fig, ax = plt.subplots()
# ax.pcolor(lon,lat,np.log(np.nansum(np.nansum(Fdust[:, :, :, :], axis=0), axis=0)))
# borderShape[borderShape['NM_MUN']=='Sul'].boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
