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

wrfoutPath='/media/leohoinaski/HDD/SC_2019/wrfout_d02_2019-01-03_18:00:00'
rootFolder =  os.path.dirname(os.path.dirname(os.getcwd()))
#wrfoutPath='/mnt/sdb1/SC_2019/wrfout_d02_2019-01-01'
GRDNAM = 'SC_2019'
RESET_GRID = True
inputFolder = os.path.dirname(os.getcwd())+'/inputs'
tablePath = os.path.dirname(os.getcwd())+'/inputs/tables'
outfolder = os.path.dirname(os.getcwd())+'/outputs'
year = 2021
idSoils = [23,30,25] #4.1. Praia, Duna e Areal  4.3. Mineração 4.4. Outras Áreas não Vegetadas
EmisD = [1,2.5,10]

for D in EmisD:
    av,al,alarea,lat,lon,domainShp = regMap.main(wrfoutPath,GRDNAM,inputFolder,outfolder,year,idSoils,RESET_GRID)
    clayRegrid,sRef = sp.main(inputFolder,outfolder,domainShp,GRDNAM,lat,lon,D,RESET_GRID)
    ustar,ustarT,ustarTd,avWRF = mp.main(wrfoutPath,tablePath,av,al,alarea,D,clayRegrid)
    Fdust = wbd.wbdFlux(avWRF,alarea,sRef,ustar,ustarT,ustarTd)
    #ncCreate.createNETCDFtemporal(rootPath,folder,name,data,mcipMETCRO3Dpath,D)
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
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.pcolor(lon,lat,np.log(np.nansum(np.nansum(Fdust[:, :, :, :], axis=0), axis=0)))
borderShape[borderShape['NM_MUN']=='Sul'].boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
