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
import os
import numpy as np
import geopandas as gpd

wrfoutPath='/media/leohoinaski/HDD/SC_2019/wrfout_d02_2019-01-03_18:00:00'
rootFolder =  os.path.dirname(os.path.dirname(os.getcwd()))
#wrfoutPath='/mnt/sdb1/SC_2019/wrfout_d02_2019-01-01'
GRDNAM = 'SC_2019'
inputFolder = os.path.dirname(os.getcwd())+'/inputs'
tablePath = os.path.dirname(os.getcwd())+'/inputs/tables'
outfolder = os.path.dirname(os.getcwd())+'/outputs'
year = 2021
idSoils = [23,24,30,25] #4.1. Praia, Duna e Areal  4.2. Área Urbanizada  4.3. Mineração 4.4. Outras Áreas não Vegetadas
D = 10
av,al,alarea,lat,lon,domainShp = regMap.main(wrfoutPath,GRDNAM,inputFolder,outfolder,year,idSoils)
clayRegrid,sRef = sp.main(inputFolder,outfolder,domainShp,GRDNAM,lat,lon,D)
ustar,ustarT,ustarTd,avWRF = mp.main(wrfoutPath,tablePath,av,al,alarea,D,clayRegrid)
Fdust = wbd.wbdFlux(avWRF,alarea,sRef,ustar,ustarT,ustarTd)

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
fig, ax = plt.subplots()
ax.pcolor(lon,lat,np.log(np.nansum(np.nansum(Fdust[:, :, :, :], axis=0), axis=0)))
borderShape[borderShape['NM_MUN']=='Sul'].boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)
