#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 12:58:58 2024

This class is used to prepare the soil properties for windBlowDuscCalc

Inputs : http://geoinfo.cnps.embrapa.br/documents/3295
        https://geoftp.ibge.gov.br/informacoes_ambientais/pedologia/vetores/brasil_5000_mil/
        https://geo.anm.gov.br/portal/apps/webappviewer/index.html?id=6a8f5ccc4b6a4c2bba79759aa952d908
        
@author: leohoinaski
"""


import pandas as pd
import numpy as np
import os
import netCDF4 as nc
import rioxarray as riox
import regridMAPBIOMAS as regMap
from scipy import optimize,stats
from rasterio.enums import Resampling



def rasterLatLon(raster):
    """

    Parameters
    ----------

    raster : xarray
        rioxarray variable

    Returns
    -------
    x : longitude
        values depend on coordinate system
    y : latitude
        values depend on coordinate system.

    """
    
    # Reprojetando para o EPSG 4326
    raster = raster.rio.reproject("EPSG:4326")
    
    # Extraindo matriz de x
    x = raster.x.values
    
    # Extraindo matriz de y
    y = raster.y.values
    
    return x, y



def cutSoil(domainShp,inputFolder,outfolder,GRDNAM):
    """
    
    Parameters
    ----------
    domainShp : TYPE
        DESCRIPTION.
    inputFolder : TYPE
        DESCRIPTION.
    outfolder : TYPE
        DESCRIPTION.
    GRDNAM : TYPE
        DESCRIPTION.

    Returns
    -------
    raster : TYPE
        DESCRIPTION.

    """
    
    print('Starting cutSoil function - windBlowDust')
    
    # Abrindo arquivo com o teor de argila
    raster = riox.open_rasterio(inputFolder+'/br_clay_content_30-60cm_pred_g_kg/br_clay_content_30-60cm_pred_g_kg.tif')
    
    # Reduzindo a dimensão do raster 1/5
    downscale_factor = 1/5
    
    # nova largura e altura
    new_width = raster.rio.width * downscale_factor
    new_height = raster.rio.height * downscale_factor
    
    # faz o downscaling
    raster = raster.rio.reproject(raster.rio.crs, shape=(int(new_height), int(new_width)), resampling=Resampling.bilinear)
    
    # VERIFICAR!!! conversão da unidade de g/kg para % 
    # https://angeo.copernicus.org/articles/17/149/1999/angeo-17-149-1999.pdf
    # equação 4 :https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016MS000823
    raster = raster/1000 
    
    raster = raster.where(raster>0)
    
    return raster



def rasterInGrid(domainShp,raster,x,y,lat,lon):
    
    
    lonsIdx = []
    for i in x:
        idx = (np.abs(i - lon[0,:])).argmin()
        lonsIdx.append(idx)
    latsIdx = []
    for i in y:
        idx = (np.abs(i - lat[:,0])).argmin()
        latsIdx.append(idx)
    matRegrid=np.empty((lat.shape[0],lat.shape[1]))
    matRegrid[:,:] = np.nan
    raster = raster.rio.reproject("EPSG:4326")
    for ii in range(0,lat.shape[0]):
        for jj in range(0,lon.shape[1]):
            #print(matArr[idr,idc].sum())
            print(str(ii)+' '+str(jj))
            if (np.size(np.where(np.array(latsIdx)==ii)[0])>0) and \
                (np.size(np.where(np.array(lonsIdx)==jj)[0])>0):
                matArr = raster[0,np.where(np.array(latsIdx)==ii)[0],np.where(np.array(lonsIdx)==jj)[0]].data
                matArr[np.array(matArr==raster._FillValue)]=np.nan
                matRegrid[ii,jj]=np.nanmedian(matArr)
                print('------------>Contain Clay')
            else:
                matRegrid[ii,jj]=0
    matRegrid[np.isnan(matRegrid)] = 0
    return matRegrid

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def regridSoilTexture(outfolder,inputFolder,lat,lon,GDNAM):
    raster = riox.open_rasterio(inputFolder+'/Solos_5000mil/SolosTextureRaster.tif', masked=True).squeeze()
    x = raster.x.values
    y = raster.y.values
    raster = raster.rio.reproject("EPSG:4326")
    lonsIdx = []
    for i in x:
        idx = (np.abs(i - lon[0,:])).argmin()
        lonsIdx.append(idx)
    latsIdx = []
    for i in y:
        idx = (np.abs(i - lat[:,0])).argmin()
        latsIdx.append(idx)
    matRegrid=np.empty((1,lat.shape[0],lat.shape[1]))
    matRegrid[:,:,:] = 0
    raster = raster.rio.reproject("EPSG:4326")
    for ii in range(0,lat.shape[0]):
        for jj in range(0,lon.shape[1]):
            #print(matArr[idr,idc].sum())
            print(str(ii)+' '+str(jj))
            if (np.size(np.where(np.array(latsIdx)==ii)[0])>0) and \
                (np.size(np.where(np.array(lonsIdx)==jj)[0])>0):
                matArr = raster[np.where(np.array(latsIdx)==ii)[0],np.where(np.array(lonsIdx)==jj)[0]].data
                counts = np.bincount(matArr[~np.isnan(matArr)].astype(int))
                if len(counts)>0:
                    matRegrid[0,ii,jj]=np.argmax(counts)
                    print('------------>Texture type '+str(np.argmax(counts))+' '+str(counts))
            else:
                matRegrid[0,ii,jj]=0
    matRegrid[np.isnan(matRegrid)] = 0
    regMap.createNETCDF(outfolder,'regridedSoilTexture_'+GDNAM,matRegrid,lon,lat)
    return matRegrid

def soilType(inputFolder,outfolder,lat,lon,D,GDNAM):

    raster = nc.Dataset(outfolder+'/regridedSoilTexture_'+GDNAM+'.nc', masked=True)['MAT'][:]
    #x = raster.x.values
    #y = raster.y.values
    #raster = raster.rio.reproject("EPSG:4326")
    sRef=np.empty((lat.shape[0],lat.shape[1]))
    sRef[:,:] = 0
    soilNames=['Sand', 'Silt', 'Clay']
    for kk,soiln in enumerate(soilNames):
        soilDist = pd.read_csv(inputFolder+'/tables/particleDist/'+soiln+'.csv')
        xs=soilDist['D']*1000
        ys=soilDist['P']
        f = lambda x,mu,sigma: stats.norm(mu,sigma).cdf(x)
        mu,sigma = optimize.curve_fit(f,xs,ys/100)[0]
        xx = np.arange(0,800,0.01)
        deriv = np.append(np.nan,np.diff(f(xx,mu,sigma)*100))
        idx = find_nearest(xx, D)
        sRef[raster[0,:,:].astype(int)==kk+1]=deriv[idx]
    
    
    # lonsIdx = []
    # for i in x:
    #     idx = (np.abs(i - lon[0,:])).argmin()
    #     lonsIdx.append(idx)
    # latsIdx = []
    # for i in y:
    #     idx = (np.abs(i - lat[:,0])).argmin()
    #     latsIdx.append(idx)
    
    # sRef=np.empty((1,lat.shape[0],lat.shape[1]))
    # sRef[:,:,:] = np.nan
    # for ii in range(0,lat.shape[0]):
    #     for jj in range(0,lon.shape[1]):
    #         #print(matArr[idr,idc].sum())
    #         print(str(ii)+' '+str(jj))
    #         if (np.size(np.where(np.array(latsIdx)==ii)[0])>0) and \
    #             (np.size(np.where(np.array(lonsIdx)==jj)[0])>0):
    #             matArr = raster[np.where(np.array(latsIdx)==ii)[0],np.where(np.array(lonsIdx)==jj)[0]].data
    #             #matArr[np.array(matArr==raster._FillValue)]=np.nan
    #             print(np.nanmedian(matArr))
    #             sRef[0,ii,jj]=np.nanmedian(matArr)
    #         else:
    #             sRef[0,ii,jj]=np.nan

    return sRef
    
def main(inputFolder,outfolder,domainShp,GDNAM,lat,lon,D,RESET_GRID):
    if (os.path.exists(outfolder+'/regridClay_'+GDNAM+'.nc')) and (os.path.exists(
            outfolder+'/regridedSoilTexture_'+GDNAM+'.nc')):
        if RESET_GRID==False:
            print ('You already have the regridClay_'+GDNAM+'.nc file')
            ds = nc.Dataset(outfolder+'/regridClay_'+GDNAM+'.nc')
            clayRegrid = ds['MAT'][:]
            print ('You already have the regridedSoilTexture_'+GDNAM+'.nc file')
            ds = nc.Dataset(outfolder+'/regridedSoilTexture_'+GDNAM+'.nc')
            sRef = soilType(inputFolder,outfolder,lat,lon,D,GDNAM)
        else:
          # inputFolder = os.path.dirname(os.getcwd())+'/inputs'
          # outfolder = os.path.dirname(os.getcwd())+'/outputs'
          raster = cutSoil(domainShp,inputFolder,outfolder,GDNAM)
          x, y = rasterLatLon(raster)
          clayRegrid = rasterInGrid(domainShp,raster,x,y,lat,lon)
          print('Creating netCDF')
          regMap.createNETCDF(outfolder,'regridClay_'+GDNAM,clayRegrid,lon,lat)
          regridSoilTexture(outfolder,inputFolder,lat,lon,GDNAM)
          sRef = soilType(inputFolder,outfolder,lat,lon,D,GDNAM)
          sRef[np.isnan(sRef)] = 0  
    else:
        # inputFolder = os.path.dirname(os.getcwd())+'/inputs'
        # outfolder = os.path.dirname(os.getcwd())+'/outputs'
        raster = cutSoil(domainShp,inputFolder,outfolder,GDNAM)
        x, y = rasterLatLon(raster)
        clayRegrid = rasterInGrid(domainShp,raster,x,y,lat,lon)
        print('Creating netCDF')
        regMap.createNETCDF(outfolder,'regridClay_'+GDNAM,clayRegrid,lon,lat)
        regridSoilTexture(outfolder,inputFolder,lat,lon,GDNAM)
        sRef = soilType(inputFolder,outfolder,lat,lon,D,GDNAM)
        sRef[np.isnan(sRef)] = 0
    return clayRegrid,sRef
