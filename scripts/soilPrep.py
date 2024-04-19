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

#import rasterio as rs
#import rasterio.mask
import pandas as pd
import numpy as np
import os
import netCDF4 as nc
#import geopandas as gpd
#from shapely.geometry import Point
#import nctoolkit as nctools
import rioxarray as riox
#from shapely.geometry import mapping
import regridMAPBIOMAS as regMap
import scipy

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
    raster = raster.rio.reproject("EPSG:4326")
    x = raster.x.values
    y = raster.y.values
    return x, y

def cutSoil(domainShp,inputFolder,outfolder,GRDNAM):
    from rasterio.enums import Resampling
    print('Using original raster')
    raster = riox.open_rasterio(inputFolder+'/br_clay_content_30-60cm_pred_g_kg/br_clay_content_30-60cm_pred_g_kg.tif')
    downscale_factor = 1/5
    new_width = raster.rio.width * downscale_factor
    new_height = raster.rio.height * downscale_factor
    raster = raster.rio.reproject(raster.rio.crs, shape=(int(new_height), int(new_width)), resampling=Resampling.bilinear)
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
        for jj in range(lon.shape[1]):
            #print(matArr[idr,idc].sum())
            print(str(ii)+' '+str(jj))
            if (np.size(np.where(np.array(latsIdx)==ii)[0])>0) and \
                (np.size(np.where(np.array(lonsIdx)==jj)[0])>0):
                matArr = raster[0,np.where(np.array(latsIdx)==ii)[0],np.where(np.array(lonsIdx)==jj)[0]].data
                matRegrid[ii,jj]=np.nanmedian(matArr)
                print('------------>Contain Clay')
            else:
                matRegrid[ii,jj]=0
    return matRegrid

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def soilType(inputFolder,lat,lon,D):

    raster = riox.open_rasterio(inputFolder+'/Solos_5000mil/SolosTextureRaster.tif', masked=True).squeeze()
    x = raster.x.values
    y = raster.y.values
   
    soilNames=['Sand', 'Silt', 'Clay']
    for kk,soiln in enumerate(soilNames):
        soilDist = pd.read_csv(inputFolder+'/tables/particleDist/'+soiln+'.csv')
        xs=soilDist['D']*1000
        ys=soilDist['P']
        f = lambda x,mu,sigma: scipy.stats.norm(mu,sigma).cdf(x)
        mu,sigma = scipy.optimize.curve_fit(f,xs,ys/100)[0]
        xx = np.arange(0,800,0.01)
        deriv = np.append(np.nan,np.diff(f(xx,mu,sigma)*100))
        idx = find_nearest(xx, D)
        raster.data[raster==kk+1]=deriv[idx]
    lonsIdx = []
    for i in x:
        idx = (np.abs(i - lon[0,:])).argmin()
        lonsIdx.append(idx)
    latsIdx = []
    for i in y:
        idx = (np.abs(i - lat[:,0])).argmin()
        latsIdx.append(idx)
    
    sRef=np.empty((1,lat.shape[0],lat.shape[1]))
    sRef[:,:,:] = np.nan
    for ii in range(0,lat.shape[0]):
        for jj in range(lon.shape[1]):
            #print(matArr[idr,idc].sum())
            print(str(ii)+' '+str(jj))
            if (np.size(np.where(np.array(latsIdx)==ii)[0])>0) and \
                (np.size(np.where(np.array(lonsIdx)==jj)[0])>0):
                matArr = raster[np.where(np.array(latsIdx)==ii)[0],np.where(np.array(lonsIdx)==jj)[0]].data
                sRef[0,ii,jj]=np.nanmedian(matArr)
            else:
                sRef[0,ii,jj]=np.nan

    return sRef
    
def main(inputFolder,outfolder,domainShp,GRDNAM,lat,lon,D,RESET_GRID):
    if os.path.exists(outfolder+'/regridClay_'+GRDNAM+'.nc'):
        if RESET_GRID==False:
            print ('You already have the regridClay_'+GRDNAM+'.nc file')
            ds = nc.Dataset(outfolder+'/regridClay_'+GRDNAM+'.nc')
            clayRegrid = ds['MAT'][:]
            if os.path.exists(outfolder+'/sRef_'+GRDNAM+'.nc') :
                print ('You already have the sRef_'+GRDNAM+'.nc file')
                ds = nc.Dataset(outfolder+'/sRef_'+GRDNAM+'.nc')
                sRef = ds['MAT'][:]
            else:
                sRef = soilType(inputFolder,lat,lon,D)
                regMap.createNETCDF(outfolder,'sRef_'+GRDNAM,sRef,lon,lat)
        else:
            # inputFolder = os.path.dirname(os.getcwd())+'/inputs'
            # outfolder = os.path.dirname(os.getcwd())+'/outputs'
            raster = cutSoil(domainShp,inputFolder,outfolder,GRDNAM)
            x, y = rasterLatLon(raster)
            clayRegrid = rasterInGrid(domainShp,raster,x,y,lat,lon)
            print('Creating netCDF')
            regMap.createNETCDF(outfolder,'regridClay_'+GRDNAM,clayRegrid,lon,lat)
            sRef = soilType(inputFolder,lat,lon,D)
            regMap.createNETCDF(outfolder,'sRef_'+GRDNAM,sRef,lon,lat)
    else:
        # inputFolder = os.path.dirname(os.getcwd())+'/inputs'
        # outfolder = os.path.dirname(os.getcwd())+'/outputs'
        raster = cutSoil(domainShp,inputFolder,outfolder,GRDNAM)
        x, y = rasterLatLon(raster)
        clayRegrid = rasterInGrid(domainShp,raster,x,y,lat,lon)
        print('Creating netCDF')
        regMap.createNETCDF(outfolder,'regridClay_'+GRDNAM,clayRegrid,lon,lat)
        sRef = soilType(inputFolder,lat,lon,D)
        regMap.createNETCDF(outfolder,'sRef_'+GRDNAM,sRef,lon,lat)
    return clayRegrid,sRef