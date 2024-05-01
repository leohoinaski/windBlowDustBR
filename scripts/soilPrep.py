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
from scipy import optimize,stats
from shapely.geometry import Polygon
import rasterio as rs
import rasterio.mask

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

def soilType(inputFolder,lat,lon,D):

    raster = riox.open_rasterio(inputFolder+'/Solos_5000mil/SolosTextureRaster.tif', masked=True).squeeze()
    x = raster.x.values
    y = raster.y.values
    raster = raster.rio.reproject("EPSG:4326")
    rasterOriginal = raster.copy()
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
        raster.data[rasterOriginal==kk+1]=deriv[idx]
    
    lonCorner = np.append(np.append(lon[0,:-1]- np.diff(lon[0,:])/2,lon[0,-1]),
                          lon[0,-1]+np.diff(lon[0,-3:-1])/2)
    
    latCorner = np.append(np.append(lat[:-1,0]- np.diff(lat[:,0])/2,lat[-1,0]),
                          lat[-1,0]+np.diff(lat[-3:-1,0])/2)
    
    grids=[]
    for ii in range(1,lonCorner.shape[0]):
        #Loop over each cel in y direction
        for jj in range(1,latCorner.shape[0]):
            #roadClip=[]
            lat_point_list = [latCorner[jj-1], latCorner[jj], latCorner[jj], latCorner[jj-1]]
            lon_point_list = [lonCorner[ii-1], lonCorner[ii-1], lonCorner[ii], lonCorner[ii]]
            cel = Polygon(zip(lon_point_list, lat_point_list))
            grids.append(cel)
    
    pixelInRaster=[]
    sRef=np.empty((1,lat.shape[0],lat.shape[1]))
    sRef[:,:,:] = np.nan
    for i, shape in enumerate(grids):
        print(str(i) +' from ' + str(len(grids)))
        try:
            clipped = raster.rio.clip([shape])
            if np.nanmedian(clipped)>0:
                print('------------>Texture in grid')
                pixelInRaster.append(np.nanmedian(clipped))
                
        except:
            print('Pixel outside raster')
            pixelInRaster.append(np.nan)

    sRef[0,:,:] = np.array(pixelInRaster).reshape((lon.shape[1],lon.shape[0])).transpose() 
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
    if os.path.exists(outfolder+'/regridClay_'+GDNAM+'.nc'):
        if RESET_GRID==False:
            print ('You already have the regridClay_'+GDNAM+'.nc file')
            ds = nc.Dataset(outfolder+'/regridClay_'+GDNAM+'.nc')
            clayRegrid = ds['MAT'][:]
            if os.path.exists(outfolder+'/sRef_'+GDNAM+'.nc') :
                print ('You already have the sRef_'+GDNAM+'.nc file')
                ds = nc.Dataset(outfolder+'/sRef_'+GDNAM+'.nc')
                sRef = ds['MAT'][:]
            else:
                sRef = soilType(inputFolder,lat,lon,D)
                regMap.createNETCDF(outfolder,'sRef_'+GDNAM,sRef,lon,lat)
        else:
            # inputFolder = os.path.dirname(os.getcwd())+'/inputs'
            # outfolder = os.path.dirname(os.getcwd())+'/outputs'
            raster = cutSoil(domainShp,inputFolder,outfolder,GDNAM)
            x, y = rasterLatLon(raster)
            clayRegrid = rasterInGrid(domainShp,raster,x,y,lat,lon)
            print('Creating netCDF')
            regMap.createNETCDF(outfolder,'regridClay_'+GDNAM,clayRegrid,lon,lat)
            sRef = soilType(inputFolder,lat,lon,D)
            sRef[np.isnan(sRef)] = 0
            regMap.createNETCDF(outfolder,'sRef_'+GDNAM,sRef,lon,lat)
    else:
        # inputFolder = os.path.dirname(os.getcwd())+'/inputs'
        # outfolder = os.path.dirname(os.getcwd())+'/outputs'
        raster = cutSoil(domainShp,inputFolder,outfolder,GDNAM)
        x, y = rasterLatLon(raster)
        clayRegrid = rasterInGrid(domainShp,raster,x,y,lat,lon)
        print('Creating netCDF')
        regMap.createNETCDF(outfolder,'regridClay_'+GDNAM,clayRegrid,lon,lat)
        sRef = soilType(inputFolder,lat,lon,D)
        sRef[np.isnan(sRef)] = 0
        regMap.createNETCDF(outfolder,'sRef_'+GDNAM,sRef,lon,lat)
    return clayRegrid,sRef
