#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 12:58:58 2024

@author: leohoinaski
"""

import rasterio as rs
import rasterio.mask
#import pandas as pd
import numpy as np
import os
import netCDF4 as nc
#import geopandas as gpd
#from shapely.geometry import Polygon
#import nctoolkit as nctools
import rioxarray as riox
#from shapely.geometry import mapping
import regridMAPBIOMAS as regMap

def rasterLatLon(outfolder,GRDNAM):
    raster = riox.open_rasterio(outfolder+'/clay_'+GRDNAM+'.tif')
    x = raster.x.values
    y = raster.y.values
    return x, y

def cutSoil(domainShp,inputFolder,outfolder,GRDNAM):
    # raster = riox.open_rasterio(inputFolder+'/br_clay_content_30-60cm_pred_g_kg/br_clay_content_30-60cm_pred_g_kg.tif', masked=True).squeeze()
    # raster = raster.rio.reproject('EPSG:4326')
    #raster = riox.open_rasterio(inputFolder+'/br_clay_content_30-60cm_pred_g_kg/br_clay_content_30-60cm_pred_g_kg.tif', masked=True)

    with rs.open(inputFolder+'/br_clay_content_30-60cm_pred_g_kg/br_clay_content_30-60cm_pred_g_kg.tif') as src:
        print('crop extent crs: ', domainShp.crs)
        print('raster crs: ', src.crs)
        domainShp5880  = domainShp.to_crs({'init':  src.crs})
        print('crop extent crs: ', domainShp5880.crs)
        try:
            out_image, out_transform = rasterio.mask.mask(src,domainShp5880.geometry,
                                                          crop=True)
            out_meta = src.meta
            out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})
        except:
            out_meta=None
            raster=[]
        if out_meta:   
            with rs.open(outfolder+'/clay_'+GRDNAM+'.tif', "w", **out_meta) as dest:
                dest.write(out_image)
            raster = riox.open_rasterio(outfolder+'/'+GRDNAM+'.tif', masked=True).squeeze()
            raster = raster.rio.reproject('EPSG:4326')    
            
    return out_meta,raster

def rasterInGrid(raster,x,y,lat,lon):
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
    matArr = raster.values.copy()
    for ii in range(0,lat.shape[0]):
        for jj in range(lon.shape[1]):
            idr,idc = np.meshgrid(np.where(np.array(latsIdx)==ii),np.where(np.array(lonsIdx)==jj))
            #print(matArr[idr,idc].sum())
            print(str(ii)+' '+str(jj))
            matRegrid[ii,jj]=np.nanmedian(matArr[idr,idc])
    return matRegrid

def main(inputFolder,outfolder,domainShp,GRDNAM,lat,lon):
    if os.path.exists(outfolder+'/regridClay_'+GRDNAM+'.nc'):
        print ('You already have the regridClay_'+GRDNAM+'.nc file')
        ds = nc.Dataset(outfolder+'/regridClay_'+GRDNAM+'.nc')
        clayRegrid = ds['MAT'][:]
    else:
        inputFolder = os.path.dirname(os.getcwd())+'/inputs'
        outfolder = os.path.dirname(os.getcwd())+'/outputs'
        GRDNAM = 'SC_2019'
        out_meta,raster = cutSoil(domainShp,inputFolder,outfolder,GRDNAM)
        x, y = rasterLatLon(outfolder,GRDNAM)
        clayRegrid = rasterInGrid(raster,x,y,lat,lon)
        regMap.createNETCDF(outfolder,'regridClay_'+GRDNAM,clayRegrid,lon,lat)
    return clayRegrid