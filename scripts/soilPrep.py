#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 12:58:58 2024

Inputs : http://geoinfo.cnps.embrapa.br/documents/3295
        https://geoftp.ibge.gov.br/informacoes_ambientais/pedologia/vetores/brasil_5000_mil/

@author: leohoinaski
"""

import rasterio as rs
import rasterio.mask
import pandas as pd
import numpy as np
import os
import netCDF4 as nc
import geopandas as gpd
from shapely.geometry import Point
#import nctoolkit as nctools
import rioxarray as riox
#from shapely.geometry import mapping
import regridMAPBIOMAS as regMap
import scipy

def rasterLatLon(raster):
    #raster = riox.open_rasterio(outfolder+'/clay_'+GRDNAM+'.tif')
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
            raster = riox.open_rasterio(outfolder+'/clay_'+GRDNAM+'.tif', masked=True).squeeze()
            raster = raster.rio.reproject('EPSG:4326')    
    del out_image, out_transform, out_meta
    raster = raster/1000
    return raster

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
    del matArr
    return matRegrid

def log_interp1d(xx, yy, kind='linear'):
    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = scipy.interpolate.interp1d(logx, logy, 
                                            kind=kind,
                                            fill_value='extrapolate')
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp

def soilType(inputFolder,lat,lon,D):
    points = gpd.GeoSeries(map(Point, zip(lon.flatten(), lat.flatten()))).set_crs(4326, allow_override=True)
    points = gpd.GeoDataFrame(geometry=points)
    soil = gpd.read_file(inputFolder+'/Solos_5000mil/Solos_5000.shp').set_crs(4326, allow_override=True)
    soil['SoilType'] = np.nan
    soil['SoilType'][soil.DSC_TEXTUR.values.astype(str)=='arenosa'] = 'Sand'
    soil['SoilType'][soil.DSC_TEXTUR.values.astype(str)=='media'] = 'Silt'
    soil['SoilType'][soil.DSC_TEXTUR.values.astype(str)=='argilosa  ou muito argilosa'] = 'Clay'
    pointInPolys = gpd.tools.sjoin(points, soil, how='left')
    soilDist = pd.read_csv(inputFolder+'/tables/soil_particleDist.csv')
    x=soilDist.columns[1:].values.astype(float)[::-1]
    soilNames=['Sand', 'Silt', 'Clay']
    pointInPolys['Sref'] = np.nan
    for soiln in soilNames:
        soilSect = soilDist[soilDist.Type==soiln]
        y=soilSect.iloc[0,1:][::-1].cumsum()
        cs = log_interp1d(x.astype(float),y.astype(float))
        pointInPolys['Sref'][np.array(pointInPolys['SoilType']).astype(str)==soiln]=cs(D)/100
    sRef =  np.reshape(pointInPolys['Sref'],lat.shape)
    #import matplotlib.pyplot as plt 
    # fig,ax = plt.subplots()
    # soil.plot(ax=ax)
    # points.plot(ax=ax)
    # fig, ax = plt.subplots(figsize=(6.5, 4))
    # ax.plot(x,y, 'o', label='data')
    # ax.plot(np.arange(1,1000,1),
    #         cs(np.arange(1,1000,1)), label="S")
    # fig, ax = plt.subplots(figsize=(6.5, 4))
    # ax.pcolor(sRef)
    return sRef
    
def main(inputFolder,outfolder,domainShp,GRDNAM,lat,lon,D,RESET_GRID):
    if os.path.exists(outfolder+'/regridClay_'+GRDNAM+'.nc'):
        if RESET_GRID==False:
            print ('You already have the regridClay_'+GRDNAM+'.nc file')
            ds = nc.Dataset(outfolder+'/regridClay_'+GRDNAM+'.nc')
            clayRegrid = ds['MAT'][:]
            sRef = soilType(inputFolder,lat,lon,D)
        else:
            inputFolder = os.path.dirname(os.getcwd())+'/inputs'
            outfolder = os.path.dirname(os.getcwd())+'/outputs'
            GRDNAM = 'SC_2019'
            raster = cutSoil(domainShp,inputFolder,outfolder,GRDNAM)
            x, y = rasterLatLon(raster)
            clayRegrid = rasterInGrid(raster,x,y,lat,lon)
            regMap.createNETCDF(outfolder,'regridClay_'+GRDNAM,clayRegrid,lon,lat)
            sRef = soilType(inputFolder,lat,lon,D)
    else:
        inputFolder = os.path.dirname(os.getcwd())+'/inputs'
        outfolder = os.path.dirname(os.getcwd())+'/outputs'
        GRDNAM = 'SC_2019'
        raster = cutSoil(domainShp,inputFolder,outfolder,GRDNAM)
        x, y = rasterLatLon(raster)
        clayRegrid = rasterInGrid(raster,x,y,lat,lon)
        regMap.createNETCDF(outfolder,'regridClay_'+GRDNAM,clayRegrid,lon,lat)
        sRef = soilType(inputFolder,lat,lon,D)
    return clayRegrid,sRef