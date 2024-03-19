#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 19:45:10 2024

@author: leohoinaski
"""



import rasterio as rs
import rasterio.mask
import pandas as pd
import numpy as np
import os
import netCDF4 as nc
import geopandas as gpd
from shapely.geometry import Polygon
import nctoolkit as nctools
import rioxarray as riox

def createDomainShp(wrfoutPath):
    ds = nc.Dataset(wrfoutPath)
    lat = ds['XLAT'][0,:,:]
    lon = ds['XLONG'][0,:,:]
    lat_point_list = [lat.min(), lat.max(), lat.max(), lat.min(), lat.min()]
    lon_point_list = [lon.min(), lon.min(), lon.max(), lon.max(), lon.min()]
    polygon_geom = Polygon(zip(lon_point_list, lat_point_list))
    domainShp = gpd.GeoDataFrame(index=[0], crs='epsg:4326', geometry=[polygon_geom])   
    return  domainShp,lat,lon
        
def cutMapbiomas(domainShp,inputFolder,outfolder,year,GRDNAM):
    with rs.open(inputFolder+'/mapbiomas/brasil_coverage_'+str(year)+'.tif') as src:
        try:
            out_image, out_transform = rasterio.mask.mask(src,domainShp.geometry,
                                                          crop=True)
            out_meta = src.meta
            out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})
            out_image[out_image==0]=26
        except:
            out_meta=None
            arr=[]
        if out_meta:   
            with rs.open(outfolder+'/'+GRDNAM+'.tif', "w", **out_meta) as dest:
                dest.write(out_image)
            with rs.open(outfolder+'/'+GRDNAM+'.tif', 'r') as ds:
                arr = ds.read()[0,:,:] 
    return out_meta,arr

def rasterLatLon(outfolder,GRDNAM):
    raster = riox.open_rasterio(outfolder+'/'+GRDNAM+'.tif')
    x = raster.x.values
    y = raster.y.values
    return x, y

def rasterInGrid(arr,x,y,lat,lon,idSoils):
    lonsIdx = []
    for i in x:
        idx = (np.abs(i - lon[0,:])).argmin()
        lonsIdx.append(idx)
    latsIdx = []
    for i in y:
        idx = (np.abs(i - lat[:,0])).argmin()
        latsIdx.append(idx)
      
    matRegrid=np.empty((len(idSoils),lat.shape[0],lat.shape[1]))
    pixelsIn=np.empty((lat.shape[0],lat.shape[1]))
    matRegrid[:,:,:] = np.nan
    for kk, soilid in enumerate(idSoils):
        matArr = arr.copy()
        matArr[matArr!=soilid]=0
        matArr[matArr==soilid]=1
        for ii in range(0,lat.shape[0]):
            for jj in range(lon.shape[1]):
                idr,idc = np.meshgrid(np.where(np.array(latsIdx)==ii),np.where(np.array(lonsIdx)==jj))
                #print(matArr[idr,idc].sum())
                print(str(ii)+' '+str(jj))
                matRegrid[kk,ii,jj]=matArr[idr,idc].sum()
                pixelsIn[ii,jj] = np.size(matArr[idr,idc])
    av = (pixelsIn-np.sum(matRegrid, axis=0))/pixelsIn
    al = matRegrid/pixelsIn
    return matRegrid,pixelsIn,av,al

def main(wrfoutPath,GRDNAM,inputFolder,outfolder,year,idSoils):  
    domainShp,lat,lon =  createDomainShp(wrfoutPath)
    out_meta,arr = cutMapbiomas(domainShp,inputFolder,outfolder,year,GRDNAM)
    x, y = rasterLatLon(outfolder,GRDNAM)
    matRegrid,pixelsIn,av,al= rasterInGrid(arr,x,y,lat,lon,idSoils)
    return matRegrid,av,al,lat,lon,domainShp


# wrfoutPath='/media/leohoinaski/HDD/SC_2019/wrfout_d02_2019-01-03_18:00:00'
# #wrfoutPath='/mnt/sdb1/SC_2019/wrfout_d02_2019-01-01'
# GRDNAM = 'SC_2019'
# inputFolder = os.path.dirname(os.getcwd())+'/inputs'
# outfolder = os.path.dirname(os.getcwd())+'/outputs'
# year = 2021
# idSoils = [23,24,30,25] #4.1. Praia, Duna e Areal  4.2. Área Urbanizada  4.3. Mineração 4.4. Outras Áreas não Vegetadas
# matRegrid,av,lat, lon = main(wrfoutPath,GRDNAM,inputFolder,outfolder,year,idSoils)
# import matplotlib.pyplot as plt 
# rootFolder =  os.path.dirname(os.path.dirname(os.getcwd()))
# shape_path= rootFolder+'/shapefiles/BR_regions.shp'   
# borderShape = gpd.read_file(shape_path)
# fig, ax = plt.subplots()
# plt.pcolor(lon,lat,av)
# borderShape[borderShape['NM_MUN']=='Sul'].boundary.plot(edgecolor='black',linewidth=0.5,ax=ax)