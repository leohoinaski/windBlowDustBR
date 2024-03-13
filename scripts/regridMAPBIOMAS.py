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
                arr = ds.read()[0,:,:]  # read all raster values
                no_data=ds.nodata
                v = [arr[x,y] for x,y in np.ndindex(arr.shape)]
                lonMapbio = [i[0] for i in data]
                latMapbio = [i[1] for i in data]
              
    return out_meta,arr,lonMapbio,latMapbio

wrfoutPath='/media/leohoinaski/HDD/SC_2019/wrfout_d02_2019-01-03_18:00:00'
GRDNAM = 'SC_2019'
inputFolder = os.path.dirname(os.getcwd())+'/inputs'
outfolder = os.path.dirname(os.getcwd())+'/outputs'
year = 2021
idSoil = [23,30,25]

domainShp,lat,lon =  createDomainShp(wrfoutPath)
out_meta,arr,lonMapbio,latMapbio = cutMapbiomas(domainShp,inputFolder,outfolder,year,GRDNAM)

