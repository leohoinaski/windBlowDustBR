#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 12:58:58 2024

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

def rasterLatLon(outfolder,GRDNAM):
    raster = riox.open_rasterio(outfolder+'/'+GRDNAM+'.tif')
    x = raster.x.values
    y = raster.y.values
    return x, y

def cutSoil(domainShp,inputFolder,outfolder,GRDNAM):
    raster = riox.open_rasterio(inputFolder+'/br_clay_content_30-60cm_pred_g_kg/br_clay_content_30-60cm_pred_g_kg.tif', masked=True).squeeze()
    raster_new = raster.rio.reproject(CRS.from_string('EPSG:4326'))
    
    with rs.open(inputFolder+'/br_clay_content_30-60cm_pred_g_kg/br_clay_content_30-60cm_pred_g_kg.tif') as src:
        try:
            out_image, out_transform = rasterio.mask.mask(src,domainShp.geometry,
                                                          crop=True)
            out_meta = src.meta
            out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})
        except:
            out_meta=None
            arr=[]
        if out_meta:   
            with rs.open(outfolder+'/'+GRDNAM+'.tif', "w", **out_meta) as dest:
                dest.write(out_image)
            with rs.open(outfolder+'/'+GRDNAM+'.tif', 'r') as ds:
                arr = ds.read()[0,:,:] 
    return out_meta,arr

inputFolder = os.path.dirname(os.getcwd())+'/inputs'
outfolder = os.path.dirname(os.getcwd())+'/outputs'
out_meta,arr = cutSoil(domainShp,inputFolder,outfolder,GRDNAM)
