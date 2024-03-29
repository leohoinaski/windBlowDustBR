#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 19:45:10 2024

inputs:  https://storage.googleapis.com/mapbiomas-public/initiatives/brasil/collection_8/lclu/coverage/brasil_coverage_2022.tif


https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2008JD011236


@author: leohoinaski
"""



import rasterio as rs
import rasterio.mask
#import pandas as pd
import numpy as np
import os
import netCDF4 as nc
import geopandas as gpd
from shapely.geometry import Polygon
#import nctoolkit as nctools
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
            
            if out_meta:   
                with rs.open(outfolder+'/mapbiomas_'+GRDNAM+'.tif', "w", **out_meta) as dest:
                    dest.write(out_image)
                with rs.open(outfolder+'/mapbiomas_'+GRDNAM+'.tif', 'r') as ds:
                    arr = ds.read()[0,:,:] 
        except:
            out_meta=src.meta
            arr = src 
    return out_meta,arr

def rasterLatLon(outfolder,GRDNAM,inputFolder,year):
    try:
        raster = riox.open_rasterio(outfolder+'/mapbiomas_'+GRDNAM+'.tif')
    except:
        raster = riox.open_rasterio(inputFolder+'/mapbiomas/brasil_coverage_'+str(year)+'.tif')
        
    x = raster.x.values
    y = raster.y.values
    return x, y

def rasterInGrid(arr,x,y,lat,lon,idSoils,year,inputFolder):
    if type(arr)==rasterio.io.DatasetReader:
        arr = riox.open_rasterio(inputFolder+'/mapbiomas/brasil_coverage_'+str(year)+'.tif')
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
        
        try:
            matArr = arr.copy()
            matArr[matArr!=soilid]=0
            matArr[matArr==soilid]=1
            #my_rast.where(my_rast != 1, 10)
            for ii in range(0,lat.shape[0]):
                for jj in range(lon.shape[1]):
                    idr,idc = np.meshgrid(np.where(np.array(latsIdx)==ii),np.where(np.array(lonsIdx)==jj))
                    #print(matArr[idr,idc].sum())
                    print(str(ii)+' '+str(jj))
                    matRegrid[kk,ii,jj]=matArr[idr,idc].sum()
                    pixelsIn[ii,jj] = np.size(matArr[idr,idc])
        except:
            #my_rast.where(my_rast != 1, 10)
            for ii in range(0,lat.shape[0]):
                matArr = arr[0,ii].data
                matArr[matArr!=soilid]=0
                matArr[matArr==soilid]=1
                for jj in range(lon.shape[1]):
                    idr,idc = np.meshgrid(np.where(np.array(latsIdx)==ii),np.where(np.array(lonsIdx)==jj))
                    #print(matArr[idr,idc].sum())
                    print(str(ii)+' '+str(jj))
                    matRegrid[kk,ii,jj]=matArr[idc].sum()
                    pixelsIn[ii,jj] = np.size(matArr[idc])
    av = (pixelsIn-np.nansum(matRegrid, axis=0))/pixelsIn
    al = np.nansum(matRegrid,axis=0)/pixelsIn
    alarea = matRegrid*30*30
    return matRegrid,pixelsIn,av,al,alarea


def createNETCDF(outfolder,name,data,xlon,ylat):
    print('===================STARTING createNETCDF.py=======================')
          
    f2 = nc.Dataset(outfolder+'/'+name+'.nc','w') #'w' stands for write 
    # # Specifying dimensions
    #tempgrp = f.createGroup('vehicularEmissions_data')
    if len(data.shape)>2:
        f2.createDimension('VAR', data.shape[0])
        f2.createDimension('ROW', data.shape[1])
        f2.createDimension('COL', data.shape[2])
    else:
        f2.createDimension('VAR', 1)
        f2.createDimension('ROW', data.shape[0])
        f2.createDimension('COL', data.shape[1])
    
    # Building variables
    # Passing data into variables
    LON = f2.createVariable('LON', 'f4', ( 'VAR','ROW','COL'))
    LAT = f2.createVariable('LAT', 'f4', ( 'VAR','ROW','COL'))
    LAT[:,:] =  ylat
    LON[:,:] = xlon
    LON.units = 'degrees '
    LAT.units = 'degrees '
    MAT = f2.createVariable('MAT', np.float32, ('VAR','ROW','COL'))
    if len(data.shape)>2:
        for ii in range(0,data.shape[0]):
            MAT[ii,:,:] = data[ii,:,:]
    else:
        MAT[:,:,:] = data[:,:]
    f2.close()
    return f2
    
def main(wrfoutPath,GRDNAM,inputFolder,outfolder,year,idSoils,RESET_GRID): 
    if os.path.exists(outfolder+'/regridMAPBIOMAS_'+GRDNAM+'.nc'):
        if RESET_GRID==False:
            print ('You already have the regridMAPBIOMAS_'+GRDNAM+'.nc file')
            domainShp,lat,lon =  createDomainShp(wrfoutPath)
            ds = nc.Dataset(outfolder+'/regridMAPBIOMAS_'+GRDNAM+'.nc')
            #mapbioRegrid = ds['MAT'][0:len(idSoils)-1,:,:]
            al= ds['MAT'][0,:,:] 
            av= ds['MAT'][1,:,:] 
            alarea= ds['MAT'][(len(idSoils)+1):,:,:] 
        else:
            domainShp,lat,lon =  createDomainShp(wrfoutPath)
            out_meta,arr = cutMapbiomas(domainShp,inputFolder,outfolder,year,GRDNAM)
            x, y = rasterLatLon(outfolder,GRDNAM,inputFolder,year)
            mapbioRegrid,pixelsIn,av,al,alarea= rasterInGrid(arr,x,y,lat,lon,idSoils,year,inputFolder)
            matRegrid2=np.empty((2+alarea.shape[0],lat.shape[0],lat.shape[1]))
            matRegrid2[:,:,:] = np.nan
            matRegrid2[0,:,:] = al
            matRegrid2[1,:,:] = av
            matRegrid2[2:(2+alarea.shape[0]),:,:] = alarea
            #name = 'regridMAPBIOMAS_'+GRDNAM
            #data = matRegrid2
            createNETCDF(outfolder,'regridMAPBIOMAS_'+GRDNAM,matRegrid2,lon,lat)
            
    else:
        domainShp,lat,lon =  createDomainShp(wrfoutPath)
        out_meta,arr = cutMapbiomas(domainShp,inputFolder,outfolder,year,GRDNAM)
        x, y = rasterLatLon(outfolder,GRDNAM,inputFolder,year)
        mapbioRegrid,pixelsIn,av,al,alarea= rasterInGrid(arr,x,y,lat,lon,idSoils,year,inputFolder)
        matRegrid2=np.empty((2+alarea.shape[0],lat.shape[0],lat.shape[1]))
        matRegrid2[:,:,:] = np.nan
        matRegrid2[0,:,:] = al
        matRegrid2[1,:,:] = av
        matRegrid2[2:(2+alarea.shape[0]),:,:] = alarea
        #name = 'regridMAPBIOMAS_'+GRDNAM
        #data = matRegrid2
        createNETCDF(outfolder,'regridMAPBIOMAS_'+GRDNAM,matRegrid2,lon,lat)
    return av,al,alarea,lat,lon,domainShp


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