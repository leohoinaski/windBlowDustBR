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
    # x = lon[0,:]
    # y = lat[:,0]
    # hlines = [((x1, yi), (x2, yi)) for x1, x2 in zip(x[:-1], x[1:]) for yi in y]
    # vlines = [((xi, y1), (xi, y2)) for y1, y2 in zip(y[:-1], y[1:]) for xi in x]
    # from shapely.ops import polygonize
    # from shapely.geometry import MultiLineString
    # grids = list(polygonize(MultiLineString(hlines + vlines)))
    pixelInRaster=[]
    matRegrid=np.empty((len(idSoils),lat.shape[0],lat.shape[1]))
    pixelsIn =[] 
    for kk, soilid in enumerate(idSoils):
        for i, shape in enumerate(grids):
            print(str(i) +' from ' + str(len(grids)))
            with rasterio.open(inputFolder+'/mapbiomas/brasil_coverage_'+str(year)+'.tif') as src:
                try:
                    out_image, out_transform = rasterio.mask.mask(src, [shape], crop=True)
                    #out_meta = src.meta
                    out_image = np.array(out_image)
                    out_image[out_image!=soilid]=0
                    out_image[out_image==soilid]=1
                    if out_image.sum()>0:
                            print('------------>Soil in grid')
                    pixelInRaster.append(out_image.sum())
                    pixelsIn.append(np.size(out_image))
                except:
                    print('Pixel outside raster')
                    pixelInRaster.append(np.nan)
                    pixelsIn.append(np.nan)
        matRegrid[kk,:,:] = np.array(pixelInRaster).reshape(lon.shape) 
        pixelsIn = np.array(pixelsIn).reshape(lon.shape) 
    
    matRegrid[np.isnan(matRegrid)]=0
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
    
def main(wrfoutPath,GDNAM,inputFolder,outfolder,year,idSoils,RESET_GRID): 
    if os.path.exists(outfolder+'/regridMAPBIOMAS_'+str(year)+'_'+GDNAM+'.nc'):
        if RESET_GRID==False:
            print ('You already have the regridMAPBIOMAS_'+str(year)+'_'+GDNAM+'.nc file')
            domainShp,lat,lon =  createDomainShp(wrfoutPath)
            ds = nc.Dataset(outfolder+'/regridMAPBIOMAS_'+str(year)+'_'+GDNAM+'.nc')
            #mapbioRegrid = ds['MAT'][0:len(idSoils)-1,:,:]
            al= ds['MAT'][0,:,:] 
            av= ds['MAT'][1,:,:] 
            alarea= ds['MAT'][(len(idSoils)+1):,:,:] 
        else:
            domainShp,lat,lon =  createDomainShp(wrfoutPath)
            out_meta,arr = cutMapbiomas(domainShp,inputFolder,outfolder,year,GDNAM)
            x, y = rasterLatLon(outfolder,GDNAM,inputFolder,year)
            mapbioRegrid,pixelsIn,av,al,alarea= rasterInGrid(arr,x,y,lat,lon,idSoils,year,inputFolder)
            matRegrid2=np.empty((2+alarea.shape[0],lat.shape[0],lat.shape[1]))
            matRegrid2[:,:,:] = np.nan
            matRegrid2[0,:,:] = al
            matRegrid2[1,:,:] = av
            matRegrid2[2:(2+alarea.shape[0]),:,:] = alarea
            #name = 'regridMAPBIOMAS_'+GRDNAM
            #data = matRegrid2
            createNETCDF(outfolder,'regridMAPBIOMAS_'+str(year)+'_'+GDNAM,matRegrid2,lon,lat)
            
    else:
        domainShp,lat,lon =  createDomainShp(wrfoutPath)
        out_meta,arr = cutMapbiomas(domainShp,inputFolder,outfolder,year,GDNAM)
        x, y = rasterLatLon(outfolder,GDNAM,inputFolder,year)
        mapbioRegrid,pixelsIn,av,al,alarea= rasterInGrid(arr,x,y,lat,lon,idSoils,year,inputFolder)
        matRegrid2=np.empty((2+alarea.shape[0],lat.shape[0],lat.shape[1]))
        matRegrid2[:,:,:] = np.nan
        matRegrid2[0,:,:] = al
        matRegrid2[1,:,:] = av
        matRegrid2[2:(2+alarea.shape[0]),:,:] = alarea
        #name = 'regridMAPBIOMAS_'+GRDNAM
        #data = matRegrid2
        createNETCDF(outfolder,'regridMAPBIOMAS_'+str(year)+'_'+GDNAM,matRegrid2,lon,lat)
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
