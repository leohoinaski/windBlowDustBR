#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 12:58:58 2024

Inputs : http://geoinfo.cnps.embrapa.br/documents/3295
        https://geoftp.ibge.gov.br/informacoes_ambientais/pedologia/vetores/brasil_5000_mil/
        https://geo.anm.gov.br/portal/apps/webappviewer/index.html?id=6a8f5ccc4b6a4c2bba79759aa952d908
        
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

def rasterLatLon(domainShp,raster):
    df_mercator = domainShp.to_crs("epsg:3857")
    areaDomain = df_mercator.area/10**6
    #raster = riox.open_rasterio(outfolder+'/clay_'+GRDNAM+'.tif')
    if areaDomain[0]<2*10**6:
        x = raster.x.values
        y = raster.y.values
    else:
        raster = raster.rio.reproject("EPSG:4326")
        x = raster.x.values
        y = raster.y.values
        #xx,yy = np.meshgrid(x,y)
        # xnew=[]
        # ynew=[]
        # print('CONVERTING LATLON')
        # for ii,xi in enumerate(x):
        #     print(str(ii)+' of '+ str(np.size(x)))
        #     df1 = pd.DataFrame({'X':np.repeat(xi,np.size(y)), 'Y':y})
        #     df1['coords'] = list(zip(df1['X'], df1['Y']))
        #     df1['coords'] = df1['coords'].apply(Point)
        #     gdf1 = gpd.GeoDataFrame(df1, geometry='coords')
        #     gdf1.crs = "EPSG:5880"
        #     gdf1 = gdf1.to_crs(4326)
        #     xnew.append(gdf1.geometry.x)
        #     ynew.append(gdf1.geometry.y)
        # x = np.unique(np.array(xnew))
        # y = np.unique(np.array(ynew))
        
    return x, y

def cutSoil(domainShp,inputFolder,outfolder,GRDNAM):
    # raster = riox.open_rasterio(inputFolder+'/br_clay_content_30-60cm_pred_g_kg/br_clay_content_30-60cm_pred_g_kg.tif', masked=True).squeeze()
    # raster = raster.rio.reproject('EPSG:4326')
    #raster = riox.open_rasterio(inputFolder+'/br_clay_content_30-60cm_pred_g_kg/br_clay_content_30-60cm_pred_g_kg.tif', masked=True)
    df_mercator = domainShp.to_crs("epsg:3857")
    areaDomain = df_mercator.area/10**6
    if areaDomain[0]<2*10**6:
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
    else:
        from rasterio.enums import Resampling
        print('Using original raster')
        raster = riox.open_rasterio(inputFolder+'/br_clay_content_30-60cm_pred_g_kg/br_clay_content_30-60cm_pred_g_kg.tif')
        downscale_factor = 1/5
        #Caluculate new height and width using downscale_factor
        new_width = raster.rio.width * downscale_factor
        new_height = raster.rio.height * downscale_factor
        #downsample raster
        raster = raster.rio.reproject(raster.rio.crs, shape=(int(new_height), int(new_width)), resampling=Resampling.bilinear)
        # print(raster.rio.resolution(), down_sampled.rio.resolution())
        # # ((500.0, -500.0), (1000.4340277777777, -1000.0))
        # print(raster.shape, down_sampled.shape)
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
    df_mercator = domainShp.to_crs("epsg:3857")
    areaDomain = df_mercator.area/10**6
    matRegrid=np.empty((lat.shape[0],lat.shape[1]))
    matRegrid[:,:] = np.nan
    if areaDomain[0]<2*10**6:
        matArr = raster.values.copy()
        for ii in range(0,lat.shape[0]):
            for jj in range(lon.shape[1]):
                #idr,idc = np.meshgrid(np.where(np.array(latsIdx)==ii),np.where(np.array(lonsIdx)==jj))
                #print(matArr[idr,idc].sum())
                print(str(ii)+' '+str(jj))
                matRegrid[ii,jj]=np.nanmedian(matArr[np.where(np.array(latsIdx)==ii),
                                                     np.where(np.array(lonsIdx)==jj)])
                #matRegrid[ii,jj]=np.nanmedian(matArr[idr,idc])
    else:
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
    #del matArr,raster
    return matRegrid


# def log_interp1d(xx, yy, kind='linear'):
#     # from scipy.interpolate import UnivariateSpline
#     # spl = UnivariateSpline(xx, yy)
#     from scipy.interpolate import interp1d
#     spl = interp1d(xx, yy, kind='cubic',bounds_error=False)
#     # logx = np.log10(xx)
#     # logy = np.log10(yy)
#     # lin_interp = scipy.interpolate.interp1d(logx, logy, 
#     #                                         kind=kind,
#     #                                         fill_value='extrapolate')
#     # log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
#     return spl
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
        #ff = f(D,mu,sigma)
        deriv = np.append(np.nan,np.diff(f(xx,mu,sigma)*100))
        # plt.plot(xx,deriv,'-r')
        # plt.plot(xx,f(xx,mu,sigma)*100,'-b')
        # plt.plot(x,y,'og')
        # #cs = log_interp1d(x.astype(float),y.astype(float))
        # #fd1 = cs.derivative()
        idx = find_nearest(xx, D)
        raster.data[raster==kk+1]=deriv[idx]
        #pointInPolys['Sref'][np.array(pointInPolys['SoilType']).astype(str)==soiln]=ff
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
    # soil = gpd.read_file(inputFolder+'/Solos_5000mil/Solos_5000.shp').set_crs(4326, allow_override=True)
    # soil['SoilType'] = np.nan
    # soil['SoilType'][soil.DSC_TEXTUR.values.astype(str)=='arenosa'] = 'Sand'
    # soil['SoilType'][soil.DSC_TEXTUR.values.astype(str)=='media'] = 'Silt'
    # soil['SoilType'][soil.DSC_TEXTUR.values.astype(str)=='argilosa  ou muito argilosa'] = 'Clay'
    # soil['SoilTypeNumber'] = np.nan
    # soil['SoilTypeNumber'][soil.DSC_TEXTUR.values.astype(str)=='arenosa'] = 1
    # soil['SoilTypeNumber'][soil.DSC_TEXTUR.values.astype(str)=='media'] = 2
    # soil['SoilTypeNumber'][soil.DSC_TEXTUR.values.astype(str)=='argilosa  ou muito argilosa'] = 3
    # soil.to_file(inputFolder+'/Solos_5000mil/SolosTextures.shp', driver='ESRI Shapefile')
    # import matplotlib.pyplot as plt 
    # fig, ax = plt.subplots(figsize=(6.5, 4))
    # ax.pcolor(np.log(sRef[0,:,:]))
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
            inputFolder = os.path.dirname(os.getcwd())+'/inputs'
            outfolder = os.path.dirname(os.getcwd())+'/outputs'
            #GRDNAM = 'SC_2019'
            raster = cutSoil(domainShp,inputFolder,outfolder,GRDNAM)
            x, y = rasterLatLon(domainShp,raster)
            clayRegrid = rasterInGrid(domainShp,raster,x,y,lat,lon)
            print('Creating netCDF')
            regMap.createNETCDF(outfolder,'regridClay_'+GRDNAM,clayRegrid,lon,lat)
            sRef = soilType(inputFolder,lat,lon,D)
            regMap.createNETCDF(outfolder,'sRef_'+GRDNAM,sRef,lon,lat)
    else:
        inputFolder = os.path.dirname(os.getcwd())+'/inputs'
        outfolder = os.path.dirname(os.getcwd())+'/outputs'
        #GRDNAM = 'SC_2019'
        raster = cutSoil(domainShp,inputFolder,outfolder,GRDNAM)
        x, y = rasterLatLon(domainShp,raster)
        clayRegrid = rasterInGrid(domainShp,raster,x,y,lat,lon)
        print('Creating netCDF')
        regMap.createNETCDF(outfolder,'regridClay_'+GRDNAM,clayRegrid,lon,lat)
        sRef = soilType(inputFolder,lat,lon,D)
        regMap.createNETCDF(outfolder,'sRef_'+GRDNAM,sRef,lon,lat)
    return clayRegrid,sRef