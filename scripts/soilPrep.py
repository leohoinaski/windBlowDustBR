#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Created on Tue Mar 19 12:58:58 2024

This class is used to prepare the soil properties for windBlowDuscCalc

Inputs : http://geoinfo.cnps.embrapa.br/documents/3295
        https://geoftp.ibge.gov.br/informacoes_ambientais/pedologia/vetores/brasil_5000_mil/
        https://geo.anm.gov.br/portal/apps/webappviewer/index.html?id=6a8f5ccc4b6a4c2bba79759aa952d908
        
        Vasques, G.M., Coelho, M.R., Dart, R.O., Cintra, L.C., Baca, J.F.M. 
        (2021). Soil Clay, Silt and Sand Content Maps for Brazil at 0-5, 5-15,
        15-30, 30-60, 60-100 and 100-200 cm Depth Intervals with 90 m Spatial 
        Resolution. Version 2021. Embrapa Solos, Rio de Janeiro, Brazil.
        
@author: leohoinaski


"""

import pandas as pd
import numpy as np
import os
import netCDF4 as nc
import rioxarray as riox
import regridMAPBIOMAS as regMap
from scipy import optimize,stats
from rasterio.enums import Resampling
from shapely.geometry import Polygon
import rasterio
import geopandas as gpd


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
    
    # Reprojetando para o EPSG 4326
    raster = raster.rio.reproject("EPSG:4326")
    
    # Extraindo matriz de x
    x = raster.x.values
    
    # Extraindo matriz de y
    y = raster.y.values
    
    return x, y



def cutSoil(domainShp,inputFolder,outfolder,GRDNAM):
    """
    
    Esta função é utilizada para fazer reduzir a resolução do dado de clay 
    content. 
    
    Parameters
    ----------
    domainShp : TYPE
        DESCRIPTION.
    inputFolder : TYPE
        DESCRIPTION.
    outfolder : TYPE
        DESCRIPTION.
    GRDNAM : TYPE
        DESCRIPTION.

    Returns
    -------
    raster : TYPE
        DESCRIPTION.

    """
    
   # print('Starting cutSoil function - windBlowDust')
    
    # Abrindo arquivo com o teor de argila
    raster = riox.open_rasterio(inputFolder+'/br_clay_content_30-60cm_pred_g_kg/br_clay_content_30-60cm_pred_g_kg.tif')
    
    # Reduzindo a dimensão do raster 1/5
    downscale_factor = 1/5
    
    # nova largura e altura
    new_width = raster.rio.width * downscale_factor
    new_height = raster.rio.height * downscale_factor
    
    # faz o downscaling
    raster = raster.rio.reproject(raster.rio.crs, shape=(int(new_height), 
                                                         int(new_width)), 
                                  resampling=Resampling.bilinear)
    
    # VERIFICAR!!! conversão da unidade de g/kg para % 
    # https://angeo.copernicus.org/articles/17/149/1999/angeo-17-149-1999.pdf
    # equação 4 :https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016MS000823
    # dividido por 1000 para transformar g em kg e vezes 100 para colocar em %
    raster = (raster/1000)*100  
    
    # remove os valores faltantes iguais a -999
    raster = raster.where(raster>0)
    
    return raster


def rasterInGrid(domainShp,raster,x,y,lat,lon):
    """
    
    Esta função faz o regrid da matriz de clay content para o domínio de 
    modelagem.

    Parameters
    ----------
    domainShp : geodataframe
        geodataframe com o poligono do dominio de modelagem.
    raster : rioxarray
        matriz com os dados do raster clay content.
    x : np.array
        matriz de x do raster clay content.
    y : np.array
        matriz de y do raster clay content.
    lat : np.array
        matriz de latitudes do domínio.
    lon : np.array
        matriz de longitudes do domínio.

    Returns
    -------
    matRegrid : np.array
        matriz com o regrid do clay content para o domínio de modelagem.

    """
    
    # Inicializando matriz de longitudes dentro da cada célula do domínio
    lonsIdx = []
    
    # loop para cada valor de longitude do clay content
    for i in x:
        
        # indice que possui da celula que possui a longitude do raster
        idx = (np.abs(i - lon[0,:])).argmin()
        lonsIdx.append(idx)
        
    # Inicializando matriz de latitudes dentro da cada célula do domínio    
    latsIdx = []
    
    # loop para cada valor de longitude do clay content
    for i in y:
        
        # indice que possui da celula que possui a latitude do raster
        idx = (np.abs(i - lat[:,0])).argmin()
        latsIdx.append(idx)
    
    # Incializando a matriz de regrid do clay content
    matRegrid=np.empty((lat.shape[0],lat.shape[1]))
    # Todos os valores como nan
    matRegrid[:,:] = np.nan
    # Definindo o EPSG
    raster = raster.rio.reproject("EPSG:4326")
    
    # loop em cada latitude do dominio
    for ii in range(0,lat.shape[0]):
        
        # loop em cada longitude do domínio
        for jj in range(0,lon.shape[1]):
            
            #print(matArr[idr,idc].sum())
            # acha os indices da matriz do raster que estão dentro da célula do 
            #domínio
            print(str(ii)+' '+str(jj))
            if (np.size(np.where(np.array(latsIdx)==ii)[0])>0) and \
                (np.size(np.where(np.array(lonsIdx)==jj)[0])>0):
                    
                # define uma matriz de pixels do raster que estão dentro
                # da respectiva célula
                matArr = raster[0,np.where(np.array(latsIdx)==ii)[0],
                                np.where(np.array(lonsIdx)==jj)[0]].data
                
                # nan para quando a matriz tiver valores faltantes
                matArr[np.array(matArr==raster._FillValue)]=np.nan
               
                # preenche a nova matriz (dimensão igual ao domínio) com valores
                # da médiana do clay content na célula
                matRegrid[ii,jj]=np.nanmedian(matArr)
                print('------------>Contain Clay')
            
            # se não tiver pixels dentro da célula...    
            else:
                matRegrid[ii,jj]=0
                
    # substitui nan por 0
    matRegrid[np.isnan(matRegrid)] = 0
    
    return matRegrid

def find_nearest(array, value):
    """
    Encontra o indice da matriz que possui o valor mais próximo.

    Parameters
    ----------
    array : np.array
        array com matriz de valores.
    value : float
        valor para ser encontrado na matriz.

    Returns
    -------
    idx : int
        índice da matriz.

    """
    # transforma em np.array
    array = np.asarray(array)
    
    # Acha o indice
    idx = (np.abs(array - value)).argmin()
    
    return idx

def regridSoilTexture(outfolder,inputFolder,lat,lon,GDNAM):
    """
    

    Parameters
    ----------
    outfolder : path
        caminho para a pasta de outputs.
    inputFolder : path
        caminho para a pasta de inputs.
    lat : np.array
        matriz de latitudes do domínio.
    lon : np.array
        matriz de longitudes do domínio.
    GDNAM : str
        nome do domínio conforme o MCIP.

    Returns
    -------
    matRegrid : np.array
        matriz com o regrid da textura do solo.

    """
    
  
    # Extraindo os cantos das longitudes e latitudes
    lonCorner = np.append(np.append(lon[0,:-1]- np.diff(lon[0,:])/2,lon[0,-1]),
                          lon[0,-1]+np.diff(lon[0,-3:-1])/2)
    
    latCorner = np.append(np.append(lat[:-1,0]- np.diff(lat[:,0])/2,lat[-1,0]),
                          lat[-1,0]+np.diff(lat[-3:-1,0])/2)
    
    # Inicializando a grid
    grids=[]
    
    # Loop para cada longitude
    for ii in range(1,lonCorner.shape[0]):
        
        #Loop over each cel in y direction
        for jj in range(1,latCorner.shape[0]):
            
            #Criando retângulo de de cada célula
            lat_point_list = [latCorner[jj-1], latCorner[jj], latCorner[jj], latCorner[jj-1]]
            lon_point_list = [lonCorner[ii-1], lonCorner[ii-1], lonCorner[ii], lonCorner[ii]]
            
            # Criando um polígono para cada celula
            cel = Polygon(zip(lon_point_list, lat_point_list))
            grids.append(cel)
            
    
    

    # Abrindo shapefile com a textura do solo
    shapeSolos = gpd.read_file(inputFolder+'/Solos_5000mil/Solos_5000.shp')
    shapeSolos = shapeSolos.set_crs('epsg:4326')
    df = pd.DataFrame({'geometry':grids})
    gdf = gpd.GeoDataFrame(df, crs="EPSG:4326")
    gdf.set_geometry('geometry', inplace=True)
 
    # Recortando o tipo de solo para cada célula do domínio.
    soilIdx=[]
    for index, row in gdf.iterrows():
        clipped = gpd.clip(shapeSolos, row.geometry)
        if clipped['DSC_TEXTUR'].iloc[0] == 'argilosa  ou muito argilosa':
            soilIdx.append(3)
        elif clipped['DSC_TEXTUR'].iloc[0] == 'media':
            soilIdx.append(2)
        elif clipped['DSC_TEXTUR'].iloc[0] == 'arenosa':
            soilIdx.append(1)
            
        else:
            soilIdx.append(np.nan)
    
    # # Inicializando a matriz de pixels de cada idSoil no domínio
    matRegrid = np.empty((lat.shape[0], lat.shape[1]))
    
    # fazendo o rashape
    matRegrid[:,:] = np.array(soilIdx).reshape((lon.shape[1],lon.shape[0])).transpose() 

    # substitui nan por 0
    matRegrid[np.isnan(matRegrid)] = 0
    
    # escreve o arquivo netCDF com a textura do solo para não precisar fazer 2 vezes
    regMap.createNETCDF(outfolder,'regridedSoilTexture_'+GDNAM,matRegrid,lon,lat)
    
    return matRegrid



def soilType(inputFolder,outfolder,lat,lon,D,GDNAM):
    """
    função utilizada para determinar a porcentagem de particulas de um determinado
    diâmetro em cada celula do dominio. 
    
    ESTA FUNÇÃO PRECISA SER VERIFICADA 

    Parameters
    ----------
    inputFolder : path
        caminho para a pasta de inputs.
    outfolder : pat
        caminho para a pasta de outputs.
    lat : np.array
        matriz de lat do dominio.
    lon : np.array
        matriz de lon do dominio.
    D : float
        diametro da particula.
    GDNAM : str
        nome da grade de acordo com o MCIP.

    Returns
    -------
    sRef : np.array
        matriz com os valores da porcentagem de particulas com um determinado 
        diametro.

    """
    
    # abre o raster com o regrid da soiltexture
    raster = nc.Dataset(outfolder+'/regridedSoilTexture_'+GDNAM+'.nc', 
                        masked=True)['MAT'][:]

    # inicializa a matriz que conterá os valores de porcentagem
    sRef=np.empty((lat.shape[0],lat.shape[1]))
    sRef[:,:] = 0
    
    # lista com os tipos de soilTextures
    soilNames=['Sand', 'Silt', 'Clay']
    
    # loop para cara soilTextures
    for kk,soiln in enumerate(soilNames):
        
        # abre csv com a distribuição de particulas para cada uso do solo
        # https://www.slideshare.net/slideshow/classificac3a7c3a3o-dossolosaashtosucs/49327763
        # VERIFICAR
        soilDist = pd.read_csv(inputFolder+'/tables/particleDist/'+soiln+'.csv')
        
        # Converte para micrometros os diâmetros
        xs=soilDist['D']*1000
        
        # proporção acumulada de particulas em cada diametro
        ys=soilDist['P']
        
        # função para fitar a curva de granulometrica acumulada
        f = lambda x,mu,sigma: stats.norm(mu,sigma).cdf(x)
        
        # encontrando o mu e sigma da curva
        mu,sigma = optimize.curve_fit(f,xs,ys/100)[0]
        
        # valores de diâmetros de 0 a 800 micrometros para usar na funçaõ fitada
        xx = np.arange(0,800,0.01)
        
        # derivada da curva acumulada, ou seja, o valor de porcentagem de um 
        # determinado diametro. 
        deriv = np.append(np.nan,np.diff(f(xx,mu,sigma)*100))
        
        # indice da matriz que possui o determinado diâmetro
        idx = find_nearest(xx, D)
        
        # estabelece o valor de porcentagem para um determinado diametro na matriz
        # com o mesmo tamanho do dominio
        sRef[raster[0,:,:].astype(int)==kk+1]=deriv[idx]

    return sRef


    
def main(inputFolder,outfolder,domainShp,GDNAM,lat,lon,D,RESET_GRID):
    """
    Esta função controla a geração de arquivos de solo - claycontent, soiltexture,
    porcentagem de particulas de um determinado diametro.

    Parameters
    ----------
    inputFolder : path
        DESCRIPTION.
    outfolder : path
        DESCRIPTION.
    domainShp : geodataframe
        geodataframe com o shape do domínio.
    GDNAM : str
        nome da grade de acordo com o MCIP.
    lat : np.array
        matriz de latitudes do dominio.
    lon : np.array
        matriz de longitudes do dominio.
    D : float
        diametro da particula.
    RESET_GRID : boolean
        True ou False para resetar a matriz e reescrever o netCDF com o clay content.

    Returns
    -------
    clayRegrid : np.array
        matriz com o clay content.
    sRef : np.array
        matriz de porcentagem de particulas com um determinado diâmetro.

    """
    print('=====STARTING soilPrep.py=====' )
    
    # se existir o arquivo de regridClay para o domínio
    if (os.path.exists(outfolder+'/regridClay_'+GDNAM+'.nc')) and (os.path.exists(
            outfolder+'/regridedSoilTexture_'+GDNAM+'.nc')):
        
        # se não quiser resetar a grade = usa o regrid que já tem
        if RESET_GRID==False:
            print ('You already have the regridClay_'+GDNAM+'.nc file')
            
            # abre o netCDF com o regridClay já feito
            ds = nc.Dataset(outfolder+'/regridClay_'+GDNAM+'.nc')
            clayRegrid = ds['MAT'][:]
            
            # abre o arquivo de regridedSoilTexture já feito
            print ('You already have the regridedSoilTexture_'+GDNAM+'.nc file')
            ds = nc.Dataset(outfolder+'/regridedSoilTexture_'+GDNAM+'.nc')
            
            # roda a função de soilType que precisa ser executada 
            # para cada diâmetro
            sRef = soilType(inputFolder,outfolder,lat,lon,D,GDNAM)
        
        # se quiser resetar a grade    
        else:
            
          # executa a função cutSoil para cortar o arquivo original
          raster = cutSoil(domainShp,inputFolder,outfolder,GDNAM)
          
          # extraindo x e y dos pixels do raster 
          x, y = rasterLatLon(raster)
          
          # executa a função para fazer o regrid do clayContent
          clayRegrid = rasterInGrid(domainShp,raster,x,y,lat,lon)
          
          # cria o netCDF com o regrid do clayCOntent
          print('Creating netCDF')
          regMap.createNETCDF(outfolder,'regridClay_'+GDNAM,clayRegrid,lon,lat)
          
          # executa a função para fazer o regrid do soilTexture
          regridSoilTexture(outfolder,inputFolder,lat,lon,GDNAM)
          
          # executa a função para estimar a porcentagem de particulas na grade
          sRef = soilType(inputFolder,outfolder,lat,lon,D,GDNAM)
          sRef[np.isnan(sRef)] = 0  
    
    # se não existir o arquivo de regrid      
    else:
        
        # executa a função cutSoil para cortar o arquivo original
        raster = cutSoil(domainShp,inputFolder,outfolder,GDNAM)
        
        # extraindo x e y dos pixels do raster 
        x, y = rasterLatLon(raster)
        
        # executa a função para fazer o regrid do clayContent
        clayRegrid = rasterInGrid(domainShp,raster,x,y,lat,lon)
        
        # cria o netCDF com o regrid do clayCOntent        
        print('Creating netCDF')
        regMap.createNETCDF(outfolder,'regridClay_'+GDNAM,clayRegrid,lon,lat)
        
        # executa a função para fazer o regrid do soilTexture
        regridSoilTexture(outfolder,inputFolder,lat,lon,GDNAM)
        
        # executa a função para estimar a porcentagem de particulas na grade
        sRef = soilType(inputFolder,outfolder,lat,lon,D,GDNAM)
        sRef[np.isnan(sRef)] = 0
        
    return clayRegrid,sRef
