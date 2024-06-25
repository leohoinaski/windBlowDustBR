#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------shRunnerWindBlowDust.py-------------------------------

Este script em python é utilizado para rodar o módulo windBlowDust e gerar os 
inputs para o CMAQ. O script usa argumentos de entrada na linha de comando do 
terminal. Os argumentos são:

   windBlowDustFolder = caminho da pasta master do módulo
   mcipPath = caminho para os arquivos do MCIP
   wrfoutFolder = caminho para os arquivos do WRF
   domain = número do domínio do WRF
   GDNAM = nome da grade conforme o MCIP
   YEAR = ano de referência da simulação
   RESET_GRID = True para resetar os arquivos intermediários gerados pelo módulo 

Created on Fri Mar 15 16:04:37 2024

Referências:
    https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016MS000823
    https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2010JD014649

@author: leohoinaski

"""

import regridMAPBIOMAS as regMap
import soilPrep as sp
import metPrep as mp
import windBlowDustCalc as wbd
import netCDFcreator as ncCreate
import os
import numpy as np
import netCDF4 as nc
import wrf
import pandas as pd
from datetime import timedelta
import ismember
import argparse
import windBlowDustSpeciation as wbds



if __name__ == '__main__':
    
    # Argumentos de entrada
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', default=0, action='count')
    parser.add_argument('windBlowDustFolder')
    parser.add_argument('mcipPath')
    parser.add_argument('wrfoutFolder')
    parser.add_argument('domain')
    parser.add_argument('GDNAM')
    parser.add_argument('YEAR', type=int)
    parser.add_argument('RESET_GRID', type=int) # 0 OU 1
    args = parser.parse_args()
    
    # passando os argumentos para variáveis
    windBlowDustFolder = args.windBlowDustFolder
    mcipPath = args.mcipPath
    wrfoutFolder = args.wrfoutFolder
    domain = args.domain
    domain = 'd0'+domain
    GDNAM  = args.GDNAM
    YEAR = args.YEAR
    RESET_GRID = args.RESET_GRID
    
    # definindo o caminho para o arquivo METCRO3D do MCIP
    mcipMETCRO3Dpath = mcipPath+'/METCRO3D_'+GDNAM+'.nc'
    
    # definição dos caminhos das pastas de input, output e table
    inputFolder = windBlowDustFolder+'/inputs'
    tablePath = windBlowDustFolder+'/inputs/tables'
    outfolder = windBlowDustFolder+'/Outputs/'+GDNAM
    
    
    # condição para restar ou não os arquivos intermediários.
    if RESET_GRID == 0:
        RESET_GRID=False
    else:
        RESET_GRID=True
    
    # Dicionários de poluentes
    PM25 = {
      "Unit": '$\g.S^{-1}$',
      "tag":'PMFINE',
      "range":[0,2.5] # micrometers
    }
    
    PMC = {
      "Unit": '$\g.S^{-1}$',
      "tag":'PMC',
      "range":[2.5,10] # micrometers
    }
    
    PM10 = {
      "Unit": '$\g.S^{-1}$',
      "tag":'PM10',
      "range":[0,10] # micrometers
    }
    
    PM1 = {
      "Pollutant": "$NO_{2}$",
      "Unit": '$\g.S^{-1}$',
      "tag":'PMULTRAFINE',
      "range":[0,1] # micrometers
    }
    
    ALL = {
      "Unit": '$\g.S^{-1}$',
      "tag":'AllFractions',
      "fractions":['PMFINE','PMC','PM10'] # micrometers
      #"fractions":['PMFINE','PMC'] # micrometers
    }

    
    # Definição dos ids do MAPBIOMAS que serão utilizados na estimativa 
    # das emissões no windblowdust
    idSoils = [23,30] #4.1. Praia, Duna e Areal  4.3. Mineração 4.4. 25 Outras Áreas não Vegetadas
    
    # espaçamento entre diâmetros para integração dos valores
    dx = 0.1 # dx para integração das emissões das frações.
    
    # frações que serão calculadas
    Fractions = [PM25,PMC] # Lista com tipo de emissão por diâmetro. 
                            #Não precisa incluir o PM10 se já tiver PM25 e PM10
       
    
    # condição para verificar se a pasta de output existe
    if os.path.isdir(outfolder):
        print('You have the outputs folder')
    else:
        # se não existir, cria a pasta.
        os.makedirs(outfolder, exist_ok=True)
    
    # abre o arquivo METCROD3D
    ds = nc.Dataset(mcipMETCRO3Dpath)
    
    # Extrai as datas do arquivo do MCIP
    datesTimeMCIP = ncCreate.datePrepCMAQ(ds)
    
    # Deinindo os arquivos do WRF que serão abertos. Devem ser compatíveis com 
    # as datas do respectivo arquivo do MCIP. 
    # file = [i for i in os.listdir(wrfoutFolder) if os.path.isfile(os.path.join(wrfoutFolder,i)) and \
    #          'wrfout_'+domain+'_'+str((datesTimeMCIP.datetime[0]- timedelta(days=1)).year).zfill(4)+'-'+\
    #              str((datesTimeMCIP.datetime[0]- timedelta(days=1)).month).zfill(2)+'-'+\
    #                  str((datesTimeMCIP.datetime[0]- timedelta(days=1)).day).zfill(2) in i]
    files=['wrfout_'+domain+'_'+str((datesTimeMCIP.datetime[0]).year).zfill(4)+'-'+\
                 str((datesTimeMCIP.datetime[0]).month).zfill(2)+'-'+\
                     str((datesTimeMCIP.datetime[0]).day).zfill(2)+'_00:00:00',\
           'wrfout_'+domain+'_'+str((datesTimeMCIP.datetime[0]+ timedelta(days=1)).year).zfill(4)+'-'+\
                 str((datesTimeMCIP.datetime[0]+ timedelta(days=1)).month).zfill(2)+'-'+\
                     str((datesTimeMCIP.datetime[0]+ timedelta(days=1)).day).zfill(2)+'_00:00:00']
    
    # pega o caminho da pasta atual
    cwd = os.getcwd()
    
    # move para a pasta do WRF
    os.chdir(wrfoutFolder)
    #wrfoutPath = wrfoutFolder+'/'+file[0]
    wrfoutPath = wrfoutFolder+'/'+files[0]
    
    # abre os arquivos do WRF
    ds = nc.MFDataset(files)
    
    # Volta para a pasta root do windblow/scripts
    os.chdir(cwd)
    
    # extrai as datas dos arquivos do WRF que foram abertos
    datesTime = ncCreate.datePrepWRF(pd.to_datetime(
        wrf.extract_times(ds,wrf.ALL_TIMES)))
    
    # identifica datas coincidentes no MCIP e WRF
    lia, loc = ismember.ismember(np.array(datesTime.datetime),
                                 np.array(datesTimeMCIP.datetime))
    print(datesTime.iloc[lia,:])
    
    # executa a função de regridMAPBIOMAS
    av,al,alarea,lat,lon,domainShp = regMap.main(wrfoutPath,GDNAM,inputFolder,
                                                 outfolder,YEAR,idSoils,RESET_GRID)
    
    # loop para cada fração do PM
    for EmisD  in Fractions:
        
        # determina os ranges das particulas - min e max
        Dmax = np.max(EmisD['range'])
        Dmin = np.min(EmisD['range'])
        
        # cria um arranjo de diâmetros de particulas para integração
        diam = np.arange(np.min(EmisD['range']),np.max(EmisD['range'])+dx,dx)
        
        # inicializa a variável FdustTotal que acumulará as estimativas para 
        # cada diâmetro
        FdustTotal=[]
        
        # seleciona os valores dentro da faixa da fração ex-0 a 2.5 micrometros
        diamSelect = diam[(Dmin<=diam) & (Dmax>=diam)]
        
        # loop em cada diâmetro
        for jj,diameters in enumerate(diamSelect):
            
            print(diameters)
            
            # executa a função soilPrep
            clayRegrid,sRef = sp.main(inputFolder,outfolder,domainShp,GDNAM,
                                      lat,lon,diameters,RESET_GRID)
            
            # executa a função metPrep
            ustar,ustarT,ustarTd,avWRF,ustarWRF = mp.main(ds,tablePath,av,
                                                          al,diameters,
                                                          clayRegrid,lia)
            
            # executa a função windBlowDustCalc
            Fdust,Fhd,Fhtot,Fvtot = wbd.wbdFlux(avWRF,alarea,sRef,clayRegrid,
                                                ustar,ustarT,ustarTd)
            
            # já rodou uma vez, logo, não precisa resetar os arquivos 
            # intermediários
            RESET_GRID = False
            
            # acumula os valores em cada diâmetro
            FdustTotal.append(Fdust)
        
        # empilha os valores em um array numpy
        FdustTotal = np.stack(FdustTotal)
        
        # estima a massa total de particulas dentro da faixa da fração
        # faz a integral dos dados estimados
        #FdustD = np.trapz(FdustTotal,dx=dx, axis=0)
        
        
        # faz a média do fluxo para cada diâmetro
        FdustD = np.nanmedian(FdustTotal, axis=0)   
        
        # cria o netCDF com a estimativa das particulas
        ncCreate.createNETCDFtemporal(outfolder,'windBlowDust_',
                                      FdustD,datesTime.iloc[lia,:],
                                      mcipMETCRO3Dpath,EmisD)
        
        # faz a especiação química das particulas
        if EmisD==PM25:
            FdustFINE = FdustD
            FdustFINESpec = wbds.speciate(windBlowDustFolder, FdustFINE)
        elif EmisD==PMC:
            FdustCOARSE = FdustD
            FdustCOARSEpec = wbds.speciate(windBlowDustFolder, FdustCOARSE)
        else:
            print('You have selected an awkward fraction')
        
        # se existir as parcelas PMFINE e PMC, calcula o PM10
        try:
            FdustPM10 = FdustFINE+FdustCOARSE
            ncCreate.createNETCDFtemporal(outfolder,'windBlowDust_',
                                          FdustPM10,datesTime[lia],
                                          mcipMETCRO3Dpath,PM10)
        except:
            print('You do not have the fractions required for PM10')   
            
    
    # Acumula todas as estimativas de particulas sem especiação
    FdustALL = [FdustFINE,FdustCOARSE,FdustPM10]
    #FdustALL = [FdustFINE,FdustCOARSE]
    
    # empilha em uma array numpy
    FdustALL = np.stack(FdustALL,axis=0)
    
    # cria o netCDF com todas as especies de particulas/frações
    ncCreate.createNETCDFtemporal(outfolder,'windBlowDust_',FdustALL,
                                  datesTime[lia],mcipMETCRO3Dpath,ALL)
    
    # soma as emissões de cada especie no PM25 e PMC
    FdustSpeciated = FdustFINESpec + FdustCOARSEpec
    
    # cria o netCDF especiado
    ncCreate.createNETCDFtemporalSpeciated(windBlowDustFolder,
                                           outfolder,'windBlowDust_',
                                           FdustSpeciated,
                                           datesTime[lia],
                                           mcipMETCRO3Dpath)
    
