#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
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
import os

wrfoutPath='/media/leohoinaski/HDD/SC_2019/wrfout_d02_2019-01-03_18:00:00'
#wrfoutPath='/mnt/sdb1/SC_2019/wrfout_d02_2019-01-01'
GRDNAM = 'SC_2019'
inputFolder = os.path.dirname(os.getcwd())+'/inputs'
tablePath = os.path.dirname(os.getcwd())+'/inputs/tables'
outfolder = os.path.dirname(os.getcwd())+'/outputs'
year = 2021
idSoils = [23,24,30,25] #4.1. Praia, Duna e Areal  4.2. Área Urbanizada  4.3. Mineração 4.4. Outras Áreas não Vegetadas
D = 10
av,al,alarea,lat,lon,domainShp = regMap.main(wrfoutPath,GRDNAM,inputFolder,outfolder,year,idSoils)
clayRegrid,sRef = sp.main(inputFolder,outfolder,domainShp,GRDNAM,lat,lon,D)
ustar,ustarT,ustarTd = mp.main(wrfoutPath,tablePath,av,al,alarea,D,clayRegrid)
#Fdust = wbd.wbdFlux(av,al,SrefD,ustar,ustarT,ustarTd)