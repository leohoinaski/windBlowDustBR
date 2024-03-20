#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 16:04:37 2024

@author: leohoinaski
"""

import regridMAPBIOMAS as regMap
import os

#wrfoutPath='/media/leohoinaski/HDD/SC_2019/wrfout_d02_2019-01-03_18:00:00'
wrfoutPath='/mnt/sdb1/SC_2019/wrfout_d02_2019-01-01'
GRDNAM = 'SC_2019'
inputFolder = os.path.dirname(os.getcwd())+'/inputs'
outfolder = os.path.dirname(os.getcwd())+'/outputs'
year = 2021
idSoils = [23,24,30,25] #4.1. Praia, Duna e Areal  4.2. Área Urbanizada  4.3. Mineração 4.4. Outras Áreas não Vegetadas

matRegrid,av,al,lat,lon,domainShp = regMap.main(wrfoutPath,GRDNAM,inputFolder,outfolder,year,idSoils)