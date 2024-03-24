#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------------------------------------------------
                            netCDFcreator_v1.py

This function creates the netCDF files ready to use in CMAQ.

Inputs:
    
    folder: folter to output files
    
    name: output names
    
    data: matrix with data ready to convert in netCDF
    
    xv, yv: meshgrid outputs - grid definition
    
    center: pixel centroids 
    
    year: respective years of emission inventories
    
    month: respective month of emission inventories

Outputs:
        
    netdCDF files
      

Last update = 29/09/2021

Author: Leonardo Hoinaski - leonardo.hoinaski@ufsc.br
-------------------------------------------------------------------------------
"""
import netCDF4 as nc4
#from pyproj import Proj, transform
import numpy as np
#import datetime
#import pandas as pd
#import tarfile
#import os

def createNETCDFtemporal(rootPath,folder,name,data,mcipMETCRO3Dpath,D):
    print('===================STARTING netCDFcreator_v1.py=======================')
          
    f2 = nc4.Dataset(folder+'/'+name+'PM'+str(D).replace('.',''),'w', format='NETCDF3_CLASSIC') #'w' stands for write   
    #Add global attributes
    f2.IOAPI_VERSION ='$Id: @(#) ioapi library version 3.1 $'
    f2.EXEC_ID = '???????????????'
    f2.FTYPE =  1
    ds3 = nc4.Dataset(mcipMETCRO3Dpath)
    for gatr in ds3.ncattrs() :
        setattr(f2, gatr, ds3.getncattr(gatr))
    # f2.VGTYP= -9999
    #f2.VGTOP= 0.0
    #f2.VGLVLS= [0,0]
    #f2.NLAYS = 1 
    f2.NVARS= data.shape[1]
    # f2.NCOLS = f2.NCOLS
    # f2.NROWS = f2.NROWS
    f2.FILEDESC= 'File produced by WindBlowDustBr.py'
    f2.HISTORY =''
    # # Specifying dimensions
    #tempgrp = f.createGroup('vehicularEmissions_data')
    f2.createDimension('TSTEP', None )
    f2.createDimension('DATE-TIME', 2)
    f2.createDimension('LAY', ds3.NLAYS)
    f2.createDimension('VAR', data.shape[1])
    f2.createDimension('ROW', ds3.NROWS)
    f2.createDimension('COL', ds3.NCOLS)
    
    # Building variables
    TFLAG = f2.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
    
    strVAR = ''
    polids = ['PM'+str(D).replace('.','')]
    for ids in polids:
        strVAR = strVAR + ids.ljust(16)
    #strVAR ='ACET            ACROLEIN        ALD2            ALD2_PRIMARY    ALDX            BENZ            BUTADIENE13     CH4             CH4_INV         CL2             CO              CO2_INV         ETH             ETHA            ETHY            ETOH            FORM            FORM_PRIMARY    HCL             HONO            IOLE            ISOP            KET             MEOH            N2O_INV         NAPH            NH3             NH3_FERT        NO              NO2             NVOL            OLE             PAL             PAR             PCA             PCL             PEC             PFE             PH2O            PK              PMC             PMG             PMN             PMOTHR          PNA             PNCOM           PNH4            PNO3            POC             PRPA            PSI             PSO4            PTI             SO2             SOAALK          SULF            TERP            TOL             UNK             UNR             VOC_INV         XYLMN           '
    setattr(f2, 'VAR-LIST', strVAR)
    
    for ii, pid in enumerate(polids):
        globals()[pid] = f2.createVariable(pid, np.float32, ('TSTEP', 'LAY', 'ROW','COL'))
        globals()[pid][:,:,:,:] = data[ii,:,:,:,:]
        globals()[pid].var_desc = pid+'[1]'
        globals()[pid].units = 'VERIFICAR!!'


    # Passing data into variables
    TFLAG[:,:,:] = np.repeat(ds3['TFLAG'][:,0,:][:, np.newaxis,:], 
                             data.shape[1], axis=1)
    TFLAG.units = '<YYYYDDD,HHMMSS>'
    TFLAG.var_desc = 'Timestep-valid flags:  (1) YYYYDDD or (2) HHMMSS'
   
    f2.close()
    print('WindBlowDust emissions in netCDF format is ready for CMAQ!')
    return TFLAG