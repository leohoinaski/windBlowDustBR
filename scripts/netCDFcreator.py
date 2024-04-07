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
import netCDF4 as nc
#from pyproj import Proj, transform
import numpy as np
import datetime
import pandas as pd
#import tarfile
#import os

def datePrepCMAQ(ds):
    tf = np.array(ds['TFLAG'][:][:,1,:])
    date=[]
    for ii in range(0,tf.shape[0]):
        date.append(datetime.datetime.strptime(tf[:,0].astype(str)[ii] + (tf[:,1]/10000).astype(int).astype(str)[ii], '%Y%j%H').strftime('%Y-%m-%d %H:00:00'))
    
    date = np.array(date,dtype='datetime64[s]')
    dates = pd.DatetimeIndex(date)
    datesTime=pd.DataFrame()
    datesTime['year'] = dates.year
    datesTime['month'] = dates.month
    datesTime['day'] = dates.day
    datesTime['hour'] = dates.hour
    datesTime['datetime']=dates
    return datesTime

def datePrepWRF(date):
    date = np.array(date,dtype='datetime64[s]')
    dates = pd.DatetimeIndex(date)
    datesTime=pd.DataFrame()
    datesTime['year'] = dates.year
    datesTime['month'] = dates.month
    datesTime['day'] = dates.day
    datesTime['hour'] = dates.hour
    datesTime['datetime']=dates
    return datesTime

def createNETCDFtemporal(folder,name,data,datesTime,mcipMETCRO3Dpath,EmisD):
    data[np.isnan(data)]=0
    print('===================STARTING netCDFcreator_v1.py=======================')
    print(data.shape)
    ds = nc.Dataset(mcipMETCRO3Dpath)  
    #datesTime = datePrepCMAQ(ds)
    print('Initial date: '+str(datesTime.iloc[0,-1]))
    print('Final date: '+str(datesTime.iloc[-1,-1]))
    f2 = nc.Dataset(folder+'/'+name+EmisD['tag']+'_'+\
                    str(datesTime.iloc[0,-1]).replace(' ','-')+'_'+\
                    str(datesTime.iloc[-1,-1]).replace(' ','-')+\
                    '.nc','w', format='NETCDF3_CLASSIC') #'w' stands for write   
    #Add global attributes
    f2.IOAPI_VERSION ='$Id: @(#) ioapi library version 3.1 $'
    f2.EXEC_ID = '???????????????'
    f2.FTYPE =  1
    for gatr in ds.ncattrs() :
        setattr(f2, gatr, ds.getncattr(gatr))
    # f2.VGTYP= -9999
    #f2.VGTOP= 0.0
    #f2.VGLVLS= [0,0]
    #f2.NLAYS = 1 
    # f2.NCOLS = f2.NCOLS
    # f2.NROWS = f2.NROWS
    f2.FILEDESC= 'File produced by WindBlowDustBr.py'
    f2.HISTORY =''
    # # Specifying dimensions
    #tempgrp = f.createGroup('vehicularEmissions_data')
    f2.createDimension('TSTEP', None )
    f2.createDimension('DATE-TIME', 2)
    f2.createDimension('LAY', 1)
    
    if len(data.shape)==3:
        f2.NVARS= 1
        f2.createDimension('VAR', 1)
        f2.createDimension('ROW', data.shape[1])
        f2.createDimension('COL', data.shape[2])
        tflag = np.empty([data.shape[0],1,2],dtype='i4')
        for ii in range(0,data.shape[0]):
            tflag[ii,:,0]=int(datesTime['year'].iloc[0]*1000 + datesTime.datetime.iloc[ii].timetuple().tm_yday)
            tflag[ii,:,1]=int(str(datesTime['hour'].iloc[ii])+'0000')
        # Building variables
        TFLAG = f2.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        TFLAG[:,:,:] = tflag
        TFLAG.units = '<YYYYDDD,HHMMSS>'
        
        strVAR = ''
        polids = [EmisD['tag']]
        for ids in polids:
            strVAR = strVAR + ids.ljust(16)
        #strVAR ='ACET            ACROLEIN        ALD2            ALD2_PRIMARY    ALDX            BENZ            BUTADIENE13     CH4             CH4_INV         CL2             CO              CO2_INV         ETH             ETHA            ETHY            ETOH            FORM            FORM_PRIMARY    HCL             HONO            IOLE            ISOP            KET             MEOH            N2O_INV         NAPH            NH3             NH3_FERT        NO              NO2             NVOL            OLE             PAL             PAR             PCA             PCL             PEC             PFE             PH2O            PK              PMC             PMG             PMN             PMOTHR          PNA             PNCOM           PNH4            PNO3            POC             PRPA            PSI             PSO4            PTI             SO2             SOAALK          SULF            TERP            TOL             UNK             UNR             VOC_INV         XYLMN           '
        setattr(f2, 'VAR-LIST', strVAR)
        
        for ii, pid in enumerate(polids):
            globals()[pid] = f2.createVariable(pid, np.float32, ('TSTEP', 'ROW','COL'))
            globals()[pid][:,:,:] = data[:,:,:]
            globals()[pid].var_desc = pid+'[1]'
            globals()[pid].units = EmisD['Unit']
    else:
        f2.NVARS= data.shape[0]
        f2.createDimension('VAR', data.shape[0])
        f2.createDimension('ROW', data.shape[2])
        f2.createDimension('COL', data.shape[3])
        tflag = np.empty([data.shape[1],1,2],dtype='i4')
        for ii in range(0,datesTime.shape[0]):
            tflag[ii,:,0]=int(datesTime['year'].iloc[0]*1000 + datesTime.datetime.iloc[ii].timetuple().tm_yday)
            tflag[ii,:,1]=int(str(datesTime['hour'].iloc[ii])+'0000')
        # Building variables
        TFLAG = f2.createVariable('TFLAG', 'i4', ('TSTEP', 'VAR', 'DATE-TIME'))
        TFLAG[:,:,:] = tflag
        TFLAG.units = '<YYYYDDD,HHMMSS>'
        
        strVAR = ''
        #polids = EmisD['fractions']
        polid = pd.DataFrame()
        polid['ID'] = EmisD['fractions']
        for ids in polid.ID:
            strVAR = strVAR + ids.ljust(16)
        #strVAR ='ACET            ACROLEIN        ALD2            ALD2_PRIMARY    ALDX            BENZ            BUTADIENE13     CH4             CH4_INV         CL2             CO              CO2_INV         ETH             ETHA            ETHY            ETOH            FORM            FORM_PRIMARY    HCL             HONO            IOLE            ISOP            KET             MEOH            N2O_INV         NAPH            NH3             NH3_FERT        NO              NO2             NVOL            OLE             PAL             PAR             PCA             PCL             PEC             PFE             PH2O            PK              PMC             PMG             PMN             PMOTHR          PNA             PNCOM           PNH4            PNO3            POC             PRPA            PSI             PSO4            PTI             SO2             SOAALK          SULF            TERP            TOL             UNK             UNR             VOC_INV         XYLMN           '
        setattr(f2, 'VAR-LIST', strVAR)
        print(strVAR)
        
        for ii in range(0,polid.shape[0]):
            globals()[polid.ID[ii]] = f2.createVariable(polid.ID[ii], np.float32, ('TSTEP', 'ROW','COL'))
            globals()[polid.ID[ii]][:,:,:,:] = data[ii,:,:,:]
            #globals()[polid.ID[ii]].var_desc = polid.ID[ii]+'[1]'
            globals()[polid.ID[ii]].units = EmisD['Unit']
            

    
        



   
    f2.close()
    print('WindBlowDust emissions in netCDF format is ready for CMAQ!')
    return TFLAG
