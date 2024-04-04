#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 13:04:17 2024

@author: leohoinaski
"""
import numpy as np
import pandas as pd
import netCDF4 as nc
import wrf

def ustarCalc(uz,z0):
    """
    Esta função calcula a velocidade frictiva

    Parameters
    ----------
    uz : numpy array
        velocidade do vento.
        
    z0 : numpy array
        rugosidade.
    Returns
    -------
    ustar : numpy array
        Velocidade frictiva.

    """
    k= 0.4 # von karman constant
    z= 10 # altura de referência
    ustar = k*uz/np.log(z/z0)
    ustar[np.isnan(z0)]=np.nan
    return ustar

def roughness(tablePath,av,al):
    alpS = pd.read_csv(tablePath + '/hs.csv')
    hv = pd.read_csv(tablePath + '/hv.csv')
    alphaS = av*alpS['alphaS'][0] + np.nansum(al,axis=0)*alpS['alphaS'][2] 
    alphaV = -0.35*np.log(1-av)
    alphaV[alphaV == np.inf] = np.nanmax(alphaV)
    alpha = alphaS + alphaV
    hS= av*alpS['hs'][0] + al*alpS['hs'][2]
    hV = av*hv['mean'][0] + al*hv['mean'][2]
    h = (hV*alphaV + hS*alphaS)/(alphaS + alphaV)
    z0 = alpha.copy() 
    # z0[alpha<0.2] = h[alpha<0.2]*0.96*alpha[alpha<0.2]**(1.07)
    # z0[alpha>=0.2] = h[alpha>=0.2]*0.083*alpha[alpha>=0.2]**(-0.46)
    z0[alpha<0.2] = (0.96*alpha[alpha<0.2]**(1.07))*h[alpha<0.2]
    z0[alpha>=0.2] = (0.083*alpha[alpha>=0.2]**(-0.46))*h[alpha>=0.2]
    z0[:,np.isnan(al)]=np.nan
    alphaV[:,np.isnan(al)]=np.nan
    alphaS[:,np.isnan(al)]=np.nan
    return z0,alphaV,alphaS

def ustarThreshold(D,clayRegrid,w,alphaV,alphaS,av):
    An = 0.0123
    tao = 1.65*10**(-4)
    roa = 1.227# g/cm³
    rop = 2.650 # g/cm³ assumi - montar matriz de densidades
    g = 9.81
    ustarTd = np.sqrt(An*((rop*g*(D*10**(-6))/roa)+(tao/roa*D*10**(-6))))
    clayRegrid = clayRegrid*100
    wl = 0.0014*(clayRegrid**2)+0.17*clayRegrid
    fm = np.empty(w.shape)
    fm[:,:,:]=np.nan
    fm=fm.astype(float)
    fm = (1+1.21*(w-wl)**(0.68))**(0.5)
    fm[w<wl] = 1 
    # sigmaV = 1.45
    # mV = 0.16
    # betaV = 202
    sigmaS = 1
    mS = 0.5
    betaS = 90
    #fr = ((1-sigmaV*mV*alphaV)*(1+betaV*mV*alphaV)*(1-sigmaS*mS*(alphaS/(1-av)))*(1+betaS*mS*(alphaS/(1-av))))**0.5
    # Referência do fr = https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2010JD014649
    fr = ((1-sigmaS*mS*(alphaS/(1-av)))*(1+betaS*mS*(alphaS/(1-av))))**0.5 # ref = )
    ustarT=ustarTd*fm.data*fr
    return ustarT,ustarTd

def main(wrfoutPath,tablePath,av,al,D,clayRegrid,lia):
    #tablePath = '/mnt/sdb1/windBlowDustBR/inputs/tables'
    ds = nc.Dataset(wrfoutPath)
    av = ds['VEGFRA'][lia,:,:]/100
    ustarWRF = ds['UST'][lia,:,:]
    z0,alphaV,alphaS = roughness(tablePath,av,al)
    #metPath = '/mnt/sdb1/SC_2019/wrfout_d02_2019-01-01'
    uz = np.array(wrf.g_wind.get_destag_wspd_wdir10(ds,timeidx=wrf.ALL_TIMES)[0,lia,:,:])
    ustar = ustarCalc(uz,z0)
    w = ds['SMOIS'][lia,0,:,:]
    #ust = ds['UST'][:]
    ustarT,ustarTd = ustarThreshold(D,clayRegrid,w,alphaV,alphaS,av)
    return ustar,ustarT,ustarTd,av,ustarWRF