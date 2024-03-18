#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 13:04:17 2024

@author: leohoinaski
"""
import numpy as np
import pandas as pd

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
    ustar = k*uz/np.nanlog(z/z0)
    return ustar

def roughness(tablePath,av,al):
    
    alpS = pd.read_csv(tablePath + '/hs.csv')
    hv = pd.read_csv(tablePath + '/hv.csv')
    alphaS = av*alpS['alphaS'][0] + al*alpS['alphaS'][2] 
    alphaV = -0.35*np.nanlog(1-av)
    alpha = alphaS + alphaV
    hS= av*alpS['hs'][0] + al*alpS['hs'][2] 
    hV = av*hv['mean'][0] + al*hv['mean'][2] 
    h = (hV*alphaV + hS*alphaS)/(alphaS + alphaV)
    
    
    if alpha<0.2:
        z0 = h*0.96*alpha**(1.07)
    else:
        z0 = h*0.083*alpha**(-0.46)
    
    return z0,alphaV,alphaS

def ustarThreshold(D,clay,w,alphaV,alphaS,av):
    An = 0.0123
    tao = 1.65*10**(-4)
    roa = 1.2923 # g/cm³
    rop = 2.6 # g/cm³ assumi - montar matriz de densidades
    g = 9.81
    ustarTd = np.sqrt(An*((rop*g*D/roa)+(tao/roa*D)))
    wl = 0.0014*(clay**2)+0.17*clay
    if w<wl:
        fm = 1 
    else:
        fm = (1+1.21*(w-wl)**(0.68))**(0.5)
     
    sigmaV = 1.45
    mV = 0.16
    betaV = 202
    sigmaS = 1
    mS = 0.5
    betaS = 90
    fr = ((1-sigmaV*mV*alphaV)*(1-betaV*mV*alphaV)*(1-sigmaS*mS*(alphaS/(1-av)))*(1+betaS*mS*(alphaS/(1-av))))**0.5
    ustarT = ustarTd*fm*fr
    return ustarT