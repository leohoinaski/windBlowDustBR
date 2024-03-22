#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 15:20:44 2024

@author: leohoinaski
"""
import matplotlib.pyplot as plt
import numpy as np

def wbdFlux(avWRF,alarea,sRef,ustar,ustarT,ustarTd):
    """
    Parameters
    ----------
    av : numpy array 2D
        porcentagem com vegetação em cada pixel.
    al : numpy array 3D
        porcentagem de cada uso do solo descoberto.
    sRef : numpy array 2D
        is the relative surface area covered with particles with diameter 
        https://agupubs.onlinelibrary.wiley.com/doi/10.1002/jgrd.50313
    ustar : numpy array 2D
        Friction velocity .
    ustarT : numpy array 2D
        threshold friction velocity
    ustarTd : numpy array 2D
        threshold friction velocity based on the particle size.

    Returns
    -------
    Fdust : numpy array 3D
        total dust emission in each computational grid cell .
    
    
    References:
        ρb and ρp are the bulk-soil and soil particle density
        A common range of bulk density of fine-textured surface soils is 1.0–1.3 g/cm3. 
        Coarse-textured surface soils usually have bulk density in the range of 1.3–1.8 g/cm3 (Millar et al., 1965, p. 51–52).
        https://www.sciencedirect.com/topics/agricultural-and-biological-sciences/bulk-density-of-soil#:~:text=A%20common%20range%20of%20bulk,than%201%20g%2Fcm3.
        particle densities for soils range from 2.60 to 2.75 g/cm3 for mineral particles. 
        However, they can be as high as 3.0 g/cm3 for very dense particles and as 
        low as 0.9 g/ cm3 for organic particles.
        Plastic pressure, p (N/m2)	Sand: 5000 Loam: 10000 Sandy clay loam: 10000 
        Clay: 30000
        # ρb and ρp are the bulk-soil and soil particle density
        # Cα and Cβ are constant values
        # cα	0.0002∼0.001 = assume = 0.0006
        # cβ	0.63∼2.09 = assume = 1.36
        # f is the fraction of dust contained in the volume - matriz com proporções no volume
        #Plastic pressure, p (N/m2)	Sand: 5000 Loam: 10000 Sandy clay loam: 10000 Clay: 30000
        
    """
    Ca = 0.0006
    Cb = 1.36
    g = 9.81
    # ===================CUIDADO!!!!
    p = 0.5 # asumi - montar matriz de Plastic pressure com base no solo
    rob = 1.3# g/cm³ - assumi - montar matriz de densidades
    rop = 2.6 # g/cm³ assumi - montar matriz de densidades
    f = 0.2 # Assumi - montar matriz de fração de poeira de um determinado diâmetro para cada tipo de solo
    roa = 1.227 # g/cm³
    c = 1
    Fhd = ((c*roa*(ustar**3))/g)*(1-(ustarTd/ustar))*((1+(ustarTd/ustar))**2)
    Fhd[ustarT>ustar] = 0
    Fhtot = Fhd*sRef
    alpha = (Ca*g*f*rob/(2*p))*(0.24+Cb*ustar*np.sqrt(rop/p))
    Fvtot = alpha*Fhtot
    Fdu = np.empty(Fvtot[:,:,:].shape)
    Fdu[:,:,:] = np.nan
    Fdust = []
    for ii in range(0,alarea.shape[0]):
        for jj in range(0,ustar.shape[0]):
            Fdu[jj,:,:] = Fvtot[jj,:,:]*alarea[ii,:,:]*(1-avWRF[ii,:,:])
        Fdust.append(Fdu)
    Fdust = np.array(Fdust)
    return Fdust



