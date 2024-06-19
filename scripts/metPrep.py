#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

---------------------------metPrep.py------------------------------------------

Esta script é utilizado para determinar as condições meteorlógicas (velocidade
de fricção) e rugosidade da superfície para o domínio de modelagem. Estas condi
ções são utilizadas na estimativa da ressuspensão em cada pixel. 

Utilizei o artigo https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016MS000823
como referência. 

Ref

https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2000JD900304

Created on Mon Mar 18 13:04:17 2024

@author: leohoinaski

"""



import numpy as np
import pandas as pd
import wrf



def ustarCalc(uz,z0):
    """
    Esta função calcula a velocidade frictiva. Não estamos usando esta função
    pois os resultados não são bons. 

    Parameters
    ----------
    uz : numpy array
        velocidade do vento.
        
    z0 : numpy array
        rugosidade. 
        Recebe o dado da função roughness neste mesmo script.
        
    Returns
    -------
    ustar : numpy array
        Velocidade frictiva

    """
    
    k= 0.4 # von karman constant
    
    z= 10.0 # altura de referência
    
    z0[z0<=0]=np.nan # Remove os nans 
    
    print('z0 shape: '+ str(z0.shape))
    
    print('uz shape: '+ str(uz.shape))
    
    #print('z0 shape: '+ str(z0.shape))
    
    ustar = k*uz/np.log(z/z0) # Equação que estima a friction velocity
    
    ustar[np.isnan(z0)]=np.nan # NaN para quando z0 for NaN
    
    ustar[z0==0]=np.nan # NaN para quando z0 for 0
    
    print('ustar shape: '+ str(ustar.shape))

    return ustar

def roughness(tablePath,av,al):
    """
    
    Esta função é utilizada para calcular a rugosidade da superfície de acordo
    com o artigo https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016MS000823

    Parameters
    ----------
    tablePath : path
        Caminho para o arquivo com hs, hv e  alphas de acordo com o artigo 
        https://agupubs.onlinelibrary.wiley.com/doi/10.1002/2016MS000823
        Table 1. Monthly Vegetation Geometric Height 
        and Geometric Height and Density of Nonvegetation (Solid) 
        Elements Based on the Land Type
        
    av : numpy array
        porcentagem de área vegetada.
    al : TYPE
        porcentagem de área com o soilId.

    Returns
    -------
    z0 : np.array
        matriz com rugosidades no domínio.
    alphaV : np.array
        matriz com alphaV no domínio..
    alphaS : np.array
        matriz com alphaS no domínio..

    """
    
    # lendo csv com os dados de hs de acordo com o artigo
    alpS = pd.read_csv(tablePath + '/hs.csv')
    
    # lendo csv com os dados de hv de acordo com o artigo
    hv = pd.read_csv(tablePath + '/hv.csv')
    
    # estimando o alphaS com base na proporção de área vegetada e não vegetada
    alphaS = av*alpS['alphaS'][0] + np.nansum(al,axis=0)*alpS['alphaS'][2] 
    
    # estimando alphaV de acordo com o artigo
    alphaV = -0.35*np.log(1-av)
    
    # removing inf values 
    alphaV[alphaV == np.inf] = np.nanmax(alphaV)
    
    # calculando alpha total
    alpha = alphaS + alphaV
    
    # estimando hs pela ponderação entre area vegetada e não vegetada
    hS= av*alpS['hs'][0] + al*alpS['hs'][2]
    
    # estimando hv pela ponderação entre area vegetada e não vegetada
    # usa a média dos meses 
    hV = av*hv['mean'][0] + al*hv['mean'][2]
    
    # estimando h conforme artigo
    h = (hV*alphaV + hS*alphaS)/(alphaS + alphaV)
    
    # inicializando z0
    z0 = alpha.copy() 
    # z0[alpha<0.2] = h[alpha<0.2]*0.96*alpha[alpha<0.2]**(1.07)
    # z0[alpha>=0.2] = h[alpha>=0.2]*0.083*alpha[alpha>=0.2]**(-0.46)
    
    # calculando o z0 de acordo com o artigo para as duas condições de alpha
    z0[alpha<0.2] = (0.96*alpha[alpha<0.2]**(1.07))*h[alpha<0.2]
    z0[alpha>=0.2] = (0.083*alpha[alpha>=0.2]**(-0.46))*h[alpha>=0.2]
    
    # remove nans
    z0[:,np.isnan(al)]=np.nan
    alphaV[:,np.isnan(al)]=np.nan
    alphaS[:,np.isnan(al)]=np.nan
    
    return z0,alphaV,alphaS

def ustarThreshold(D,clayRegrid,w,alphaV,alphaS,av):
    """
    Esta função calcula o limiar que define quando vai ocorrer a ressuspensão, 
    com base na velocidade frictiva ustarThreshold. 

    Parameters
    ----------
    D : int
        diâmetro da partícula.
    clayRegrid : np.array
        porporção de argila - matriz de clayContent .
    w : np.array
        matriz de umidade do solo de acordo com o output do WRF.
    alphaV : np.array
        roughness densities vegetation elements.
    alphaS : np.array
        roughness densities based on solid (nonvegetation)
    av : np.array
        matriz de proporção de área vegetada.

    Returns
    -------
    ustarT : np.array
        matriz da velocidade frictiva threshold.
    ustarTd : np.array
        matriz de ideal threshold friction velocity based on the particle size.

    """
    
    # constante (accounting for the magnitude of the interparticle cohesive forces) 
    An = 0.0123
    
    # constante (accounting for the magnitude of the interparticle cohesive forces) 
    tao = 1.65*(10**(-4)) #  kg s-2
    
    # densidade do ar
    roa = 1.227# g/cm³
    
    # densidade da particula
    # mantive esta unidade após comparar com o gráfico da figura 3 do artigo
    # https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2000JD900304
    rop = 2665 # kg/m³ assumi - montar matriz de densidades
    
    # gravidade
    g = 9.81 # m/s2
        
    # separando a equação em dois termos
    t1 = rop*g*(D*10**(-6))/roa
    t2 = tao/(roa*D*10**(-6))
    
    # equação para estimar ustarTd, conforme o artigo
    ustarTd = np.sqrt(An*(t1+t2))
    
    # Usar para teste e comparação com artigo
    #D = np.arange(0.1,1000)
    # plt.plot(D,ustarTd)
    # plt.xscale('log')
        
    # calculo da umidade do solo, conforme o artigo
    wl = 0.0014*(clayRegrid**2)+0.17*clayRegrid
    
    # Soil moisture increases the threshold friction velocity through 
    # increasing the cohesive forces suppressing the mobilization of particles.
    # inicializando a matriz de fm
    fm = np.empty(w.shape)
    fm[:,:,:]=np.nan
    fm=fm.astype(float)
    
    # estimando conforme o artigo
    fm = (1+1.21*(w-wl)**(0.68))**(0.5)
    fm[w<wl] = 1.0 
    
    # constantes
    sigmaV = 1.45
    mV = 0.16
    betaV = 202
    sigmaS = 1.0
    mS = 0.5
    betaS = 90
    
    # para teste
    # av = np.random.rand(1000)
    # al = 1-av
    # alphaV = -0.35*np.log(1-av)
    # tablePath = '/mnt/sdb1/windBlowDustBR/inputs/tables'
    # alpS = pd.read_csv(tablePath + '/hs.csv')
    # hv = pd.read_csv(tablePath + '/hv.csv')
    # alphaS = av*alpS['alphaS'][0] + al*alpS['alphaS'][2] 
    
    #===========================VERIFICAR!!!
    fr = ((1-sigmaV*mV*alphaV)*(1+betaV*mV*alphaV)*(1-sigmaS*mS*(alphaS/(1-av)))*(1+betaS*mS*(alphaS/(1-av))))**0.5
    # Referência do fr = https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2010JD014649
    #fr2 = ((1-sigmaS*mS*(alphaS/(1-av)))*(1+betaS*mS*(alphaS/(1-av))))**0.5 # ref = )
    
    # Estimativa do threshold
    ustarT=ustarTd*fm*fr
    
    return ustarT,ustarTd

def main(ds,tablePath,av,al,D,clayRegrid,lia):
    """
    Esta função controla todo o script de metPrep.

    Parameters
    ----------
    ds : netCDF object 
        objeto com o netCDF do WRF
    tablePath : path
        caminho para a pasta com as tabelas do artigo.
    av : np.array
        proporção de área vegetada.
    al : np.array
        proporção de área de cada soilId.
    D : int
        diâmetro da partícula.
    clayRegrid : np.array
        matriz com os dados de claycontent.
    lia : np.array
        vetor com datas compatíveis entre MCIP e WRF.

    Returns
    -------
    ustar : np.array
        matriz com valores de ustar para cada tempo.
    ustarT : np.array
        matriz de ustar threshold.
    ustarTd : np.array
        matriz de ustar threshold.
    av : np.array
        matriz de proporção de área vegetada do WRF.
    ustarWRF : np.array
        matriz de ustar do WRF.

    """
    
    # Extraindo dado de fração de vegetação no domínio usado no WRF
    avWRF = ds['VEGFRA'][lia,:,:]/100
    
    # extraindo dados de ustar do WRF
    ustarWRF = ds['UST'][lia,:,:]
    
    # estimando rugosidade
    z0,alphaV,alphaS = roughness(tablePath,avWRF,al)
    
    # extraindo velocidade do vento do WRF
    uz = np.array(wrf.g_wind.get_destag_wspd_wdir10(ds,timeidx=wrf.ALL_TIMES)[0,lia,:,:])
    
    # calculando ustar de acordo com artigo
    ustar = ustarCalc(uz,z0)
    
    # extraindo a umidade do solo do WRF
    w = ds['SMOIS'][lia,0,:,:]
    
    # estimando ustarT e ustarTd de acordo com o artigo
    ustarT,ustarTd = ustarThreshold(D,clayRegrid,w,alphaV,alphaS,avWRF)
    
    return ustar,ustarT,ustarTd,avWRF,ustarWRF
