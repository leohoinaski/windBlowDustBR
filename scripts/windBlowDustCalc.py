#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 15:20:44 2024

@author: leohoinaski
"""
import numpy as np

def wbdCalc(av,al,Fhd,Sd):
    Fhtot = Fhd*Sd
    alpha = (Calpha*g*rob/(2*p))*(0.24+Cb*ustar*np.sqrt(rop/p))
    Fvtot = alpha*Fhtot
    Fdust = (1-av)*np.nansum(al*Fvtot)
    return Fdust