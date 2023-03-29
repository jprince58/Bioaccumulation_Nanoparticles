# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 13:43:48 2023

@author: joshu
"""

import numpy as np
import scipy as sp

def lambda_calculator():
    
    return lam

def const_calc(n,lamn):
    int1=
    int2=
    return int1/int2



sumt=10
tv=np.arange(0,5,0.1)
xv=np.arange(0,1,0.05)
y=np.zeros((len(t),len(x)))
for t in tv:
    for x in xv:
        yw=0
        for n in np.arange(1,sumt):
            lamn= #Calculate bessel zero
            An=const_calc(n,lamn) #Calculate the An for this n (uses numerical integrator)
            yn=An*np.exp(-lamn**2*t)*sp.special.jv(0,lamn*x)
            yw=yw+yn
            
    