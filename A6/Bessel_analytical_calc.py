# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 13:43:48 2023

@author: joshu
"""

import numpy as np
import scipy as sp

def getRJ(x,y):
    #calcualte residual and jacobian
    
    return

def Newton_Rhapson(epislon):
    
    
    #Newton Rhapson Structure
    
    return (alpha, lam)

def const_calc(epsilon):
    
    #For loop over a series of initial guesses
        #(lam, alpha)=Newton-Rhpason function call (, initial guess iepsilon)
    #Sorting function to not include in final output if these are redundant
    return #Spit out all the lambda and corresponding alpha's I was able to find in increasing order.



epsilon=
[lamn,alpha]=const_calc(epsilon) #Calculate bessel zero

sumt=10
tv=np.arange(0,5,0.1)
xv=np.arange(0,1,0.05)
y=np.zeros((len(t),len(x)))
for t in tv:
    for x in xv:
        yw=0
        for n in np.arange(1,sumt):
            [lamn,alpha]=const_calc(epsilon,) #Calculate bessel zero
            An=const_calc(n,lamn) #Calculate the An for this n (uses numerical integrator)
            yn=An*np.exp(-lamn**2*t)*sp.special.jv(0,lamn*x)
            yw=yw+yn
            
    