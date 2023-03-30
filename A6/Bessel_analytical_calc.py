# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 13:43:48 2023

@author: joshu
"""

import numpy as np
import scipy
import scipy.special as sp
import matplotlib.pyplot as plt

sumt=5
epsilon=0.2
tv=np.arange(0,2,0.02)
xv=np.arange(epsilon,epsilon+1,0.05)
const_coeff=np.zeros((6,3))
const_coeff[0,0]=3.25093 
const_coeff[0,1]=27.4144
const_coeff[0,2]=0.165236
const_coeff[1,0]=5.8317
const_coeff[1,1]=-1.51014
const_coeff[1,2]=-1.52909
const_coeff[2,0]=8.32035
const_coeff[2,1]=-0.159299
const_coeff[2,2]=-2.70904
const_coeff[3,0]=10.7796
const_coeff[3,1]=0.905449
const_coeff[3,2]=-1.38796
const_coeff[4,0]=13.2256
const_coeff[4,1]=-11.2964
const_coeff[4,2]=0.180238
const_coeff[5,0]=15.6645
const_coeff[5,1]=-0.604154
const_coeff[5,2]=1.60086



y=np.zeros((len(tv),len(xv)))
ti=0
for t in tv:
    xi=0
    for x in xv:
        yw=0
        for n in np.arange(0,sumt):
            [lamn,alphan,An]=const_coeff[n,:] #Calculate bessel zero
            yn=An*np.exp(-lamn**2*t)*x*(sp.j0(2*lamn*x**(1/2))-alphan*sp.y0(2*lamn*x**(1/2)))
            yw=yw+yn
        y[ti,xi]=yw
        xi=xi+1
    ti=ti+1
    
plt.figure(0) 
tp_u=10 #number of time points to plot
nt=len(tv)
lognt_u=np.log10(nt) #Logarthmic timepoints
logspace_u=round((lognt_u)/tp_u,10) #Logarthmic timepoints
logtindex_u=np.arange(0,lognt_u,logspace_u) #Logarthmic timepoints
for logi_u in logtindex_u:
    i_u=int(10**logi_u)
    plt.plot(xv,y[i_u,:],label=f't={tv[i_u]}')
    plt.legend()
    
            
    