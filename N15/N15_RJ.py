#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_RJ=1.1

import numpy as np
from N15_BCcalc import *

def RJ(x,y,p):
#def RJ(x,y,p,t):
    nx=len(x)-1 #Grab the mesh size for position
    ny=len(y)-2 #Grab number of y-points
    dx=1/nx #Calculate the distance between nodes (assumes domain is from 0 to 1)
    gam=p[0] ##Define effective diffusivity  for this iteration
    mu=p[1] #Define dimensionless and porosity adjusted binding rate constant for this iteration
    nu=p[2] ##Define dimensionless and porosity adjusted maximum binding site density for this iteration
    eps=p[3] #Define dimensionless minimum interstital porosity for this iteration
    rho=p[4] #Define dimensionless minimum traditional porosity for this iteration
    # phi_min=p[4] #Define liquid-interface interstitial porosity
    kappa=p[5]  #Define dimensionless and porosity adjusted equilibrium constant  for this iteration
    a=p[6] #Define shape paramter for binding site profile  for this iteration
    b=p[7]  #Define shape paramter for intersitital porosity profile  for this iteration
    # rho=phi_min+(np.exp(b)-1)*eps #Define phi_max
    beta=p[8]  #Define shape paramter for traditional proosity profile  for this iteration
    Kp=p[9] #Define partition coefecient 

    R=np.zeros(ny+2) #Initialize R
    J=np.zeros((ny+2,ny+2)) #Initialize J
    index=np.arange(0,ny+2) #Creater index vector
    for i in index:
        if i==0 :
            l=int(i/2)
            R[i]=2*gam/dx**2*(y[i+2]-y[i])-y[i]*(mu*(nu*x[l]**a+beta-y[i+1]))+mu*kappa*y[i+1]
            J[i,i]=-2*gam/dx**2-mu*(nu*x[l]**a+beta-y[i+1])
            J[i,i+2]=2*gam/dx**2
            J[i,i+1]=mu*(y[i]+kappa)
        elif i%2==1 :
            l=int((i-1)/2)
            R[i]=mu*(y[i-1]*(nu*x[l]**a+beta-y[i])-kappa*y[i])
            J[i,i]=-mu*(y[i-1]+kappa)
            J[i,i-1]=mu*(nu*x[l]**a+beta-y[i])
        elif i==ny:
            l=int(i/2)
            R[i]=y[i]-1
            J[i,i]=1
            # phi0=rho-eps*(np.exp(b*x[l])-1)
            # R[i]=(Kp*y[i]+y[i+1])-BCcalc(t,fit_coeff,tmax,kconv)/phi0
            # test=BCcalc(t,fit_coeff,tmax,kconv)
            # J[i,i]=Kp
            # J[i,i+1]=1
        elif i%2==0 and i!=0:
            l=int(i/2)
            phi0=rho-eps*(np.exp(b*x[l])-1)
            phi1=-b*eps*np.exp(b*x[l])
            R[i]=gam/dx**2*(y[i+2]-2*y[i]+y[i-2])+gam*phi1/phi0/(2*dx)*(y[i+2]-y[i-2])-y[i]*(mu*(nu*x[l]**a+beta-y[i+1]))+mu*kappa*y[i+1]
            J[i,i]=-2*gam/dx**2-mu*(nu*x[l]**a+beta-y[i+1])
            J[i,i+2]=gam/dx**2+gam*phi1/phi0/(2*dx)
            J[i,i-2]=gam/dx**2-gam*phi1/phi0/(2*dx)
            J[i,i+1]=mu*(y[i]+kappa)
        else:
            print('Uh oh')
    return (R,J,vn_RJ)

"""
Version 1.1

Purpose: Code to return the values of the Residual vector and Jacobian matrix for use in method of lines. Equations can be found in N6.0 documentation

Created on Thu Oct 29 21:55:18 2020

@author: joshuaprince
"""