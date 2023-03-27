#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_RJ=1.1

import numpy as np

def RJ(x,y,p):
    nx=len(x)-1 #Grab the mesh size for position
    ny=len(y)-2 #Grab number of y-points
    dx=1/nx #Calculate the distance between nodes (assumes domain is from 0 to 1)
    gam=p[0] ##Define effective diffusivity  for this iteration
    alp=p[1] #Define dimensionless and porosity adjusted binding rate constant for this iteration
    xi=p[2] ##Define dimensionless and porosity adjusted maximum binding site density for this iteration
    eps=p[3] #Define dimensionless minimum interstital porosity for this iteration
    rho=p[4] #Define dimensionless minimum traditional porosity for this iteration
    K=p[5]  #Define dimensionless and porosity adjusted equilibrium constant  for this iteration
    a=p[6] #Define shape paramter for binding site profile  for this iteration
    b=p[7]  #Define shape paramter for intersitital porosity profile  for this iteration
    c=p[8]  #Define shape paramter for traditional proosity profile  for this iteration

    R=np.zeros(ny+2) #Initialize R
    J=np.zeros((ny+2,ny+2)) #Initialize J
    index=np.arange(0,ny+2) #Creater index vector
    for i in index:
        if i==0 :
            l=int(i/2)
            por_p_0=np.tanh(c*x[l])+rho #Pore-volume based porosity- 0th derivative
            por_p_1=c*(1/np.cosh(c*x[l]))**2 #Pore-volume based porosity- 1st derivative
            por_p_2=-2*c**2*np.tanh(c*x[l])*(1/np.cosh(c*x[l]))**2 #Pore-volume based porosity- 2nd derivative
            por_i_0=np.tanh(b*x[l])+eps #Interstitial-volume based porosity- 0th derivative
            por_i_1=b*(1/np.cosh(b*x[l]))**2 #Interstitial-volume based porosity- 1st derivative
            beta1=por_i_0
            beta2=por_i_1-2*por_i_0*por_p_1/por_p_0
            beta3=2*por_i_0*por_p_1**2/por_p_0**2-(por_i_1*por_p_1+por_i_0*por_p_2)/por_p_0
            R[i]=2*gam*beta1*(y[i+2]-y[i])/dx**2+y[i]*(gam*beta3-alp*(xi-y[i+1]))+alp*K*y[i+1]
            J[i,i]=-2*gam*beta1/dx**2+gam*beta3-alp*(xi-y[i+1])
            J[i,i+2]=2*gam*beta1/dx**2
            J[i,i+1]=alp*(y[i]+K)
        elif i%2==1 :
            l=int((i-1)/2)
            R[i]=alp*(y[i-1]*(xi*(1-x[l]**a)-y[i])-K*y[i])
            J[i,i]=-alp*(y[i-1]+K)
            J[i,i-1]=alp*(xi*(1-x[l]**a)-y[i])
        elif i==ny:
            l=int(i/2)
            R[i]=y[i]-1
            J[i,i]=1;
        elif i%2==0 and i!=0:
            l=int(i/2)
            por_p_0=np.tanh(c*x[l])+rho #Pore-volume based porosity- 0th derivative
            por_p_1=c*(1/np.cosh(c*x[l]))**2 #Pore-volume based porosity- 1st derivative
            por_p_2=-2*c**2*np.tanh(c*x[l])*(1/np.cosh(c*x[l]))**2 #Pore-volume based porosity- 2nd derivative
            por_i_0=np.tanh(b*x[l])+eps #Interstitial-volume based porosity- 0th derivative
            por_i_1=b*(1/np.cosh(b*x[l]))**2 #Interstitial-volume based porosity- 1st derivative
            beta1=por_i_0
            beta2=por_i_1-2*por_i_0*por_p_1/por_p_0
            beta3=2*por_i_0*por_p_1**2/por_p_0**2-(por_i_1*por_p_1+por_i_0*por_p_2)/por_p_0
            R[i]=gam*beta1*(y[i+2]-2*y[i]+y[i-2])/dx**2+gam*beta2*(y[i+2]-y[i-2])/(2*dx)+y[i]*(gam*beta3-alp*(xi*(1-x[l]**a)-y[i+1]))+alp*K*y[i+1]
            J[i,i]=-2*gam*beta1/dx**2+gam*beta3-alp*(xi*(1-x[l]**a)-y[i+1])
            J[i,i+2]=gam*beta1/dx**2+gam*beta2/(2*dx)
            J[i,i-2]=gam*beta1/dx**2-gam*beta2/(2*dx)   
            J[i,i+1]=alp*(y[i]+K)
        else:
            print('Uh oh')
    return (R,J,vn_RJ)

"""
Version 1.1

Purpose: Code to return the values of the Residual vector and Jacobian matrix for use in method of lines. Equations can be found in N6.0 documentation

Created on Thu Oct 29 21:55:18 2020

@author: joshuaprince
"""