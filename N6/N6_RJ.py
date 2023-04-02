#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_RJ=1.1

import numpy as np
from N6_BCcalc import *

def RJ(x,y,p,fit_coeff,t,tmax,kconv):
#def RJ(x,y,p,t):
    nx=len(x)-1 #Grab the mesh size for position
    ny=len(y)-2 #Grab number of y-points
    dx=1/nx #Calculate the distance between nodes (assumes domain is from 0 to 1)
    omega=p[0] ##Define effective diffusivity  for this iteration
    mu=p[1] #Define dimensionless and porosity adjusted binding rate constant for this iteration
    nu=p[2] ##Define dimensionless and porosity adjusted maximum binding site density for this iteration
    eps=p[3] #Define dimensionless minimum interstital porosity for this iteration
    rho=p[4] #Define dimensionless minimum traditional porosity for this iteration
    kappa=p[5]  #Define dimensionless and porosity adjusted equilibrium constant  for this iteration
    a=p[6] #Define shape paramter for binding site profile  for this iteration
    b=p[7]  #Define shape paramter for intersitital porosity profile  for this iteration
    c=p[8]  #Define shape paramter for traditional proosity profile  for this iteration
    Kp=p[9] #Define partition coefecient 

    R=np.zeros(ny+2) #Initialize R
    J=np.zeros((ny+2,ny+2)) #Initialize J
    index=np.arange(0,ny+2) #Creater index vector
    for i in index:
        if i==0 :
            l=int(i/2)
            #Calculate intermediate values because these are getting too long of lines of code
            intm1=omega*(x[l]**b+eps)/dx**2 
            intm2=omega/(2*dx)*(b*x[l]**(b-1)-2*c*x[l]**(c-1)*(x[l]**b+eps)/(x[l]**c+rho)) #Calcualte intermediate value 2 
            intm3=omega*((c*b*x[l]**(c+b-2)+c*(c-1)*x[l]**(c-2)*(x[l]**b+eps))/(x[l]**c+rho)-2*c**2*x[l]**(2*c-2)*(x[l]**b+eps)/(x[l]**c+rho)**2)
            intm4=mu*(nu*(1-x[l]**a)-y[i+1])
            intm5=mu*kappa
            R[i]=intm1*(2*y[i+2]-2*y[i])+intm2*(y[i+2]-y[i-2])-(intm3+intm4)*y[i]+intm5*y[i+1]
            J[i,i]=-2*intm1-intm3-intm4
            J[i,i+2]=intm1+intm2
            J[i,i-2]=intm1-intm2
            J[i,i+1]=mu*(y[i]+kappa)
        elif i%2==1 :
            l=int((i-1)/2)
            R[i]=mu*(y[i-1]*(nu*(1-x[l]**a)-y[i])-kappa*y[i])
            J[i,i]=-mu*(y[i-1]+kappa)
            J[i,i-1]=mu*(nu*(1-x[l]**a)-y[i])
        elif i==ny:
            l=int(i/2)
            #R[i]=y[i]-1-rho
            R[i]=y[i]-BCcalc(t,fit_coeff,tmax,kconv)*(1+rho)/Kp
            # test=BCcalc(t,fit_coeff,tmax)+rho
            J[i,i]=1;
        elif i%2==0 and i!=0:
            l=int(i/2)
            #Calculate intermediate values because these are getting too long of lines of code
            intm1=omega*(x[l]**b+eps)/dx**2 
            intm2=omega/(2*dx)*(b*x[l]**(b-1)-2*c*x[l]**(c-1)*(x[l]**b+eps)/(x[l]**c+rho)) #Calcualte intermediate value 2 
            intm3=omega*((c*b*x[l]**(c+b-2)+c*(c-1)*x[l]**(c-2)*(x[l]**b+eps))/(x[l]**c+rho)-2*c**2*x[l]**(2*c-2)*(x[l]**b+eps)/(x[l]**c+rho)**2)
            intm4=mu*(nu*(1-x[l]**a)-y[i+1])
            intm5=mu*kappa
            R[i]=intm1*(y[i+2]-2*y[i]+y[i-2])+intm2*(y[i+2]-y[i-2])-(intm3+intm4)*y[i]+intm5*y[i+1]
            J[i,i]=-2*intm1-intm3-intm4
            J[i,i+2]=intm1+intm2
            J[i,i-2]=intm1-intm2
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