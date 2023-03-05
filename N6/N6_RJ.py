#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_RJ=1.1

import numpy as np

def RJ(x,y,p):
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

    R=np.zeros(ny+2) #Initialize R
    J=np.zeros((ny+2,ny+2)) #Initialize J
    index=np.arange(0,ny+2) #Creater index vector
    for i in index:
        if i==0 :
            l=int(i/2)
            R[i]=2*omega*(x[l]**b+eps)/dx**2*(y[i+2]-y[0])-y[i]*mu/(x[l]**b+eps)*(nu*(1-x[l]**a)-y[i+1])+mu*y[i+1]*kappa/(x[l]**b+eps)
            J[i,i]=-2*omega*(x[l]**b+eps)/dx**2-mu/(x[l]**b+eps)*(nu*(1-x[l]**a)-y[i+1])
            J[i,i+1]=mu*(y[i-1]+kappa)/(x[l]**b+eps)
            J[i,i+2]=2*omega*(x[l]**b+eps)/dx**2
        elif i%2==1 :
            l=int((i-1)/2)
            R[i]=mu/(x[l]**b+eps)*(y[i-1]*(nu*(1-x[l]**a)-y[i])-kappa*y[i])
            J[i,i]=mu/(x[l]**b+eps)*(nu*(1-x[l]**a)-y[i])
            J[i,i-1]=-mu*(y[i-1]+kappa)/(x[l]**b+eps)
        elif i==ny:
            l=int(i/2)
            R[i]=y[i]-(1+rho)
            J[i,i]=1;
        elif i%2==0 and i!=0:
            l=int(i/2)
            #Calculate intermediate values because these are getting too long of lines of code
            intm1=omega*(x[l]**b+eps)/dx**2 
            intm2=omega/(2*dx)*(b*x[l]**(b-1)-2*c*x[l]**(c-1)*(x[l]**b+eps)/(x[l]**c+rho)) #Calcualte intermediate value 2 
            intm3=omega*((c*b*x[l]**((c-1)*(b-1))+c*(c-1)*x[l]**(c-2)*(x[l]**b+eps))/(x[l]**c+rho)-2*c**2*x[l]**(2*c-2)*(x[l]**b+eps)/(x[l]**c+rho)**2)
            intm4=mu/(x[l]**b+eps)*(nu*(1-x[l]**a)-y[i+1])
            intm5=mu*kappa/(x[l]**b+eps)
            R[i]=intm1*(y[i+2]-2*y[i]+y[i-2])+intm2*(y[i+2]-y[i-2])-(intm3+intm4)*y[i]+intm5*y[i+1]
            J[i,i]=-2*intm1-intm3-intm4
            J[i,i+2]=intm1+intm2
            J[i,i-2]=intm1-intm2
            J[i,i+1]=mu*(y[i-1]+kappa)/(x[l]**b+eps)
        else:
            print('Uh oh')
    return (R,J,vn_RJ)

"""
Version 1.1

Purpose: Code to return the values of the Residual vector and Jacobian matrix for use in method of lines. Equations can be found in N2.0 documentation

Created on Thu Oct 29 21:55:18 2020

@author: joshuaprince
"""