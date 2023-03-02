#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_RJ=1.1

import numpy as np

def RJ(x,y,p):
    nx=len(x)-1 #Grab the mesh size for position
    ny=len(y)-2 #Grab number of y-points
    dx=1/nx #Calculate the distance between nodes (assumes domain is from 0 to 1)
    gam=p[0] #Grab the dimenionless ratio of diffusivities
    beta=p[1] #Grab the dimensionless ratio of potentials
    F=p[2] #Grab the dimensionless forward reaction rate constant
    Re=p[3] #Grab the dimensionless reverse reaction rate constant
    n=p[4] #Grab the hill coeffecient
    R=np.zeros(ny+2) #Initialize R
    J=np.zeros((ny+2,ny+2)) #Initialize J
    index=np.arange(0,ny+2) #Creater index vector
    for i in index:
        if i==0 :
            R[i]=y[i]-1
            J[i,i]=1;
        elif i%2==1 :
            R[i]=F*y[i-1]**n*(1-y[i])-Re*y[i]
            J[i,i]=-(F*y[i-1]**n+Re)
            J[i,i-1]=n*F*y[i-1]**(n-1)*(1-y[i])
        elif i==ny:
            l=int(i/2)
            R[i]=2*(y[i-2]-y[i])/dx**2*(1-x[l]**2+gam)-y[i]*(F*y[i]**(n-1)*(1-y[i+1])-2*beta)+Re*y[i+1]
            J[i,i]=(-2/dx**2)*(1-x[l]**2+gam)-(n*F*y[i]**(n-1)*(1-y[i-1]-2*beta))
            J[i,i-2]=2/dx**2*(1-x[l]**2+gam)
            J[i,i+1]=F*y[i]**n+Re
        elif i%2==0 and i!=0:
            l=int(i/2)
            R[i]=(y[i+2]-2*y[i]+y[i-2])/dx**2*(1-x[l]**2+gam)-(y[i+2]-y[i-2])/2/dx*(2*x[l]*(1-beta))-y[i]*(F*y[i]**(n-1)*(1-y[i+1])-2*beta)+Re*y[i+1]
            J[i,i]=(-2/dx**2)*(1-x[l]**2+gam)-(n*F*y[i]**(n-1)*(1-y[i-1]-2*beta))
            J[i,i+2]=1/dx**2*(1-x[l]**2+gam)-1/2/dx*(2*x[l]*(1-beta))
            J[i,i-2]=1/dx**2*(1-x[l]**2+gam)+1/2/dx*(2*x[l]*(1-beta))
            J[i,i+1]=F*y[i]**n+Re
        else:
            print('Uh oh')
    return (R,J,vn_RJ)

"""
Version 1.1

Purpose: Code to return the values of the Residual vector and Jacobian matrix for use in method of lines. Equations can be found in N2.0 documentation

Created on Thu Oct 29 21:55:18 2020

@author: joshuaprince
"""