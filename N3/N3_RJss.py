#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_RJss=1.0

import numpy as np

def RJss(x,y,yold,p):
    nx=len(x)-1 #Grab the mesh size for position
    ny=len(y) #Grab number of y-points
    dx=1/nx #Calculate the distance between nodes (assumes domain is from 0 to 1)
    gam=p[0] #Grab the dimenionless ratio of diffusivities
    F=p[1] #Grab the dimensionless forward reaction rate constant
    K=p[2] #Grab the dimensionless Eqilbrium constant for NP binding
    eps=p[3] #Grab the ratio of total NP binding sites to supernatant NP concentration
    omega=p[4] #Grab ratio of NP contribution to electrical potential profile to their electrokinetic mobility
    ups=p[5] #Grab ratio of Biofilm cotridubtion to electrical potential profile to NP electrokinetic mobility
    Kp=p[6] #Grab partition coeffecient of NP into biofilm at water-biofilm interface
    beta=p[7] #Grab ratio of medium mobility from diffusion to medium mobility ffrom eletrokinesis
    R=np.zeros(ny) #Initialize R
    J=np.zeros((ny,ny)) #Initialize J
    index=np.arange(0,ny) #Creater index vector
    for i in index:
        if i==0 : #solid-biofilm interface node for unbound NP balance
            R[i]=2*gam*(y[2]-y[0])/dx**2+2*beta*y[0]*(y[3]-y[1])/dx**2
            J[i,0]=-2*gam/dx**2+2*beta*(y[3]-y[1])/dx**2
            J[i,1]=-2*beta*y[0]/dx**2
            J[i,2]=2*gam/dx**2
            J[i,3]=2*beta*y[0]/dx**2
        elif i==1 : #solid-biofilm interface node for Gauss' Law
            R[i]=2*(y[3]-y[1])/dx**2-1+omega*y[0]*(Kp+eps/(y[0]+K))
            J[i,0]=omega*(Kp+eps*K/(y[0]+K)**2)
            J[i,1]=-2/dx**2
            J[i,3]=2/dx**2
        elif i==ny-2: #water-biofilm interface node for unbound NP balance
            R[i]=y[i]-1
            J[i,i]=1
        elif i==ny-1: #water-biofilm interface node for Gauss' Law
            R[i]=y[i]
            J[i,i]=1
        elif i%2==0 : #interior biofilm nodes for unbound NP balance
            l=int(i/2)
            R[i]=(y[i+2]-2*y[i]+y[i-2])*(x[l]+gam)/dx**2+(y[i+2]-y[i-2])*(1+beta*(y[i+3]-y[i-1])/(2*dx))/(2*dx)+y[i]*beta*(y[i+3]-2*y[i+1]+y[i-1])/dx**2
            J[i,i-2]=(x[l]+gam)/dx**2-(1+beta*(y[i+3]-y[i-1])/(2*dx))/(2*dx)
            J[i,i-1]=beta*y[i]/dx**2-(y[i+2]-y[i-2])/(4*dx**2)
            J[i,i]=-2*(x[l]+gam)/dx**2+beta*(y[i+3]-2*y[i+1]+y[i-1])/dx**2
            J[i,i+1]=-beta*2*y[i]/dx**2
            J[i,i+2]=(x[l]+gam)/dx**2+(1+beta*(y[i+3]-y[i-1])/(2*dx))/(2*dx)
            J[i,i+3]=beta*y[i]/dx**2+(y[i+2]-y[i-2])/(4*dx**2)
        elif i%2==1 : #interior biofilm nodes for Gauss' Law
            l=int((i-1)/2)
            R[i]=(y[i+2]-2*y[i]+y[i-2])/dx**2-(1-x[l])+omega*y[i-1]*(Kp+eps/(y[i-1]+K))
            J[i,i-2]=1/dx**2
            J[i,i-1]=omega*(Kp+eps*K/(y[i-1]+K)**2)
            J[i,i]=-2/dx**2
            J[i,i+2]=1/dx**2

        else:
            print('Uh oh')
    return (R,J,vn_RJss)

"""
Version 1.0

Purpose: Code to return the values of the Residual vector and Jacobian matrix for use in method of lines. Equations can be found in N2.0 documentation

Changes from Version 0.1 to 1.0 (2/20/2022 1:30 am)
    Got initial ocde to run. Had to comment out most functionalities
    

Created on Thu Oct 29 21:55:18 2020

@author: joshuaprince
"""