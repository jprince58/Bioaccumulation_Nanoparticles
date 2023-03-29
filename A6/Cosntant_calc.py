# -*- coding: utf-8 -*-
"""
Created on Sun Mar 26 13:43:48 2023

@author: joshu
"""

import numpy as np
import scipy
from scipy import sparse
from scipy.sparse import linalg
import scipy.special as sp

def RJ(y,epsilon):
    #calcualte residual and jacobian
    eps=epsilon**(1/2)
    eps_1=(1+epsilon)**(1/2)
    R=np.zeros(2)
    J=np.zeros((2,2))
    R[0]=y[1]*eps*(y[0]*sp.y1(2*y[1]*eps)-sp.j1(2*y[1]*eps))+sp.j0(2*y[1]*eps)-y[0]*sp.y0(2*y[1]*eps)
    R[1]=sp.j0(2*y[1]*eps_1)-y[0]*sp.y0(2*y[1]*eps_1)
    J[0,0]=eps*y[1]*sp.y1(2*y[1]*eps)-sp.y0(2*y[1]*eps)
    term1=-2*eps*sp.j1(2*y[1]*eps)
    term2=2*eps*y[0]*sp.y1(2*y[1]*eps)
    term3=eps*(y[0]*sp.y1(2*y[1]*eps)-sp.j1(2*y[1]*eps))
    term4=y[1]*epsilon*y[0]*(sp.y0(2*y[1]*eps)-sp.yn(2,2*y[1]*eps))
    term5=-y[1]*epsilon*(sp.j0(2*y[1]*eps)-sp.jn(2,2*y[1]*eps))
    J[0,1]=term1+term2+term3+term4+term5
    J[1,0]=sp.y0(2*y[1]*eps_1)
    J[1,1]=-2*eps_1*sp.j1(2*y[1]*eps_1)+2*y[0]*eps_1*sp.y1(2*y[1]*eps_1)
    
    return [R,J]

def Newton_Rhapson(y,epislon):
    yw=y #Initalize the working concentration vector
    whoops=0 #initialize error function
    [R,J]=RJ(yw,epsilon); #Calculate Residual and Jacobian from new y value
    k=0
    while np.linalg.norm(R)>10**-8 :
        k=k+1
        current_tol=np.linalg.norm(R)
        dif=-np.linalg.solve(J,R) #Apply built in sparse Linear solver to find delta from J and R
        yw=yw+dif ; #Update y
        [R,J]=RJ(yw,epsilon); #Calculate Residual and Jacobian from new y value
        if k>100:
            print('Whoops')
            whoops=whoops+1
            break
    y=yw
    return (y,whoops)

epsilon=0.2
y0=[1,1]
first_guess=Newton_Rhapson(y0,epsilon)

#Going to start by calculating a value using alpha and lamba equal one,

#define inital values to test
#Initialize the solution vector

#loop over initial values to test, a while loop until 4 smallest values are found

 


# [lamn,alpha]=const_calc(epsilon,) #Calculate bessel zero

# #Export values as a .csv file so I don't need to recalculate them everytime.

# epsilon=
# [lamn,alpha]=const_calc(epsilon) #Calculate bessel zero

# sumt=10
# tv=np.arange(0,5,0.1)
# xv=np.arange(0,1,0.05)
# y=np.zeros((len(t),len(x)))
# for t in tv:
#     for x in xv:
#         yw=0
#         for n in np.arange(1,sumt):
#             An=const_calc(n,lamn) #Calculate the An for this n (uses numerical integrator)
#             yn=An*np.exp(-lamn**2*t)*sp.special.jv(0,lamn*x)
#             yw=yw+yn
            
    