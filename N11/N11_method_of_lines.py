#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_method_of_lines=1.1

import numpy as np
import scipy as sp
from scipy import sparse
from scipy.sparse import *

from N11_RJ import *
from N11_BCcalc import *

def method_of_lines(t,x,y,h,p,tol,fit_coeff,kconv):
#def method_of_lines(t,x,y,h,p,tol):
    yw=y[:,0] #Initalize the working concentration vector
    index=np.arange(0,len(t)) #Creater index vector
    whoops=0 #initialize error function
    tmax=t[-1]
    rho=p[4] #Define dimensionless minimum traditional porosity for this iteration
    Kp=p[9] #Define partition coefecient 
    for i in index: #Begin for loop which iterates over all the entries in the time vector, assigning them the value td
        yold=yw #update the old y-value
        #yw[-2]=BCcalc(t[i],fit_coeff,tmax,kconv)*(1+rho)/Kp #Pass along boudnary condition as initial guess
        #think I need to hard-code in the BC as the 
        [R,J,vn_RJ]=RJ(x,yw,p,fit_coeff,t[i],tmax,kconv); #Calculate Residual and Jacobian from new y value
        #[R,J,vn_RJ]=RJ(x,yw,p,t[i]); #Calculate Residual and Jacobian from new y value
        R=yw-yold-h*R
        k=0
        while np.linalg.norm(R)>tol :
            k=k+1
            J=np.eye(len(yw))-h*J; #Calculate new Jacobian from new y
            J=sp.sparse.csc_matrix(J)
            dif=-sp.sparse.linalg.spsolve(J,R) #Apply built in sparse Linear solver to find delta from J and R
            yw=yw+dif ; #Update y
            [R,J,vn_RJ]=RJ(x,yw,p,fit_coeff,t[i],tmax,kconv); #Calculate Residual and Jacobian from new y value
            #[R,J,nv_RJ]=RJ(x,yw,p,t[i]); #Calculate Residual and Jacobian from new y value
            R=yw-yold-h*R ; #Update Residual
            if k>100:
                print('Whoops')
                whoops=whoops+1
                break
        y[:,i]=yw
    return (y,whoops,vn_method_of_lines,vn_RJ)

"""
Purpose: Script for running method of lines to solve the equations for the model

Version 1.1

Created on Sat Oct 31 17:25:29 2020

@author: joshuaprince
"""