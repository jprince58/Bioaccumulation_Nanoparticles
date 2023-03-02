#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_method_of_lines=1.0

import numpy as np
import scipy as sp
from scipy import sparse
from scipy.sparse import linalg

from N3_RJ import *

def method_of_lines(t,x,y,h,p,tol,fit_coeff):
    yw=y[:,0] #Initalize the working concentration vector
    index=np.arange(0,len(t)) #Creater index vector
    whoops=0 #initialize error function
    for i in index: #Begin for loop which iterates over all the entries in the time vector, assigning them the value td
        if i==0: continue
        else:
            yold=yw #update the old y-value
            [R,J,vn_RJ]=RJ(x,yw,yold,p,h,t[i],fit_coeff); #Calculate Residual and Jacobian from new y value
            k=0
            while np.linalg.norm(R)>tol :
                k=k+1
                current_tol=np.linalg.norm(R)
                J=sp.sparse.csc_matrix(J)
                dif=-sp.sparse.linalg.spsolve(J,R) #Apply built in sparse Linear solver to find delta from J and R
                yw=yw+dif ; #Update y
                [R,J,nv_RJ]=RJ(x,yw,yold,p,h,t[i],fit_coeff); #Calculate Residual and Jacobian from new y value
                if k>100:
                    print('Whoops')
                    whoops=whoops+1
                    break
            y[:,i]=yw
    return (y,whoops,vn_method_of_lines,vn_RJ)

"""
Purpose: Script for running method of lines to solve the equations for the model

Version 1.0

Changes from Version 0.1 to 1.0 (2/20/2022 1:30 am)
    Got initial ocde to run. 
    

Created on Sat Oct 31 17:25:29 2020

@author: joshuaprince
"""