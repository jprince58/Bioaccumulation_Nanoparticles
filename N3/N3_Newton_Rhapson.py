# -*- coding: utf-8 -*-

vn_Newton_Rhapson=1.0

import numpy as np
import scipy as sp
from scipy import sparse
from scipy.sparse import linalg

from N3_RJss import *

def Newton_Rhapson(x,y,p,tol):
    yw=y #Initalize the working concentration vector
    whoops=0 #initialize error function
    yold=yw #update the old y-value
    [R,J,vn_RJss]=RJss(x,yw,yold,p); #Calculate Residual and Jacobian from new y value
    k=0
    while np.linalg.norm(R)>tol :
        k=k+1
        current_tol=np.linalg.norm(R)
        J=sp.sparse.csc_matrix(J)
        dif=-sp.sparse.linalg.spsolve(J,R) #Apply built in sparse Linear solver to find delta from J and R
        yw=yw+dif ; #Update y
        [R,J,nv_RJ]=RJss(x,yw,yold,p); #Calculate Residual and Jacobian from new y value
        if k>100:
            print('Whoops')
            whoops=whoops+1
            break
    y=yw
    return (y,whoops,vn_Newton_Rhapson,vn_RJss)

"""
Purpose: Script for running Newton-Rhapson to solve steady-state version of N3

Version 1.0

Changes from Version 0.1 to 1.0 (2/20/2022 1:30 am)
    Got initial ocde to run. 

Created on Sun Feb 27 22:39:28 2022

@author: joshu
"""

