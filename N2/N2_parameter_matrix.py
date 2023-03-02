#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_parameter_matrix_generator=1.1

import numpy as np

def parameter_matrix_generator(h,tol,t1,t2,nx,gam,beta,F,Re,n):
    parameter_combos_count=len(h)*len(tol)*len(t1)*len(t2)*len(nx)*len(gam)*len(beta)*len(F)*len(Re)*len(n) #Number of combination of parameters needed to test
    parameter_matrix = np.zeros((parameter_combos_count,10)) #Initialize matrix with new combo of paramters in each row. 9 is because there are 9 paramteres, can change if model changes
    k=0 #Begin counter for parameter matrix loop
    for h_i in h:  #Loop for timestep
        for tol_i in tol: #Loop for tolerances
            for t1_i in t1:  #Loop for inital time
                for t2_i in t2:  #Loop for final time
                    for nx_i in nx:  #Loop for mesh-size
                        for gam_i in gam:  #Loop for dimenionless ratio of diffusivities
                            for beta_i in beta:  #Loop for dimenionless ratio of potentials
                                for F_i in F:  #Loop for dimensionless forward reaction rate constants
                                    for Re_i in Re:  #Loop for dimensionless reverse reaction rate constants
                                        for n_i in n:  #Loop for hill coeffecient
                                            parameter_matrix[k,:]=[h_i,tol_i,t1_i,t2_i,nx_i,gam_i,beta_i,F_i,Re_i,n_i]
                                            k=k+1
    
    return [parameter_matrix,parameter_combos_count,vn_parameter_matrix_generator]

"""
Version 1.1

Purpose: Script for taking the desired vectors of parameters which want to be tested and converting them to usable variables for method of lines script

Created on Sat Oct 31 17:45:01 2020

@author: joshuaprince
"""