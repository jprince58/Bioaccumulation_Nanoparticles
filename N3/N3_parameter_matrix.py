#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_parameter_matrix_generator=0.1

import numpy as np

def parameter_matrix_generator(h,tol,t1,t2,nx,gam,F,K,eps,omega,ups,Kp,beta):
    parameter_combos_count=len(h)*len(tol)*len(t1)*len(t2)*len(nx)*len(gam)*len(F)*len(K)*len(eps)*len(omega)*len(ups)*len(Kp)*len(beta) #Number of combination of parameters needed to test
    parameter_matrix = np.zeros((parameter_combos_count,13)) #Initialize matrix with new combo of paramters in each row. 9 is because there are 9 paramteres, can change if model changes
    k=0 #Begin counter for parameter matrix loop
    for h_i in h:  #Loop for timestep
        for tol_i in tol: #Loop for tolerances
            for t1_i in t1:  #Loop for inital time
                for t2_i in t2:  #Loop for final time
                    for nx_i in nx:  #Loop for mesh-size
                        for gam_i in gam:  #Loop for dimenionless ratio of diffusivities
                            for F_i in F:  #Loop for dimensionless forward reaction rate constants
                                for K_i in K:  #Loop for dimesnionless NP binding constant
                                    for eps_i in eps:  #Loop for NP binding sites to NP supernatant concentration ratio
                                        for omega_i in omega: #Loop for ratio of electrokinetic nonsense
                                            for ups_i in ups: #Loop for other electrokinetic nonsense
                                                for Kp_i in Kp: #Loop for equilibrium parition coeffecient of NP between water and biofilm
                                                    for beta_i in beta: #Loop for ratio of medium mobility from diffusion to medium mobility ffrom eletrokinesis
                                                        parameter_matrix[k,:]=[h_i,tol_i,t1_i,t2_i,nx_i,gam_i,F_i,K_i,eps_i,omega_i,ups_i,Kp_i,beta_i]
                                                        k=k+1
    
    return [parameter_matrix,parameter_combos_count,vn_parameter_matrix_generator]

"""
Version 1.0

Purpose: Script for taking the desired vectors of parameters which want to be tested and converting them to usable variables for method of lines script

Changes from Version 0.1 to 1.0 (2/20/2022 1:30 am)
    Got initial ocde to run. 

Created on Sat Oct 31 17:45:01 2020

@author: joshuaprince
"""