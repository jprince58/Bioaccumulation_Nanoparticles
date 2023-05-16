#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_parameter_matrix_generator=1.1

import numpy as np

def parameter_matrix_generator(h,tol,t1,t2,nx,omega,mu,nu,eps,rho,kappa,a,b,c,Kp):
    parameter_combos_count=len(h)*len(tol)*len(t1)*len(t2)*len(nx)*len(omega)*len(mu)*len(nu)*len(eps)*len(rho)*len(kappa)*len(a)*len(b)*len(c)*len(Kp) #Number of combination of parameters needed to test
    parameter_matrix = np.zeros((parameter_combos_count,15)) #Initialize matrix with new combo of paramters in each row. 9 is because there are 9 paramteres, can change if model changes
    k=0 #Begin counter for parameter matrix loop
    for h_i in h:  #Loop for timestep
        for tol_i in tol: #Loop for tolerances
            for t1_i in t1:  #Loop for inital time
                for t2_i in t2:  #Loop for final time
                    for nx_i in nx:  #Loop for mesh-size
                        for omega_i in omega:  
                            for mu_i in mu:  
                                for nu_i in nu:  
                                    for eps_i in eps:  
                                        for rho_i in rho:  
                                            for kappa_i in kappa:
                                                for a_i in a:  
                                                    for b_i in b:
                                                        for c_i in c:
                                                            for Kp_i in Kp:
                                                                parameter_matrix[k,:]=[h_i,tol_i,t1_i,t2_i,nx_i,omega_i,mu_i,nu_i,eps_i,rho_i,kappa_i,a_i,b_i,c_i,Kp_i]
                                                                k=k+1
    # parameter_combos_count=len(h)*len(tol)*len(t1)*len(t2)*len(nx)*len(omega)*len(mu)*len(nu)*len(eps)*len(rho)*len(kappa)*len(a)*len(c)*len(Kp) #Number of combination of parameters needed to test
    # parameter_matrix = np.zeros((parameter_combos_count,15)) #Initialize matrix with new combo of paramters in each row. 9 is because there are 9 paramteres, can change if model changes
    # k=0 #Begin counter for parameter matrix loop
    # b_i=-1 #Coutner for b-value to pull
    # for h_i in h:  #Loop for timestep
    #     for tol_i in tol: #Loop for tolerances
    #         for t1_i in t1:  #Loop for inital time
    #             for t2_i in t2:  #Loop for final time
    #                 for nx_i in nx:  #Loop for mesh-size
    #                     for omega_i in omega:  
    #                         for mu_i in mu:  
    #                             for nu_i in nu:  
    #                                 for eps_i in eps:
    #                                     b_i=b_i+1
    #                                     for rho_i in rho:  
    #                                         for kappa_i in kappa:
    #                                             for a_i in a:  
    #                                                 for c_i in c:
    #                                                     for Kp_i in Kp:
    #                                                         parameter_matrix[k,:]=[h_i,tol_i,t1_i,t2_i,nx_i,omega_i,mu_i,nu_i,eps_i,rho_i,kappa_i,a_i,b[b_i],c_i,Kp_i]
    #                                                         k=k+1
    return [parameter_matrix,parameter_combos_count,vn_parameter_matrix_generator]

"""
Version 1.1

Purpose: Script for taking the desired vectors of parameters which want to be tested and converting them to usable variables for method of lines script

Created on Sat Oct 31 17:45:01 2020

@author: joshuaprince
"""