# -*- coding: utf-8 -*-
"""
Created on Wed Dec  8 12:56:18 2021

@author: joshu
"""
from scipy.integrate import odeint
import numpy as np
from N11_exp_data_fitter import *
from N11_experimental_data_extractor import *
from N11_parameter_matrix import *
from N11_parameter_tester import *


def residual_calc(paras):  
    direct_export_path=r'C:\Users\joshu\Box\Quantum Biofilms\Raw Data\N6_results\Direct Exports'
    internal_export_path=r'C:\Users\joshu\Box\Quantum Biofilms\Raw Data\N6_results\Internal Exports' 

    
    #Grab parameter values
    try:
        h=np.array([paras['h'].value])
        tol=np.array([paras['tol'].value])
        t1=np.array([paras['t1'].value])
        t2=np.array([paras['t2'].value])
        nx=np.array([paras['nx'].value])
        omega=np.array([paras['omega'].value])
        mu=np.array([paras['mu'].value])
        nu=np.array([paras['nu'].value])
        eps=np.array([paras['eps'].value])
        rho=np.array([paras['rho'].value])
        kappa=np.array([paras['kappa'].value])
        a=np.array([paras['a'].value])
        b=np.array([paras['b'].value])
        c=np.array([paras['c'].value])
        Kp=np.array([paras['Kp'].value])
        kconv=np.array([paras['kconv'].value])
        
    except KeyError:
        h,tol,t1,t2,nx,omega,mu,nu,eps,rho,kappa,a,b,c,Kp,kconv=paras

    ci=10**(-10) #Define the inital concentration in the biofilm (Can't be zero, if one wants to be zero, set it to a very small number instead)
    
    #Grab Experimental Results
    experimental_data_file=r'C:\Users\joshu\Box\Quantum Biofilms\Processed Data\Extracted data from literature\tseng_fits_Fig2B_Cy5_incubation_bump.csv'
    [experimental_results,fit_coeff] = experimental_data_extractor(experimental_data_file)
    
    #Calculate Model Results
    [parameter_matrix,parameter_combos_count,vn_parameter_matrix_generator]=parameter_matrix_generator(h,tol,t1,t2,nx,omega,mu,nu,eps,rho,kappa,a,b,c,Kp)
    [c_set,vn_parameter_checker,vn_method_of_lines,vn_RJ] = parameter_checker(parameter_matrix,ci,fit_coeff,kconv) #output the set of concentration over time and space results for each set of parameters tested
    
    #Use Experimental Data Fitter to take model results and align with literature values
    fitting_results=exp_data_fitter(c_set,experimental_results,parameter_combos_count,internal_export_path,kconv,t2)
    
    err=fitting_results[:,2]-fitting_results[:,3] #Calculate error between model and experimental results
    return err