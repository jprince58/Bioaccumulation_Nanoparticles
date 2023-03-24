#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_parameter_tester=1.1

import numpy as np
from N6_method_of_lines import *
from N6_Newton_Rhapson import *


def parameter_checker(parameter_matrix,ci): #unpack paramteres and test
    
    #Calculate other internal paramters to model
    parameter_combos_count=np.shape(parameter_matrix) [0]
    c_set = [[0 for i in range(12)] for j in range(parameter_combos_count)]
    for i in np.arange(0,parameter_combos_count,1): #Begin for loop to test the different model paramters using MOL #Check if you got he upperbound right
        h=parameter_matrix[i,0] #Define timesteps to test
        tol=parameter_matrix[i,1] #Define tolerance to test
        t1=parameter_matrix[i,2] #Define initial time for this iteration
        t2=parameter_matrix[i,3] #Define final time for this iteration
        nx=int(parameter_matrix[i,4]) #Define mesh size for this iteration
        omega=parameter_matrix[i,5] ##Define effective diffusivity  for this iteration
        mu=parameter_matrix[i,6] #Define dimensionless and porosity adjusted binding rate constant for this iteration
        nu=parameter_matrix[i,7] ##Define dimensionless and porosity adjusted maximum binding site density for this iteration
        eps=parameter_matrix[i,8] #Define dimensionless minimum interstital porosity for this iteration
        rho=parameter_matrix[i,9] #Define dimensionless minimum traditional porosity for this iteration
        kappa=parameter_matrix[i,10]  #Define dimensionless and porosity adjusted equilibrium constant  for this iteration
        a=parameter_matrix[i,11]  #Define shape paramter for binding site profile  for this iteration
        b=parameter_matrix[i,12]  #Define shape paramter for intersitital porosity profile  for this iteration
        c=parameter_matrix[i,13]  #Define shape paramter for traditional proosity profile  for this iteration
        Kp=parameter_matrix[i,14] #define partition coeffecient for this iteration
        
        #Calculate internal paramters to model
        x=np.linspace(0,1,nx+1) #Define x (note, if x doesn't range from 0 to 1, should edit this)
        c_set[i][0]=x #Pass along x-vector for this parameter set for plotting
        t=np.arange(t1,t2+h,h) #Define t
        c_set[i][1]=t #Pass along t-vector for this parameter set for plotting
        ny=2*nx #Solution vector size
        nt=len(t) #Define the number of timepoints
        c_set[i][2]=nt #Pass along number of time-points used for this parameter set for plotting
        y=np.zeros((ny+2,nt)) #Initialize y
        y=y+10**(-8) #Make starting values not exactly equal to zero (divide by zero erros pop up)
        y[ny,:]=(1+rho)
        p=[omega, mu, nu, eps, rho, kappa, a, b, c, Kp] #dimensionless parameter matrix

        #Run calculation for parameters of interest
        [c,whoops,vn_method_of_lines,vn_RJ]=method_of_lines(t,x,y,h,p,tol) #Find the concntration profiles in space and time using Method of Lines (MOL)
        print('you whoopsed {} many times'.format(whoops))
        
        #Unpack the data
        cm=np.zeros((nx+1,nt)) #Initalize new concentration array where bound and unbound NP concentrations are "unpacked" such that they occupy two different matrices in the same 3-D array 
        ca=np.zeros((nx+1,nt))
        ct=np.zeros((nx+1,nt))
        xindex=np.arange(0,nx+1)
        yss_guess=np.zeros(2*nx+2)
        for x_i in xindex:
            j=2*x_i #secondary index (position of unbound NP concentration in original concentration matrix)
            k=2*x_i+1 #Secondary index (position of bound NP concentration in original concentration matrix)
            cm[x_i,:]=c[j,:]
            ca[x_i,:]=c[k,:]
            ct[x_i,:]=Kp*c[j,:]+c[k,:]
            yss_guess[j]=c[j,nt-1]
            yss_guess[k]=c[k,nt-1]
        c_set[i][3]=cm
        c_set[i][4]=ca
        c_set[i][5]=ct
        
        # %% Calcualte Steady-state profile
        [css,whoops2,vn_Newton_Rhapson,vn_RJss]=Newton_Rhapson(x,yss_guess,p,tol) #Find concentration profile at steady-state using Newton-Rhapson method
        print('you were bamboozled on ss calc {} many times'.format(whoops2))
        
        #Unpack SS profile and pass along to c_set
        cm_ss=np.zeros((nx+1,nt)) #Initalize new concentration array where bound and unbound NP concentrations are "unpacked" such that they occupy two different matrices in the same 3-D array 
        ca_ss=np.zeros((nx+1,nt))
        ct_ss=np.zeros((nx+1,nt))
        
        for x_i in xindex:
            j=2*x_i #secondary index (position of unbound NP concentration in original concentration matrix)
            k=2*x_i+1 #Secondary index (position of bound NP concentration in original concentration matrix)
            for x_i in xindex:
                j=2*x_i #secondary index (position of unbound NP concentration in original concentration matrix)
                k=2*x_i+1 #Secondary index (position of bound NP concentration in original concentration matrix)
                cm_ss[x_i]=css[j]
                ca_ss[x_i]=css[k]
                ct_ss[x_i]=Kp*css[j]+css[k]
        c_set[i][9]=cm_ss #Pass along steady-state mobile nanoparticle concentration
        c_set[i][10]=ca_ss #Pass along steady-state attached nanoparticle concentration
        c_set[i][11]=ct_ss #Pass along steady-state total nanoparticle concentration
        
    
        #Find average NP Concentration Overtime
        average_mobile_conc_overtime=np.zeros(nt)
        average_attached_conc_overtime=np.zeros(nt)
        average_total_conc_overtime=np.zeros(nt)
        tindex=np.arange(0,nt)
        for t_i in tindex:
            average_mobile_conc_overtime[t_i]=np.average(cm[:,t_i])
            average_attached_conc_overtime[t_i]=np.average(ca[:,t_i])
            average_total_conc_overtime[t_i]=Kp*average_mobile_conc_overtime[t_i]+average_attached_conc_overtime[t_i]
        c_set[i][6]=average_mobile_conc_overtime
        c_set[i][7]=average_attached_conc_overtime
        c_set[i][8]=average_total_conc_overtime
        
            
    return [c_set,vn_parameter_tester,vn_method_of_lines,vn_RJ]


"""
Purpose: Script for running through each parameter combination and passing to method of lines vector. Then unpacks the results into usable numpy arrays for unbound and bound concentrations

Version 1.1

Created on Sat Oct 31 17:32:29 2020

@author: joshuaprince
"""