#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_parameter_tester=1.0

import numpy as np
from N3_method_of_lines import *
from N3_Newton_Rhapson import *


def parameter_checker(parameter_matrix,fit_coeff): #unpack paramteres and test

    #Calculate other internal paramters to model
    parameter_combos_count=np.shape(parameter_matrix) [0]
    c_set = [[0 for i in range(15)] for j in range(parameter_combos_count)] #Number in the i-range bracket is number of matrices being sotred in c_set. Can always makde more than used
    for i in np.arange(0,parameter_combos_count,1): #Begin for loop to test the different model paramters using Newton-Rhapson with Implicit Euler
        h=parameter_matrix[i,0] #Define timesteps to test
        tol=parameter_matrix[i,1] #Define tolerance to test
        t1=parameter_matrix[i,2] #Define initial time for this iteration
        t2=parameter_matrix[i,3] #Define final time for this iteration
        nx=int(parameter_matrix[i,4]) #Define mesh size for this iteration
        gam=parameter_matrix[i,5] #Define dimensionless ratio of diffusivities for this iteration
        F=parameter_matrix[i,6] #Define dimensionless forward rate constant for this iteration
        K=parameter_matrix[i,7] #Define dimensionless Eqilbrium constant for NP binding
        eps=parameter_matrix[i,8] #Define ratio of total NP binding sites to supernatant NP concentration
        omega=parameter_matrix[i,9] #Define ratio of NP contribution to electrical potential profile to their electrokinetic mobility
        ups=parameter_matrix[i,10] #Define ratio of Biofilm cotridubtion to electrical potential profile to NP electrokinetic mobility
        Kp=parameter_matrix[i,11] #Define partition coeffecient of NP into biofilm at water-biofilm interface
        beta=parameter_matrix[i,12] #Define ratio of medium mobility from diffusion to medium mobility ffrom eletrokinesis
        #Calculate internal paramters to model
        x=np.linspace(0,1,nx+1) #Define x (note, if x doesn't range from 0 to 1, should edit this)
        c_set[i][9]=x #Pass along x-vector for this parameter set for plotting
        t=np.arange(t1,t2+h,h) #Define t
        c_set[i][7]=t #Pass along t-vector for this parameter set for plotting
        ny=3*nx #Solution vector size
        nt=len(t) #Define the number of timepoints
        c_set[i][8]=nt #Pass along number of time-points used for this parameter set for plotting
        y=np.zeros((ny+3,nt)) #Initialize y
        y=y+10**(-8) #Make starting values not exactly equal to zero (divide by zero erros pop up)
        yss=np.zeros(2*nx+2)
        yss=yss+10**(-8)
        #Define initial condition for electrical potential
        for j in np.arange(0,len(x),1): #Loop over all nodes
            l=j*3+2 #Convert node # to corresponding global variables for potential
            y[l,0]=ups/6*(3*x[j]**2-x[j]**3-2)
        p=[gam,F,K,eps,omega,ups,Kp,beta] #dimensionless parameter matrix
        #Run calculation for parameters of interest
        [c,whoops,vn_method_of_lines,vn_RJ]=method_of_lines(t,x,y,h,p,tol,fit_coeff) #Find the concntration profiles in space and time using Implicit Euler 
        print('you whoopsed {} many times'.format(whoops))
        
        #Unpack the time-data
        cb=np.zeros((nx+1,nt)) #Initalize new concentration array where bound and unbound NP concentrations are "unpacked" such that they occupy two different matrices in the same 3-D array 
        cu=np.zeros((nx+1,nt))
        pot=np.zeros((nx+1,nt))
        ct=np.zeros((nx+1,nt))
        xindex=np.arange(0,nx+1)
        yss_guess=np.zeros(2*nx+2)
        for x_i in xindex:
            j=3*x_i #secondary index (position of unbound NP concentration in original concentration matrix)
            k=3*x_i+1 #Secondary index (position of bound NP concentration in original concentration matrix)
            l=3*x_i+2 #Secondary index (position of potential value in original concentration matrix)
            cu[x_i,:]=c[j,:]
            cb[x_i,:]=c[k,:]
            pot[x_i,:]=c[l,:]
            ct[x_i,:]=Kp*c[j,:]+c[k,:]
            j_guess=2*x_i
            l_guess=2*x_i+1
            yss_guess[j_guess]=c[j,nt-1]
            yss_guess[l_guess]=c[l,nt-1]
        c_set[i][0]=cb
        c_set[i][1]=cu
        c_set[i][2]=pot
        c_set[i][14]=ct

        # solve for steady-state starting with final value from time-solution
        [css,whoops2,vn_Newton_Rhapson,vn_RJss]=Newton_Rhapson(x,yss_guess,p,tol) #Find concentration profile at steady-state using Newton-Rhapson method
        print('you were bamboozled on ss calc {} many times'.format(whoops2))
        #Unpack the steady-state data and do basic data analysis
        cu_ss=np.zeros(nx+1)
        cb_ss=np.zeros(nx+1)
        pot_ss=np.zeros(nx+1)
        xindex=np.arange(0,nx+1)
        for x_i in xindex:
            j=2*x_i #Secondary index (posotion of unbound NP in original cocnentration matrix)
            l=2*x_i+1 #Secondary index (position of potential value in original concentration matrix)
            cu_ss[x_i]=css[j] 
            cb_ss[x_i]=cu_ss[x_i]*(eps-cu_ss[x_i])/K #Calculate bound NP cocnentration from unbound
            pot_ss[x_i]=css[l]
        c_set[i][10]=cu_ss
        c_set[i][11]=cb_ss
        c_set[i][12]=pot_ss
        average_uconc_ss=np.average(cu_ss)
        average_bconc_ss=np.average(cb_ss)
        average_tconc_ss=Kp*average_uconc_ss+average_bconc_ss
        c_set[i][13]=average_tconc_ss
    
        #Find NP Concentrations Overtime
        average_uconc_overtime=np.zeros(nt)
        average_bconc_overtime=np.zeros(nt)
        average_tconc_overtime=np.zeros(nt)
        tindex=np.arange(0,nt)
        for t_i in tindex:
            average_uconc_overtime[t_i]=np.average(cu[:,t_i])
            average_bconc_overtime[t_i]=np.average(cb[:,t_i])
            average_tconc_overtime[t_i]=Kp*average_uconc_overtime[t_i]+average_bconc_overtime[t_i]
        c_set[i][3]=average_uconc_overtime
        c_set[i][4]=average_bconc_overtime
        c_set[i][5]=average_tconc_overtime
        
        
        #Find NP Concentrations Overtime
        average_tconc_ss=average_tconc_overtime[nt-1]
        lognorm_tconc_overtime=np.zeros(nt)
        tindex=np.arange(0,nt)
        for t_i in tindex:
            average_uconc_overtime[t_i]=np.average(cu[:,t_i])
            average_bconc_overtime[t_i]=np.average(cb[:,t_i])
            average_tconc_overtime[t_i]=Kp*average_uconc_overtime[t_i]+average_bconc_overtime[t_i]
            lognorm_tconc_overtime[t_i]=np.log(round(average_tconc_ss,3)+0.00000001-round(average_tconc_overtime[t_i],3))
        c_set[i][6]=lognorm_tconc_overtime
        
        
        #Find Average total NP concentration Overtime (This was always a flawed metric)
        # taverage_conc_overtime=np.zeros(nt)
        # t4index=np.arange(0,nt)
        # for t4_i in t4index:
        #     taverage_conc_overtime[t4_i]=np.average(cu[:,t4_i])+np.average(cb[:,t4_i])
        # c_set[i][5]=taverage_conc_overtime
        
        #Find Change in Total NP Concentration overtime
        # tchange_in_concentration=np.zeros(nt)
        # t5index=t4index[:-1].copy()
        # for t5_i in t5index:
        #     tchange_in_concentration[t5_i]=(taverage_conc_overtime[t5_i+1]-taverage_conc_overtime[t5_i])/h
        # c_set[i][6]=tchange_in_concentration
        
            
    return [c_set,vn_parameter_tester,vn_method_of_lines,vn_RJ]


"""
Purpose: Script for running through each parameter combination and passing to method of lines vector. Then unpacks the results into usable numpy arrays for unbound and bound concentrations

Version 1.0

Changes from Version 0.1 to 1.0 (2/20/2022 1:30 am)
    Got initial ocde to run. Had to comment out some functionalities
    

Created on Sat Oct 31 17:32:29 2020

@author: joshuaprince
"""