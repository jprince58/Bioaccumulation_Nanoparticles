#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_parameter_tester=1.1

import numpy as np
from N2_method_of_lines import *


def parameter_checker(parameter_matrix,ci): #unpack paramteres and test
    
    #Calculate other internal paramters to model
    parameter_combos_count=np.shape(parameter_matrix) [0]
    c_set = [[0 for i in range(9)] for j in range(parameter_combos_count)]
    for i in np.arange(0,parameter_combos_count,1): #Begin for loop to test the different model paramters using MOL #Check if you got he upperbound right
        h=parameter_matrix[i,0] #Define timesteps to test
        tol=parameter_matrix[i,1] #Define tolerance to test
        t1=parameter_matrix[i,2] #Define initial time for this iteration
        t2=parameter_matrix[i,3] #Define final time for this iteration
        nx=int(parameter_matrix[i,4]) #Define mesh size for this iteration
        gam=parameter_matrix[i,5] #Define dimensionless ratio of diffusivities for this iteration
        beta=parameter_matrix[i,6] #Define dimensionless ratio of potentials for this iteration
        F=parameter_matrix[i,7] #Define dimensionless forward rate constant for this iteration
        Re=parameter_matrix[i,8] #Define dimensionless reverse rate constant for this iteration
        n=parameter_matrix[i,9] #Define hill coeffecient for this rest iteration
        #Calculate internal paramters to model
        x=np.linspace(0,1,nx+1) #Define x (note, if x doesn't range from 0 to 1, should edit this)
        c_set[i][8]=x #Pass along x-vector for this parameter set for plotting
        t=np.arange(t1,t2+h,h) #Define t
        c_set[i][6]=t #Pass along t-vector for this parameter set for plotting
        ny=2*nx #Solution vector size
        nt=len(t) #Define the number of timepoints
        c_set[i][7]=nt #Pass along number of time-points used for this parameter set for plotting
        y=np.zeros((ny+2,nt)) #Initialize y
        y=y+10**(-8) #Make starting values not exactly equal to zero (divide by zero erros pop up)
        p=[gam,beta,F,Re,n] #dimensionless parameter matrix
        Kp=10 #PArtition coeffecient
        #Run calculation for parameters of interest
        [c,whoops,vn_method_of_lines,vn_RJ]=method_of_lines(t,x,y,h,p,tol) #Find the concntration profiles in space and time using Method of Lines (MOL)
        print('you whoopsed {} many times'.format(whoops))
        
        #Unpack the data
        cb=np.zeros((nx+1,nt)) #Initalize new concentration array where bound and unbound NP concentrations are "unpacked" such that they occupy two different matrices in the same 3-D array 
        cu=np.zeros((nx+1,nt))
        xindex=np.arange(0,nx+1)
        for x_i in xindex:
            j=2*x_i #secondary index (position of unbound NP concentration in original concentration matrix)
            k=2*x_i+1 #Secondary index (position of bound NP concentration in original concentration matrix)
            cu[x_i,:]=c[j,:]
            cb[x_i,:]=c[k,:]
        c_set[i][0]=cb
        c_set[i][1]=cu
    
        #Find Average Unbound Concentration Overtime
        average_conc_overtime=np.zeros(nt)
        t2index=np.arange(0,nt)
        for t2_i in t2index:
            average_conc_overtime[t2_i]=np.average(cu[:,t2_i])
        c_set[i][2]=average_conc_overtime
        
        #Find Change in Concentration overtime
        change_in_concentration=np.zeros(nt)
        t3index=t2index[:-1].copy()
        for t3_i in t3index:
            change_in_concentration[t3_i]=(average_conc_overtime[t3_i+1]-average_conc_overtime[t3_i])/h
        c_set[i][3]=change_in_concentration
        
        #Find Average total NP concentration Overtime
        taverage_conc_overtime=np.zeros(nt)
        t4index=np.arange(0,nt)
        for t4_i in t4index:
            taverage_conc_overtime[t4_i]=np.average(Kp*cu[:,t4_i])+np.average(cb[:,t4_i])
        c_set[i][4]=taverage_conc_overtime
        
        #Find Change in Total NP Concentration overtime
        tchange_in_concentration=np.zeros(nt)
        t5index=t4index[:-1].copy()
        for t5_i in t5index:
            tchange_in_concentration[t5_i]=(taverage_conc_overtime[t5_i+1]-taverage_conc_overtime[t5_i])/h
        c_set[i][5]=tchange_in_concentration
        
            
    return [c_set,vn_parameter_tester,vn_method_of_lines,vn_RJ]


"""
Purpose: Script for running through each parameter combination and passing to method of lines vector. Then unpacks the results into usable numpy arrays for unbound and bound concentrations

Version 1.1

Created on Sat Oct 31 17:32:29 2020

@author: joshuaprince
"""