#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_linear_fitting=1.2

import numpy as np
import matplotlib.pyplot as plt
import os

def linear_fit(c_set,parameter_combos_count,parameter_matrix,internal_export_path):
    
    # %% Calculate linear best fit, slope and intercept for each parameter set
    lin_fit=np.zeros((parameter_combos_count,3))#initialize matrix to store linear best fit parameters
    perc_acc_matrix= [[0 for i in range(1)] for j in range(parameter_combos_count)]
    for pc_i in np.arange(0,parameter_combos_count,1): #Begin for loop over different model paramter sets 
        F= parameter_matrix[pc_i,7] #Grab the Dimensionless forward rate constant for parameter set
        Re= parameter_matrix[pc_i,8] #Grab the Dimensionless reverse rate constant for parameter set
        Kp=10
        Eq=F/(F+Re)+Kp #Calculate Equilibrium total concentration value (assumes thetaequilibrates to one, which it is defined to)
        tavg_conc=c_set[pc_i][4] #Grab current change in concentration data to plot (total NP)
        perc_acc_model=tavg_conc/Eq #convert concentration vectors to percent accumulated vectors
        norm_tavg_conc=Eq-tavg_conc #Normalize average concentration by equilibrium concentration
        t=c_set[pc_i][6] #Grab time-vector for this parameter set
        # %% Find the 99% accumulation time and cutoff perc_acc_model
        cutoff= 0.95 #cutoff percentage to "reach equilibrium" (used to determine where to start fit)
        j=0 #reset counter for time-search loop
        for t_i in t:
            if perc_acc_model[j]>cutoff:
                mod_cutoff=j #Index neccesary for model to reach evaluated percent accumulated
                break
            j=j+1 #Update counter in time search loop
        t_cutoff=t[:mod_cutoff]
        norm_tavg_conc_cutoff=norm_tavg_conc[:mod_cutoff]
        lognorm_tavg_conc_cutoff=np.zeros(mod_cutoff) #Initialize log of normalized concentration vector
        j=0
        for t_i in t_cutoff:
            lognorm_tavg_conc_cutoff[j]=np.log(norm_tavg_conc_cutoff[j])
            j=j+1
        lognorm_tavg_conc=np.zeros(len(t)) #Initialize log of normalized concentration vector
        j=0
        for t_i in t:
            lognorm_tavg_conc[j]=np.log(norm_tavg_conc[j])
            j=j+1
        #%% Run linear fit
        [m,b]=np.polyfit(t_cutoff,lognorm_tavg_conc_cutoff,1) #Find linear fit for plot
        lin_fit[pc_i,0]=m
        lin_fit[pc_i,1]=b #Store fit into matrix
        
    # %% Use paramteres to generate best-fit data
        t=c_set[pc_i][6] #Grab time-vector for this parameter set
        fit_conc=np.zeros(len(t))#Initialize concentration vector
        fit_lognorm_avg=np.zeros(len(t)) #Initialize longnormal concentration vactor
        count=0 #Begin counter
        for t_i in t:
            fit_conc[count]=Eq-Eq*np.exp(m*t_i)
            fit_lognorm_avg[count]=m*t_i+b
            count=count+1
        r_sq=np.corrcoef(lognorm_tavg_conc_cutoff,t_cutoff)
        lin_fit[pc_i,2]=-r_sq[0,1]
            
        
    # %% Generate Plot and save    
        plt.figure(1000+pc_i)
        plt.plot(t,tavg_conc,label='Model Results')
        plt.plot(t,fit_conc,label='First-Order Approximation')
        upper_1 = np.amax(tavg_conc)*1.1 #Upper Bound on Average total Concentration
        upper_2 = np.amax(fit_conc)*1.1 #Upper bound on fit average total concentration overtime
        if upper_1>upper_2:
            uplimit=upper_1
        else: uplimit=upper_2
        plt.xlim(left=parameter_matrix[pc_i,2],right=parameter_matrix[pc_i,3]) 
        plt.ylim(bottom=0,top=uplimit)
        plt.xlabel('Time',fontsize=14)
        plt.ylabel('Average Dimensionless Concentration',fontsize=14)
        plt.title('Model vs First-order Mass Transfer',fontsize=16)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(loc=(0.5,0.3))
        plt.text(0.4,0.5,f'm={m}, b={b}')
        linear_filename_partial=f'Linearplot{pc_i}.png'
        linear_filename_full=os.path.join(internal_export_path,linear_filename_partial)
        plt.savefig(linear_filename_full)
        plt.close()
        
        plt.figure(1001+pc_i)
        plt.plot(t_cutoff,lognorm_tavg_conc[:mod_cutoff],label='Model Results')
        plt.plot(t_cutoff,fit_lognorm_avg[:mod_cutoff],label='First-Order Approximation')
        upper_1 = np.amax(lognorm_tavg_conc[:mod_cutoff])*1.1 #Upper bound on fit average total concentration overtime
        lower_1 = np.amin(lognorm_tavg_conc[:mod_cutoff])*1.1  #Lower Bound on fit average total concentration overtime
        plt.xlim(left=0,right=np.amax(t_cutoff))  
        plt.ylim(bottom=lower_1,top=upper_1)
        plt.xlabel('Time',fontsize=14)
        plt.ylabel('log(Normalized Average Concentration',fontsize=14)
        plt.title('Model vs First-order Mass Transfer',fontsize=16)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(loc=(0.5,0.75))
        plt.text(0.3,-3.7,f'$R^2$={-r_sq[0,1]}')
        log_filename_partial=f'Logplot{pc_i}.png'
        log_filename_full=os.path.join(internal_export_path,log_filename_partial)
        plt.savefig(log_filename_full)
        plt.close()
        
    # %% Determine fit of first-order approximation
        percents=[0.1,0.5,0.9,0.95] #The percent bioaccumulated that will be checked 
        perc_acc_approx=fit_conc/Eq 
        perc_acc_table=np.zeros((len(percents),4))
        k=0 #counter for percents for loop
        for p_i in percents:
            j=0 #reset counter for time-search loop
            for t_i in t:
                if perc_acc_model[j]>p_i:
                    mod_time=t_i #Time neccesary for model to reach evaluated percent accumulated
                    break
                j=j+1 #Update counter in time search loop
            j=0 #reset counter for time-search loop
            for t_i in t:
                if perc_acc_approx[j]>p_i:
                    approx_time=t_i #time neccesary for first-order approximation to reach evaluated percent accumulated
                    break
                j=j+1 #Update counter in time search loop
            perc_acc_table[k,0]=p_i
            perc_acc_table[k,1]=mod_time
            perc_acc_table[k,2]=approx_time
            perc_acc_table[k,3]=(approx_time-mod_time)/mod_time*100
            k=k+1 #Update counter in percent loop  
        perc_acc_matrix[pc_i][0]=perc_acc_table
    return [perc_acc_matrix,vn_linear_fitting]



"""
Created on Wed Nov 11 10:19:35 2020

@author: joshuaprince (prince.j@northeastern.edu)

Purpose: Script to fit concentration data to first-order mass-transfer relation, plot this compared to the model results, and create a percent accumulated table to compare the two

Version 1.3

changes From Version 1.2 (2/13/2021 7:40 pm):
    -added R^2 value comparing 1st order fit and data on log(normalized concentration) plot.

changes From Version 1.1 to 1.2 (11/13/2020 1:12 am):
    -changing method of fitting. going from fitting the dC vs C plot to fitting log(normalized concentration) plot. 
"""

