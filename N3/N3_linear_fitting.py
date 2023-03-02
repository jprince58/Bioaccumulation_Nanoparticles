#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_linear_fitting=1.0

import numpy as np
import matplotlib.pyplot as plt
import os
import numpy as np
import pandas as pd
import sklearn
from sklearn import linear_model

def linear_fit(c_set,parameter_combos_count,parameter_matrix,internal_export_path,gam):
    
    
    # %%Initialize Sherwood Relations plots
    
    
    # %% Calculate linear best fit, slope and intercept for each parameter set
    lin_fit=np.zeros((parameter_combos_count,11))#initialize matrix to store linear best fit parameters
    lin_fit[:,3:]=parameter_matrix[:,5:]
    perc_acc_matrix= [[0 for i in range(6)] for j in range(parameter_combos_count)]
    for pc_i in np.arange(0,parameter_combos_count,1): #Begin for loop over different model paramter sets 
        #Eq=c_set[pc_i][13] #Grab Equilibrium total concentration value (assumes thetaequilibrates to one, which it is defined to)
        tavg_conc=c_set[pc_i][5] #Grab current total average concentration data to plot 
        Eq=np.amax(tavg_conc) #Calculate steady-state value
        perc_acc_model=tavg_conc/Eq #convert concentration vectors to percent accumulated vectors
        norm_tavg_conc=Eq-tavg_conc #Normalize average concentration by equilibrium concentration
        t=c_set[pc_i][7] #Grab time-vector for this parameter set
        # %% Find the 99% accumulation time and cutoff perc_acc_model
        cutoff= 0.99 #cutoff percentage to "reach equilibrium" (used to determine where to start fit)
        j=0 #reset counter for time-search loop
        for t_i in t:
            if perc_acc_model[j]>cutoff:
                mod_cutoff=j #Index neccesary for model to reach evaluated percent accumulated
                break
            else: mod_cutoff=len(t)
            j=j+1 #Update counter in time search loop
        t_cutoff=t[:mod_cutoff]
        norm_tavg_conc_cutoff=norm_tavg_conc[:mod_cutoff]
        lognorm_tavg_conc_cutoff=np.zeros(mod_cutoff) #Initialize log of normalized concentration vector
        j=0
        for t_i in t_cutoff:
            lognorm_tavg_conc_cutoff[j]=np.log10(norm_tavg_conc_cutoff[j])
            j=j+1
        lognorm_tavg_conc=np.zeros(len(t)) #Initialize log of normalized concentration vector
        j=0
        for t_i in t:
            lognorm_tavg_conc[j]=np.log10(norm_tavg_conc[j])
            j=j+1
        #%% Run linear fit
        [m,b]=np.polyfit(t_cutoff,lognorm_tavg_conc_cutoff,1) #Find linear fit for plot
        lin_fit[pc_i,0]=m
        lin_fit[pc_i,1]=b #Store fit into matrix
        
    # %% Use paramteres to generate best-fit data
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
        plt.figure(pc_i)
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
        linear_filename_partial=f'Linearplot{pc_i}.png'
        linear_filename_full=os.path.join(internal_export_path,linear_filename_partial)
        plt.savefig(linear_filename_full, bbox_inches='tight')
        plt.close()
        
        plt.figure(2*pc_i+1)
        plt.plot(t,lognorm_tavg_conc,label='Model Results')
        plt.plot(t,fit_lognorm_avg,label='First-Order Approximation')
        upper_1 = np.amax(lognorm_tavg_conc_cutoff)*1.1 #Upper bound on fit average total concentration overtime
        lower_1 = np.amin(lognorm_tavg_conc_cutoff)*0.9  #Lower Bound on fit average total concentration overtime
        plt.xlim(left=parameter_matrix[pc_i,2],right=np.amax(t_cutoff))  
        plt.ylim(bottom=lower_1,top=upper_1)
        plt.xlabel('Time',fontsize=14)
        plt.ylabel('log(Normalized Average Concentration)',fontsize=14)
        plt.title('Model vs First-order Mass Transfer',fontsize=16)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(loc=(0.5,0.75))
        plt.figtext(0.2,0.2,f'$R^2$={-r_sq[0,1]}')
        log_filename_partial=f'Logplot{pc_i}.png'
        log_filename_full=os.path.join(internal_export_path,log_filename_partial)
        plt.savefig(log_filename_full, bbox_inches='tight')
        plt.close()
        
    # %% Determine fit of first-order approximation
        percents=[0.5,0.8,0.9,0.95,0.99] #The percent bioaccumulated that will be checked 
        perc_acc_approx=fit_conc/Eq 
        perc_acc_table=np.zeros((len(percents),5))
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
    
    
    
    # plt.figure(1000)
    # #gam=parameter_matrix[:,5] #Dimensionless numbers checked (getting gamma directly right now)
    # Sh1=np.zeros(int(parameter_combos_count/4))
    # Sh2=np.zeros(int(parameter_combos_count/4))
    # Sh3=np.zeros(int(parameter_combos_count/4))
    # Sh4=np.zeros(int(parameter_combos_count/4))
    # for i in np.arange(0,parameter_combos_count):
    #     # if i%4==0:
    #     #     Sh1[int(i/4)]=-lin_fit[i,0] #Sherwood #'s from calculations
    #     # elif i%4==1:
    #     #     Sh2[int((i-1)/4)]=-lin_fit[i,0] #Sherwood #'s from calculations
    #     # elif i%4==2:
    #     #     Sh3[int((i-2)/4)]=-lin_fit[i,0] #Sherwood #'s from calculations
    #     # else:
    #     #     Sh4[int((i-3)/4)]=-lin_fit[i,0] #Sherwood #'s from calculations

    #     if i<8:
    #         Sh1[int(i)]=-lin_fit[i,0] #Sherwood #'s from calculations
    #     elif i<16:
    #         Sh2[int(i-8)]=-lin_fit[i,0] #Sherwood #'s from calculations
    #     elif i<24:
    #         Sh3[int(i-16)]=-lin_fit[i,0] #Sherwood #'s from calculations
    #     else:
    #         Sh4[int(i-24)]=-lin_fit[i,0] #Sherwood #'s from calculations
    # [m1,b1]=np.polyfit(gam,Sh1,1) #Find linear fit for plot
    # [m2,b2]=np.polyfit(gam,Sh2,1) #Find linear fit for plot
    # [m3,b3]=np.polyfit(gam,Sh3,1) #Find linear fit for plot
    # [m4,b4]=np.polyfit(gam,Sh4,1) #Find linear fit for plot
    # r_sq1=np.corrcoef(gam,Sh1)
    # r_sq2=np.corrcoef(gam,Sh2)
    # r_sq3=np.corrcoef(gam,Sh3)
    # r_sq4=np.corrcoef(gam,Sh4)
    # perc_acc_matrix[0][4]=[[m1,b1,r_sq1],[m2,b2,r_sq2],[m3,b3,r_sq3],[m4,b4,r_sq4]]
    # plt.scatter(gam,Sh1, label='epsilon=10 , omega=0.1') #Plot Varied Dimensionless number vs Sh. 
    # plt.scatter(gam,Sh2, label='epsilon=10 , omega=10') #Plot Varied Dimensionless number vs Sh.
    # plt.scatter(gam,Sh3, label='epsilon=100 , omega=0.1') #Plot Varied Dimensionless number vs Sh.
    # plt.scatter(gam,Sh4, label='epsilon=100 , omega=10') #Plot Varied Dimensionless number vs Sh.
    # plt.plot(gam,m1*gam+b1) #Plot best fit line
    # plt.plot(gam,m2*gam+b2) #Plot best fit line
    # plt.plot(gam,m3*gam+b3) #Plot best fit line
    # plt.plot(gam,m4*gam+b4) #Plot best fit line
    # #plt.loglog(basex=10, basey=10)
    # upper_2 = np.amax(gam)*1.1 #Upper bound on fit average total concentration overtime
    # lower_2 = np.amin(gam)*0.9  #Lower Bound on fit average total concentration overtime
    # upper_3 = np.amax(np.concatenate((Sh1,Sh2,Sh3,Sh4)))*1.1 #Upper bound on fit average total concentration overtime
    # lower_3 = np.amin(np.concatenate((Sh1,Sh2,Sh3,Sh4)))*0.9  #Lower Bound on fit average total concentration overtime
    # plt.xlim(left=lower_2,right=upper_2)  
    # plt.ylim(bottom=lower_3,top=upper_3)
    # plt.xlabel('Dimenionless Ratio of Electrophoresis to Diffusion',fontsize=14)
    # plt.ylabel('Modelled Sherwood #',fontsize=14)
    # #plt.title('Sherwood # vs Diffusion Constant',fontsize=16)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.legend()
    # sherwood_filename_partial=f'Sherwdoodplot{pc_i}.png'
    # sherwood_filename_full=os.path.join(internal_export_path,sherwood_filename_partial)
    # plt.savefig(sherwood_filename_full, bbox_inches='tight')
    # plt.close()
    
    # plt.figure(1001)
    # loggam=np.log10(gam) #Dimensionless numbers checked (getting gamma directly right now)
    # logSh1=np.zeros(int(parameter_combos_count/4))
    # logSh2=np.zeros(int(parameter_combos_count/4))
    # logSh3=np.zeros(int(parameter_combos_count/4))
    # logSh4=np.zeros(int(parameter_combos_count/4))
    # for i in np.arange(0,parameter_combos_count):
    #     # if i%4==0:
    #     #     logSh1[int((i-0)/4)]=np.log10(-lin_fit[i,0]) #Sherwood #'s from calculations
    #     # elif i%4==1:
    #     #     logSh2[int((i-1)/4)]=np.log10(-lin_fit[i,0]) #Sherwood #'s from calculations
    #     # elif i%4==2:
    #     #     logSh3[int((i-2)/4)]=np.log(-lin_fit[i,0]) #Sherwood #'s from calculations
    #     # else:
    #     #     logSh4[int((i-3)/4)]=np.log10(-lin_fit[i,0]) #Sherwood #'s from calculations
   
    #     if i<8:
    #         logSh1[int(i)]=np.log10(-lin_fit[i,0]) #Sherwood #'s from calculations
    #     elif i<16:
    #         logSh2[int(i-8)]=np.log10(-lin_fit[i,0]) #Sherwood #'s from calculations
    #     elif i<24:
    #         logSh3[int(i-16)]=np.log(-lin_fit[i,0]) #Sherwood #'s from calculations
    #     else:
    #         logSh4[int(i-24)]=np.log10(-lin_fit[i,0]) #Sherwood #'s from calculations
    # [m1,b1]=np.polyfit(loggam,logSh1,1) #Find linear fit for plot
    # [m2,b2]=np.polyfit(loggam,logSh2,1) #Find linear fit for plot
    # [m3,b3]=np.polyfit(loggam,logSh3,1) #Find linear fit for plot
    # [m4,b4]=np.polyfit(loggam,logSh4,1) #Find linear fit for plot
    # r_sq1=np.corrcoef(loggam,logSh1)
    # r_sq2=np.corrcoef(loggam,logSh2)
    # r_sq3=np.corrcoef(loggam,logSh3)
    # r_sq4=np.corrcoef(loggam,logSh4)
    # perc_acc_matrix[0][5]=[[m1,b1,r_sq1],[m2,b2,r_sq2],[m3,b3,r_sq3],[m4,b4,r_sq4]]
    # plt.scatter(loggam,logSh1,label='epsilon=10 , omega=0.1') #Plot Varied Dimensionless number vs Sh. 
    # plt.scatter(loggam,logSh2, label='epsilon=10 , omega=10') #Plot Varied Dimensionless number vs Sh.
    # plt.scatter(loggam,logSh3, label='epsilon=100 , omega=0.1') #Plot Varied Dimensionless number vs Sh.
    # plt.scatter(loggam,logSh4, label='epsilon=100 , omega=10') #Plot Varied Dimensionless number vs Sh.
    # plt.plot(loggam,m1*loggam+b1) #Plot best fit line
    # plt.plot(loggam,m2*loggam+b2) #Plot best fit line
    # plt.plot(loggam,m3*loggam+b3) #Plot best fit line
    # plt.plot(loggam,m4*loggam+b4) #Plot best fit line
    # upper_2 = np.amax(loggam)*1.1 #Upper bound on fit average total concentration overtime
    # lower_2 = np.amin(loggam)*1.1  #Lower Bound on fit average total concentration overtime
    # upper_3 = np.amax(np.concatenate((logSh1,logSh2,logSh3,logSh4)))*1.1 #Upper bound on fit average total concentration overtime
    # lower_3 = np.amin(np.concatenate((logSh1,logSh2,logSh3,logSh4)))*0.9  #Lower Bound on fit average total concentration overtime
    # plt.xlim(left=lower_2,right=upper_2)  
    # plt.ylim(bottom=lower_3,top=upper_3)
    # plt.xlabel('Log of Dimenionless Ratio of Electrophoresis to Diffusion',fontsize=14)
    # plt.ylabel('Log Modelled Sherwood #',fontsize=14)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # plt.legend()
    # logsherwood_filename_partial=f'LogSherwdoodplot{pc_i}.png'
    # logsherwood_filename_full=os.path.join(internal_export_path,logsherwood_filename_partial)
    # plt.savefig(logsherwood_filename_full, bbox_inches='tight')
    # plt.close()
    
    """Multiple linear regression approach not working. going to just check indivdual variables to make it more sensible"""
    #%% run Multiple linear regression on Sherwood vs log10 of dimensionless numbers 
    # X=np.log10(lin_fit[:,3:])
    # y=np.log10(-lin_fit[:,0])
    # model_ols = linear_model.Lasso(normalize=True, alpha=0.55)
    # model_ols.fit(X,y) 
    # R2=model_ols.score(X,y)
    # coef = model_ols.coef_
    # perc_acc_matrix[0][1]=coef
    # intercept = model_ols.intercept_
    # perc_acc_matrix[0][2]=intercept
    # perc_acc_matrix[0][3]=R2
    
    # %% Calculate predicted sherwood number based on inputs
    # logSh_Pred=np.zeros(parameter_combos_count)
    # Sh_Pred=np.zeros(parameter_combos_count)
    # Sh_Act=-lin_fit[:,0]
    # for pc_i in np.arange(0,parameter_combos_count,1):
    #     logSh_Pred[pc_i]=intercept
    #     for i in np.arange(0,len(coef)):
    #         logSh_Pred[pc_i]=coef[i]*np.log(lin_fit[pc_i,i+3])+logSh_Pred[pc_i]
    # Sh_Pred=10**logSh_Pred
    
    # plt.figure(1000)
    # plt.scatter(Sh_Pred,Sh_Act)
    # upper_2 = np.amax(Sh_Pred)*1.1 #Upper bound on fit average total concentration overtime
    # lower_2 = np.amin(Sh_Pred)*0.9  #Lower Bound on fit average total concentration overtime
    # upper_3 = np.amax(Sh_Act)*1.1 #Upper bound on fit average total concentration overtime
    # lower_3 = np.amin(Sh_Act)*0.9  #Lower Bound on fit average total concentration overtime
    # plt.xlim(left=lower_2,right=upper_2)  
    # plt.ylim(bottom=lower_3,top=upper_3)
    # plt.xlabel('Linear Regression Sherwood #',fontsize=14)
    # plt.ylabel('Modelled Sherwood #',fontsize=14)
    # plt.title('Linear Regression Sherwood vs Modelled',fontsize=16)
    # plt.xticks(fontsize=12)
    # plt.yticks(fontsize=12)
    # sherwood_filename_partial=f'Sherwdoodplot{pc_i}.png'
    # sherwood_filename_full=os.path.join(internal_export_path,sherwood_filename_partial)
    # plt.savefig(sherwood_filename_full)
    # plt.close()
        
    
    return [perc_acc_matrix,vn_linear_fitting]



"""
Created on Wed Nov 11 10:19:35 2020

@author: joshuaprince (prince.j@northeastern.edu)

Purpose: Script to fit concentration data to first-order mass-transfer relation, plot this compared to the model results, and create a percent accumulated table to compare the two

Version 0.1


"""

