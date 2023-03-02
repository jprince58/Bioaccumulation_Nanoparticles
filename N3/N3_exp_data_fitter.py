#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_experimental_data_extractor=1.0

import os
import numpy as np
import matplotlib.pyplot as plt
import docx
from docx.shared import Pt
from datetime import datetime
import matplotlib.animation as anim
from matplotlib.animation import FuncAnimation
import pandas as pd


def exp_data_fitter(c_set,exp_data,parameter_combos_count,internal_export_path):
    
    for pc_i in np.arange(0,parameter_combos_count,1): #Begin for loop to plot the different model paramters tested
    
        # %% Grab relevant data from model     
        ct=c_set[pc_i][14] #Grab total dimensionless concentration plot
        x=c_set[pc_i][9] #Grab the position-vector for this parameter set for plotting
        t=c_set[pc_i][7] #Grab time-vector for this parameter set for plotting
        
        # %% Convert dimensionless mdoel results inot dimenional model results
        to= 0.5 #guess at dimensionless time [min]
        kconv=0.75 #guess at absorbance units-NP concentration conversion factor
        t_d=t*to #convert dimensionless time into dimensional time [min]
        ct_d=kconv*ct #convert dimensionless total concentration to absorbance units
        
        # %% Get model results into array similar to experimental results
        mod_data=np.zeros((len(t)*len(x),3)) #Initialize array
        for xi in np.arange(0,len(x)):
            position=x[xi]
            for ti in np.arange(0,len(t_d)):
                time=t_d[ti]
                index=xi+len(x)*ti #sent index in the mod_data matrix
                mod_data[index,0]=time
                mod_data[index,1]=position
                mod_data[index,2]=ct_d[xi,ti]
            
        
        # %% Replace zero positions in literature data with arbitary small number
        for i in np.arange(0,len(exp_data)):
            if exp_data[i,1]==0:
                exp_data[i,1]=10**-5
        
        # %% Initialize collated results array
        exp_data_shape=np.shape(exp_data) #Shape of experimental data
        collated_results=np.append(exp_data,np.zeros((exp_data_shape[0],1)),axis=1) #Initialize vector which ahs expeirmental and model results. First three columns ar eexperimental data, last column is model fit
        
        # %% Calculate model values to put in collated vector

        ti_old=-1
        for i in np.arange(0,len(collated_results)): #loop over all expeirmental values
            ti_new=collated_results[i,0] #grab time-point to use
            if ti_old==ti_new:
                pass
            else: #when we have a new time to extract
                ti_old=ti_new
                t_exp_subset=collated_results[collated_results[:, 0] == ti_old, :] #grab only experimental values at that time-point
                t_mod_subset=mod_data[mod_data[:, 0] == ti_old, :] #grab model values at this timepoint
                if ti_old==14:
                    t_mod_subset=mod_data[71400:71450,:]
                elif ti_old==28:
                    t_mod_subset=mod_data[142800:142850,:]
                for j in np.arange(0,len(t_exp_subset)):
                    x_inter=t_exp_subset[j,1] #position value which needs its corresponding model value interpolated
                    for k in np.arange(0,len(t_mod_subset)):
                        if x_inter>t_mod_subset[k,1] and x_inter<t_mod_subset[k+1,1]: #check if x-value is greater current model-x and less than next (where to interpolate)
                            collated_results[i+j,3] = t_mod_subset[k,2]+(t_mod_subset[k+1,2]-t_mod_subset[k,2])*(x_inter-t_mod_subset[k,1])/(t_mod_subset[k+1,1]-t_mod_subset[k,1]) #standard interpolation formula
                        else: pass
        
        # %% Plot Collated results for visualizing "fit"
        ti_old=-1
        for i in np.arange(0,len(collated_results)): #loop over all expeirmental values
            ti_new=collated_results[i,0] #grab time-point to use
            if ti_old==ti_new:
                pass
            else: #when we have a new time to extract
                ti_old=ti_new
                collated_subset=collated_results[collated_results[:, 0] == ti_old, :] #grab only experimental values at that time-point
                plt.figure(0, figsize=[6, 6])
                plt.plot(collated_subset[:,1],collated_subset[:,2],'o',label=f'Literature value for t={ti_old}')
                plt.plot(collated_subset[:,1],collated_subset[:,3],label=f'Model Value for t={ti_old}')
        plt.legend(loc=[0.75,0],fontsize=7)
        plt.xlim(left=0,right=1.4)
        plt.rcParams['figure.dpi'] = 300
        modelfit_filename_partial=f'Modelfitplot{pc_i}.png'
        modelfit_filename_full=os.path.join(internal_export_path,modelfit_filename_partial)
        plt.savefig(modelfit_filename_full)
        plt.close()
    return 0

# %%
"""
Created on 11/2/2022

@author: joshuaprince

Purpose; Script to extract data from Tseng fits. Hmm, once I did the edited on excel this is actually a pretty simple code. Oh well. 
Version 1.0

    
"""