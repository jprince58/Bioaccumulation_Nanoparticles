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


def experimental_data_extractor(experimental_data_file):
    data = pd.read_csv(experimental_data_file, header=None) #Import data into pandas dataframe
    data=data.to_numpy() #Immediately convert to a numpy array
    
    # %%Converting Data into Dimensionless Fomr
    #First, convert distance into dimensionless form
    min_x=np.min(data[:,1]) #Calculate minimum distance
    data[:,1]=data[:,1]-min_x #Substract minimum baseline
    max_x=np.max(data[:,1]) #Calculate max distance
    data[:,1]=data[:,1]/max_x #Divide each distance by max-distance
    data[:,1]=1-data[:,1]#since in the model the biofilm-water interface is x=1, will need to convert distances to reflect that
     
    #No need to convert absorbance values. that will be handled later. 
    
    # %% Remove baseline from experimental AU values
    baseline=8.20 #This is apprently the baseline AU value, when the experimental results should really be zero (eye-balled it)
    data[:,2]=data[:,2]-baseline
    
    # %%Fit x=0 values to polynomial fit 
    
    #Grab data to fit
    fit_data=np.zeros((1,2))
    fit_data[0,0]=data[0,0]
    fit_data[0,1]=data[0,2]
    ti_old=0
    for i in np.arange(0,len(data)): #loop over all expeirmental values
        ti_new=data[i,0] #grab time-point to use
        if ti_old==ti_new:
            pass
        else: #when we have a new time found
            ti_old=ti_new
            add_data=np.zeros((1,2))
            add_data[0,0]=data[i,0]
            add_data[0,1]=data[i,2]
            fit_data=np.append(fit_data,add_data, axis=0)
    
    #Fit data to polynomial
    fit_res=np.polyfit(fit_data[:,0],fit_data[:,1],6,full=True)
    
    #Plot fit to check how good it is
    appender=np.zeros((len(fit_data),1))
    fit_data=np.append(fit_data,appender,axis=1)
    for i in np.arange(0,len(fit_data)):
        t=fit_data[i,0]
        fit_coeff=fit_res[0]
        fit_data[i,2]=fit_coeff[0]*t**6+fit_coeff[1]*t**5+fit_coeff[2]*t**4+fit_coeff[3]*t**3+fit_coeff[4]*t**2+fit_coeff[5]*t**1+fit_coeff[6]*t**0
        
    plt.plot(fit_data[:,0],fit_data[:,1],label='Literature Values')   
    plt.plot(fit_data[:,0],fit_data[:,2],label='Fit Values')  
    plt.legend()
    plt.pause(1)
    
    return [data,fit_coeff]

# %%
"""
Created on 11/2/2022

@author: joshuaprince

Purpose; Script to extract data from Tseng fits. Hmm, once I did the edited on excel this is actually a pretty simple code. Oh well. 
Version 1.0

    
"""