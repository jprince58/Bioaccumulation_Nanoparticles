#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_csv_generator=1.0

import numpy as np
import pandas as pd
import os

def csv_generator(c_set,parameter_combos_count,parameter_matrix,direct_export_path,new_count_number,machine_number):
    for pc_i in np.arange(0,parameter_combos_count,1): #Begin for loop to plot the different model paramters using MOL 
        cb=c_set[pc_i][0] #Grab current bound concentration data to plot
        cu=c_set[pc_i][1] #Grab current unbound concentration data to plot
        pot=c_set[pc_i][2] #Grab current potential data to plot
        t=c_set[pc_i][7] #Grab time-vector for this parameter set for plotting
        t_oom=np.floor(np.log10(np.abs(t[1]))).astype(int) #order of magninutde of time-step
        t=np.around(t,-t_oom) #Remove some rounding errors from time data-set
        x=c_set[pc_i][9] #Grab the position-vector for this parameter set for plotting
        x=np.around(x,3) #Remove some rounding errors from position data-set
        cb_data=pd.DataFrame(data=cb, index=x, columns=t) #Create pandas dataframe for bound concentration
        cu_data=pd.DataFrame(data=cu, index=x, columns=t) #Create pandas dataframe for unbound concentration
        pot_data=pd.DataFrame(data=pot, index=x, columns=t) #Create pandas dataframe for potential
        cb_data_filename_partial=f'bounddata{pc_i}_{new_count_number}-{machine_number}.csv'
        cb_data_filename_full=os.path.join(direct_export_path,cb_data_filename_partial)
        cu_data_filename_partial=f'unbounddata{pc_i}_{new_count_number}-{machine_number}.csv'
        cu_data_filename_full=os.path.join(direct_export_path,cu_data_filename_partial)
        pot_data_filename_partial=f'potentialdata{pc_i}_{new_count_number}-{machine_number}.csv'
        pot_data_filename_full=os.path.join(direct_export_path,pot_data_filename_partial)
        cb_data.to_csv(cb_data_filename_full) #Export bound dataframe as csv
        cu_data.to_csv(cu_data_filename_full) #Export unbound dataframe as csv
        pot_data.to_csv(pot_data_filename_full) #Export unbound dataframe as csv
    return vn_csv_generator


"""
Created on Mon Nov  2 20:02:34 2020

@author: joshuaprince

Purpose of Code: Take in concentration data, convert numpy arrays to pandas dataframes as .csv files, then export to direct exports

Version 1.0

Changes from Version 0.1 to 1.0 (2/20/2022 1:30 am)
    Got initial ocde to run. 
    
"""