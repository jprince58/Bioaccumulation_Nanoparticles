#!/usr/bin/env python3
# -*- coding: utf-8 -*-

vn_csv_generator=1.2

import numpy as np
import pandas as pd
import os

def csv_generator(c_set,parameter_combos_count,parameter_matrix,direct_export_path,new_count_number,machine_number):
    for pc_i in np.arange(0,parameter_combos_count,1): #Begin for loop to plot the different model paramters using MOL 
        cm=c_set[pc_i][3] #Grab current Mobile concentration data to plot
        ca=c_set[pc_i][4] #Grab current Attached concentration data to plot
        t=c_set[pc_i][1] #Grab time-vector for this parameter set for plotting
        x=c_set[pc_i][0] #Grab the position-vector for this parameter set for plotting
        t_oof=np.floor(np.log10(np.abs(t[1]))).astype(int) #order of magninutde of time-step
        t=np.around(t,-t_oof) #Remove some rounding errors from time data-set
        x=np.around(x,3) #Remove some rounding errors from position data-set
        ca_data=pd.DataFrame(data=ca, index=x, columns=t) #Create pandas dataframe for bound concentration
        cm_data=pd.DataFrame(data=cm, index=x, columns=t) #Create pandas dataframe for unbound concentration
        ca_data_filename_partial=f'bounddata{pc_i}_{new_count_number}-{machine_number}.csv'
        ca_data_filename_full=os.path.join(direct_export_path,ca_data_filename_partial)
        cm_data_filename_partial=f'unbounddata{pc_i}_{new_count_number}-{machine_number}.csv'
        cm_data_filename_full=os.path.join(direct_export_path,cm_data_filename_partial)
        ca_data.to_csv(ca_data_filename_full) #Export bound dataframe as csv
        cm_data.to_csv(cm_data_filename_full) #Export unbound dataframe as csv
    return vn_csv_generator


"""
Created on Mon Nov  2 20:02:34 2020

@author: joshuaprince

Purpose of Code: Take in concentration data, convert numpy arrays to pandas dataframes as .csv files, then export to direct exports

Version 1.3

changes from Version 1.2 (2/13/2021 4:00 pm)
    Added machine number to names of csv files as part of effort to make this code run on multiple machines

Changes from Version 1.2 to Version 1.1 (11/3/2020 7:45 am)
    Added run number labeling for csv-data exports so data isn't so easily lost
    Lines 29 & 31: '_{new_count_number}' added to file names
"""