 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Code Version Numbers
vn_N3=0.1 #See meta-data at bottom for details
vn_Main_Code=0.1 #See meta-data fro details


"""Machine code feature is commented out as we aren't running the code on multiple machines for now, just my local machine Replaced with simple direct_export_path and internal_export_path variables"""
# %% Current Machine Running Code with neccesary adjustments
machine_number=1 #Input the machine you are running this code on

"""Machine codes:
    Josh's Thinkpad- 1

"""

if machine_number == 1: #Code to assign the right export paths for each machine, since each is unique
   direct_export_path=r'C:\Users\joshu\Box Sync\Quantum Biofilms\Raw Data\N3_results\Direct Exports'
   internal_export_path=r'C:\Users\joshu\Box Sync\Quantum Biofilms\Raw Data\N3_results\Internal Exports' 

   

# %% I don't remember what this is for but I'm scared to delete it 
from IPython import get_ipython
get_ipython().magic('reset -sf')

# %% Imports
import time
import os
import numpy as np
from N3_RJ import *
from N3_method_of_lines import *
from N3_parameter_checker import *
from N3_parameter_matrix import *
from N3_report_generator import *
from N3_csv_generator import *
from N3_linear_fitting import *

# %% Start Timer
t_start=time.time()

# %% Determine what run number this is then update run number
counter_file = open(f"counter_file_{machine_number}.txt",'r+')
old_count_number=int(counter_file.read())
new_count_number=str(old_count_number+1)
counter_file = open(f"counter_file_{machine_number}.txt",'w')
counter_file.write(new_count_number)
counter_file.close()

# %%Inputs Code Block
h=np.array([0.01]) #Define timesteps to test
tol=np.array([10**(-8)])  #Define the tolerance the code will run with when running Newton-Rhapson
t1=np.array([0]) #Define initialtime vector of values to test
t2=np.array([2.5]) #Final Time
nx=np.array([200]) #Mesh size
gam=np.array([0.1]) #Define dimenionless ratio of diffusivities to test
F=np.array([1]) #Define the dimensionless forward reaction rate constant to test
K=np.array([1]) #Define the Eqilbrium constant for NP binding
eps=np.array([1.5]) #Define ratio of total NP binding sites to supernatant NP concentration
omega=np.array([1]) #Define ratio of NP contribution to electrical potential profile to their electrokinetic mobility
ups= np.array([1])#Define ratio of Biofilm cotridubtion to electrical potential profile to NP electrokinetic mobility
Kp= np.array([10]) #Define partition coeffecient of NP into biofilm at water-biofilm interface
ci=10**(-10) #Define the inital concentration in the biofilm (Can't be zero, if one wants to be zero, set it to a very small number instead)


# %% Generate Parameter Matrix for Testing
[parameter_matrix,parameter_combos_count,vn_parameter_matrix_generator]=parameter_matrix_generator(h,tol,t1,t2,nx,gam,F,K,eps,omega,ups,Kp)
                    
# %% Run parameters through numerical model (Heart of the Code)               
[c_set,vn_parameter_checker,vn_method_of_lines,vn_RJ] = parameter_checker(parameter_matrix,ci) #output the set of concentration over time and space results for each set of parameters tested

# %% Export results to csv files
vn_csv_generator = csv_generator(c_set,parameter_combos_count,parameter_matrix,direct_export_path,new_count_number,machine_number)

# %% Fit model to first order approximation, plot approximation, and determine fit of approximation
[perc_acc_matrix,vn_linear_fitting]=linear_fit(c_set,parameter_combos_count,parameter_matrix,internal_export_path)

# %% Report Generator: Exports Plots as Word Document to Seperate Directory (see file N2_report_generator.py)
report=plot_generator(c_set,parameter_combos_count,parameter_matrix,new_count_number,vn_N2,vn_Main_Code,vn_parameter_matrix_generator,vn_parameter_checker,vn_csv_generator,vn_method_of_lines,vn_RJ,perc_acc_matrix,vn_linear_fitting,machine_number,internal_export_path)

# %% Stop Timer
#End timer
t_end=time.time()
total_time=t_end-t_start
print('Total time is {} sec'.format(total_time))

#%% To export report, turn on this code block
#Finish Report
para5=report.add_paragraph(f'Time to Run (sec): {total_time}     ')
report_filename_partial=f'N2_report{new_count_number}-{machine_number}.docx'
report_filename_full=os.path.join(direct_export_path,report_filename_partial)
report.save(report_filename_full)

"""
Created on Tue Jun 16 16:13:07 2020

@author: joshuaprince

Overall N3.1 Computer Model Meta-data
Version 0.1

This version of our model code improves from N2 by cinorporating electrodynamics via Gauss Law and direct incorporation of electrical potential. 
Dropped the hill coeffecient binding kinetics because couldn't find examples in enough outside sources. Also made diffusion coeffecient and electrical potential profiles from biofilm linear instead of quadratic because previous literatre indicated they were pretty close to linear. 


Main_Code File Meta-data
Version 0.1

This main code is what is the script which is actually run to "run the model". Has lots of sub-function calls and such to the core components to the model.
This is where model parameters are entered. They can be entered as vectors to run the code over mutiple parameter sets. 
"""



