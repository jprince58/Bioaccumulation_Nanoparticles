 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Code Version Numbers
vn_N6=1.0 #See meta-data at bottom for details
vn_Main_Code=1.0 #See meta-data fro details

# %% Current Machine Running Code with neccesary adjustments
"""Machine code feature is commented out as we aren't running the code on multiple machines for now, just my local machine Replaced with simple direct_export_path and internal_export_path variables"""
# %% Current Machine Running Code with neccesary adjustments
machine_number=1 #Input the machine you are running this code on

"""Machine codes:
    Josh's Thinkpad- 1

"""

if machine_number == 1: #Code to assign the right export paths for each machine, since each is unique
   direct_export_path=r'C:\Users\joshu\Box\Quantum Biofilms\Raw Data\N14_results\Direct Exports'
   internal_export_path=r'C:\Users\joshu\Box\Quantum Biofilms\Raw Data\N14_results\Internal Exports' 

   
# %% I don't remember what this is for but I'm scared to delete it 
from IPython import get_ipython
get_ipython().magic('reset -sf')

# %% Imports
import time
import os
import numpy as np
from N14_RJ import *
from N14_method_of_lines import *
from N14_parameter_tester import *
from N14_parameter_matrix import *
from N14_report_generator import *
from N14_csv_generator import *
from N14_linear_fitting import *
from N14_experimental_data_extractor import *
from N14_exp_data_fitter import *

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
h=np.array([0.01*2/3]) #Define timesteps to test
tol=np.array([10**(-8)])  #Define the tolerance the code will run with when running Newton-Rhapson
t1=np.array([0]) #Define initialtime vector of values to test
t2=np.array([20]) #Final Time
nx=np.array([100]) #Mesh size
omega=np.array([1]) #Define effective diffusivity 
mu=np.array([0.05]) #Define dimensionless and porosity adjusted binding rate constant
nu=np.array([0]) #Define dimnesionless binding site density difference
kappa=np.array([5]) #Define dimensionless and porosity adjusted equilibrium constant
a=np.array([7.5]) #Define shape parameter for binding site profile
b=np.array([30]) #Define shape parameter for intersitital porosity profile
c=np.array([1]) #Define dimensionless minimum binding site concentration
phi_min=np.array([0.25]) #Define minimum interstitial porosity
phi_max=np.array([0.8]) #Define maximum interstitial porosity
rho=np.zeros(len(phi_min)) #initialize rho vector
eps=np.zeros(len(phi_min)) #initialize epsilon vector
for i in np.arange(0,len(rho)):
    rho[i]=phi_max
    eps[i]=0 #hard-coding in beta equal to zero to simualte no porosity profile
    # eps[i]=(phi_max[i]-phi_min[i])/(np.exp(b)+1)
Kp=np.array([1]) #Define partition coeffecient
kconv=130 #guess a AU to particle conversion factor
ci=10**(-10) #Define the inital concentration in the biofilm (Can't be zero, if one wants to be zero, set it to a very small number instead)


#%% Grab Experimental Results to fit to model
# experimental_data_file=r'C:\Users\joshu\Box\Quantum Biofilms\Processed Data\Extracted data from literature\tseng_fits_Fig2B_Cyt5tob_incubation.csv'
experimental_data_file=r'C:\Users\joshu\Box\Quantum Biofilms\Processed Data\Extracted data from literature\tseng_fits_Fig2B_Cy5cipro_incubation_bump.csv'
[experimental_results,fit_coeff] = experimental_data_extractor(experimental_data_file)

# %% Generate Parameter Matrix for Testing
[parameter_matrix,parameter_combos_count,vn_parameter_matrix_generator]=parameter_matrix_generator(h,tol,t1,t2,nx,omega,mu,nu,eps,rho,kappa,a,b,c,Kp)
                    
# %% Run parameters through numerical model (Heart of the Code)               
[c_set,vn_parameter_checker,vn_method_of_lines,vn_RJ] = parameter_checker(parameter_matrix,ci,fit_coeff,kconv) #output the set of concentration over time and space results for each set of parameters tested
#[c_set,vn_parameter_checker,vn_method_of_lines,vn_RJ] = parameter_checker(parameter_matrix,ci) #output the set of concentration over time and space results for each set of parameters tested

# %% Compare Model to Experimental Data
fitting_results=exp_data_fitter(c_set,experimental_results,parameter_combos_count,internal_export_path,kconv,t2,direct_export_path,new_count_number,machine_number)

# %% 

# %% Export results to csv files
vn_csv_generator = csv_generator(c_set,parameter_combos_count,parameter_matrix,direct_export_path,new_count_number,machine_number)

# %% Fit model to first order approximation, plot approximation, and determine fit of approximation
#[perc_acc_matrix,vn_linear_fitting]=linear_fit(c_set,parameter_combos_count,parameter_matrix,internal_export_path)

# %% Report Generator: Exports Plots as Word Document to Seperate Directory (see file N2_report_generator.py)
#report=plot_generator(c_set,parameter_combos_count,parameter_matrix,new_count_number,vn_N6,vn_Main_Code,vn_parameter_matrix_generator,vn_parameter_checker,vn_csv_generator,vn_method_of_lines,vn_RJ,perc_acc_matrix,vn_linear_fitting,machine_number,internal_export_path)
report=plot_generator(c_set,parameter_combos_count,parameter_matrix,new_count_number,vn_N6,vn_Main_Code,vn_parameter_matrix_generator,vn_parameter_checker,vn_csv_generator,vn_method_of_lines,vn_RJ,machine_number,internal_export_path,phi_min,phi_max)

# %% Stop Timer
#End timer
t_end=time.time()
total_time=t_end-t_start
print('Total time is {} sec'.format(total_time))

#%% To export report, turn on this code block
#Finish Report
para5=report.add_paragraph(f'Time to Run (sec): {total_time}     ')
report_filename_partial=f'N14_report{new_count_number}-{machine_number}.docx'
report_filename_full=os.path.join(direct_export_path,report_filename_partial)
report.save(report_filename_full)

"""
Created on Tue Jun 16 16:13:07 2020

@author: joshuaprince

Overall N2.1 Computer Model Meta-data
Version 0.2

The code runs again and exports the report.

"""

