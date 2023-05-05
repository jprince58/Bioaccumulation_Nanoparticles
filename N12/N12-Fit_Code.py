 
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
   direct_export_path=r'C:\Users\joshu\Box\Quantum Biofilms\Raw Data\N12_results\Direct Exports'
   internal_export_path=r'C:\Users\joshu\Box\Quantum Biofilms\Raw Data\N12_results\Internal Exports' 

   
# %% I don't remember what this is for but I'm scared to delete it 
from IPython import get_ipython
get_ipython().magic('reset -sf')

# %% Imports
import time
import os
import numpy as np
from N12_RJ import *
from N12_method_of_lines import *
from N12_parameter_tester import *
from N12_parameter_matrix import *
from N12_report_generator import *
from N12_csv_generator import *
from N12_linear_fitting import *
from N12_experimental_data_extractor import *
from N12_exp_data_fitter import *
from N12_residual_calc import *
from lmfit import minimize, Parameters, Parameter, report_fit, Model

#Way to keep iterations from going on forever
def check(params, iter, resid):
    if iter>5:
        x=True
    else:
        x=False
    return x


# %% Start Timer
t_start=time.time()

# %% Determine what run number this is then update run number
counter_file = open(f"counter_file_{machine_number}.txt",'r+')
old_count_number=int(counter_file.read())
new_count_number=str(old_count_number+1)
counter_file = open(f"counter_file_{machine_number}.txt",'w')
counter_file.write(new_count_number)
counter_file.close()

#%% Grab Experimental Results to fit to model
experimental_data_file=r'C:\Users\joshu\Box\Quantum Biofilms\Processed Data\Extracted data from literature\tseng_fits_Fig2B_Cy5cipro_incubation_bump.csv'
[experimental_results,fit_coeff] = experimental_data_extractor(experimental_data_file)

# # %% Compare Model to Experimental Data
# fitting_results=exp_data_fitter(c_set,experimental_results,parameter_combos_count,internal_export_path,kconv,t2)

# %%Inputs Code Block
# set parameters including bounds for model; you can also fix parameters (use vary=False)
params = Parameters()
params.add('h', value=0.05/2, vary=False)
params.add('tol', value=10**(-6), vary=False)
params.add('t1', value=0, vary=False)
params.add('t2', value=15, vary=False)
params.add('nx', value=50, vary=False)
# params.add('omega', value=0.8, min=0.5, max=1.)
params.add('omega', value=1, vary=False)
params.add('mu', value=2, min=0.01, max=10.)
# params.add('mu', value=5, vary=False)
params.add('nu', value=10, min=0.01, max=100.)
# params.add('nu', value=10, vary=False)
# params.add('phi_min', value=0.3, min=0.1, max=1)
params.add('phi_min', value=0.3, vary=False)
# params.add('phi_max', value=0.8, min=0.1, max=1)
params.add('phi_max', value=0.8, vary=False)
params.add('kappa', value=0.75, min=0.1, max=10.)
# params.add('kappa', value=0.75, vary=False)
# params.add('a', value=1.5, min=0.2, max=50.) 
params.add('a', value=1.5, vary=False)  
params.add('b', value=25, min=0, max=100.)
# params.add('b', value=25, vary=False)
# params.add('c', value=0.42, min=0.01, max=10.)
params.add('c', value=0.42, vary=False)
# params.add('Kp', value=5, min=0.01, max=10.)
params.add('Kp', value=5.33, vary=False)
params.add('kconv', value=188.8, vary=False)

# fit models
# fit_results = minimize(residual_calc, params, method='leastsq', iter_cb=check(params, iter, resid))  # fitting for model    
fit_results = minimize(residual_calc, params, method='leastsq')  # fitting for model    

# %% Grab parameters from model fit
h=np.array([fit_results.params['h'].value]) #Define timesteps 
tol=np.array([fit_results.params['tol'].value])  #Define the tolerance the code will run with when running Newton-Rhapson
t1=np.array([fit_results.params['t1'].value]) #Define initialtime vector of values 
t2=np.array([fit_results.params['t2'].value]) #Final Time
nx=np.array([fit_results.params['nx'].value]) #Mesh size
omega=np.array([fit_results.params['omega'].value]) #Define effective diffusivity 
mu=np.array([fit_results.params['mu'].value]) #Define dimensionless and porosity adjusted binding rate constant
nu=np.array([fit_results.params['nu'].value]) #Define dimensionless and porosity adjusted maximum binding site density
kappa=np.array([fit_results.params['kappa'].value]) #Define dimensionless and porosity adjusted equilibrium constant
a=np.array([fit_results.params['a'].value]) #Define shape paramter for binding site profile
b=np.array([fit_results.params['b'].value]) #Define shape paramter for intersitital porosity profile
phi_min=np.array([fit_results.params['phi_min'].value])
phi_max=np.array([fit_results.params['phi_max'].value])   
rho=np.zeros(len(phi_min)) #initialize rho vector
eps=np.zeros(len(phi_min)) #initialize beta vector     
for i in np.arange(0,len(rho)):
    rho[i]=1+phi_max[i]/(phi_max[i]-phi_min[i])*(np.exp(b)+1)
    eps[i]=(phi_max[i]-phi_min[i])/(np.exp(b)+1)
c=np.array([fit_results.params['c'].value]) #Define shape paramter for traditional proosity profile
Kp=np.array([fit_results.params['Kp'].value]) #Define partition coeffecient  
kconv=188.8 #guess a AU to particle conversion factor
ci=10**(-10) #Define the inital concentration in the biofilm (Can't be zero, if one wants to be zero, set it to a very small number instead)


#%% Grab Experimental Results to fit to model
experimental_data_file=r'C:\Users\joshu\Box\Quantum Biofilms\Processed Data\Extracted data from literature\tseng_fits_Fig2B_Cy5cipro_incubation_bump.csv'
[experimental_results,fit_coeff] = experimental_data_extractor(experimental_data_file)

# %% Generate Parameter Matrix for Testing
[parameter_matrix,parameter_combos_count,vn_parameter_matrix_generator]=parameter_matrix_generator(h,tol,t1,t2,nx,omega,mu,nu,eps,rho,kappa,a,b,c,Kp)
                    
# %% Run parameters through numerical model (Heart of the Code)               
[c_set,vn_parameter_checker,vn_method_of_lines,vn_RJ] = parameter_checker(parameter_matrix,ci,fit_coeff,kconv) #output the set of concentration over time and space results for each set of parameters tested
#[c_set,vn_parameter_checker,vn_method_of_lines,vn_RJ] = parameter_checker(parameter_matrix,ci) #output the set of concentration over time and space results for each set of parameters tested

# %% Compare Model to Experimental Data
fitting_results=exp_data_fitter(c_set,experimental_results,parameter_combos_count,internal_export_path,kconv,t2)

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

# %% To export report, turn on this code block
# Finish Report
para5=report.add_paragraph(f'Time to Run (sec): {total_time}     ')
report_filename_partial=f'N12_report{new_count_number}-{machine_number}.docx'
report_filename_full=os.path.join(direct_export_path,report_filename_partial)
report.save(report_filename_full)

"""
Created on Tue Jun 16 16:13:07 2020

@author: joshuaprince

Overall N2.1 Computer Model Meta-data
Version 0.2

The code runs again and exports the report.

"""

