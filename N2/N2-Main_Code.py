 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# %% Code Version Numbers
vn_N2=1.8 #See meta-data at bottom for details
vn_Main_Code=1.5 #See meta-data fro details

# %% Current Machine Running Code with neccesary adjustments
"""Machine code feature is commented out as we aren't running the code on multiple machines for now, just my local machine Replaced with simple direct_export_path and internal_export_path variables"""
# %% Current Machine Running Code with neccesary adjustments
machine_number=1 #Input the machine you are running this code on

"""Machine codes:
    Josh's Thinkpad- 1

"""

if machine_number == 1: #Code to assign the right export paths for each machine, since each is unique
   direct_export_path=r'C:\Users\joshu\Box\Quantum Biofilms\Raw Data\N2_results\Direct Exports'
   internal_export_path=r'C:\Users\joshu\Box\Quantum Biofilms\Raw Data\N2_results\Internal Exports' 

   
# %% I don't remember what this is for but I'm scared to delete it 
from IPython import get_ipython
get_ipython().magic('reset -sf')

# %% Imports
import time
import os
import numpy as np
from N2_RJ import *
from N2_method_of_lines import *
from N2_parameter_tester import *
from N2_parameter_matrix import *
from N2_report_generator import *
from N2_csv_generator import *
from N2_linear_fitting import *

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
t2=np.array([10]) #Final Time
nx=np.array([200]) #Mesh size
gam=np.array([1]) #Define dimenionless ratio of diffusivities to test
beta=np.array([0]) #Define the dimensionless ratio of potentials to test
F=np.array([5]) #Define the dimensionless forward reaction rate constant to test
Re=np.array([0.01]) #Define the dimensionless reverse reaction rate constant to test
n=np.array([1]) #Define the hill coeffecient to test
ci=10**(-10) #Define the inital concentration in the biofilm (Can't be zero, if one wants to be zero, set it to a very small number instead)


# %% Generate Parameter Matrix for Testing
[parameter_matrix,parameter_combos_count,vn_parameter_matrix_generator]=parameter_matrix_generator(h,tol,t1,t2,nx,gam,beta,F,Re,n)
                    
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

Overall N2.1 Computer Model Meta-data
Version 1.8

Changes from version 1.7 (2/13/2021 4:00 pm):
    Added new functionalities to run on different machines, including
    -machine_number variable
    -new internal/direct export paths for a given machine number
    -adding machine numbers to the end of all direct exports

Changes from version 1.6 to 1.7 (11/11/2020 9;00 pm):
    Added linear fits functions, along with neccesary report plotting

Changes from version 1.5 to 1.6 (11/3/2020 4:35 pm):
    Added version number tracker for all scripts

Changes from  version 1.4 and 1.5 (11/3/2020 8:15 am):
    In report_generator, added plot close functionalities for all plots (python was complaining about holding so many plots in memory)

Changes from version 1.4 to 1.3: (11/3/2020 8:00 am)
    Main_code and report_generator were changed so that the unbound animation isn't a direct export, but an internal export, since it is now included in the report.
    
Changes from version 1.3 to 1.2: (11/3/2020 7:50 am)
    Main_Code.py and csv_generator were changed so that when csv files with data were exported, they are labelled with run number so data isn't overridden everytime model is run



Main_Code File Meta-data
Version 1.6

Changes from version 1.5 to 1.6 (2/18/2022 10:44 pm):
    
changed machine and export paths for new laptop in Duke system

Changes from Version 1.4 to 1.5 (11/11/2020 9:00 pm):
    Added linear fit function

Changes from Version 1.3 to 1.4 (11/3/2020 9:00 am):
    Added Version number tracker for all scripts
    Line 5, 6: Added new manual entry variable which tracks version number for overall N2 model and Main_Code
    

Changes from Version 1.2 to Version 1.3 (11/3/2020 8:00 am):
    Changed animated gif from direct exporting to internal exporting
    Line 90: direct_export_path variable removed from function call

Changes from version 1.2 to version 1.1 (11/3/2020 7:45 am):
    Added run number labeling for csv-data exports so data isn't so easily lost
    Line 88: variable "new_count_number" was added to function pass
"""



