#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 14 16:13:07 2020

@author: joshuaprince
"""


#Inputs Code Block
"""
This cell take in the inputs to this code. This does so by:
1a) takes in the user-specified dimensional values for the code, or 
1b) takes in user-specified dimensionless numbers then
2) defines the relevant inputs and outputs which go to the heart of the code 

This is the primary code one should be modifying when testing different physical systems for the model
"""


import time
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from numpy import linalg
from numpy.linalg import solve
from scipy import sparse
from scipy.sparse import linalg
from scipy.sparse import csc_matrix

#Editing Code

#Start Timer
t_start=time.time()

def RJ(x,y,p):
    nx=len(x)-1 #Grab the mesh size for position
    ny=len(y)-2 #Grab number of y-points
    dx=1/nx #Calculate the distance between nodes (assumes domain is from 0 to 1)
    gam=p[0] #Grab the dimenionless ratio of diffusivities
    beta=p[1] #Grab the dimensionless ratio of potentials
    F=p[2] #Grab the dimensionless forward reaction rate constant
    Re=p[3] #Grab the dimensionless reverse reaction rate constant
    n=p[4] #Grab the hill coeffecient
    R=np.zeros(ny+2) #Initialize R
    J=np.zeros((ny+2,ny+2)) #Initialize J
    index=np.arange(0,ny+2) #Creater index vector
    for i in index:
        if i==0 :
            R[i]=y[i]-1
            J[i,i]=1;
        elif i%2==1 :
            R[i]=F*y[i-1]**n*(1-y[i])-Re*y[i]
            J[i,i]=-(F*y[i-1]**n+Re)
            J[i,i-1]=n*F*y[i-1]**(n-1)*(1-y[i])
        elif i==ny:
            l=int(i/2)
            R[i]=2*(y[i-2]-y[i])/dx**2*(1-x[l]**2+gam)-y[i]*(F*y[i]**(n-1)*(1-y[i+1])-2*beta)+Re*y[i+1]
            J[i,i]=(-2/dx**2)*(1-x[l]**2+gam)-(n*F*y[i]**(n-1)*(1-y[i-1]-2*beta))
            J[i,i-2]=2/dx**2*(1-x[l]**2+gam)
            J[i,i+1]=F*y[i]**n+Re
        elif i%2==0 and i!=0:
            l=int(i/2)
            R[i]=(y[i+2]-2*y[i]+y[i-2])/dx**2*(1-x[l]**2+gam)-(y[i+2]-y[i-2])/2/dx*(2*x[l]*(1-beta))-y[i]*(F*y[i]**(n-1)*(1-y[i+1])-2*beta)+Re*y[i+1]
            J[i,i]=(-2/dx**2)*(1-x[l]**2+gam)-(n*F*y[i]**(n-1)*(1-y[i-1]-2*beta))
            J[i,i+2]=1/dx**2*(1-x[l]**2+gam)-1/2/dx*(2*x[l]*(1-beta))
            J[i,i-2]=1/dx**2*(1-x[l]**2+gam)+1/2/dx*(2*x[l]*(1-beta))
            J[i,i+1]=F*y[i]**n+Re
        else:
            print('Uh oh')
    return (R,J)



def MOL(t,x,y,h,p,tol):
    yw=y[:,0] #Initalize the working concentration vector
    index=np.arange(0,len(t)) #Creater index vector
    whoops=0 #initialize error function
    for i in index: #Begin for loop which iterates over all the entries in the time vector, assigning them the value td
        if i==0: continue
        else:
            yold=yw #update the old y-value
            yw[0]=1 #hardcode in boundary condition
            out=RJ(x,yw,p); #Calculate Residual and Jacobian from new y value
            R=out[0] #Grab the Residual
            J=out[1] #Grab the Jacobian
            #print('This is the residual')
            #print(R)
            R=yw-yold-h*R
            k=0
            while np.linalg.norm(R)>tol :
                k=k+1
                J=np.eye(len(yw))-h*J; #Calculate new Jacobian from new y
                #print('This is the Jacobian and the determinant')
                #print(J)
                #print(np.linalg.det(J))
                J=sp.sparse.csc_matrix(J)
                dif=-sp.sparse.linalg.spsolve(J,R) #Apply built in sparse Linear solver to find delta from J and R
                #dif=-np.linalg.solve(J,R)  #Regular Solver. Just keeping it in there in case the sparse solver isn't working for some reason
                #print('This is delta')
                #print(dif)
                #This section is to determine if the differecne will make the working y outside the domain, and apply a different change
                #ywc=yw+dif
                
                #This is the working area
                yw=yw+dif ; #Update y
                #print('This is the new y')
                #print(yw)
                out=RJ(x,yw,p); #Calculate Residual and Jacobian from new y value
                R=out[0] #Grab the Residual
                #print('This is the residual')
                #print(R)
                J=out[1] #Grab the Jacobian
                R=yw-yold-h*R ; #Update Residual
                if k>100:
                    print('Whoops: Did not converge')
                    whoops=whoops+1
                    break
            y[:,i]=yw
    return (y,whoops)
    
#Working code to test MOL, turn off when testing with Inputs Code Block
h=0.02 #Define timestep
t1=0 #Initial time
t2=4 #Final Time
nx=100 #Mesh size
ny=2*nx #Solution vector size
x=np.linspace(0,1,nx+1) #Define x (note, if x doesn't range from 0 to 1, should edit this)
t=np.arange(t1,t2+h,h) #Define t
nt=len(t) #Define the number of timepoints
y=np.zeros((ny+2,nt)) #Initialize y
y=y+10**(-8) #Make starting values not exactly equal to zero (divide by zero erros pop up)
p=[0.1,1,10,0.1,1]
tol=10**(-8)

[c,whoops]=MOL(t,x,y,h,p,tol) #Find the concntration profiles in space and time using Method of Lines (MOL)
print('you whoopsed {} many times'.format(whoops))
cb=np.zeros((nx+1,nt)) #Initalize new concentration array where bound and unbound NP concentrations are "unpacked" such that they occupy two different matrices in the same 3-D array 
cu=np.zeros((nx+1,nt))
xindex=np.arange(0,nx+1)
for i in xindex:
    j=2*i #secondary index (position of unbound NP concentration in original concentration matrix)
    k=2*i+1 #Secondary index (position of bound NP concentration in original concentration matrix)
    cu[i,:]=c[j,:]
    cb[i,:]=c[k,:]
    
#Find Average Unbound Concentration Overtime
average_conc_overtime=np.zeros(nt)
t2index=np.arange(0,nt)
for i in t2index:
    average_conc_overtime[i]=np.average(cu[:,i])

#Find Change in Concentration overtime
change_in_concentration=np.zeros(nt)
t3index=t2index[:-1].copy()
for i in t3index:
    change_in_concentration[i]=(average_conc_overtime[i+1]-average_conc_overtime[i])/h

"""Plotting"""

#Find relelvant maximums and minimums
upper_1 = np.amax(cu)*1.1 #Upper bound on Unbound Concentration
upper_2 = np.amax(average_conc_overtime)*1.1 #Upper Bound on Average Unbound Concentration
upper_3 = np.amax(change_in_concentration)*1.1 #Upper Bound on Change in Average Concentration

#Unbound
tp=20 #number of time points to plot
space=int((nt-1)/tp)
tindex=np.arange(0,nt,space)
plt.figure(0)
for i in tindex:
    cc=cu[:,i]
    ti=round(t[i],1)
    plt.plot(x,cc,label='t={}'.format(ti))
plt.xlim(left=0,right=1)
plt.ylim(bottom=0,top=upper_1)
plt.xlabel('Position',fontsize=14)
plt.ylabel('Dimensionless Concentration',fontsize=14)
plt.title('Dimensionless Unbound Concentration plot',fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc=(1.02,0))

#Bound
tp=5 #number of time points to plot
space=int((nt-1)/tp)
tindex=np.arange(0,nt,space)
plt.figure(1)
for i in tindex:
    cc=cb[:,i]
    ti=round(t[i],1)
    plt.plot(x,cc,label='t={}'.format(ti))

plt.xlim(left=0,right=1)
plt.ylim(bottom=0,top=1)
plt.xlabel('Position',fontsize=14)
plt.ylabel('Dimensionless Concentration',fontsize=14)
plt.title('Dimensionless Bound Concentration plot',fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc=(1.02,0))

#Average Concetration Overtime
plt.figure(2)
plt.plot(t,average_conc_overtime)
plt.xlim(left=t1,right=t2)
plt.ylim(bottom=0,top=upper_2)
plt.xlabel('Time',fontsize=14)
plt.ylabel('Dimensionless Concentration',fontsize=14)
plt.title('Average Dimensionless Unbound Concentration plot',fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

#Change in Concentration vs Concentration
plt.figure(3)
plt.plot(average_conc_overtime,change_in_concentration)
plt.xlim(left=0,right=upper_2)
plt.ylim(bottom=0,top=upper_3)
plt.xlabel('Concentration',fontsize=14)
plt.ylabel('Dimensionless Change in Concentration',fontsize=14)
plt.title('dC vs C plot',fontsize=16)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)

#print(average_conc_overtime)
#print(change_in_concentration)

#End timer
t_end=time.time()
total_time=t_end-t_start
print('Total time is {} sec'.format(total_time))
