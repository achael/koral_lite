import pylab
import matplotlib.pyplot as plt
import matplotlib
import os
import numpy as np
import json
import glob

###############################################
# Input file directory ########################
###############################################
folder_path = '../diagnostics/'

#Load sorted list of data
list_of_files = sorted( filter( os.path.isfile,
                        glob.glob(folder_path + '*') ) )

###############################################
# Initialize plotting variables ###############
###############################################
zone=0 #zone to plot (incex of u/p tuples i.e. cpu_prims[zone][[cell],[prims]])?
NV=0 #number of variables, grabs at 0th file
i=0

#arrays for plotting - TODO add more as needed
# Note arrays filled as cpu_u[NV][Nfiles] so cpu_u[0] returns array of u[0](t)
t_arr = []
cpu_u = []
cpu_p = []
gpu_u = []
gpu_p = []

for filename in list_of_files:
  with open(filename, 'r') as f:
    y=json.load(f)
    tin = y["time"]
    cpu_pin = y["cpu_prims"] #can be an array
    cpu_uin = y["cpu_cons"] #can be an array
    gpu_pin = y["gpu_prims"] #can be an array
    gpu_uin = y["gpu_cons"] #can be an array
    
    #if initializing, create empty arrays to append primitives over time
    if i==0:
      NV = len(cpu_uin[zone][1])
      for j in range(NV):
        cpu_p.append([])
        cpu_u.append([])
        gpu_p.append([])
        gpu_u.append([])
      
    #append new data
    t_arr.append(tin)
    for j in range(NV):
      cpu_p[j].append(cpu_pin[zone][1][j])
      cpu_u[j].append(cpu_uin[zone][1][j])
      gpu_p[j].append(gpu_pin[zone][1][j])
      gpu_u[j].append(gpu_uin[zone][1][j])

    #increment counter
    i+=1  
      
print("Json files loaded...")
    
###############################################
# Plotting variables ##########################
###############################################    
#TODO - get automated matplotlib array sizes
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(9, 9))
tfs=32 #title font size
fs=16 #font size for labels

#host primitives
for i in range(0,3):
  for j in range(0,3):
    nv = 3*i+j
    label = 'p['+str(nv)+']'
    ax[i][j].plot(t_arr,cpu_p[nv],label=r'CPU')
    ax[i][j].set_ylabel(label,fontsize=fs)
    
#device primitives
for i in range(0,3):
  for j in range(0,3):
    nv = 3*i+j
    ax[i][j].plot(t_arr,gpu_p[nv],marker='o',ls='',label=r'GPU')
    ax[i][j].legend()
    
ax[2][1].set_xlabel(r'time ($t_g$)',fontsize=fs)

plt.suptitle('Primitive Evolution',fontsize=tfs)
plt.tight_layout()
plt.show()


#clear and do conserved variables
plt.close('all')
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(9, 9))

#host primitives
for i in range(0,3):
  for j in range(0,3):
    nv = 3*i+j
    label = 'u['+str(nv)+']'
    ax[i][j].plot(t_arr,cpu_u[nv],label=r'CPU')
    ax[i][j].set_ylabel(label,fontsize=fs)
    
#device primitives
for i in range(0,3):
  for j in range(0,3):
    nv = 3*i+j
    ax[i][j].plot(t_arr,gpu_u[nv],marker='o',ls='',label=r'GPU')
    ax[i][j].legend()
    
ax[2][1].set_xlabel(r'time ($t_g$)',fontsize=fs)

plt.suptitle('Conserved Evolution',fontsize=tfs)
plt.tight_layout()
plt.show()

#I don't really find this one useful for now
#clear and do diff
plt.close('all')
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(9, 9))

#host primitives
for i in range(0,3):
  for j in range(0,3):
    nv = 3*i+j
    label = 'nv='+str(nv)
    d_prim = np.subtract(cpu_p[nv],gpu_p[nv]) 
    d_con = np.subtract(cpu_u[nv],gpu_u[nv]) 
    ax[i][j].plot(t_arr,d_prim,label=r'Primitive')
    ax[i][j].plot(t_arr,d_con,label=r'Conserved')
    ax[i][j].set_ylabel(label,fontsize=fs)
    ax[i][j].legend()
    
ax[2][1].set_xlabel(r'time ($t_g$)',fontsize=fs)

plt.suptitle('Difference',fontsize=tfs)
plt.tight_layout()
plt.show()
