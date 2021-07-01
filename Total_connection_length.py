# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 14:20:28 2020

@author: tpngu28
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import csv
from Load_and_plot_data import LoadGrid, PlotRasters

#%%
### Calculate Probability distribution function of the connectivity length


### Create 2D list of array with size (nrows,ncols)
def Array2D(nrows, ncols):
    cc = [[ np.array([]) for col in range(ncols)] for row in range(nrows)]
    return cc

### Stretch the 2D list into 1D array
def list2array(array):
    array = np.array(array,dtype=object)
    
    for j in range(len(array)):
        for i in range(len(array[0])):
            if len(array[j][i]) == 0:
                array[j][i] = np.array([0])
    
    array_cc = array.reshape(-1)
    array_cc = np.concatenate(array_cc).astype(None) 
    return array_cc
    
######## save data to txt file
        
def Writedata(data, filename, path ):
    file_path =  path + filename
    file = open(file_path, 'a+')
    np.savetxt(file, data, fmt='%.1f', delimiter = ' ')
    file.close()

#%%
########### read Qx, Qy and WD file #####################
    
#path = 'E:\\Narran3_peaktime400\\result\\res-' 
#num = 198
#file_path_wd = path + str(num).zfill(4) + '.wd'
#file_path_Qx = path + str(num).zfill(4) + '.Qx'
#file_path_Qy = path + str(num).zfill(4) + '.Qy'
    
#change path with different model 
file_path_wd = 'E:\\Flood_2013_model\\result\\res-0019.wd' 
file_path_Qx = 'E:\\Flood_2013_model\\result\\res-0019.Qx'    
file_path_Qy = 'E:\\Flood_2013_model\\result\\res-0019.Qy'


wd,wd_params = LoadGrid(file_path_wd)
Qx,Qx_params = LoadGrid(file_path_Qx)
Qy,Qy_params = LoadGrid(file_path_Qy)

    
### Calculate the connectivity number for all the cells
    
cc = Array2D(wd_params['nrows'] ,wd_params['ncols'])
save_cc = np.array([0])


### assign input water cell, change with different input gauge station
cc[70][9] = np.array([1])
 
### find the cell where water start to follow marked as 1
# check cell with water depth > 1 cm and Qx or Qy > 0.01 m^3/s 
index_wd = np.where(wd>0.01 ) 
index_wet_0 =np.array([])
index_wet_1 =np.array([])


for n in range(len(index_wd[1])):
    i = index_wd[1][n]
    j = index_wd[0][n] 
    if abs(Qx[j,i]) >= 0.01 or abs(Qx[j,i+1]) >= 0.01 or abs(Qy[j,i]) >=0.01 or abs(Qy[j+1,i]) >=0.01:
        index_wet_0 = np.append(index_wet_0,j).astype('int64')
        index_wet_1 = np.append(index_wet_1,i).astype('int64')
        
    index_wet = tuple((index_wet_0,index_wet_1))

  
for n in range(len(index_wet[1])):
    i = index_wet[1][n]
    j = index_wet[0][n] 
    if Qx[j,i] <=0.01 and Qx[j,i+1] >=0.01  and Qy[j,i] <=0.01  and Qy[j+1,i] >=0.01:
        cc[j][i] = np.array([1])
            
           
save_cc_next = np.array([1])
        

###start to number the connection, stop where no change in cc array: 
###write the loop for iteration
        
iter = 0 # count number of iteration

while (len(save_cc) == len(save_cc_next)) != True or (save_cc == save_cc_next).all() != True: 
#while iter <50: 
    save_cc = np.array(save_cc_next)
    #save = np.array(cc,dtype=object)
    iter += 1
    for n in range(len(index_wet[1])):
        i = index_wet[1][n]
        j = index_wet[0][n]        
        # check the flow into the cell (i,j)
        #
        if Qx[j,i] > 0.01:
            value1 = cc[j][i-1] + 1   
        else: value1 = cc[j][i]
        if Qx[j,i+1] < -0.01:
            value2 = cc[j][i+1] + 1
        else: value2 = cc[j][i]           
        if Qy[j,i] > 0.01:
            value3 = cc[j-1][i] + 1
        else: value3 = cc[j][i]
        if Qy[j+1,i] < -0.01:
            value4 = cc[j+1][i] + 1
        else: value4 = cc[j][i]
        
        value = np.concatenate((value1,value2,value3,value4))
        value = value[value !=0]
        value_list = list(value)
        value_save = list(dict.fromkeys(value_list)) 
        cc[j][i] = np.sort(np.array(value_save))[::-1]
    
    save_cc_next = list2array(cc)

#%%        
### plot the relative frequency distribution of the flow connection:
### plot pr distribution
connect_count = list2array(cc) 

###count frequency distribution          

plt.ion()

raster = np.array(connect_count[connect_count !=0])
### Plot the probability distribution of the data
def ProDist(raster):
    
    fig, ax = plt.subplots()
    num_bins = 20
    n, bins, patches = ax.hist(raster, num_bins, density=True)
    
    # add a 'best fit' line
    mu = np.mean(raster)  # mean of distribution
    sigma = np.std(raster)  # standard deviation of distribution
    
    y = ((1 / (np.sqrt(2 * np.pi) * sigma)) *
         np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
    ax.plot(bins, y, '--')
    ax.set_xlabel('Distance (pixels)')
    ax.set_ylabel('Probability density')
    title  = '$\mu$ = ' + str(np.round(mu,1)) + ',' + '$\sigma$ = ' + str(np.round(sigma,1))
    ax.set_title(title)
    
    # Tweak spacing to prevent clipping of ylabel
    fig.tight_layout()
    plt.show()
    
    return fig

### Plot probability of length of connection
figure = plt.figure(figsize=(2.49, 2)) 
plt.rc('font', size=10)
value, counts = np.unique(raster, return_counts=True)
prob = counts/sum(counts)

plt.plot(value, prob, label = 't = 200 hrs')
plt.xlabel('Connection length (pixels)')
plt.ylabel('Relative frequency')

max_length = max(raster)
plt.xticks(np.arange(0,max_length+1, 500))

figure.tight_layout()    
plt.legend(loc='upper right')

text1 = 'Mean: ' + str(round(np.mean(raster),1)) 
text2 = 'Max: ' + str(round(np.max(raster),1)) 
plt.text(1400, .0007, text1)
plt.text(1400, .0006, text2)
plt.show()


Writedata(raster, 'Flood2013_t19_totalCF.txt', path = 'E:\\Model_Narran\\Catchment data\\Connectivity_check\\')

#%%
# ##### option to plot the data from the saved file
# ##### open raster
# #### reload from the save data
# file_path = 'E:\\Model_Narran\\Catchment data\\Connectivity_check\\Flood2011_t190_totalCF.txt'
# raster = np.genfromtxt(file_path, unpack=True)
# #
# ###
# figure = plt.figure(figsize=(3, 2.5))  

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif')
# plt.rc('font', size=10)

# value, counts = np.unique(raster, return_counts=True)
# prob = counts/sum(counts)

# plt.plot(value, prob, color = 'Firebrick')
# plt.xlabel('Connection length (pixels)')
# plt.ylabel('Relative frequency')

# max_length = max(raster)
# plt.xticks(np.arange(0,max_length+1, 500))

# figure.tight_layout()  

# idx = np.where(counts == max(counts))
# Peak_length = value[idx[0][0]]

# Mean = 'Mean: ' + str.format('{0:.0f}', np.mean(raster))
# Max = 'Max: ' + str.format('{0:.0f}', np.max(raster))
# Peak = 'Peak: ' + str.format('{0:.0f}', Peak_length)

# plt.text(900, .0023, Mean)
# plt.text(900, .0021, Max)
# plt.text(900, .0019, Peak)

  
# plt.show()















         