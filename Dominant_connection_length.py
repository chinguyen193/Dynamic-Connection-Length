# -*- coding: utf-8 -*-
"""
Created on Wed Jan  6 17:08:56 2021

@author: tpngu28

Code to calculate the dominant connection length 
Input: 
1.Raster layer of water depth
2.Raster layer of discharge in x and y directions

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from datetime import datetime
import matplotlib.dates as mdates
from matplotlib import rc
from matplotlib import cm
from Load_and_plot_data import LoadGrid, PlotRasters
import os


#### calculate connectivity number
    
def Connect_nb(wd,wd_params,Qx,Qx_params,Qy,Qy_params):

    cc = np.zeros(( wd_params['nrows'] ,wd_params['ncols']))
    save = np.zeros(( wd_params['nrows'] ,wd_params['ncols']))
    ### assign input water cell, change with different input gauge stations
    cc[70,9] = 1   
    ###
    # find the cell where water start to follow marked as 1
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
        i = index_wd[1][n]
        j = index_wd[0][n]        
        if Qx[j,i] <=0.01 and Qx[j,i+1] >=0.01 and Qy[j,i] <=0.01 and Qy[j+1,i] >=0.01:
            cc[j,i] = 1
    
                
    # start to number the connection, stop where no change in cc array:
    iter = 0 # count number of iteration
    while (save == cc).all() != True:
        save = np.array(cc)
        iter += 1
        for n in range(len(index_wet[1])):
            i = index_wet[1][n]
            j = index_wet[0][n]        
            # check the flow into the cell (i,j)
            #
            if Qx[j,i] > 0.01:
                value1 = cc[j,i-1] + 1   
            else: value1 = cc[j,i]
            if Qx[j,i+1] < -0.01:
                value2 = cc[j,i+1] + 1
            else: value2 = cc[j,i]            
            if Qy[j,i] > 0.01:
                value3 = cc[j-1,i] + 1
            else: value3 = cc[j,i]
            if Qy[j+1,i] < -0.01:
                value4 = cc[j+1,i] + 1
            else: value4 = cc[j,i]
            
            # assign connection number depending on the amount of water flow into the cell
            Q_arround = np.array([Qx[j,i],-Qx[j,i+1],Qy[j,i],-Qy[j+1,i]])
            value = np.array([value1, value2, value3, value4])
            cc[j,i] = value[Q_arround.argmax()]

    return cc 


#%%
############# plot figure
#read the water level and discharge file
#change the path to the location of the files
file_path_wd = 'E:\\Flood_2013_model\\result\\res-0198.wd' 
file_path_Qx = 'E:\\Flood_2013_model\\result\\res-0198.Qx'    
file_path_Qy = 'E:\\Flood_2013_model\\result\\res-0198.Qy'


wd,wd_params = LoadGrid(file_path_wd)
Qx,Qx_params = LoadGrid(file_path_Qx)
Qy,Qy_params = LoadGrid(file_path_Qy)
cc = Connect_nb(wd,wd_params,Qx,Qx_params,Qy,Qy_params)



fig = PlotRasters(wd, Legend ='Water level (m)' )
#plt.savefig('Flood_2013_t198_WL.tif', dpi=500, format="tiff")

fig = PlotRasters(cc, Legend = 'Connection length (pixels)')
#plt.savefig('Flood_2013_t198_CN.tif', dpi=500, format="tiff")

###plot pr distribution
figure = plt.figure(figsize=(4.9, 4)) 
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=20)

#for i in range(len(connect_save)):
max_length = np.max(cc)
len_connect = np.arange(50,max_length+1,1)

### plot distribution of length of flow connect
dist = np.zeros(len(len_connect))
prob = np.zeros(len(len_connect))
count = 0 

#value, counts = np.unique(raster, return_counts=True)
#prob = counts/sum(counts)

for k in len_connect:
    dist[count]  =  np.count_nonzero(cc==k)
    count += 1
#
prob = dist/sum(dist) #relative frequency

#plt.plot(len_connect,prob, label = 't = ' + str(i*4*10) + 'hr')
plt.plot(len_connect,prob, label = 't = 4260 h')
plt.xlabel('Connection length (pixels)')
plt.ylabel('Relative frequency')
plt.xticks(np.arange(0,max_length+1, 200))

figure.tight_layout()    
#plt.legend()
plt.show()


#figure.savefig('Flood_2013_t198_RF.tif', dpi=500, format="tiff")