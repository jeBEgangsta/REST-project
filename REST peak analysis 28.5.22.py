# -*- coding: utf-8 -*-
"""
Created on Sat May 28 22:55:11 2022

@author: ftead
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import trapz
from scipy.signal import find_peaks

import seaborn as sns
#%%

#testing fofr calculation of a single nucleus
df1 = pd.read_excel('wt_cells.xls', 'WT 4')


x = list(df1['distance'])
y = list(df1['red'])


#%%
'''
def get_peaks(data):
    peaks_x = []
    for i in range(1, len(data)-1):
        if (data[i] > data[i-1] and data[i] > data[i+1]):
            peaks_x.append(i)
    
    return peaks_x
'''
#%%

#scipy_peaks = list(find_peaks(y, prominence = )[0])



#peaks_x = get_peaks(y)

from peakdetect import peakdetect

peaks = peakdetect(y, lookahead=25) 
higherPeaks = np.array(peaks[0])
lowerPeaks = np.array(peaks[1])

a= np.full(len(y),((peaks[0][0][1]-(peaks[0][0][1] - peaks[1][0][1])/2)))
idx = np.argwhere(np.diff(np.sign(a - y))).flatten()
#%%




#%%


plt.plot(y)
plt.plot(higherPeaks[:,0], higherPeaks[:,1], 'ro')
plt.plot(lowerPeaks[:,0], lowerPeaks[:,1], 'ko')
plt.plot(idx,[0]*len(idx), 'ro')
plt.show()

#%%

#Generating a dunction that will quantify peak thickness for all the cells

#%%


def peak_is_thicc(data):
    results = pd.DataFrame(columns = ['cell','peak_width'])
    for i in data:
        distance = list(data[i]['distance'])#List of all the distance values across the measured cell.
        rest_signal = list(data[i]['red'])#List of the the intensity values of REST across the measured cell.
        peaks = peakdetect(rest_signal, lookahead=20)#Identifies peaks in array
        if len(peaks[0]) == 2:#If there is more then one peak in the signal
           middle_of_peak = np.full(len(rest_signal),((peaks[0][0][1]-(peaks[0][0][1] - peaks[1][0][1])/2)))#Calculates the middle of the peak
           idx = np.argwhere(np.diff(np.sign(middle_of_peak - rest_signal))).flatten()#Identifying the index of intersextion between the signal and the middle of the peak
           peak_width = abs(distance[idx[0]] - distance[idx[1]])
           values = pd.DataFrame([[i,peak_width]], columns = ['cell','peak_width'])
           results = results.append(values)
        else:
            print(i)
            '''
            middle_of_peak = np.full(len(rest_signal),((max(list(peaks[0][j][1] for j in range(len(peaks[0]))))-min(rest_signal))/2))
            idx = np.argwhere(np.diff(np.sign(middle_of_peak - rest_signal))).flatten()
            peak_width = abs(distance[idx[0]] - distance[idx[1]])
            values = pd.DataFrame([[i,peak_width]], columns = ['cell','peak_width'])
            results = results.append(values)
             '''                     
           
    return results

#%%
wt_cells_dict = pd.read_excel('wt_cells.xls', sheet_name=None)
ko_cells_dict = pd.read_excel('ko_cells.xls', sheet_name=None)




results_wt = peak_is_thicc(wt_cells_dict)
results_ko = peak_is_thicc(ko_cells_dict)



results_wt.to_csv('wt_peak_width.csv')
results_ko.to_csv('ko_peak_width.csv')



