# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 10:28:54 2021

@author: ftead
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import trapz

import seaborn as sns
#%%

#testing fofr calculation of a single nucleus
df1 = pd.read_excel('ko_cells.xls', 'KO 21')

#normalization of every channle values to its maximum vlaue. puts all channels on the same scale.
df1['normal_red'] = df1['red']/df1['red'].max()
df1['normal_green'] = df1['green']/df1['green'].max()
df1['normal_blue'] = df1['blue']/df1['blue'].max()
df1['normal_distance'] = df1['distance']/df1['distance'].max()

#dx means the distance interval along the x axis
dx = df1.at[1,'distance'] - df1.at[0,'distance']

#clculation of are under the curve
red_area = trapz(list(df1['normal_red']),x = df1['normal_distance'])
green_area = trapz(list(df1['normal_green']),x = df1['distance'])
blue_area = trapz(list(df1['normal_blue']),dx = dx)

#value represents the measure of similarity of REST between LAMINB and DAPI.
#the higher the difference the less the sigle resembles the standard
#the lower the difference the curve of REST looks more like the standard
# The ratio between the diffeerences is estimating the distribution of REST."
# The higher the value (red to blue - high, red to green - low) the more closer REST is to the LAMINA
# The lower the value (red to blue - low, red to green - high) the more REST is uniform.
dist_value = (red_area - blue_area)/(red_area - green_area)

#%%

#making the calculation for n cells, saved in xls file. each sheet is one cells.




def rest_dist_calculation(cells):
    results = pd.DataFrame(columns = ['cell','dist_value'])
    for i in cells:
        cells[i]['normal_red'] = cells[i]['red']/cells[i]['red'].max()
        cells[i]['normal_green'] = cells[i]['green']/cells[i]['green'].max()
        cells[i]['normal_blue'] = cells[i]['blue']/cells[i]['blue'].max()
        cells[i]['normal_distance'] = cells[i]['distance']/cells[i]['distance'].max()
        dx = cells[i].at[1,'normal_distance'] - cells[i].at[0,'normal_distance']
        red_area = trapz(list(cells[i]['normal_red']),dx = dx)
        green_area = trapz(list(cells[i]['normal_green']),dx = dx)
        blue_area = trapz(list(cells[i]['normal_blue']),dx = dx)
        dist_value = abs((red_area - blue_area)/(red_area - green_area))
        values = pd.DataFrame([[i,dist_value]], columns =  ['cell','dist_value'])
        results = results.append(values)
    return results


#%%
wt_cells_dict = pd.read_excel('wt_cells.xls', sheet_name=None)
ko_cells_dict = pd.read_excel('ko_cells.xls', sheet_name=None)

#%%

results_WT = rest_dist_calculation(wt_cells_dict)
results_KO = rest_dist_calculation(ko_cells_dict)
results_WT.to_csv('distribution_ratio_WT.csv')
results_KO.to_csv('distribution_ratio_KO.csv')


#%%

def only_area_dif(cells):
    results = pd.DataFrame(columns = ['cell','red-green','red-blue','green-blue'])
    for i in cells:
        cells[i]['normal_red'] = cells[i]['red']/cells[i]['red'].max()
        cells[i]['normal_green'] = cells[i]['green']/cells[i]['green'].max()
        cells[i]['normal_blue'] = cells[i]['blue']/cells[i]['blue'].max()
        cells[i]['normal_distance'] = cells[i]['distance']/cells[i]['distance'].max()
        dx = cells[i].at[1,'normal_distance'] - cells[i].at[0,'normal_distance']
        red_area = trapz(list(cells[i]['normal_red']),dx = dx)
        green_area = trapz(list(cells[i]['normal_green']),dx = dx)
        blue_area = trapz(list(cells[i]['normal_blue']),dx = dx)
        red_green = abs(red_area - green_area)
        red_blue = abs(red_area - blue_area)
        green_blue = green_area - blue_area
        values = pd.DataFrame([[i,red_green,red_blue,green_blue]], columns = ['cell','red-green','red-blue','green-blue'])
        results = results.append(values)
    return results

#%%
results_WT_raw = only_area_dif(wt_cells_dict)
results_KO_raw = only_area_dif(ko_cells_dict)

#%%
def where_is_REST(row):
    if row['red-blue'] > row['red-green']:
        return 'Laminar REST'
    else:
        return 'Uniform REST'

results_WT_raw['group'] = results_WT_raw.apply(where_is_REST, axis=1)
results_KO_raw['group'] = results_KO_raw.apply(where_is_REST, axis=1)


results_WT_raw.to_csv('results_wt_raw.csv')
results_KO_raw.to_csv('results_ko_raw.csv')

#%%

fig, ax = plt.subplots(figsize = (5,5))
plt.axis('equal')
plt.xlim(-0.01,0.4)
plt.ylim(-0.01,0.4)
p1 = sns.scatterplot(data = results_KO_raw, x = 'red-green', y = 'red-blue', ax = ax, hue = 'group')

