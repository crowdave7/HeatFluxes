#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

list_of_variables = ['tran', 'pr', 'mrsor', 'lai', 'par', 'vpd', 'tas', 'hurs', 'mf', 'huss', 'ws', 'co2']

data = np.zeros((len(list_of_variables), 13))

# data[0] = [1.94819953,	2.1493272,	2.19513115,	2.35791256,	1.500080955,	0.638841645,	0.36752354,	0.35287752,	0.49867262,	0.909735635,	1.40357493,	1.756619575,	1.94819953]
#
# data[1] = [6.649,	6.789,	6.959,	3.408,	0.324,	0.033,	0.013,	0.234,	1.103,	3.155,	6.054,	7.532,	6.649]
#
# data[2] = [0.37687988,	0.38976528,	0.39666251,	0.383457135,	0.32644878,	0.27914338,	0.249866185, 0.229208885,	0.217630545, 0.229112645,	0.283159155, 0.345875725,	0.37687988]
#
# data[3] = [2.8497,	2.8715,	2.6876,	2.2882,	1.6765,	1.1042,	0.81838, 0.82693, 1.1209, 1.7111,	2.3708,	2.712,	2.8497]
#
# data[4] = [97.09134293,	100.8722763,	102.4113808,	105.9181099,	110.1431046,	105.928997,	107.5653801,	110.5575371,	110.540554,	106.0336571,	100.408493,	96.84759522,	97.09134293]
#
# data[5] = [0.901694512,	0.974296735,	0.976730806,	1.218595656,	1.760236229, 1.821631004,	2.117426873, 3.013363005,	3.660581398,	2.340649705,	1.304926218,	0.989352221,	0.901694512]
#
# data[6] = [295.2956122,	295.5202553,	295.7871262,	295.7379782,	294.586,	292.5310334,	292.6782578,	295.9414506,	298.8023808,	298.0068696,	296.414635,	295.6360789,	295.2956122]
#
# data[7] = [80.04609628,	78.73881162,	79.03614123,	73.76501465,	59.28236249,	52.02374313,	44.75315421,	35.93940082,	34.7000543,	56.17579598,	73.05918305,	78.56514345,	80.04609628]
#
# data[8] = [213.1213644,	139.3968614,	41.704169,	951.2744365,	1395.929773,	1564.732088,	1589.877798,	1229.07752,	535.3658373,	73.22984423,	29.01345546,	86.89787363,	213.1213644]
#
# data[9] = [14.72767376,	14.66866434,	14.981791,	13.86026895,	10.21557909,	7.74082118,	6.64392124,	6.5558969,	7.62846059,	11.78094482,	14.21370618,	14.71785657,	14.72767376]
#
# data[10] = [1.55324495,	1.48713599,	0.92761431,	3.11453786,	4.45818525,	5.47661756,	6.11650264,	5.64280842,	4.21162664,	1.44924184,	0.84046844,	1.05042139,	1.55324495]
#
# data[11] = [373.31,	374.53,	374.43,	376.1,	378.01,	379.39,	380.01,	380.59,	379.81,	379.31,	378.92,	379.02,	373.31]

data[0] = [2.32811036,	2.171634235,	2.273459215,	2.290101945,	2.331474745,	2.31643507,	2.110419565,	2.130110465,	2.29835421,	2.212871055,	2.21672765,	2.362752095,	2.32811036]

data[1] = [1.948,	2.699,	4.518,	5.850,	5.459,	4.991,	4.755,	5.790,	6.304,	6.510,	5.063,	2.463,	1.948]

data[2] = [0.3671788,	0.34454368,	0.346626395,	0.366943755,	0.380901795,	0.38381721,	0.384104035,	0.38927002,	0.396929795,	0.40314178,	0.403800205,	0.392346755,	0.3671788]

data[3] = [5.7964,	5.9926,	6.136,	6.2034,	6.1927,	6.0956,	6.0234,	6.0442,	6.0887,	6.074,	5.9588,	5.7888,	5.7964]

data[4] = [107.9403877,	111.5307655,	108.023243,	101.7953491,	95.89775467,	88.23717118, 86.74423218,	90.03637696,	95.88925553,	91.51346588,	92.05929947,	99.94840241,	107.9403877]

data[5] = [1.630956957,	1.710325402,	1.305728455,	1.017991434,	0.944497447,	0.958558339,	0.958514993,	0.919193214,	0.922746196,	0.863016495,	0.8459612,	1.162627042,	1.630956957]

data[6] = [296.2533237,	296.9013537,	296.4315649,	295.7212756,	295.265293,	294.7745404,	294.534256,	294.5059092,	294.6107898,	294.5152984,	294.744703,	295.4823427,	296.2533237]

data[7] = [65.99078243,	65.73322564,	73.07079005,	78.06105789,	79.05937804,	78.08647416,	77.7553567,	78.6299846,	78.68787812,	79.94782126,	80.62441996,	74.56915704,	65.99078243]

data[8] = [6.176143637,	55.69631143,	287.8593142,	165.4332641,	169.4372477,	234.7378433,	144.7581508,	141.0313827,	220.3906831, 106.3807789,	27.22164086,	19.71080077,	6.176143637]

data[9] = [12.80215223,	13.23524722,	14.37356504,	14.82146827,	14.59617497,	13.99387167,	13.74366807,	13.90460846,	14.0053093,	14.15374009,	14.48119423,	13.93771504,	12.80215223]

data[10] = [1.58371287,	1.84880766,	2.50276831,	2.47619024,	1.98977067,	1.94899937,	1.73172407,	1.67464755,	1.70770312,	1.35808496,	1.20177748,	1.29254905,	1.58371287]

data[11] = [380.498,	381.088,	381.558,	379.754,	378.138,	376.874,	376.23,	376.446,	376.404,	376.33,	377.954,	380.108,	380.498]


subplot_columns = 4
subplot_rows = 3

if subplot_columns == 3 and subplot_rows == 3:
    fig = plt.figure(figsize=(10,10))

if subplot_columns == 4 and subplot_rows == 2:
    fig = plt.figure(figsize=(8,6))

if subplot_columns == 4 and subplot_rows == 3:
    fig = plt.figure(figsize=(10,8))

if subplot_columns == 4 and subplot_rows == 4:
    fig = plt.figure(figsize=(10, 10))

if subplot_columns == 4 and subplot_rows == 5:
    fig = plt.figure(figsize=(10, 12))

if subplot_columns == 1 and subplot_rows == 2:
    fig = plt.figure(figsize=(10, 12))

variable_number = 0

for i in range(len(list_of_variables)):

    ax = plt.subplot(subplot_rows, subplot_columns, variable_number+1)
    ax.set_xlim([0, 12])
    ax.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on')
    ax.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on')
    ax.patch.set_visible(False)

    objects = ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', 'J')
    x_pos = np.arange(len(objects))
    plt.xticks(x_pos, objects, fontsize=7)

    if list_of_variables[i] == 'tran':
        lower_y_lim = 0.0
        upper_y_lim = 3.0
        y_tick_interval = 0.5
        plt.ylabel('Transpiration (mm $\mathregular{day^{-1}}$)', fontsize=7)
        plt.title('a) AGCM', fontsize=7)
        line_colour = 'dodgerblue'

    if list_of_variables[i] == 'pr':
        lower_y_lim = 0.0
        upper_y_lim = 10.0
        y_tick_interval = 1.0
        plt.ylabel('Precipitation (mm $\mathregular{day^{-1}}$)', fontsize=7)
        plt.title('b) GPCC', fontsize=7)
        line_colour = 'darkblue'

    if list_of_variables[i] == 'mrsor':
        lower_y_lim = 0.18
        upper_y_lim = 0.42
        y_tick_interval = 0.06
        plt.ylabel('Root-Zone Soil Moisture ($\mathregular{m^{3}}$ $\mathregular{m^{-3}}$)', fontsize=7)
        plt.title('c) GLEAM', fontsize=7)
        line_colour = 'saddlebrown'

    if list_of_variables[i] == 'lai':
        lower_y_lim = 0.0
        upper_y_lim = 7.0
        y_tick_interval = 1.0
        plt.ylabel('Leaf Area Index', fontsize=7)
        plt.title('d) MODIS', fontsize=7)
        line_colour = 'darkgreen'

    if list_of_variables[i] == 'par':
        lower_y_lim = 85.0
        upper_y_lim = 115.0
        y_tick_interval = 5.0
        plt.ylabel(' PAR (W $\mathregular{m^{-2}}$)', fontsize=7)
        plt.title('e) GEWEX', fontsize=7)
        line_colour = 'firebrick'

    if list_of_variables[i] == 'vpd':
        lower_y_lim = 0.0
        upper_y_lim = 4.0
        y_tick_interval = 0.5
        plt.ylabel('Vapour Pressure Deficit (kPa)', fontsize=7)
        plt.title('f) MERRA-2', fontsize=7)
        line_colour = 'steelblue'

    if list_of_variables[i] == 'tas':
        lower_y_lim = 275.0
        upper_y_lim = 310.0
        y_tick_interval = 5.0
        plt.ylabel('Surface Air Temperature (K)', fontsize=7)
        plt.title('g) MERRA-2', fontsize=7)
        line_colour = 'darkorange'

    if list_of_variables[i] == 'hurs':
        lower_y_lim = 30.0
        upper_y_lim = 90.0
        y_tick_interval = 10.0
        plt.ylabel('Relative Humidity (%)', fontsize=7)
        plt.title('h) MERRA-2', fontsize=7)
        line_colour = 'olive'

    if list_of_variables[i] == 'mf':
        lower_y_lim = 0
        upper_y_lim = 1800
        y_tick_interval = 200
        plt.ylabel('Moisture Flux (g kg $\mathregular{^{-1}}$ m s $\mathregular{^{-1}}$)', fontsize=7)
        plt.title('i) MERRA-2', fontsize=7)
        line_colour = 'blueviolet'

    if list_of_variables[i] == 'huss':
        lower_y_lim = 4.0
        upper_y_lim = 16.0
        y_tick_interval = 2.0
        plt.ylabel('Specific Humidity (g kg $\mathregular{^{-1}}$)', fontsize=7)
        plt.title('j) MERRA-2', fontsize=7)
        line_colour = 'orange'

    if list_of_variables[i] == 'ws':
        lower_y_lim = 0.0
        upper_y_lim = 7.0
        y_tick_interval = 1.0
        plt.ylabel('Wind Speed (m s $\mathregular{^{-1}}$)', fontsize=7)
        plt.title('k) MERRA-2', fontsize=7)
        line_colour = 'orchid'

    if list_of_variables[i] == 'co2':
        lower_y_lim = 360
        upper_y_lim = 390
        y_tick_interval = 3.0
        plt.ylabel('CO$_{2}$ (ppm)', fontsize=7)
        plt.title('l) ESA GHG-CCI', fontsize=7)
        line_colour = 'lightseagreen'

    plt.ylim(lower_y_lim, upper_y_lim)
    plt.yticks(np.arange(lower_y_lim, upper_y_lim+y_tick_interval, y_tick_interval), fontsize=7)

    line = ax.plot(x_pos, data[i], zorder=2, linestyle='-', linewidth=2.0, color=line_colour)

    variable_number +=1

fig.tight_layout()

if subplot_columns == 4 and subplot_rows == 2:

    #plt.subplots_adjust(left=0.05, right=0.95, wspace=0.8)
    plt.subplots_adjust(left=0.05, right=0.95, wspace=0.6, hspace=0.2)

if subplot_columns == 4 and subplot_rows == 3:

    plt.subplots_adjust(left=0.05, right=0.95, wspace=0.6, hspace=0.2)

if subplot_columns == 3 and subplot_rows == 3:

    plt.subplots_adjust(left=0.05, right=0.95, wspace=0.6, hspace=0.2)

if subplot_columns == 4 and subplot_rows == 4:

    plt.subplots_adjust(left=0.05, right=0.95, wspace=0.6, hspace=0.2)

if subplot_columns == 4 and subplot_rows == 5:

    plt.subplots_adjust(left=0.05, right=0.95, wspace=0.6, hspace=0.2)

if subplot_columns == 1 and subplot_rows == 2:

    plt.subplots_adjust(left=0.05, right=0.95, wspace=0.6, hspace=0.2)

print "saving final figure"
fig.savefig("test_multiple_models.png", bbox_inches='tight')
plt.close()
print "plot done"
