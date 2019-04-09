#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np

list_of_variables = ['tran', 'lai', 'mrsor', 'par', 'vpd', 'tas', 'mf', 'ws', 'co2']

data = np.zeros((len(list_of_variables), 13))

data[0] = [1.94819953,	2.1493272,	2.19513115,	2.35791256,	1.500080955,	0.638841645,	0.36752354,	0.35287752,	0.49867262,	0.909735635,	1.40357493,	1.756619575,	1.94819953]

data[1] = [2.8497,	2.8715,	2.6876,	2.2882,	1.6765,	1.1042,	0.81838, 0.82693, 1.1209, 1.7111,	2.3708,	2.712,	2.8497]

data[2] = [0.37687988,	0.38976528,	0.39666251,	0.383457135,	0.32644878,	0.27914338,	0.249866185, 0.229208885,	0.217630545, 0.229112645,	0.283159155, 0.345875725,	0.37687988]

data[3] = [97.09134293,	100.8722763,	102.4113808,	105.9181099,	110.1431046,	105.928997,	107.5653801,	110.5575371,	110.540554,	106.0336571,	100.408493,	96.84759522,	97.09134293]

data[4] = [0.901694512,	0.974296735,	0.976730806,	1.218595656,	1.760236229, 1.821631004,	2.117426873, 3.013363005,	3.660581398,	2.340649705,	1.304926218,	0.989352221,	0.901694512]

data[5] = [295.2956122,	295.5202553,	295.7871262,	295.7379782,	294.586,	292.5310334,	292.6782578,	295.9414506,	298.8023808,	298.0068696,	296.414635,	295.6360789,	295.2956122]

data[6] = [213.1213644,	139.3968614,	41.704169,	951.2744365,	1395.929773,	1564.732088,	1589.877798,	1229.07752,	535.3658373,	73.22984423,	29.01345546,	86.89787363,	213.1213644]

data[7] = [1.55324495,	1.48713599,	0.92761431,	3.11453786,	4.45818525,	5.47661756,	6.11650264,	5.64280842,	4.21162664,	1.44924184,	0.84046844,	1.05042139,	1.55324495]

data[8] = [373.31,	374.53,	374.43,	376.1,	378.01,	379.39,	380.01,	380.59,	379.81,	379.31,	378.92,	379.02,	373.31]


subplot_columns = 3
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
        upper_y_lim = 2.5
        y_tick_interval = 0.5
        plt.title('a) AGCM Transpiration (mm $\mathregular{day^{-1}}$)', fontsize=7)
        line_colour = 'dodgerblue'

    if list_of_variables[i] == 'lai':
        lower_y_lim = 0.0
        upper_y_lim = 3.5
        y_tick_interval = 0.5
        plt.title('b) MODIS Leaf Area Index', fontsize=7)
        line_colour = 'darkgreen'

    if list_of_variables[i] == 'mrsor':
        lower_y_lim = 0.18
        upper_y_lim = 0.42
        y_tick_interval = 0.06
        plt.title('c) GLEAM Root-Zone Soil Moisture ($\mathregular{m^{3}}$ $\mathregular{m^{-3}}$)', fontsize=7)
        line_colour = 'saddlebrown'

    if list_of_variables[i] == 'par':
        lower_y_lim = 90.0
        upper_y_lim = 115.0
        y_tick_interval = 5.0
        plt.title('d) GEWEX PAR (W $\mathregular{m^{-2}}$)', fontsize=7)
        line_colour = 'firebrick'

    if list_of_variables[i] == 'vpd':
        lower_y_lim = 0.0
        upper_y_lim = 4.0
        y_tick_interval = 0.5
        plt.title('e) MERRA-2 Vapour Pressure Deficit (kPa)', fontsize=7)
        line_colour = 'steelblue'

    if list_of_variables[i] == 'tas':
        lower_y_lim = 275.0
        upper_y_lim = 310.0
        y_tick_interval = 5.0
        plt.title('f) MERRA-2 Surface Air Temperature (K)', fontsize=7)
        line_colour = 'darkorange'

    if list_of_variables[i] == 'mf':
        lower_y_lim = 0
        upper_y_lim = 1800
        y_tick_interval = 200
        plt.title('g) MERRA-2 Moisture Flux (g kg $\mathregular{^{-1}}$ m s $\mathregular{^{-1}}$)', fontsize=7)
        line_colour = 'blueviolet'

    if list_of_variables[i] == 'ws':
        lower_y_lim = 0.0
        upper_y_lim = 7.0
        y_tick_interval = 1.0
        plt.title('h) MERRA-2 Wind Speed (m s $\mathregular{^{-1}}$)', fontsize=7)
        line_colour = 'peru'

    if list_of_variables[i] == 'co2':
        lower_y_lim = 360
        upper_y_lim = 390
        y_tick_interval = 3.0
        plt.title('i) ESA GHG-CCI CO$_{2}$ (ppm)', fontsize=7)
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
