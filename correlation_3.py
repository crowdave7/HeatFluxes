#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Import necessary modules for this code."""

import iris
import iris.analysis
import iris.analysis.cartography
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mpl_toolkits.basemap import Basemap, maskoceans, interp
import numpy as np
from scipy.stats.stats import pearsonr
matplotlib.use('Agg')

fig = plt.figure()
ax1 = fig.add_subplot(111)

ax1.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on', labelsize=7)
ax1.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on', labelsize=7)

# #Evaporation v LandFlux-EVAL rmse
# variable_y = [0.52513046, 1.0515319, 0.30871224, 0.24810952, 0.77596768, 0.46037268, 1.23658859, 0.50148155, 1.02871562, 0.71197554, 0.61166957, 0.62293586, 0.56419109, 0.07051655]
#
# # Precip v Mean (CHIRPS/GPCC) RMSE
# variable_x = [1.05976145, 3.41556894, 1.45201831, 0.63344791, 0.82590173, 0.99572356, 0.87939634, 3.01854602, 2.81287744, 1.42418161, 2.92457451, 1.49827281, 0.31040471, 0.1935775]

#
# variable_x = [0.442,	0.497,	1.317,	-0.366,	1.416,	0.609,	-0.187,	-0.618,	1.049,	0.811,	0.321,	0.433,	0.292,	0.345]
#
# variable_y = [0.024,	-0.012,	-0.395,	-0.575,	-0.031,	-0.488,	-0.600,	-0.314,	-0.009,	-0.381,	-0.152,	-0.287,	-0.249,	-0.430]

# NR - Tran

#NCB Mar-Jul
# variable_x = [-6.558, 3.315, -3.642, 	-15.733, -17.975,	-44.677, -31.426,	-20.181,	-44.709,	-11.759,	-12.425,	-18.131,	-35.817]
#
# variable_y = [-1.104,	0.505,	0.895,	-0.419,	-0.392,	-0.888,	-0.994,	0.085,	-0.136,	0.423,	0.273,	-0.157,	-0.471]

#NCB Jul-Nov
# variable_x = [4.508,	-12.613,	-21.339,	14.625,	12.469,	33.114,	13.702,	14.248,	39.240,	-8.891,	-0.255,	7.158,	15.308]
#
# variable_y = [0.950,	-0.667,	-1.210,	0.438,	0.547,	0.452,	0.732,	-0.414,	-0.087,	-0.833,	-0.690,	-0.095,	0.110]

#SCB Mar-Jul
# variable_x = [-43.837,	-15.67,	-24.373,	-49.193,	-51.053,	-54.392,	-43.647,	-65.56,	-54.841,	-47.903,	-37.75,	-44.429,	-42.918]
#
# variable_y = [-1.796,	0.023,	0.813,	-2.127,	-1.854,	-1.182,	-1.882,	-1.638,	-0.505,	-0.745,	-0.727,	-1.082,	-1.743]

#SCB Jul-Nov
# variable_x = [45.552,	21.262,	21.591,	41.241,	47.585,	50.282,	50.154,	73.402, 9.273,	49.748,	51.452,	47.183,	33.576]
#
# variable_y = [1.402,	-0.001,	-1.205,	0.918,	1.454,	0.376,	1.219,	0.284,	0.455,	0.131,	0.657,	0.538,	1.201]


# SM - TRAN

#NCB Mar-Jul
# variable_x = [-7.897,	-7.358,	-9.863,	-2.801,	-1.782,	-14.204,	-13.865,	-7.412,	-7.875,	-5.482,	-12.012,	-8.867,	-4.214]
# variable_y = [-1.104,	0.505,	0.895,	-0.419,	-0.392,	-0.888,	-0.994,	0.085,	-0.136,	0.423,	0.273,	-0.157,	-0.471]

#NCB Jul-Nov
# variable_x = [7.967,	7.759,	13.452,	2.201,	6.141,	13.192,	17.786,	7.329,	10.070,	8.248,	12.456,	10.734,	6.235]
# variable_y = [0.950,	-0.667,	-1.210,	0.438,	0.547,	0.452,	0.732,	-0.414,	-0.087,	-0.833,	-0.690,	-0.095,	0.110]

#SCB Mar-Jul
# variable_x = [-9.656,	-15.915, -21.174,	-8.813,	-9.531,	-11.602,	-21.094,	-12.589,	-13.340,	-19.021,	-16.022,	-14.647,	-19.003]
# variable_y = [-1.796,	0.023,	0.813,	-2.127,	-1.854,	-1.182,	-1.882,	-1.638,	-0.505,	-0.745,	-0.727,	-1.082,	-1.743]

#SCB Jul-Nov
variable_x = [8.458,	15.478,	21.381,	5.518,	9.795,	9.515,	10.597,	12.140,	10.512,	13.095,	12.108,	11.317,	9.881]
variable_y = [1.402,	-0.001,	-1.205,	0.918,	1.454,	0.376,	1.219,	0.284,	0.455,	0.131,	0.657,	0.538,	1.201]


list_of_models = ['bcc-csm1-1-m', 'BNU-ESM', 'CanAM4', 'CNRM-CM5', 'GFDL-HIRAM-C360', 'GISS-E2-R', 'inmcm4', 'IPSL-CM5A-MR', 'MIROC5', 'MRI-AGCM3-2S', 'NorESM1-M', 'Ensemble', 'GLEAM']

cmap = plt.get_cmap('rainbow')
colours = [cmap(i) for i in np.linspace(0, 1, len(variable_x))]

"""For each model,"""
for i in np.arange(0, len(variable_x), 1):

    """Select the line colour and add one to the line count variable."""
    dot_colour = colours[i]

    """Add the scatter plot. Select dot colour, marker and label."""
    ax1.scatter(variable_x[i], variable_y[i], color=dot_colour, marker='o', s=16, label=list_of_models[i])

    """Add a legend for each dot."""
    legend = plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.525), fontsize=7)

plt.ylabel(u'Δ Transpiration Jul - Nov (mm $\mathregular{day^{-1}}$)', fontsize=7)
#plt.ylim((-0.7, 0.1))
#plt.yticks([-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1])

plt.ylim((-2.5, 2.5))
#label = ax1.set_xlabel

plt.xlabel(u'Δ Upper Layer Soil Moisture Jul-Nov (mm $\mathregular{day^{-1}}$)', fontsize=7)
#ax1.xaxis.set_label_coords(0.7, 0.96)
plt.xlim((-25, 25))
plt.xticks([-25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25])

pearson = pearsonr(variable_x, variable_y)
pearsoncoeff = pearson[0]
pearsoncoeff_round = round(pearsoncoeff,3)

plt.title ("SCB, Jul-Nov, r = "+str(pearsoncoeff_round))

plt.plot(variable_x, np.poly1d(np.polyfit(variable_x, variable_y, 1))(variable_x))

# x = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]
#
# y = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]
#
# ax1.plot(x, y, color = 'k', linestyle = '--', linewidth=1)

# Move left y-axis and bottim x-axis to centre, passing through (0,0)
#ax1.spines['left'].set_position('zero')
#ax1.spines['bottom'].set_position('zero')

# Eliminate upper and right axes
#ax1.spines['right'].set_color('none')
#ax1.spines['top'].set_color('none')

# Show ticks in the left and lower axes only
ax1.xaxis.set_ticks_position('bottom')
ax1.yaxis.set_ticks_position('left')

plt.axvline(0, color = 'black')
plt.axhline(0, color = 'black')

#plt.title("Congo Basin (14$^\circ$S - 4$^\circ$N, 18$^\circ$E - 29$^\circ$E)")

"""Save figure."""
fig.savefig("corr_SM_tran.png", bbox_extra_artists=(legend,), bbox_inches='tight', dpi=600)
