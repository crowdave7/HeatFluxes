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

#Evaporation v LandFlux-EVAL rmse
#variable_y = [0.52513046, 1.0515319, 0.30871224, 0.24810952, 0.77596768, 0.46037268, 1.23658859, 0.50148155, 1.02871562, 0.71197554, 0.61166957, 0.62293586, 0.56419109, 0.07051655]

# Precip v Mean (CHIRPS/GPCC) RMSE
#variable_x = [1.05976145, 3.41556894, 1.45201831, 0.63344791, 0.82590173, 0.99572356, 0.87939634, 3.01854602, 2.81287744, 1.42418161, 2.92457451, 1.49827281, 0.31040471, 0.1935775]

#Change in Precipitation Nov - Mar
variable_x = [0.442, 0.497, 1.317, -0.366, 1.416, 0.609, -0.187, -0.618, 1.049, 0.811, 0.321, 0.433, 0.697, 0.132, 0.292]

# Change in Evaporation Nov-Mar
variable_y = [0.024, -0.012, -0.395, -0.575, -0.031, -0.488, -0.6, -0.314, -0.009, -0.381, -0.152, -0.287, -0.025, -0.545, -0.430]

list_of_models = ['bcc-csm1-1-m', 'BNU-ESM', 'CanAM4', 'CNRM-CM5', 'GFDL-HIRAM-C360', 'GISS-E2-R', 'inmcm4', 'IPSL-CM5A-MR', 'MIROC5', 'MRI-AGCM3-2S', 'NorESM1-M', 'AGCM Ensemble', 'Weaker AGCMs', 'Better AGCMs', 'MSWEP/GLEAM']

print len(variable_x)
print len(variable_y)
print len(list_of_models)

cmap_models = plt.get_cmap('rainbow')
colours_models = [cmap_models(i) for i in np.linspace(0, 1, 11)]

cmap_reanalysis = plt.get_cmap('brg')
colours_reanalysis = [cmap_reanalysis(i) for i in np.linspace(0, 1, len(variable_x))]

"""For each model,"""
for i in np.arange(0, len(variable_x), 1):

    """Select the line colour and add one to the line count variable."""

    if list_of_models[i] == "AGCM Ensemble":
        ax1.scatter(variable_x[i], variable_y[i], color='black', marker='o', s=16, label=list_of_models[i])

    elif list_of_models[i] == "Better AGCMs":
        ax1.scatter(variable_x[i], variable_y[i], color='forestgreen', marker='o', s=16, label=list_of_models[i])

    elif list_of_models[i] == "Weaker AGCMs":
        ax1.scatter(variable_x[i], variable_y[i], color='saddlebrown', marker='o', s=16, label=list_of_models[i])

    elif list_of_models[i] == "MSWEP/GLEAM":
        ax1.scatter(variable_x[i], variable_y[i], color=colours_reanalysis[0], marker='o', s=16, label=list_of_models[i])

    else:

        """Add the scatter plot. Select dot colour, marker and label."""
        ax1.scatter(variable_x[i], variable_y[i], color=colours_models[i], marker='o', s=16, label=list_of_models[i])

    """Add a legend for each dot."""
    legend = plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.525), fontsize=7)

plt.ylabel(u'Δ Evaporation between March/November (mm $\mathregular{day^{-1}}$)', fontsize=7)
plt.ylim((-0.7, 0.1))
plt.yticks([-0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1])

plt.xlabel(u'Δ Rainfall between March/November (mm $\mathregular{day^{-1}}$)', fontsize=7)
plt.xlim((-1.0, 2.0))
plt.xticks([-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0])


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

plt.axvline(0, color = 'black', linewidth=0.9)
plt.axhline(0, color = 'black', linewidth=0.9)

#plt.title("Congo Basin (14$^\circ$S - 4$^\circ$N, 18$^\circ$E - 29$^\circ$E)")

"""Save figure."""
fig.savefig("Change_precip_evap_Congo.png", bbox_extra_artists=(legend,), bbox_inches='tight', dpi=600)
