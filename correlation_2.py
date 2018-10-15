#!/usr/bin/env python

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

ax1.tick_params(axis='x', direction='in', which='both', labelbottom='on', labeltop='off', bottom='on', top='on')
ax1.tick_params(axis='y', direction='in', which='both', labelleft='on', labelright='off', left='on', right='on')

#Evaporation v LandFlux-EVAL rmse
variable_x = [0.52513046, 1.0515319, 0.30871224, 0.24810952, 0.77596768, 0.46037268, 1.23658859, 0.50148155, 1.02871562, 0.71197554, 0.61166957, 0.62293586, 0.56419109, 0.07051655]

# Precip v Mean (CHIRPS/GPCC) RMSE
variable_y = [1.05976145, 3.41556894, 1.45201831, 0.63344791, 0.82590173, 0.99572356, 0.87939634, 3.01854602, 2.81287744, 1.42418161, 2.92457451, 1.49827281, 0.31040471, 0.1935775]

list_of_models = ['bcc-csm1-1-m', 'BNU-ESM', 'CanAM4', 'CNRM-CM5', 'GFDL-HIRAM-C360', 'GISS-E2-R', 'inmcm4', 'IPSL-CM5A-MR', 'MIROC5', 'MRI-AGCM3-2S', 'NorESM1-M', 'Ensemble', 'ERA-5', 'MSWEP/GLEAM']

cmap = plt.get_cmap('rainbow')
colours = [cmap(i) for i in np.linspace(0, 1, len(variable_x))]

"""For each model,"""
for i in np.arange(0, len(variable_x), 1):

    """Select the line colour and add one to the line count variable."""
    dot_colour = colours[i]

    """Add the scatter plot. Select dot colour, marker and label."""
    ax1.scatter(variable_x[i], variable_y[i], color=dot_colour, marker='o', label=list_of_models[i])

    """Add a legend for each dot."""
    legend = plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.525), fontsize=9)

plt.xlabel('Annual Evaporation RMSE v LandFlux-EVAL')
plt.xlim((0, 1.4))

plt.ylabel('Annual Precipitation RMSE v CHIRPS/GPCC')
plt.ylim((0, 4))

# x = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]
#
# y = [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]
#
# ax1.plot(x, y, color = 'k', linestyle = '--', linewidth=1)

plt.title("Congo Basin (14$^\circ$S - 4$^\circ$N, 18$^\circ$E - 29$^\circ$E)")

"""Save figure."""
fig.savefig("RMSE_precip_evap_Congo.png", bbox_extra_artists=(legend,), bbox_inches='tight', dpi=600)
