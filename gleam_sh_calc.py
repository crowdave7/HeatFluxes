import iris
import iris.analysis
from iris.util import unify_time_units
import numpy as np
import netCDF4
from netCDF4 import num2date, date2num
import h5py
import iris.coords as icoords

#DOWNLOAD DATA

# list_of_yearmonths = []
# for i in np.arange(1984, 2008, 1):
#     for j in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
#         yearmonth = str(i) + j
#         list_of_yearmonths = np.append(list_of_yearmonths, yearmonth)
# print list_of_yearmonths

# for i in list_of_yearmonths:
#     cube_generator = iris.fileformats.netcdf.load_cubes("https://opendap.larc.nasa.gov/opendap/hyrax/SRB/GLW/SRB_REL3.1_LW_MONTHLY_NC/srb_rel3.1_longwave_monthly_"+str(i)+".nc")
#     cube = list(cube_generator)
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/srb_rel3.1_longwave_monthly_"+str(i)+".nc")
#     print "done "+str(i)+""
#
# list_of_yearmonths = []
# for i in np.arange(1984, 2008, 1):
#     for j in ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']:
#         yearmonth = str(i) + j
#         list_of_yearmonths = np.append(list_of_yearmonths, yearmonth)
# print list_of_yearmonths

# for i in list_of_yearmonths:
#     cube_generator = iris.fileformats.netcdf.load_cubes("https://opendap.larc.nasa.gov/opendap/hyrax/SRB/GSW/SRB_REL3.0_SW_MONTHLY_UTC_NC/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc")
#     cube = list(cube_generator)
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc")
#     print "done "+str(i)+""
#

# SAVE SW DOWNWELLING

# for i in list_of_yearmonths:
#     cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc", "surface_downwelling_shortwave_flux")
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/cube/downwelling/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc")
#     print "done "+str(i)+""
# #

# SAVE LW DOWNWELLING

# for i in list_of_yearmonths:
#     cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/srb_rel3.1_longwave_monthly_"+str(i)+".nc", "surface_downwelling_longwave_flux")
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/cube/downwelling/srb_rel3.1_longwave_monthly_"+str(i)+".nc")
#     print "done "+str(i)+""

#SAVE SW UPWELLING

# for i in list_of_yearmonths:
#     cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc", "surface_upwelling_shortwave_flux")
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/cube/upwelling/srb_rel3.0_shortwave_monthly_utc_"+str(i)+".nc")
#     print "done "+str(i)+""
#
# #SAVE LW UPWELLING
#
# for i in list_of_yearmonths:
#     cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/srb_rel3.1_longwave_monthly_"+str(i)+".nc", "surface_upwelling_longwave_flux")
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/cube/upwelling/srb_rel3.1_longwave_monthly_"+str(i)+".nc")
#     print "done "+str(i)+""

# cube1 = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/LW/cube/upwelling/srb_rel3.1_longwave_monthly.nc")
# cube2 = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/SW/cube/upwelling/srb_rel3.0_shortwave_monthly_utc.nc")
# print cube1
# print cube2

# LOAD LWR DOWNWELLING AND SWR DOWNWELLING

# lw_cube_down = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/srb_rel3.1_longwave_monthly_downwelling.nc")
# #print lw_cube_down
# sw_cube_down = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/srb_rel3.0_shortwave_monthly_utc_downwelling.nc")
# #print sw_cube_down
#
# # LOAD LWR UPWELLING AND SWR UPWELLING
#
# lw_cube_up = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/srb_rel3.1_longwave_monthly_upwelling.nc")
# #print lw_cube_up
# sw_cube_up = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/Rad_data/srb_rel3.0_shortwave_monthly_utc_upwelling.nc")
# #print sw_cube_up
#
# # ADD LWRDOWN + SWRDOWN TO GET TOTAL DOWNWELLING
#
# downwelling_cube = iris.analysis.maths.add(lw_cube_down, sw_cube_down)
#
#
# # ADD LWRUP + SWRUP TO GET TOTAL UPWELLING
# upwelling_cube = iris.analysis.maths.add(lw_cube_up, sw_cube_up)
#
# # DOWNWELLING MINUS UPWELLING GIVES NET DOWNWARD SOLAR Radiation
#
# net_downward_rad_cube = iris.analysis.maths.subtract(downwelling_cube, upwelling_cube)
# net_downward_rad_cube.long_name = "Net Downward Solar Radiation"
# print net_downward_rad_cube
# iris.save(net_downward_rad_cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/srb_net_downward_rad.nc")




# CONSTRAIN YEARS IN HFLS CUBE
#
# cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/hfls_gleam.nc")
# time_range = iris.Constraint(time=lambda cell: 1984 <= cell.point.year <= 2007)
# with iris.FUTURE.context(cell_datetime_objects=True):
#     cube = cube.extract(time_range)
#     cube = iris.analysis.maths.multiply(cube, 28)
#     cube.long_name = "Surface Upward Latent Heat Flux"
#     cube.units = "W m-2"
#     iris.save(cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/gleam_ungridded.nc")


#
# # #LOAD NET RAD CUBE
#
# srb_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/srb_net_downward_rad.nc")
# print srb_cube
# #
#
# #LOAD HFLS CUBE
#
# gleam_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/gleam_ungridded.nc")
# print gleam_cube
#
# #REGRID HFLS CUBE ONTO NET RAD CUBE
#
# srb_cube.coord('longitude').guess_bounds()
# srb_cube.coord('latitude').guess_bounds()
#
# gleam_cube.coord('longitude').guess_bounds()
# gleam_cube.coord('latitude').guess_bounds()
#
# #
# regridded_gleam_cube = gleam_cube.regrid(srb_cube, iris.analysis.AreaWeighted())
# print regridded_gleam_cube
# iris.save(regridded_gleam_cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/hfls_gleam_regridded1.nc")
#
# gleam_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/hfls_gleam_regridded1.nc")
# print gleam_cube
# REGRID HFLS CUBE ONTO REANALYSIS CUBE


#
# # # LOAD REGRIDDED HFLS CUBE AND NET RAD CUBE
# #
# gleam_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/hfls_gleam_regridded.nc")
# print gleam_cube.coord('time')
# srb_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/srb_net_downward_rad.nc")
# print srb_cube.coord('time')
#
# #REMOVE TIME COORD IN HFLS CUBE AND REPLACE WITH COORD IN NET RAD CUBE
#
# srb_coord = srb_cube.coord('time')
# gleam_cube.remove_coord('time')
# gleam_cube.add_dim_coord(srb_coord, 0)
#
# # SAVE NEW HFLS CUBE
#
# iris.save(gleam_cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/hfls_gleam_regridded1.nc")
#
# # RELOAD HFLS AND SRB cubes

# srb_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/srb_net_downward_rad.nc")
#
# gleam_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/hfls_gleam_regridded.nc")
# print gleam_cube.coord('time')
# print srb_cube.coord('time')
#
# #HFSS = NET DOWNWARD RAD - HFLS
#
# sensible_cube = iris.analysis.maths.subtract(srb_cube, gleam_cube)
#
# sensible_cube.long_name = "Surface Upward Sensible Heat Flux"
# sensible_cube.units = "W m-2"
# #
# # SAVE HFSS CUBE (180x360 GRID)
#
# iris.save(sensible_cube, "/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/hfss_gleam_regridded.nc")

# LOAD SENSIBLE CUBE
# LOAD LATENT CUBE

hfss_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/hfss_gleam.nc")
hfls_cube = iris.load_cube("/ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/hfls_gleam.nc")

print hfss_cube
print hfls_cube


#gleam_cube.coord('time').points = gleam_times

#print gleam_cube.coord('time')

#gleam_cube.remove_coord('time')

#gleam_cube.add_dim_coord(icoords.DimCoord(gleam_times), 0)
