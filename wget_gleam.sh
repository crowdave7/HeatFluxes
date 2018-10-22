# #!/usr/bin/env python
# from ecmwfapi import ECMWFDataServer
# server = ECMWFDataServer()
# server.retrieve({
#     "class": "ei",
#     "dataset": "interim_land",
#     "date": "1979-01-01/to/1979-01-31",
#     "expver": "2",
#     'format': "netcdf",
#     "levtype": "sfc",
#     "param": "39.128",
#     "stream": "oper",
#     "time": "00:00:00/06:00:00/12:00:00/18:00:00",
#     "type": "an",
#     "target": "interim_land_test.grib",
# })


#!/bin/bash
#for i in {1980..2008..1}; do
    #cdo shifttime,-1days /ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/MSWEP_data/${i}_new.nc /ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/MSWEP_data/${i}_shiftedtime.nc
    #cdo monmean /ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/MSWEP_data/${i}_shiftedtime.nc /ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/MSWEP_data/${i}_monthly.nc

    #cdo monmean /ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/transpiration/Et_${i}_GLEAM_v3.2a.nc /ouce-home/students/kebl4396/Paper1/Paper1ReanalysisFiles/GLEAM_data/transpiration/Et_${i}_GLEAM_v3.2a_monthly.nc


#done


# for i in {1982..2013..1}; do
#
#   ncpdq -a dim2,dim0 /ouce-home/projects/land_surface_climate/modis/Global_Monthly_ET_${i}.nc /ouce-home/projects/land_surface_climate/modis/Global_Monthly_ET_${i}_swapped.nc
#   echo done
#
# done

# for i in {1982..2013..1}; do
#
#   cdo settunits,days -settaxis,${i}-01-15,00:00,1month /ouce-home/projects/land_surface_climate/modis/Global_Monthly_ET_${i}_swapped.nc /ouce-home/projects/land_surface_climate/modis/Global_Monthly_ET_${i}_swapped_time.nc
#
# done

#for i in bcc-csm1-1-m BNU-ESM CanAM4 CNRM-CM5 GFDL-HIRAM-C360 GISS-E2-R inmcm4 IPSL-CM5A-MR MIROC5 MRI-AGCM3-2S NorESM1-M; do

for i in bcc-csm1-1-m BNU-ESM CanAM4 CNRM-CM5 GFDL-HIRAM-C360 GISS-E2-R inmcm4 IPSL-CM5A-MR MIROC5 MRI-AGCM3-2S NorESM1-M; do

  #cdo sub /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_pr_minus_evap_${i}_amip_r1i1p1.nc /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_mrro_Amon_${i}_amip_r1i1p1.nc /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_Amon_${i}_amip_r1i1p1.nc

  #cdo muldpm /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_Amon_${i}_amip_r1i1p1.nc /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_month_Amon_${i}_amip_r1i1p1.nc
  #cdo selsmon,1 /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_Amon_${i}_amip_r1i1p1.nc /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_jan_Amon_${i}_amip_r1i1p1.nc
  #cdo muldpm /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_jan_Amon_${i}_amip_r1i1p1.nc /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_jan_month_Amon_${i}_amip_r1i1p1.nc

  #cdo sub /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_month_Amon_${i}_amip_r1i1p1.nc /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_jan_month_Amon_${i}_amip_r1i1p1.nc /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_swc_anom_Amon_${i}_amip_r1i1p1.nc

  #cdo ymonmean /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_month_Amon_${i}_amip_r1i1p1.nc /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_month_mean_Amon_${i}_amip_r1i1p1.nc
  #cdo ymonmean /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_jan_month_Amon_${i}_amip_r1i1p1.nc /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_jan_month_mean_Amon_${i}_amip_r1i1p1.nc

  #cdo sub /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_month_mean_Amon_${i}_amip_r1i1p1.nc /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_sa_jan_month_mean_Amon_${i}_amip_r1i1p1.nc /ouce-home/students/kebl4396/Paper1/Paper1SAModelFiles/atlas_swc_anom_Amon_${i}_amip_r1i1p1.nc
done
