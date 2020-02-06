# BWF_Synthesis
# This repository contains MATLAB code used for the following publication:
# Helbig et al. (2020): Increasing contribution of peatlands to boreal evapotranspiration in a warming climate. Nature Climate Change ??
############################################################
# The following code is available:
# gs.m - derives surface conductance by inverting the Penman-Monteith equation
# ga.m - calculates aerodynamic conductance based on wind speed and friction velocity measurements
# dryness_index_CRU.m - calculates dryness index from Climate Research Unit data
# future_ET_ratio.m - estimates the ratio of peatland and forest evapotranspiration under current and future climates 
# boundary_line_analysis.m - analysis of the upper boundary response of surface conductance and evapotranspiration to vapour pressure deficit
# VPD_CMIP5.m - derivation of daily maximum vapour pressure deficit from Earth system model output
# ET_map.m - creates a map of the ratio of peatland to forest evapotranspiration and its change under future climates
# ET_gs_responses.m - plots the upper boundary response of peatland and forest evapotranspiration to vapour pressure deficit
#############################################################
# following datasets are available
# CA-NOB.csv - example file for eddy covariance data (from the Nobel EC site, Canada)
# CIRCUM_BOREAL.shp - shapefile outlining the boreal biome [based on Olson, D. M. et al. Terrestrial ecoregions of the world: a new map of life on Earth. BioScience 51, 933 (2001)]
# in_boreal_CRU.mat - gridded Matlab file (0.5deg x 0.5deg) masking the boreal biome (boreal = 1; non-boreal = 0)
#############################################################
