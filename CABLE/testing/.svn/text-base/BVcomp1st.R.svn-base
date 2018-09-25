# Calls gen_test - general testing of CABLE for a range
# of Fluxnet sites - producing four plots per site in 
# a new directory named "CABLE[date,time]".
# Requires CABLEplots.R to be in the same directory.
# Gab Abramowitz, UNSW/CSIRO-MAR, 2007, gabsun@gmail.com

# Modified by BP (1 Nov 2007) to run internal testing for Bondville
# This file only takes 1997 met data and calls gen_test in CABLEplots_BV_LW.R
# to spin up. The gen_test there have checking turned off.

# Output file type for plots - 'pdf','ps','png' or 'jpg'
outtype='pdf' # 
# Site names are used to name output files and headings
sitenames=c('Bondville97')
# CABLE input met data files:		
infiles=c('../internal_test_met/metBV97_igbp.nc')
# Restart files for these sites (no spinup will be performed):
#restartfiles=c('./res_Bondville97.nc')
# plottypes can be and combination of:
# diurnalflux, seasflux, avwindow, timeseries
# - timeseries, called from here, will always produce a png file
plottypes=c('diurnalflux',
			'seasflux',
			'timeseries',
			'avwinflux')
# BP added the old version for comparison
oldfiles=c(
     '/home/pak007/cable_offline/cable_v1.4/testing/internal_test_sites/ver1.4b/out_Bondville97.nc'
    )
# Load gen_test function:
source('CABLEplots_BV_LW.R')
# Call gen_test	
gen_test(infiles,sitenames,outtype,plottypes,oldfiles)
