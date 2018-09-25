# Calls gen_test - general testing of CABLE for a range
# of Fluxnet sites - producing four plots per site in 
# a new directory named "CABLE[date,time]".
# Requires CABLEplots.R to be in the same directory.
# Gab Abramowitz, UNSW/CSIRO-MAR, 2007, gabsun@gmail.com

# Modified by BP (1 Nov 2007) to run internal testing for Bondville
# This file takes met data for 1998-9 and calls gen_test in CABLEplotsBV2nd.R
# with restart files. The gen_test there have checking and spin up turned off.

# Output file type for plots - 'pdf','ps','png' or 'jpg'
outtype='pdf' # 
# Site names are used to name output files and headings
sitenames=c(
            'Bondville98',
            'Bondville99'
           )
# CABLE input met data files:		
infiles=c(
          '../internal_test_met/metBV98_igbp.nc',
          '../internal_test_met/metBV99_igbp.nc'
         )
# Restart files for these sites (no spinup will be performed):
restartfiles=c(
               './restart_out.nc', # NB. if performed right after BVcomp1st.R
                              # else, need to copy restart_Bondville97.nc here.
               './restart_out.nc'
              )
# plottypes can be and combination of:
# diurnalflux, seasflux, avwindow, timeseries
# - timeseries, called from here, will always produce a png file
plottypes=c(
            'diurnalflux',
            'seasflux',
            'timeseries',
            'avwinflux'
           )
# BP added the old version for comparison
oldfiles=c(
           '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_Bondville98.nc',
           '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_Bondville99.nc'
    )
# Load gen_test function:
source('CABLEplotsBV2nd.R')
# Call gen_test	
gen_test(infiles,sitenames,outtype,restartfiles,plottypes,oldfiles)
