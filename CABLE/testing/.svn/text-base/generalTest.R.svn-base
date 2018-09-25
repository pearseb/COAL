# Calls gen_test - general testing of CABLE for a range
# of Fluxnet sites - producing three plots per site in 
# a new directory named "CABLE[date,time]".
# Requires CABLEplots.R to be in the same directory.
# Gab Abramowitz, UNSW/CSIRO-MAR, 2007, gabsun@gmail.com

# BP merged offline codes with online codes (Oct 2010)

# Usage on Vayu:
# Type 'module load R' to get the environment and paths
# Type 'R' to launch the software
# Type 'source("generalTest.R")' to compile, run CABLE and plot output in one go

################################################################################
# Output file type for plots - 'pdf','ps','png' or 'jpg'
outtype='pdf' # 
# Site names are used to name output files and headings
sitenames=c(
            'Tumbarumba'
           ,'Bondville'
           ,'Tharandt'
            )
# CABLE input met data files:		
infiles=c(
          'sample_met/Tumbarumba_igbp.nc'
         ,'sample_met/Bondville_igbp.nc'
         ,'sample_met/Tharandt_igbp.nc'
          )
# Restart files for these sites (no spinup will be performed):
restartfiles=c(
               './restart_Tumbarumba.nc'
              ,'./restart_Bondville.nc'
              ,'./restart_Tharandt.nc'
               )
# plottypes can be and combination of:
# diurnalflux, seasflux, avwindow, timeseries
# - timeseries, called from here, will always produce a png file
plottypes=c('diurnalflux'
           ,'seasflux'
           ,'timeseries'
           ,'avwinflux'
           )
# Load gen_test function:
source('CABLEplots.R')
# Call gen_test	
# gen_test(infiles,sitenames,outtype,plottypes)
gen_test(infiles,sitenames,outtype,restartfiles,plottypes)
