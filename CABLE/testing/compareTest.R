# Calls gen_test - general testing of CABLE for a range
# of Fluxnet sites - producing four plots per site in 
# a new directory named "CABLE[date,time]".
# Requires CABLEplotsComp.R to be in the same directory.
# Gab Abramowitz, UNSW/CSIRO-MAR, 2007, gabsun@gmail.com

# Modified by BP (2 Jan 2008) to run internal testing
# and output results compared to last version.

# BP merged offline codes with online codes (Oct 2010)

# Usage on Vayu:
# Type 'module load R' to get the environment and paths
# Type 'R' to launch the software
# Type 'source("compareTest.R")' to compile, run CABLE and plot output in one go

################################################################################
# Output file type for plots - 'pdf','ps','png' or 'jpg'
outtype='pdf' # 
# Site names are used to name output files and headings
sitenames=c(
            'Dinghushan',
            'WeidenBrunnen',
            'Tharandt',
            'Tumbarumba',
            'Harvard',
            'Metolius'
 # this next batch require check%ranges=.FALSE.
 #           'WalkerBranch',
 #           'Hyytiala',
 #           'LittleWashita',
 #           'Bondville97'
           )
# CABLE input met data files:		
infiles=c(
          'sample_met/metDH_igbp.nc',
          'sample_met/Weidenbrunnen_igbp.nc',
          'sample_met/metTH_igbp.nc',
          'sample_met/Tumbarumba2002-06_igbp.nc',
          'sample_met/metHV_igbp.nc',
          'sample_met/Metolius_igbp.nc'
 # this next batch require check%ranges=.FALSE.
 #         'sample_met/WalkerBranch_igbp.nc',
 #         'sample_met/Hyytiala_igbp.nc',
 #         'sample_met/Littlewashita_igbp.nc',
 #         'sample_met/metBV97_igbp.nc'
         )
# Restart files for these sites (no spinup will be performed):
#restartfiles=c('./restart_Tumbarumba.nc',
#		'./restart_Bondville.nc',
#		'./restart_Tharandt.nc')
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
           '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_Dinghushan.nc',
           '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_WeidenBrunnen.nc',
           '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_Tharandt.nc',
           '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_Tumbarumba.nc',
           '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_Harvard.nc',
           '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_Metolius.nc'
 #          '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_WalkerBranch.nc',
 #          '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_Hyytiala.nc',
 #          '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_LittleWashita.nc',
 #          '/home/599/bep599/Mk3L_CABLE20091222/core/src/model/CABLE/testing/ver1.4b/out_Bondville97.nc'
    )
# Load gen_test function:
source('CABLEplotsComp.R')
# Call gen_test	
gen_test(infiles,sitenames,outtype,plottypes,oldfiles)
