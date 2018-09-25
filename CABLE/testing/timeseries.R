# Script for plotting a simple time series of many variables simultaneously.
# Reads from single site CABLE output netcdf file.
# Gab Abramowitz, UNSW/CSIRO-MAR, gabsun@gmail.com
# modified by BP (Oct 2007)
###################################################
outtype='png' # choose 'screen','pdf','ps','png' or 'jpg'
sitename='Tumba'
#cablefile='../out_cable.nc' # if you are sure about the file's identity
# user need to modify the directory name for the following line (BP)
cablefile='CABLE2007-10-25_11.37.55/out_Tumbarumba.nc' # CABLE output file
timestep_start=1
timestep_stop=1000
##################################################
# Load CABLE plot functions:
source('CABLEplots.R')
timeseries(sitename,cablefile,outtype,timestep_start,timestep_stop)
