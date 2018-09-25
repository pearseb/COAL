# Script for plotting linear regression gradient, rsq and rmse for a range of
# temporally averaged CABLE output fluxes, reading from CABLE output
# netcdf file.
# Requires 'timeseries' function in CABLEplots.R
# Gab Abramowitz, UNSW/CSIRO-MAR, gabsun@gmail.com
#
# modified by BP (Oct 2007)
###########################################
sitename='Tumbarumba'
outtype='pdf' # choose 'screen','pdf','ps','png' or 'jpg'
#cablefile='../out_cable.nc' # if you are sure about the file's identity
# user need to modify the directory name for the following line (BP)
cablefile='CABLE2007-10-25_11.37.55/out_Tumbarumba.nc' # CABLE output file
obsfile='../sample_met/Tumbarumba.nc' # observed flux data file
##########################################
# Load CABLE plot functions:
source('CABLEplots.R')
avwinflux(outtype,sitename,cablefile,obsfile)
