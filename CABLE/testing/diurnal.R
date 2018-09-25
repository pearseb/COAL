# Script for plotting average diurnal cycle from CABLE output
# Gab Abramowitz, UNSW/CSIRO-MAR, gabsun@gmail.com
#
# This script will plot the average diurnal cycle of Rnet, 
# NEE, latent and sensible heat from netcdf CABLE output and 
# observations, seasonally, over an entire, integer-year 
# single-site data set. Dataset **MUST START AT JAN 1**
#
# modified by BP (Oct 2007)
#################################################
outtype='pdf' # choose 'screen','pdf','ps','png' or 'jpg'
sitename='Tumbarumba'
#cablefile='../out_cable.nc' # if you are sure about the file's identity
# user need to modify the directory name for the following line (BP)
cablefile='CABLE2007-10-25_11.37.55/out_Tumbarumba.nc' # CABLE output file
obsfile='../sample_met/Tumbarumba.nc' # observed flux data file
#################################################
# Load CABLE plot functions:
source('CABLEplots.R')
diurnalflux(outtype,sitename,cablefile,obsfile)
