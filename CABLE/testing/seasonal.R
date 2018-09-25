# Script for plotting monthly average (seasonal cycle) from CABLE output
# Gab Abramowitz, UNSW/CSIRO-MAR, gabsun@gmail.com
#
# This script will plot the monthly average of Rnet, 
# NEE, latent and sensible heat from netcdf CABLE output and 
# observations over an entire single-site data set.
# For this to be successful, the dataset **MUST START AT JAN 1**
# AND be a integer number of years in length.
#
# modified by BP (Oct 2007)
############################################
outtype='pdf' # choose 'screen','pdf','jpg' or 'ps'
sitename='Tumbarumba'
#cablefile='../out_cable.nc' # if you are sure about the file's identity
# user need to modify the directory name for the following line (BP)
cablefile='CABLE2007-10-25_11.37.55/out_Tumbarumba.nc' # CABLE output file
obsfile='../sample_met/Tumbarumba.nc' # observed flux data file
###########################################
# Load CABLE plot functions:
source('CABLEplots.R')
seasflux(outtype,sitename,cablefile,obsfile)
