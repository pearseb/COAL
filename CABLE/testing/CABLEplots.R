# This file contains a collection of functions for plotting CABLE output.
# They are called from a variety of other R scripts in this directory.
# Details about each function are above the function calls below.
# Gab Abramowitz, UNSW/CSIRO-MAR, 2007, gabsun@gmail.com

# BP noted that timeseries crashing on both Shine and Cherax platform
# due to the width and height specified in the png/jpg graph exceeded X11 
# memory (?!) specifications. They are now scaled down to make it work.
# (BP, Oct 2007)
# The pointsize for png/jpg graph in timeseries is changed from 24 to 12
# to make the legend smaller. This may be necessary for other plots
# if user choose to use png/jpg instead of pdf. (BP, Oct 2007)

# Gab fixed a bug (sqerr in averaging windows) in 1st Nov 2007

# BP fixed a bug (stid[k+4] and fnid[k+4] for DJF) in 2nd Jan 2008

# Jhan fixed a bug (meancab=mean(flav[,1])) in Mar 2009

# BP added casaCNP namelist variables; merged offline codes with online codes;
# and used same min, max values for plots in diurnal cycles (Oct 2010)
# Removed unit conversion for SoilMoist as output in the Mk3L version is
# already in m^3/m^3 (Oct 2010)
# Also copied the cable.nml file to archive as record of what input files
# and switches are used in that run (Oct 2010)
# Tair in output file is now in 3D instead of 4D (Oct 2010)

library(ncdf) # load netcdf library
##########################################################################################
# gen_test
# This function serves as a general evaluation test for CABLE. It calls
# most of the other functions for a variety of Fluxnet sites and puts 
# results in a new directory named by the date and time the test was executed.
gen_test=function(infiles,sitenames,outtype,restartfiles,plottypes){
# gen_test=function(infiles,sitenames,outtype,plottypes){
  #Sort out time of test and create new directory for plots
  now=Sys.time() # get system time
  newdirname=paste('CABLE',substr(now,1,10),'_',substr(now,12,13),'.',
             substr(now,15,16),'.',substr(now,18,19),sep='')
  system(paste('mkdir',newdirname)) # create new directory
  cnlist=c() # initialise CABLE namelist
  for(site in 1:length(infiles)){ # for each site
    # cd to CABLE directory:
    setwd('../')
    # Prepare CABLE namelist: (NB. switch off check%ranges for Bondville)
    #                         Do not use restartfiles if you want to spin up.
    cnlist[1] = '&cable'
    cnlist[2] = paste('  filename%met = \'',infiles[site],'\'',sep='')
    cnlist[3] = '  filename%out = \'out_cable.nc\''
    cnlist[4] = '  filename%log = \'log_cable.txt\''
    cnlist[5] = '  filename%restart_in  = \' \''
    # cnlist[5] = paste('  filename%restart_in  = \'',restartfiles[site],'\'',sep='')
    cnlist[6] = '  filename%restart_out = \'./restart_out.nc\''
    cnlist[7] = '  filename%LAI     = \'surface_data/lai_CCAMtoMk3L.nc\''
    # cnlist[7] = '  filename%LAI     = \'surface_data/LAI_Monthly_Global.nc\''
    cnlist[8] = '  filename%type    = \'surface_data/tmp.nc\''
    # cnlist[8] = '  filename%type    = \'surface_data/VegSoil_Type_Global_igbp.csv\''
    cnlist[9] = '  filename%veg     = \'surface_data/def_veg_params_igbp.txt\''
    # cnlist[9] = '  filename%veg     = \'surface_data/def_veg_params.txt\''
    cnlist[10] = '  filename%soil    = \'surface_data/def_soil_params.txt\''
    cnlist[11] = '  vegparmnew = .TRUE.  ! using new format'
    cnlist[12] = '  spinup = .TRUE.  ! do we spin up the model?'
    cnlist[13] = '  delsoilM = 0.001   ! allowed variation in soil moisture for spin up'
    cnlist[14] = '  delsoilT = 0.01    ! allowed variation in soil temperature for spin up'
    cnlist[15] = '  output%restart = .TRUE.  ! should a restart file be created?'
    cnlist[16] = '  output%met = .TRUE.  ! input met data'
    cnlist[17] = '  output%flux = .TRUE.  ! convective, runoff, NEE'
    cnlist[18] = '  output%soil = .TRUE.  ! soil states'
    cnlist[19] = '  output%snow = .TRUE.  ! snow states'
    cnlist[20] = '  output%radiation = .TRUE.  ! net rad, albedo'
    cnlist[21] = '  output%carbon    = .TRUE.  ! NEE, GPP, NPP, stores' 
    cnlist[22] = '  output%veg       = .TRUE.  ! vegetation states'
    cnlist[23] = '  output%params    = .TRUE.  ! input parameters used to produce run'
    cnlist[24] = '  output%balances  = .TRUE.  ! energy and water balances'
    cnlist[25] = '  check%ranges     = .TRUE.  ! variable ranges, input and output'
    cnlist[26] = '  check%energy_bal = .TRUE.  ! energy balance'
    cnlist[27] = '  check%mass_bal   = .TRUE.  ! water/mass balance'
    cnlist[28] = '  verbose = .TRUE. ! write details of every grid cell init and params to log?'
    cnlist[29] = '  leaps = .FALSE. ! calculate timing with leap years?'
    cnlist[30] = '  logn = 88      ! log file number - declared in input module'
    cnlist[31] = '  fixedCO2 = 350.0   ! if not found in met file, in ppmv'
    cnlist[32] = '  spincasainput = .FALSE.    ! input required to spin casacnp offline'
    cnlist[33] = '  spincasa      = .TRUE.     ! spin casa before running the model if TRUE, and should be set to FALSE if spincasainput = .TRUE.'
    cnlist[34] = '  icycle = 0   ! BP pull it out from casadimension and put here; 0 for not using casaCNP, 1 for C, 2 for C+N, 3 for C+N+P'
    cnlist[35] = '  casafile%cnpbiome=\'surface_data/pftlookup_igbp_DH.csv\'  ! biome specific BGC parameters'
    cnlist[36] = '  casafile%cnppoint=\'surface_data/siteDH_igbp.csv\'        ! point information'
    cnlist[37] = '  casafile%cnpepool=\'surface_data/poolcnpDH_igbpe.csv\'    ! end of run pool size'
    cnlist[38] = '  casafile%cnpipool=\'surface_data/poolcnpDH_igbpi.csv\'    ! initial pool size'
    cnlist[39] = '  casafile%cnpmet=\'surface_data/casametDH.csv\'            ! daily met forcing for spinning casacnp'
    cnlist[40] = '  casafile%phen=\'surface_data/modis_phenology.txt\'        ! modis phenology'
    cnlist[41] = '&end'
    write(cnlist,'cable.nml')
    # Run CABLE
    print(paste('Running CABLE at ',sitenames[site],':',sep=''))
    system('make -fMakefile_offline') # Will make and run CABLE
#    system('./cable') # just run CABLE
    setwd(paste('testing/',newdirname,sep='')) # cd to new directory
    # save a copy of the output
    file.copy('../../out_cable.nc',paste('out_',sitenames[site],'.nc',sep=''))
    file.copy('../../log_cable.txt',paste('log_',sitenames[site],'.txt',sep=''))
    file.copy('../../restart_out.nc',paste('restart_',sitenames[site],'.nc',
              sep=''))
    file.copy('../../cable.nml',paste('cable_',sitenames[site],'.nml',sep=''))
    # Produce plots:
    if(any(plottypes=='seasflux')) seasflux(outtype,sitenames[site],
    	'../../out_cable.nc',paste('../../',infiles[site],sep=''))
    if(any(plottypes=='diurnalflux')) diurnalflux(outtype,sitenames[site],
    	'../../out_cable.nc',paste('../../',infiles[site],sep=''))
    if(any(plottypes=='avwinflux')) avwinflux(outtype,sitenames[site],
    	'../../out_cable.nc',paste('../../',infiles[site],sep=''))
    if(any(plottypes=='timeseries')) timeseries(sitenames[site],
    	'../../out_cable.nc')
    setwd('../') # return to original (testing) directory
  }	
}
########################################################################################
# seasflux
# seasflux plots the monthly average of Rnet, 
# NEE, latent and sensible heat from netcdf CABLE output and 
# observations over an entire single-site data set.
# For this to be successful, the dataset **MUST START AT JAN 1**
# AND be a integer number of years in length.
seasflux=function(outtype='pdf',sitename,cablefile,obsfile){
  # Flux names to plot and their units:
  flux=c('NEE','Qle','Qh','Rnet') # outputs to plot and
  units=c('umol/m2/s','W/m2','W/m2','W/m2') # their units
  # Output file name if not writing to screen:
  outfilename=paste(sitename,'Season.',outtype,sep='')
  # Choose output file type, if not 'screen':
  if(outtype=='pdf'){
    pdf(file=outfilename,paper='a4r',width=11,height=8)
  }else if(outtype=='ps'){
    postscript(file=outfilename,paper='special',width=11,height=8)
  }else if(outtype=='png'){
    png(file=outfilename,width=2000,height=1414,pointsize=24)
  }else if(outtype=='jpg'){
    jpeg(file=outfilename,width=4000,height=2828,pointsize=24)
  }
  cnc=open.ncdf(cablefile,readunlim=FALSE) # open CABLE file
  onc=open.ncdf(obsfile,readunlim=FALSE) # open observed data file
  # Get timing details (1st 2 time steps of 'time' variable):
  time12=get.var.ncdf(cnc,'time',start=1,count=2) # read CABLE 'time' variable
  timestepsize=time12[2]-time12[1]
  tstepinday=86400/timestepsize # number of time steps in a day
  # plot layout:
  layout(matrix(1:length(flux),2,2))
  # days on which each month begin:
  month_start=c()
  month_start[1]=1    # Jan
  month_start[2]=32   # Feb
  month_start[3]=60   # Mar
  month_start[4]=91   # Apr
  month_start[5]=121  # May
  month_start[6]=152  # Jun
  month_start[7]=182  # Jul
  month_start[8]=213  # Aug
  month_start[9]=244  # Sep
  month_start[10]=274 # Oct
  month_start[11]=305 # Nov
  month_start[12]=335 # Dec
  month_start[13]=366 # i.e. beginning of next year

  cable_monthly=c() # initialise monthly averages
  obs_monthly=c()   # initialise monthly averages
  for(j in 1:length(flux)){ # for each flux variable
    # Read data:
    cable_data=get.var.ncdf(cnc,flux[j]) # read CABLE output data
    obs_data=get.var.ncdf(onc,flux[j])   # read observed data
    alllen=length(cable_data) # total number of timesteps
    # Reshape into column days:
    cable_days=matrix(cable_data,ncol=tstepinday,byrow=TRUE) 
    obs_days=matrix(obs_data,ncol=tstepinday,byrow=TRUE) 
    ndays=length(cable_days[,1]) # find # days in data set
    nyears=as.integer(ndays/365) # find # years in data set
    cavday=c() # initialise
    oavday=c() # initialise
    # transform data into daily averages:
    for(i in 1:ndays){
      cavday[i]=mean(cable_days[i,]) # calc daily average flux
      oavday[i]=mean(obs_days[i,]) # calc daily average flux
    }	
    # Transform daily means into monthly means:
    for(l in 1:12){ # for each month
      month_length=month_start[l+1]-month_start[l]
      cable_month=0 # initialise
      obs_month=0   # initialise
      for(k in 1:nyears){ # for each year of data set
        # add all daily averages for a given month
        # over all data set years:
        cable_month = cable_month + 
          sum(cavday[(month_start[l]+(k-1)*365):
          (month_start[l+1]-1 +(k-1)*365) ] )
        obs_month = obs_month + 
          sum(oavday[(month_start[l]+(k-1)*365):
          (month_start[l+1]-1 +(k-1)*365) ] )
      }
      # then divide by the total number of days added above:
      cable_monthly[l]=cable_month/(month_length*nyears)
      obs_monthly[l]=obs_month/(month_length*nyears)
    }
    xloc=c(1:12) # set location of x-coords
    # Plot CABLE output result:
    plot(xloc,cable_monthly,type="l",xaxt="n",xlab='Month',
      ylab=paste('Average',flux[j],'flux',units[j]),lwd=2,col='blue',
      ylim=c(min(cable_monthly,obs_monthly),
      max(cable_monthly,obs_monthly)))
    # Then plot obs result:
    lines(xloc,obs_monthly,lwd=2)
    axis(1,at=c(2,4,6,8,10,12),labels=c('2','4','6','8','10','12'))
    title(paste(sitename,flux[j])) # add title
  }
  legend(1,150,c('CABLE','obs'),lty=1,col=c('blue','black'),lwd=2,bty="n")
  # Close netcdf files
  close.ncdf(cnc)
  close.ncdf(onc)
  # close graphics file if used:
  if(outtype!='screen') dev.off()
  # Clear all variables in this environment
  rm(list=ls())
} # End seasflux function
########################################################################################
# diurnalflux
# This function will plot the average diurnal cycle of Rnet, 
# NEE, latent and sensible heat from netcdf CABLE output and 
# observations, seasonally, over an entire, integer-year 
# single-site data set. Dataset **MUST START AT JAN 1**
diurnalflux = function(outtype,sitename,cablefile,obsfile){
  outfilename=paste(sitename,'Diurnal.',outtype,sep='') # name of plot if not to screen
  flux=c('NEE','Qle','Qh','Rnet') # outputs to plot and
  units=c('umol/m2/s','W/m2','W/m2','W/m2') # their units
  # Set output file type, if not to screen:
  if(outtype=='pdf'){
    pdf(file=outfilename,paper='a4r',width=11,height=8)
  }else if(outtype=='ps'){
    postscript(file=outfilename,paper='special',width=11,height=8)
  }else if(outtype=='png'){
    png(file=outfilename,width=2000,height=1414,pointsize=24)
  }else if(outtype=='jpg'){
    jpeg(file=outfilename,width=4000,height=2828,pointsize=24)
  }
  labels=c('DJF','MAM','JJA','SON')
  stid=c(1,60,152,244,335) # seasonal divisions in a year
  fnid=c(59,151,243,334,365)
  cnc=open.ncdf(cablefile,readunlim=FALSE) # open CABLE file
  onc=open.ncdf(obsfile,readunlim=FALSE) # open observed data file
  # Get timing details (1st 2 time steps of 'time' variable):
  time12=get.var.ncdf(cnc,'time',start=1,count=2) # read CABLE 'time' variable
  timestepsize=time12[2]-time12[1]
  tstepinday=86400/timestepsize # number of time steps in a day
  # plot layout:
  layout(matrix(1:16,4,4))
  for(j in 1:length(flux)){ # for each flux variable
    # read data:
    cable_data=get.var.ncdf(cnc,flux[j]) # read CABLE output data
    obs_data=get.var.ncdf(onc,flux[j])   # read observed data
    # reshape into column days
    cable_days=matrix(cable_data,ncol=tstepinday,byrow=TRUE) 
    obs_days=matrix(obs_data,ncol=tstepinday,byrow=TRUE) 
    ndays=length(cable_days[,1]) # find # days in data set
    nyears=as.integer(ndays/365) # find # years in data set
    minFlux=0.0
    maxFlux=0.0
    cavday=array(0.0,c(4,tstepinday))
    oavday=array(0.0,c(4,tstepinday))
    for(k in 1:4){# for each season (DJF, MAM etc)
      Fmin=minFlux
      Fmax=maxFlux
#      cavday=c() # declare
#      oavday=c() # declare
#      cavday[1:tstepinday]=0 # initialise
#      oavday[1:tstepinday]=0 # initialise
      # sum up fluxes over each year of data set for current season
      for(l in 1:nyears){	
        for(i in 1:tstepinday){
          # calc average flux for each hour
          cavday[k,i]=cavday[k,i] + 
            sum(cable_days[(stid[k]+(l-1)*365):(fnid[k]+(l-1)*365),i])
          oavday[k,i]=oavday[k,i] + 
            sum(obs_days[(stid[k]+(l-1)*365):(fnid[k]+(l-1)*365),i]) 
        }
        if(k==1){ # i.e. DJF
          # add Dec to Jan/Feb
          for(i in 1:tstepinday){
            cavday[k,i]=cavday[k,i] + 
              sum(cable_days[(stid[k+4]+(l-1)*365):(fnid[k+4]+(l-1)*365),i])
            oavday[k,i]=oavday[k,i] + 
              sum(obs_days[(stid[k+4]+(l-1)*365):(fnid[k+4]+(l-1)*365),i])
          }
        }
      }
      # Then find the average of these fluxes:
      if(k==1){ # i.e. DJF
        cavday[k,1:tstepinday]=cavday[k,1:tstepinday]/(90*nyears)
        oavday[k,1:tstepinday]=oavday[k,1:tstepinday]/(90*nyears)
      }else{
        cavday[k,1:tstepinday]=cavday[k,1:tstepinday]/((fnid[k]-stid[k])*nyears)
        oavday[k,1:tstepinday]=oavday[k,1:tstepinday]/((fnid[k]-stid[k])*nyears)
      }
      # find the min and max values for plotting
      minFlux=min(cavday[k,1:tstepinday],oavday[k,1:tstepinday],Fmin)
      maxFlux=max(cavday[k,1:tstepinday],oavday[k,1:tstepinday],Fmax)
    }
    for(k in 1:4){# for each season (DJF, MAM etc)
      xloc=c(0:(tstepinday-1)) # set location of x-coords
      # Plot CABLE output result:
      plot(xloc,cavday[k,1:tstepinday],type="l",xaxt="n",xlab='Hour of day',
        ylab=paste(flux[j],'flux',units[j]),lwd=2,col='blue',
        ylim=c(minFlux,maxFlux))
        # ylim=c(min(cavday,oavday),max(cavday,oavday)))
      # Then plot obs result:
      lines(xloc,oavday[k,1:tstepinday],lwd=2)
      axis(1,at=c(0,6*tstepinday/24,12*tstepinday/24,18*tstepinday/24,
        23*tstepinday/24),labels=c('0','6','12','18','23'))
      title(paste(sitename,labels[k],flux[j])) # add title
    }
    rm(cavday,oavday)
  }		
  legend(-2,500,c('CABLE','obs'),lty=1,col=c('blue','black'),lwd=2,bty="n")
  # Close netcdf files
  close.ncdf(cnc)
  close.ncdf(onc)
  # close graphics file if used:
  if(outtype!='screen')dev.off()
  # Clear all current local environment variables
  rm(list=ls()) 	
} # End function diurnalflux
######################################################################################
# avwinflux
# Function for plotting linear regression gradient, rsq and rmse for a range of
# temporally averaged CABLE output fluxes, reading from CABLE output
# netcdf file.
avwinflux = function(outtype,sitename,cablefile,obsfile){
  library(boot) # load bootstrap library
  # Name of plot file if not to screen:
  outfilename=paste(sitename,'Window.',outtype,sep='')
  windowstepsize=2 # resolution of graph
  maxwindow = 30 # in days, the largest averaging window
  flux=c('Rnet','Qle','Qh','NEE')
  units=c('W/m2','W/m2','W/m2','umol/m2/s') # their units
  # Sort out plot output type, if not to screen:
  if(outtype=='pdf'){	
    pdf(file=outfilename,paper='a4r',width=11,height=8)
  }else if(outtype=='ps'){
    postscript(file=outfilename,paper='special',width=11,height=8)
  }else if(outtype=='png'){
    png(file=outfilename,width=2000,height=1414,pointsize=24)
  }else if(outtype=='jpg'){
    jpeg(file=outfilename,width=4000,height=2828,pointsize=24)
  }
  layout(matrix(1:12,3,4))
  cnc=open.ncdf(cablefile,readunlim=FALSE) # open CABLE file
  onc=open.ncdf(obsfile,readunlim=FALSE) # open observed data file
  # Get timing details (1st 2 time steps of 'time' variable):
  time12=get.var.ncdf(cnc,'time',start=1,count=2) # read CABLE 'time' variable
  timestepsize=time12[2]-time12[1]
  tstepinday=86400/timestepsize # number of time steps in a day
  numcalcs=maxwindow*tstepinday/windowstepsize
  for(l in 1:length(flux)){
    # read data:
    cable_data=get.var.ncdf(cnc,flux[l]) # read CABLE output data
    obs_data=get.var.ncdf(onc,flux[l])   # read observed data
    tsteps=length(cable_data) # total number of timesteps
    mvals = c() # initialise gradient values
    mwidth =c() # initialise 95% confidence interval for regression
    rvals = c() # initialise correlation coefficient values
    rmse = c()  # intialise RMSE
    xloc = c() # initialise x ticks
    for(i in 1:numcalcs){
      # calculate window size:
      windowsize=windowstepsize*i
      # reduce data set size if necessary for reshaping:
      numdump=tsteps %% windowsize # "%%" is modulo
      nwindows=(tsteps-numdump)/windowsize # number of windows
      sqerr=c()
      flav=matrix(0,nwindows,2) # init
      # Reshape data:
      cdat=matrix(cable_data[1:(tsteps-numdump)],windowsize,nwindows)
      odat=matrix(obs_data[1:(tsteps-numdump)],windowsize,nwindows)
      for(j in 1:nwindows){ # calculate average flux for each window
        flav[j,1]=mean(cdat[,j])	# vector of averages
        flav[j,2]=mean(odat[,j])	# vector of averages
        sqerr[j]=(flav[j,1]-flav[j,2])^2
      }
      # Perform least squares regression (b/w cable ad obs):
      rgrs=lsfit(flav[,2],flav[,1])
      mvals[i]=rgrs$coef[[2]] # store gradient values for plot
      # Get 95% confidence interval:
      linval=c() # init
      lindev=c() # init
      cabdev=c() # init
      meancab=mean(flav[,1])  # meancab=mean(flav[j,1])
      for(k in 1:nwindows){
        linval[k]=rgrs$coef[[2]]*flav[k,1]+rgrs$coef[[1]]
        lindev[k]=(linval[k]-flav[k,2])^2
        cabdev[k]=(flav[k,1]-meancab)^2
      }
      # 95% confidence interval:
      mwidth[i]=sqrt(sum(lindev)/nwindows)*1.96/sqrt(sum(cabdev))
      # get correlation cofficient:
      rvals[i]=corr(flav)^2
      # Calculate RMSE:
      rmse[i]=sqrt(mean(sqerr))	
      rm(cdat,odat,flav,sqerr) # clear variables which will have different lengths
      xloc[i]=windowsize/tstepinday
    }
    # Draw RMSE plot:
    plot(xloc,rmse,type="l",xaxt="n",xlab='Averaging window size (days)',
    ylab='RMSE',lwd=1,col='black')
    # Work out x-axis ticks:
    tix=c(0,as.integer(maxwindow/4),as.integer(2*maxwindow/4),
    as.integer(3*maxwindow/4),maxwindow)
    axis(1,at=tix,labels=as.character(tix))
    title(paste(sitename,'- RMSE for av.',flux[l])) # add title
    # Draw gradient plot:
    plot(xloc,mvals,type="l",xaxt="n",xlab='Averaging window size (days)',
      ylab='CABLE vs. obs reg. grad.',lwd=1,col='blue',
      ylim=c(min(mvals-mwidth),max(mvals+mwidth)))
    # Sort out confidence polygon:
    polyX=c(xloc,xloc[numcalcs:1])
    polyY=c(mvals+mwidth,mvals[numcalcs:1]-mwidth[numcalcs:1])
    polygon(polyX,polyY,col='grey')
    lines(xloc,mvals,lwd=1,col='blue')		
    axis(1,at=tix,labels=as.character(tix))
    title(paste(sitename,'- gradient for av.',flux[l])) # add title
    # Draw correlation plot:
    plot(xloc,rvals,type="l",xaxt="n",xlab='Averaging window size (days)',
      ylab=expression(R^2),lwd=1,col='black')
    axis(1,at=tix,labels=as.character(tix))
    title(paste(sitename,'- cor. coeff. for av.',flux[l])) # add title
  }
  # Close CABLE output netcdf file:	
  close.ncdf(cnc)
  # close graphics file if used:
  if(outtype!='screen') dev.off()
  # Clear all current local environment variables
  rm(list=ls()) 	
} # End function avwinflux
##################################################################################
# This function simply plots a timeseries of the entire dataset, for a wide
# range of variables.
# png and jpg types: width scaled down from 3564 to 2000,
#                    height from 2520 to 1414
#                    pointsize from 24 to 12 (BP, Oct 2007)
timeseries=function(sitename,cablefile,outtype='png',timestep_start=-1,timestep_stop=-1){
  outfilename=paste(sitename,'Series.',outtype,sep='') # name of plot if not to screen
  if(outtype=='pdf'){
    pdf(file=outfilename,paper='a4r',width=11,height=8)
  }else if(outtype=='ps'){
    postscript(file=outfilename,paper='special',width=11,height=8)
  }else if(outtype=='png'){
    png(file=outfilename,width=2000,height=1414,pointsize=12)
  }else if(outtype=='jpg'){
    jpeg(file=outfilename,width=2000,height=1414,pointsize=12)
  }
  # open CABLE file:
  cnc=open.ncdf(cablefile,readunlim=FALSE) # open CABLE file
  # set number of time steps in plot (data set length):
  if(timestep_start==-1){
    timestep_start=1
    timestep_stop = cnc$dim$time[[2]]
  }
  # sort out x-axis markings for plots:
  ntsteps=timestep_stop-timestep_start+1
  xax=c(timestep_start,as.integer(timestep_start+ntsteps/5),
    as.integer(timestep_start+2*ntsteps/5),
    as.integer(timestep_start+3*ntsteps/5),
    as.integer(timestep_start+4*ntsteps/5),timestep_stop)
  xloc=c(timestep_start:timestep_stop)	
  # set layout of plots (10 in 5*2 matrix):
  layout(matrix(1:12,3,4,byrow=TRUE))
  # Read data for first plot:---------------------------------------
  SoilTemp=get.var.ncdf(cnc,'SoilTemp',start=c(1,1,1,timestep_start),
    count=c(1,1,6,ntsteps)) # read soil temp
  # First plot:	
  plot(xloc,SoilTemp[1,],type="l",xaxt="n",xlab='Time steps',
    ylab='Degrees K',lwd=0.5,col='green',
    ylim=c(min(SoilTemp),max(SoilTemp)))
  lines(xloc,SoilTemp[2,],lwd=0.5,col='forestgreen')
  lines(xloc,SoilTemp[3,],lwd=0.5,col='darkslategrey')
  lines(xloc,SoilTemp[4,],lwd=0.5,col='deepskyblue4')
  lines(xloc,SoilTemp[5,],lwd=0.5,col='blue')
  lines(xloc,SoilTemp[6,],lwd=0.5,col='darkblue')
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Soil temperatures -',sitename))
  legend('topleft',c('1','2','3','4','5','6'),lty=1,
    col=c('green','forestgreen','darkslategrey','deepskyblue4',
    'blue','darkblue'),lwd=0.5)
  # Read data for second plot:--------------------------------------------
  Tair=get.var.ncdf(cnc,'Tair',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read air temp
  VegT=get.var.ncdf(cnc,'VegT',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read veg temperature
  # Second plot:	
  plot(xloc,SoilTemp[1,],type="l",xaxt="n",xlab='Time steps',
    ylab='Degrees K',lwd=0.5,col='peru',
    ylim=c(min(SoilTemp[1,],SoilTemp[2,],Tair,VegT),
    max(SoilTemp[1,],SoilTemp[2,],Tair,VegT)))
  lines(xloc,SoilTemp[2,],lwd=0.5,col='sienna4')
  lines(xloc,Tair,lwd=0.5,col='blue')
  lines(xloc,VegT,lwd=0.5,col='green')
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Temperatures -',sitename))
  legend('topleft',c('soil1','soil2','airT','vegT'),lty=1,
    col=c('peru','sienna4','blue','green'),lwd=0.5)
  # Read data for third plot:--------------------------------------------
  HSoil=get.var.ncdf(cnc,'HSoil',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read sensible heat from soil
  ESoil=get.var.ncdf(cnc,'ESoil',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps))*2.5104e6 # read latent heat from soil
  Qg=get.var.ncdf(cnc,'Qg',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read ground heat flux
  # Third plot:	
  plot(xloc,HSoil,type="l",xaxt="n",xlab='Time steps',
    ylab=expression(W/m^2),lwd=0.5,col='red',
    ylim=c(min(HSoil,ESoil,Qg),max(HSoil,ESoil,Qg)))
  lines(xloc,ESoil,lwd=0.5,col='blue')
  lines(xloc,Qg,lwd=0.5,col='peru')
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Soil fluxes -',sitename))
  legend('topleft',c('SSens','SLat','Gflux'),lty=1,
    col=c('red','blue','brown'),lwd=0.5)	
  # Read data for fourth plot:--------------------------------------------
  Evap=get.var.ncdf(cnc,'Evap',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps))*2.5104e6 # read total evapotranspiration
  ECanop=get.var.ncdf(cnc,'ECanop',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps))*2.5104e6 # read wet canopy evaporation
  TVeg=get.var.ncdf(cnc,'TVeg',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps))*2.5104e6 # read vegetation transpiration
  HVeg=get.var.ncdf(cnc,'HVeg',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read sensible heat from vegetation
  # Fourth plot:	
  plot(xloc,Evap,type="l",xaxt="n",xlab='Time steps',
    ylab=expression(W/m^2),lwd=0.5,col='black',
    ylim=c(min(Evap,ECanop,TVeg,HVeg),max(Evap,ECanop,TVeg,HVeg)))
  lines(xloc,ECanop,lwd=0.5,col='blue')
  lines(xloc,TVeg,lwd=0.5,col='green')
  lines(xloc,HVeg,lwd=0.5,col='red')
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Canopy fluxes -',sitename))
  legend('topleft',c('TotEvapotr','VEvap','VTransp','VSens'),
    lty=1,col=c	('black','blue','green','red'),lwd=0.5)
  # Read data for fifth plot:---------------------------------------
  SoilMoist=get.var.ncdf(cnc,'SoilMoist',start=c(1,1,1,timestep_start),
    count=c(1,1,6,ntsteps)) # read soil moisture
  zse=get.var.ncdf(cnc,'zse',start=c(1,1,1),count=c(1,1,6)) # read soildepth
  SMoist=SoilMoist
# removed the next 4 lines as there is no need to change units now.
#  SMoist=matrix(0,6,ntsteps)
#  for(i in 1:6){                         # Change units from kg/m^2 to m^3/m^3
#    SMoist[i,]=SoilMoist[i,]/zse[i]/1000
#  }
  # load filed capacity, saturation and wilting point values:
  sfc=vector(mode='numeric',length=ntsteps)
  ssat=vector(mode='numeric',length=ntsteps)
  swilt=vector(mode='numeric',length=ntsteps)
  sfc[]=get.var.ncdf(cnc,'sfc',start=c(1,1),count=c(1,1))
  ssat[]=get.var.ncdf(cnc,'ssat',start=c(1,1),count=c(1,1))
  swilt[]=get.var.ncdf(cnc,'swilt',start=c(1,1),count=c(1,1))
  # Fifth plot:	
  plot(xloc,SMoist[1,],type="l",xaxt="n",xlab='Time steps',
    ylab=expression(m^3/m^3),lwd=0.5,col='green',
    ylim=c(min(SMoist,sfc,swilt,ssat),max(SMoist,sfc,swilt,ssat)))
  lines(xloc,SMoist[2,],lwd=0.5,col='forestgreen')
  lines(xloc,SMoist[3,],lwd=0.5,col='darkslategrey')
  lines(xloc,SMoist[4,],lwd=0.5,col='deepskyblue4')
  lines(xloc,SMoist[5,],lwd=0.5,col='blue')
  lines(xloc,SMoist[6,],lwd=0.5,col='darkblue')
  lines(xloc,sfc,lwd=0.5,col='grey')
  lines(xloc,swilt,lwd=0.5,col='grey')
  lines(xloc,ssat,lwd=0.5,col='grey')
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Soil moisture -',sitename))
  legend('topleft',c('1','2','3','4','5','6'),lty=1,
    col=c('green','forestgreen','darkslategrey','deepskyblue4',
    'blue','darkblue'),lwd=0.5)
  # Read data for sixth plot:--------------------------------------------
  Rainf=get.var.ncdf(cnc,'Rainf',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps))*3600 # read rainfall
  CanopInt=get.var.ncdf(cnc,'CanopInt',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read canopy water storage
  Qs=get.var.ncdf(cnc,'Qs',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps))*3600 # read runoff
  Qsb=get.var.ncdf(cnc,'Qsb',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps))*3600 # read deep drainage
  # Sixth plot:	
  plot(xloc,Rainf,type="l",xaxt="n",xlab='Time steps',
    ylab='mm and mm/h',lwd=0.5,col='blue',
    ylim=c(min(Rainf,CanopInt,Qs,Qsb),max(Rainf,CanopInt,Qs,Qsb)))
  lines(xloc,CanopInt,lwd=0.5,col='green')
  lines(xloc,Qs,lwd=0.5,col='cyan')
  lines(xloc,Qsb,lwd=0.5,col='brown4')
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Precip/runoff -',sitename))
  legend('topleft',c('Precip','CanStore','Runoff','Drainage'),
    lty=1,col=c	('blue','green','cyan','brown4'),lwd=0.5)
  # Read data for seventh plot:--------------------------------------------
  SWnet=get.var.ncdf(cnc,'SWnet',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read net shortwave
  LWnet=get.var.ncdf(cnc,'LWnet',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read net longwave
  SWdown=get.var.ncdf(cnc,'SWdown',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read runoff
  LWdown=get.var.ncdf(cnc,'LWdown',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read deep drainage
  # Seventh plot:	
  plot(xloc,SWdown,type="l",xaxt="n",xlab='Time steps',
    ylab=expression(W/m^2),lwd=0.5,col='black',
    ylim=c(min(SWdown,SWnet,LWdown,LWnet),
    max(SWdown,SWnet,LWdown,LWnet)))
  lines(xloc,SWnet,lwd=0.5,col='blue')
  lines(xloc,LWdown,lwd=0.5,col='green')
  lines(xloc,LWnet,lwd=0.5,col='forestgreen')
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Radiation -',sitename))
  legend('topleft',c('SWdown','SWnet','LWdown','LWnet'),
    lty=1,col=c	('black','blue','green','forestgreen'),lwd=0.5)
  # Read data for eighth plot:--------------------------------------------
  NEE=get.var.ncdf(cnc,'NEE',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read net ecosystem exchange
  AutoResp=get.var.ncdf(cnc,'AutoResp',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read autotrophic respiration
  HeteroResp=get.var.ncdf(cnc,'HeteroResp',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read heterotrophic respiration
  GPP=get.var.ncdf(cnc,'GPP',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read net primary production
  # Eighth plot:	
  plot(xloc,NEE,type="l",xaxt="n",xlab='Time steps',
    ylab=expression(paste(mu,mol/m^2/s)),lwd=0.5,col='blue',
    ylim=c(min(NEE,AutoResp,HeteroResp,GPP),
    max(NEE,AutoResp,HeteroResp,GPP)))
  lines(xloc,AutoResp,lwd=0.5,col='green')
  lines(xloc,HeteroResp,lwd=0.5,col='cyan')
  lines(xloc,GPP,lwd=0.5,col='brown4')
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Carbon fluxes -',sitename))
  legend('topleft',c('NEE','AutoR','HeteroR','GPP'),
    lty=1,col=c	('blue','green','cyan','brown4'),lwd=0.5)
  # Read data for ninth plot:------------------------------------------
  Qle=get.var.ncdf(cnc,'Qle',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read latent heat
  Qh=get.var.ncdf(cnc,'Qh',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read sensible heat
  # Ninth plot:	
  plot(xloc,(SWnet+LWnet),type="l",xaxt="n",xlab='Time steps',
    ylab=expression(W/m^2),lwd=0.5,col='black',
    ylim=c(min(SWnet+LWnet,Qle,Qh,Qg),
    max(SWnet+LWnet,Qle,Qh,Qg)))
  lines(xloc,Qle,lwd=0.5,col='blue')
  lines(xloc,Qh,lwd=0.5,col='red')
  lines(xloc,Qg,lwd=0.5,col='brown4')
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Combined fluxes -',sitename))
  legend('topleft',c('Rnet','Lat','Sens','Gfl'),
    lty=1,col=c	('black','blue','red','brown4'),lwd=0.5)
  # Read data for tenth plot:------------------------------------------
  Albedo=get.var.ncdf(cnc,'Albedo',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read albedo
  SnowDepth=get.var.ncdf(cnc,'SnowDepth',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read snow depth
  LAI=get.var.ncdf(cnc,'LAI',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read LAI
  # Tenth plot:	
  plot(xloc,Albedo,type="l",xaxt="n",xlab='Time steps',
    ylab='- and m',lwd=0.5,col='red',
    ylim=c(min(Albedo,SnowDepth,LAI/10),max(Albedo,SnowDepth,LAI/10)))
  lines(xloc,LAI/10,lwd=0.5,col='forestgreen')
  lines(xloc,SnowDepth,lwd=0.5,col='blue')
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Albedo and snow depth -',sitename))
  legend('topleft',c('Albedo','LAI/10','Snowdepth'),
    lty=1,col=c	('red','forestgreen','blue'),lwd=0.5)
  # Read data for eleventh plot:------------------------------------------
  Ebal=get.var.ncdf(cnc,'Ebal',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read albedo
  plot(xloc,Ebal,type="l",xaxt="n",xlab='Time steps',
    ylab=expression(W/m^2),lwd=0.5,col='red',
    ylim=c(min(Ebal),max(Ebal)))
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Cumulative energy balance - ',sitename))
  # Read data for twelfth plot:------------------------------------------
  Wbal=get.var.ncdf(cnc,'Wbal',start=c(1,1,timestep_start),
    count=c(1,1,ntsteps)) # read albedo
  plot(xloc,Wbal,type="l",xaxt="n",xlab='Time steps',
    ylab='mm',lwd=0.5,col='blue',
    ylim=c(min(Wbal),max(Wbal)))
  axis(1,at=xax,labels=as.character(xax)) # add x-axis labels
  title(paste('Cumulative water balance - ',sitename))

  # Close CABLE output netcdf file:	
  close.ncdf(cnc)
  # close graphics file if used:
  if(outtype!='screen')	dev.off()
  # Clear all current local environment variables
  rm(list=ls()) 		
} # End function timeseries
