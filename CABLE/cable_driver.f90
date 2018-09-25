! cable_driver.f90
!
! Netcdf offline driver for CABLE land surface scheme, May 2009.
! Gab Abramowitz, University of New South Wales, gabsun@gmail.com
!
! Thanks to Peter Isaac (Monash) for introducing the namelist file (Oct 2007)
!
! Modified and added casaCNP, Bernard Pak, Sep 2010.
!
! Please send bug reports to Bernard.Pak@csiro.au
!
PROGRAM offline_driver
  USE cbm_module
  USE define_dimensions, ONLY:r_1,i_d,ms,mp,mvtype,mstype
  USE define_types
  USE io_variables, ONLY: logn,filename,leaps, &
       verbose, fixedCO2,output,check,patchout,patch_type
  USE input_module, ONLY: open_met_file,load_parameters, &
       get_met_data,close_met_file
  USE output_module, ONLY: create_restart,open_output_file, &
       write_output,close_output_file
  ! new modules related to casaCNP
  USE casa_cnp   ! icycle declared in casadimension which is used in casa_cnp
                 ! casafile declared in casavariable, also used in casa_cnp
  IMPLICIT NONE
  INTEGER(i_d)          :: kend ! no. of time steps in run
  TYPE (air_type)       :: air  ! air property variables
  TYPE (bgc_pool_type)  :: bgc  ! carbon pool variables
  TYPE (canopy_type)    :: canopy ! vegetation variables
  TYPE (met_type)       :: met  ! met input variables
  TYPE (balances_type)  :: bal  ! energy and water balance variables
  TYPE (radiation_type) :: rad  ! radiation variables
  TYPE (roughness_type) :: rough ! roughness varibles
  TYPE (soil_parameter_type) :: soil ! soil parameters	
  TYPE (soil_snow_type) :: ssoil ! soil and snow variables
  TYPE (sum_flux_type)  :: sum_flux ! cumulative flux variables
  TYPE (veg_parameter_type) :: veg  ! vegetation parameters	 
  TYPE (casa_biome)     :: casabiome
  TYPE (casa_pool)      :: casapool
  TYPE (casa_flux)      :: casaflux
  TYPE (casa_met)       :: casamet
  TYPE (casa_balance)   :: casabal
  TYPE (phen_variable)  :: phen 
  REAL(r_1)         :: dels ! time step size in seconds
  INTEGER(i_d)      :: kstart ! start of simulation #
  INTEGER(i_d)      :: ktau    ! index of time step = 1 ..  kend
  LOGICAL    :: vegparmnew   ! using new format input file (BP dec 2007)
  LOGICAL    :: spinup ! should the model spinup to soil state equilibrium?
  LOGICAL    :: spinConv ! has spinup converged?
  LOGICAL    :: spincasainput
                   ! TRUE: input required to spin casacnp wil be saved;
                   ! FALSE: input will be read in to spin casacnp 1000 years
  LOGICAL    :: spincasa ! TRUE: casacnp will spin mloop times,
                         ! FALSE: no spin up
  REAL(r_1)  :: delsoilM ! allowed variation in soil moisture for spin up
  REAL(r_1)  :: delsoilT ! allowed variation in soil temperature for spin up
  REAL(r_1),POINTER  :: soilMtemp(:,:) ! temporary storage for spin up
  REAL(r_1),POINTER  :: soilTtemp(:,:) ! temporary storage for spin up
  INTEGER(i_d) :: tstep  ! time step counter for spinup
  INTEGER(i_d) :: mloop  ! # spinup loops for casaCNP
  NAMELIST/CABLE/filename, &
                 vegparmnew,&
                 spinup,delsoilM,delsoilT,&
                 output,&
                 patchout,&
                 check,&
                 verbose,leaps,logn,fixedCO2, &
                 spincasainput,     &
                 spincasa,          &
                 icycle,            &
                 casafile
  !===================================================================!
  ! Open, read and close the namelist file.
  OPEN(10,FILE='cable.nml')
  READ(10,NML=CABLE)
  CLOSE(10)
  !=====================================================================!
  ! Open log file:
  OPEN(logn,FILE=filename%log)
  ! Open met data and get site information from netcdf file.
  ! This retrieves time step size, number of timesteps, starting date,
  ! latitudes, longitudes, number of sites. 
  CALL open_met_file(dels,kend,spinup)

  ! Checks where parameters and initialisations should be loaded from.
  ! If they can be found in either the met file or restart file, they will 
  ! load from there, with the met file taking precedence. Otherwise, they'll
  ! be chosen from a coarse global grid of veg and soil types, based on 
  ! the lat/lon coordinates. Allocation of CABLE's main variables also here.
  CALL load_parameters(met,air,ssoil,veg,bgc,soil,canopy, &
       rough,rad,sum_flux,bal,logn,vegparmnew)

  ! Open output file:
  CALL open_output_file(dels,soil,veg,bgc,rough)

!  print *, 'mp mstype mvtype = ',mp,mstype,mvtype
  if(icycle>0) then
    mloop=500
    kstart=1
    call alloc_casavariable(casabiome,casapool,casaflux,casamet,casabal,mp)
    call alloc_phenvariable(phen,mp)

    call casa_readpoint(veg,soil,casaflux,casamet,rad)
    call casa_readbiome(veg,soil,casabiome,casapool,casaflux,casamet,phen)
    call casa_readphen(veg,casamet,phen)
    call casa_init(casapool,casabal,veg)
!  print *, 'mp mstype mvtype = ',mp,mstype,mvtype
!    if (spincasa) then
!      call spincasacnp(dels,kstart,kend,mloop,veg,soil,casabiome,casapool, &
!                       casaflux,casamet,casabal,phen)
!    endif
  endif

  kstart = 1
  tstep = 0          ! initialise spinup time step
  spinConv = .FALSE. ! initialise spinup convergence variable
  ! spinup loop:
  DO
    ! time step loop:
    DO ktau = kstart, kend ! time step loop
      ! increment total timstep counter
      tstep = tstep + 1

      ! Get met data and LAI, set time variables.
      ! Rainfall input may be augmented for spinup purposes:
      CALL get_met_data(spinup,spinConv,ktau,met,soil,rad,veg,kend,dels) 
        
      ! CALL land surface scheme for this timestep, all grid points:
      CALL cbm(tstep, kstart, kend, dels, air, bgc, canopy, met, &
             & bal, rad, rough, soil, ssoil, sum_flux, veg)
!      CALL cbm(tstep, kstart, kend, dels, air, bgc, canopy, met, &
!             & bal, rad, rough, soil, ssoil, sum_flux, veg, mvtype, mstype)

      if(icycle >0) then
        call bgcdriver(ktau,kstart,kend,dels,met,ssoil,canopy,veg,soil, &
                  casabiome,casapool,casaflux,casamet,casabal,phen)
      endif 

      ! sumcflux is pulled out of subroutine cbm
      ! so that casaCNP can be called before adding the fluxes (Feb 2008, YP)
      CALL sumcflux(ktau, kstart, kend, dels, bgc, canopy,  &
                  & soil, ssoil, sum_flux, veg, met, casaflux)

      ! Write time step's output to file if either: we're not spinning up 
      ! or we're spinning up and the spinup has converged:
      IF((.NOT.spinup).OR.(spinup.AND.spinConv)) CALL write_output &
             & (ktau,dels,met,canopy,ssoil,rad,bal,air,soil,veg)

    END DO
    ! see if spinup (if conducting one) has converged:
    IF(spinup.AND..NOT.spinConv) THEN
      ! Write to screen and log file:
      WRITE(*,'(A18,I3,A24)') ' Spinning up: run ',INT(tstep/kend), &
            ' of data set complete...'
      WRITE(logn,'(A18,I3,A24)') ' Spinning up: run ',INT(tstep/kend), &
            ' of data set complete...'
      ! IF not 1st run through whole dataset:
      IF(INT(tstep/kend)>1) THEN 
        ! evaluate spinup
        IF(ANY(ABS(ssoil%wb-soilMtemp)>delsoilM).OR. &
           ANY(ABS(ssoil%tgg-soilTtemp)>delsoilT)) THEN
          ! No complete convergence yet
        ELSE ! spinup has converged
          spinConv = .TRUE.
          ! Write to screen and log file:
          WRITE(*,'(A33)') ' Spinup has converged - final run'
          WRITE(logn,'(A52)') &
               ' Spinup has converged - final run - writing all data'
          WRITE(logn,'(A37,F7.5,A28)') &
               ' Criteria: Change in soil moisture < ', &
               delsoilM, ' in any layer over whole run'
          WRITE(logn,'(A40,F7.5,A28)' ) & 
               '           Change in soil temperature < ', &
               delsoilT, ' in any layer over whole run'
        END IF
      ELSE ! allocate variables for storage
        ALLOCATE(soilMtemp(mp,ms), &
               & soilTtemp(mp,ms))
      END IF
      ! store soil moisture and temperature
      soilTtemp = ssoil%tgg
      soilMtemp = REAL(ssoil%wb,r_1)
    ELSE
      ! if not spinning up, or spin up has converged, exit:
      EXIT
    END IF
  END DO

  ! Write restart file if requested:
  IF(output%restart) CALL create_restart(logn,ktau,dels,&
       soil,veg,ssoil,canopy,rough,rad,bgc,bal)
!  IF(output%restart) CALL create_restart(logn,ktau,dels,&
!       soil,veg,ssoil,canopy,rough,rad,bgc,bal,mvtype,mstype)
 print *, sum_flux%sumpn, sum_flux%sumrp, sum_flux%sumrd, bal%wbal_tot, bal%ebal_tot

  ! Write final pool sizes for casaCNP
  IF (icycle>0) call casa_poolout(ktau,veg,soil,casabiome,casapool,casaflux, &
                                  casamet,casabal,phen)
  
  ! Close met data input file:
  CALL close_met_file
  ! Close output file and deallocate main variables:
  CALL close_output_file(bal, air, &
       bgc, canopy, met, rad, rough, soil, ssoil, sum_flux, veg)

  ! Close log file
  CLOSE(logn)

END PROGRAM offline_driver

