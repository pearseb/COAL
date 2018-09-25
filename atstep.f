c Fixing a bug in the handling of netCDF errors.
c SJP 2009/04/27
c
c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c (1) Fixing syntax errors identified by the g95 Fortran compiler.
c (2) Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c CABLE write output routine is called to avoid conflicting parallelism
c AJA 2009/04/02
c
c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/16
c
c Added IMPLICIT NONE statement, plus variable declarations for TDT, NDADD,
c NRUNDAY, IDAYPSAV, MSTEPSAV, NHOURS, NMINS, DTICE, IALBBAR, ICLHBAR, ICLMBAR,
c ICLLBAR and ICLBAR.
c SJP 2007/05/28
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Modified so that HMREAD1 is no longer called, as a constant mixed-layer ocean
c depth of 100m is now used when QFLUX=F.
c SJP 2003/05/29
c
c Modified so that SEMIIS, rather than SEMII, called for machine type 'ALPH'.
c SJP 2003/04/02
c
c Writes to f30.his removed for machine type ALPH in order to increase
c execution speed.
c SJP 2002/02/15
c
c Declaration added for IERR.
c SJP 2001/12/10
c
c Changes to add machine type 'ALPH'.
c SJP 2001/11/22
c
c Updated from v2.4 to v3 of netCDF and "include 'netcdf.inc'" added.
c SJP 2001/11/21
c
c $Log: atstep.f,v $
c Revision 1.29  2000/08/16 02:59:26  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.28  2000/06/20 02:08:35  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.27  1999/06/16 06:21:54  rot032
c HBG changes to V5-3
c
c Revision 1.26  1999/05/20 06:23:58  rot032
c HBG changes to V5-2
c
c Revision 1.25  1998/12/10  01:08:03  ldr
c Merge HBG changes with MRD changes.
c
c Revision 1.24  1998/07/29  07:41:16  ldr
c Correct formatting in f30.his file to allow year > 999.
c
c Revision 1.23.1.1  1998/12/10  00:56:01  ldr
c HBG changes to V5-1-21
c
c Revision 1.23  1997/12/03  05:35:42  ldr
c Put IGW's qflux code with variable MLO depth into main version of model.
c
c Revision 1.22  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.21  1997/08/13  00:02:58  mrd
c Close history files properly when nsstop is set.
c
c Revision 1.20  1997/07/24  05:59:17  ldr
c Introduce logical variable called "start".
c
c Revision 1.19  1996/10/24  01:03:25  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.18  1996/06/13  02:08:44  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.17  1996/03/21  03:19:14  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.16  1995/09/07  06:41:14  ldr
c Call semii not semiis in SGI version unless NPROC > 1.
c
c Revision 1.16  1995/09/07  06:38:39  ldr
c Call semii not semii in SGI version unless NPROC > 1.
c
c Revision 1.15  1995/02/27  01:40:34  ldr
c Bugfixes received from HBG in response to problems on Fujitsu.
c
c Revision 1.14  1994/08/08  17:23:24  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.13  94/07/11  12:31:31  ldr
c Rearrange variables in common block /atdat to avoid misalignment warnings
c at 64 bit on SGI.
c 
c Revision 1.12  94/03/23  09:27:24  mrd
c Removed the extra zenith angle calculation by avsol by moving call to
c zenith in radin and adding zenith angle as an argument to surfset.
c 
c Revision 1.11  93/12/17  12:05:38  ldr
c Merge of HBG and LDR changes.
c 
c Revision 1.10  93/12/06  16:55:38  ldr
c Changes to V4-4-45l to get T63 running on SGI (Seca).
c 
c Revision 1.9.1.1  93/12/17  11:51:55  ldr
c Changes to V4-4-45l from HBG for coupled model
c 
c Revision 1.9  93/11/29  11:38:46  ldr
c Changes to V4-4-32l from HBG for coupled model
c 
c Revision 1.8  93/11/03  11:53:58  ldr
c Merge of HBG with LDR changes since V4-4-24l.
c 
c Revision 1.7  93/10/15  14:17:14  ldr
c Changes to create special 'FUJI' version of model with fast FFTs and
c Legendre transform.
c 
c Revision 1.6.1.1  93/11/03  11:44:43  ldr
c Changes to V4-4-24l from HBG for coupled model
c
c     INPUT/OUTPUT:
c     Input:   from common/fewflags in FEWFLAGS.f
c                  idyn - if true, and semice-T, use dynamical sea-ice model
c                  ispec - if true, call routine przav at end of each day to
c                           generate spectral amplitude of fields, summed over
c                           l space
c                  ltrace - if true, calls Langrangian tracer routine, tracera
c                  nsstop - no of time steps for model run
c                  savefcor - if true, save stats for flux-corrections for
c                             coupled ocean model
c                  saveglmean - save global means
c                  sltrace - if T compute tracer transport by semi-
c                            Langrangian advection and by mixing and
c                            convection in the vertical
c
c              from common/glmean in GLMEAN.f
c                  albbar - mean albedo      clbar - mean total cloud
c                  clhbar - mean high cloud  cllbar - mean low cloud
c                  clmbar - mean mid cloud   evbar - mean evaporation
c                  hfbar - mean hflux        rnbar - mean rain
c                  tsbar - mean surf temp    tempbar - mean temp (all levels)
c                  sume - mean energy loss   sndbar - mean snow depth
c                  clforbar - mean net cloud radiative forcing
c                  sidbar - mean sea ice depth                 
c                  gsbal - mean surface energy balance               
c
c              from common/nest in NEST.f
c                  nestflag   - if true, call routine nestarc to write data for
c                               nesting the limited area model
c                  nest_start - minutes before first call to nestarc
c
c              from common/timex in TIMEX.f
c                  dt    - time step in seconds
c                  iyear - year counter for model run
c                  month - month counter
c                  mstep - timestep in minutes 
c
c              from common/worka in WORKA.f
c                  tbar - global mean temp at each model level
c
c              from arguments
c                  nato - no. of atmos time steps per ocean time step
c
c     Output:  from common/atm2oc in common ATM2OC
c
c     In/Out:  from common/atdat in this subroutine
c                  iday   - logical variable, T if end of day  
c                  tdt    - 2 X timestep (dt)
c
c              from common/fewflags in FEWFLAGS.f
c                  glmean_interval - interval in mins between prints of global
c                                    means to standard output
c                  iener   - if true, generate detailed energy diagnostics
c                  ifwds   - default F, only used when doing forward start
c                            from new restart file
c                  lcouple -flag to run model in coupled mode, default F  
c
c              from common/glmean in GLMEAN.f
c                  global_tmax - global maximum temperature
c                  global_tmin - global minimum temperature
c                  global_umax - global maximum zonal wind speed
c                  rgbar - mean net long wave heating at ground
c                  rtbar - mean long wave at top of atmosphere
c                  sgbar - mean short wave at ground
c                  sinbar - mean solar in at top of atmosphere
c                  soubar - mean solar out at top of atmosphere
c                  totkei - total kinetic energy
c                  totkzi - total kinetic energy, zonal average
c                  tot_water - global integral of moisture
c
c              from common/nest in NEST.f
c                  nest_interval - mins between calls to nestarc
c
c              from comon/timex in TIMEX.f
c                  ldays - month count of days
c                  mins  - current model time in mins
c                  ndays - day counter for model run
c                  nsteps- time step counter
c
c              from arguments
c                  ijdyn   - ice dynamics time step factor
c                  exptype - restart file header
c                  
c  
      subroutine atstep(exptyp,ijdyn,nato)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      character exptyp*50
      integer ijdyn
      integer nato

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f'
      include 'GLMEAN.f'     ! Common block glmean contains global means
      include 'HIST.f'
      include 'NEST.f'
      include 'WORKA.f'
      include 'TIMEX.f'
      include 'netcdf.inc'

      real tdt
      integer ndadd,nrunday,idaypsav,mstepsav
      logical iday,endofmonth
      common/atdat/tdt,ndadd,nrunday,idaypsav,mstepsav,iday,endofmonth

C Local work arrays and variables
      integer ierr, nhours, ialbbar, iclhbar, iclmbar, icllbar, iclbar
      real dtice

C Local data, functions etc

C Start code : ----------------------------------------------------------

c     CALL TIMER('ATSTEP  ',1)
C****
C**** PERFORMS 1 STEP OF THE Nl LEVEL DAR CLIMATE MODEL
C****
c     RESET ICE STRESSES AT START OF EACH DAY (FOR PURE RESTARTS)
      if (iday) call icestat(2)
c     Set Mixed layer depths for ice model, q-fluxes
CSJP      if (iday.and..not.qflux) call hmread1
c     SET UP SOURCES FOR SEMI-LAGRANGIAN TRACER MODELLING (LDR 12/92)
      if (sltrace ) call source
C*    REPLACE T-1 FIELDS BY T FIELDS IF REQUIRED
      if(ifwds)call transf

  111 continue

      if(ncepstrt.gt.0)then ! ncepstrt =2,1
c..  Half leapfrog timestep for startup from NCEP
        dt=0.5*mstep*60.0
        tdt=2.0*dt
      else                  ! ncepstrt =0,-1,-2...
c..  Otherwise usual leapfrog timestep
        dt=mstep*60.0
        tdt=2.0*dt
      endif

C*    SET UP U,V HARMONICS
c     CALL TIMER('UVHARM  ',3)
      call uvharm
c     CALL TIMER('UVHARM  ',4)
C*    SET ARRAYS TO ZERO BEFORE PHYSICS STEP
c     CALL TIMER('ZEROGI 1',3)
      call zerogi(1,exptyp)
c     CALL TIMER('ZEROGI 1',4)
      call ocforce(1)
C*    PERFORM PHYSICAL ADJUSTMENTS
c     CALL TIMER('PHYS    ',3)

CSJP  Former machine dependence at this point
      call phys(tdt, exptyp)

      call physgm
c     CALL TIMER('PHYS    ',4)
C*    SET ARRAYS TO ZERO & MANIPULATE DATA BEFORE DYNAMICS STEP
c     CALL TIMER('ZEROGI 2',3)
      call zerogi(2,exptyp)
c     CALL TIMER('ZEROGI 2',4)
C*    PERFORM THE DYNAMICS STEP
c     CALL TIMER('DYNM    ',3)

CSJP  Former machine dependence at this point
      call dynm(tdt)

c     CALL TIMER('DYNM    ',4)
C*    COMPUTE ENERGY AND SPECTRAL AMPLITUDE DIAGNOSTICS
      if(mod(mins+int(mstep, 8),int(glmean_interval, 8)).eq.0_8)then
        call energy(iener,1)
      endif
      if(ispec.and.mod(mins+int(mstep, 8),1440_8).eq.0_8)call specam
C*    ADD THE LINEAR SPECTRAL TENDANCIES
c     CALL TIMER('LINEAR  ',3)
      call linear
c     CALL TIMER('LINEAR  ',4)
C*    ADD THE ASSELIN FILTER FINAL COMPONENTS
C     CALL TIMER('ASSEL   ',3)
      call assel
C     CALL TIMER('ASSEL   ',4)

      call uvreal(1,global_umax)
c John McGregor's Semi-Lagrangian moisture transport
c-----------
      call jmcgslt(tdt,tot_water,exptyp)
c-----------
      ifwds=.false.
      if(ncepstrt.ge.0)call initax
C*    SET UP THE SEMI-IMPLICIT MATRICES
c     CALL TIMER('MATSET  ',3)
      call matset(tdt)
c     CALL TIMER('MATSET  ',4)
C*    PERFORM THE SEMI-IMPLICIT TIME INTEGRATION
c     CALL TIMER('SEMII   ',3)

CSJP  Former machine dependence at this point
      call semiis(tdt)

c     CALL TIMER('SEMII   ',4)
C*    ADD FORWARD IMPLICIT DIFFUSION
C     CALL TIMER('DIFFN   ',3)
      call diffn(tdt)
C     CALL TIMER('DIFFN   ',4)

      ncepstrt=ncepstrt-1
      if(ncepstrt.gt.0)go to 111

      start=.false.
      nsteps=nsteps+1
      mins=mins+mstep
      nhours=int(mins/60_8, 4)
      ndays=nhours/24
      nhours=nhours-24*ndays
c       write(6,437)ndays,nhours,pslbar,gsbal,(tbar(k),k=1,nl)

      ijdyn=1
      if(idyn)then
c do ice dynamics every hour to coincide with coupled model
       dtice=3600.
       if(mw.eq.64)dtice=1800.0 ! T63 model
       ijdyn=int((dtice+1.0e-5)/dt)
c.. check that the ice dynamics step coincides with the coupled model
c.. ocean step (to match stress calculations/averaging etc).
       if(lcouple.and.(ijdyn.ne.nato))then
         ijdyn=nato
c        print *,'ice dynamics step changed to match ocean step'
c        print *,'dtice=',dt*nato
       endif
       dtice=ijdyn*dt
       if(mod(nsteps,ijdyn).eq.0)then
         call ocicurr
         call icedrive(ldays,dtice,ijdyn)
         call icestat(3)
       endif
      endif

c     Complete the calculation of the stresses on the ocean -
c     both atmospheric (computed in Phys/Radin) and from the
c     dynamical ice model (computed in Icedrive/Dynice/Icefree)

      if(mod(nsteps,ijdyn).eq.0)then
        if(mw.eq.64)then
          call ocntauMOM2(ijdyn) ! T63 model ==> MOM2 ocean in use
        else 
          call ocntau(ijdyn)
        endif
      endif

c     Compute ocean forcing terms (heat & salinity)
c     Relocate runoff to ocean points and compute salinity forcing
c     Do after ice dynamics to allow for thin ice melt
      if((.not.lcouple.and.savefcor).or.lcouple)call ocforce(2)

c TRACER ROUTINE                                               

      if(ltrace) then
        call traceout(0,nsteps)
        call tracera
      endif

      ! calling the CABLE write output file here to avoid conflicting
      ! with parallelism
      ! sending nstepsa for ktau, and mstep for kend
      CALL CABLE_write_output(nstepsa,mstep) 
 

c PRINT VARIOUS GLOBAL MEANS AT PREDEFINED INTERVALS (DEFAULT 360 MINS)
c DATA IS ALSO WRITTEN TO UNIT 59.

      if(mod(mins,int(glmean_interval, 8)).eq.0_8)then
        call uvreal(2,global_umax)
         ialbbar=nint(100.*albbar)
         iclhbar=nint(100.*clhbar)
         iclmbar=nint(100.*clmbar)
         icllbar=nint(100.*cllbar)
         iclbar=nint(100.*clbar)
         write(6,437)ndays,nhours,pslbar,gsbal,tbar(nl),tsbar,
     &       nint(global_tmin),nint(global_tmax),nint(totkzi),
     &    nint(totkei),nint(rgbar),nint(rtbar),nint(sgbar),nint(sinbar),
     &    nint(soubar),clforbar,ialbbar,iclhbar,iclmbar,icllbar,iclbar,
     &       evbar*1440./mstep, rnbar*1440./mstep, hfbar,tot_water,
     &        sndbar,sidbar,tempbar,sume*365.,nint(global_umax)
 437    format(i9,1x,i2,1x,f6.1,f5.1,2f6.1,9i4,f5.1,5i3,2f4.1,3f5.1,
     &        f5.2,f6.1,f4.1,i4)
        if(saveglmean)then
chbg      write(59,439)ndays,nhours,pslbar,gsbal,tbar(nl),tsbar,
          write(59,439)ndays,nhours,pslbar,gsbal,tbar,tsbar,
     &       global_tmin,global_tmax,totkzi,totkei,rgbar,rtbar,sgbar,
     &       sinbar,soubar,clforbar,albbar,clhbar,clmbar,cllbar,clbar,
     &       evbar*1440./mstep, rnbar*1440./mstep, hfbar,tot_water,
     &        sndbar,sidbar,tempbar,sume*365.,global_umax
        endif
c439    format(i6,1x,i2,1x,14f8.3,5f7.4,2f6.2,4f7.3,f8.3,2f6.2)
 439    format(i6,1x,i2,1x,(7f20.12))
      endif

c If NSTEPS.NE.0, print energy diagnostics and stop after NSTEPS
      
      if(nsteps.eq.nsstop)then
        call energy(.false.,1) !To calculate totkzi and totkei
        write(6,438)totkzi,totkei
 438    format(' Total kzon = ',e17.11,' Total keddy = ',e17.11)

        if (savehist(1)) then
          ierr = nf_close(histid(1))
          if (ierr .ne. nf_noerr) stop "***  netCDF error in atstep"
        end if
        if (savehist(2)) then
          ierr = nf_close(histid(2))
          if (ierr .ne. nf_noerr) stop "***  netCDF error in atstep"
        end if

        write(6,*)'Model stopped after ',nsstop,' steps'
        stop
      endif

      iday=(mins/1440)*1440.eq.mins

c AT PRESET INTERVALS WRITE ARCHIVE FOR RUNNING NESTED MODEL

      if(nestflag.and.mod(mins-int(nest_start, 8),
     &                    int(nest_interval, 8)).eq.0_8)then
        write(0,*)'nestarc called at day ',ndays,'  hour ',nhours
        call nestarc(mins,ldays)
      endif

C****
C**** END OF ATMOSPHERIC TIMESTEP
C****
c     CALL TIMER('ATSTEP  ',2)

      return
      end
