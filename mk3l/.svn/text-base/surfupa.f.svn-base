c surface downward shortwave sgdn is passed from radin to be used by CABLE
c AJA 2009/03/13
c
c SURF1.f is made available for output variables like tgf and tgg to be
c updated by CABLE
c AJA 2009/03/06
c
c runoff argument is sent as an output to radin.f
c AJA 2009/02/04
c
c MASIV3.f is included to make the rgsav available
c AJA 2009/02/03
c
c LSMI.f is included to make the imsl array available
c RADISW.f and RDPARM.f are included to retrieve CO2 concentration for cable
c TIMEX.f is included to retreive time details for cable.
c vmod : the windspeed 
c passed as input arguments.
c AJA 2009/02/02
c
c the CABLE execution step is called when lsm_type =  "cable"
c AJA 2009/01/30 
c
C the old surface model updln1 removed
c nsib comparion is modified to lsm_type = "nsib "
c the if-else statement for lsm_type is modified to case statement
c AJA 2009/01/22
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: surfupa.f,v $
c Revision 1.49  2001/02/22 05:34:45  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.48  2001/02/12 05:39:56  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.47  2000/06/20 02:27:11  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.46  1999/05/20 06:23:57  rot032
c HBG changes to V5-2
c
c Revision 1.45  1998/12/10  00:56:00  ldr
c HBG changes to V5-1-21
c
c Revision 1.44  1997/12/23  04:09:43  ldr
c Remove redundant line for qflux case.
c
c Revision 1.43  1997/12/23  00:23:39  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.42  1997/12/17  23:23:03  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.41  1997/12/03  05:35:42  ldr
c Put IGW's qflux code with variable MLO depth into main version of model.
c
c Revision 1.40.2.1  1997/12/19  02:03:19  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.40  1996/10/24  01:03:18  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.39  1996/06/13  02:08:29  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.38  1996/03/21  03:19:08  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.37  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.33.2.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.36  1994/08/08  17:22:46  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.35  94/08/08  13:16:29  ldr
c Push the debug prints down into individual routines.
c 
c Revision 1.34  94/08/04  16:56:40  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.33  94/03/30  12:35:18  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.32  94/01/17  11:29:31  ldr
c Changes from Hal interpolate qfluxes and make them implicit. Also move
c read of qfluxes to datard so it is done every month.
c 
c Revision 1.31  94/01/04  17:18:25  ldr
c A little fix so that new qflux runs are initialized properly.
c 
c Revision 1.30  93/12/17  15:33:58  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.29  93/12/17  11:51:48  ldr
c Changes to V4-4-45l from HBG for coupled model
c 
c Revision 1.28  93/11/29  16:09:41  ldr
c Corrected length of common block hm.
c
c     INPUT/OUTPUT
c     Input:   from common/fewflags in FEWFLAGS.f
c                  debug    - flag to control column debugging
c                  insdebug - flag to control hemisphere debugging
c                  lgdebug  - flag to control latitude debugging
c                  mgdebug  - flag to control longitude debugging
c                  lsm_type - if "nsib ", use "New SIB" land surface scheme
c                  qflux    - set to T for qflux run, default F
c                  semice   - if T, use the Semtner sea-ice model
c
c              from common/freeze in FREEZE.f
c                  tfi - temperature at bottom of ice
c
c              from common/logimsl in LOGIMSL.f
c                  cice - seaice      mlo - mixed layer ocean 
c                  land - logical variable for surface type
c
c              from common/qflxdat in QFLXDAT.f
c                  gm1cur,gm2cur - qflux interpolate data
c
c              from common/timex in TIMEX.f
c                  dt     - time step in seconds
c                  mstep  - timestep in minutes 
c                  nsteps - time step counter
c                  ratlm  - fraction of month elapsed, runs 0-1
c
c              from arguments
c                  cie - internal conductivity of ice/depth
c                  gamms - density*specific heat*depth
c                  tintp - temperature at sea points
c                  z4 - grid elevations [non-spectral]
c                  vmod - surface windspeed (m/s) 
c                  zmin - height of lowest model level
c                  sgdn - surface downward shortwave to be used by CABLE
c
c     In/Out:  from common/qflxdat in QFLXDAT.f
c                  newqfluxrun - if T, use weight of zero at start of new 
c                                qflux run
c                  ochfa - ocean heat flux average
c
c              from common/relate in RELATE.f
c                  zfrm - zonal average MLO diffusive heating term
c
c              from arguments
c                  blg,hfrm
c                  cls - factor for latent heat of fusion
c                  fg    - sensible heat flux  eg - latent heat flux
c                  degdt - derivative of eg with respect to temp
c                  degdw  - derivative of eg with respect to soil wetness
c                  dfgdt - derivative of fg with respect to temp
c                  ns - hemisphere index      lg - latitude index 
c                  pg  - surface pressure(mbs)  tg - surface temperature
c                  qtg - mixing ratio      ttg - physical temperature
c                  rg - net long wave heating at ground
c                  scalev - scaling evaporation
c                  sg - solar absorbed in ground
c                  snowd - snow depth 
c                  tb2   - temperature of the second layer at the 
c                          previous time step
c                  tb3   - temperature of the third layer at the 
c                          previous time step
c                  tddd  - total wet evaporation
c                  totpev - potential evaporation in mm/timestep
c                  wb    - soil moisture of the bottom layer at the 
c                          previous time step
c                  wg    - soil moisture of the top layer at the 
c                          previous time step
c                  wgp   - soil moisture of the top layer at the 
c                          current time step
c     Output:  from argument
c                  runoff - ground runoff
c                  perc   - soil moisture percolation
c
      subroutine surfupa(lg,degdt,dfgdt,degdw,cie,gamms,sg,
     &                   tintp,rg,blg,z4,wg,wg2,snowd,rcondx,
     &                   pg,ttg,als,qtg,cls,vmod,zmin,rpreci, 
     &                   eg,fg,tg,tb2,tb3,wgp,tddd,totpev,scalev,
     &                   hfrm,he,mcmax,runoff,perc,cabTscrn)   !runoff output to radin.f

      implicit none

!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /RELATE/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'RDPARM.f' ! included to use RADISW.f

C Argument list
      integer lg
      real degdt(ln2)
      real dfgdt(ln2)
      real degdw(ln2)
      real cie(ln2)
      real gamms(ln2)
      real sg(ln2)
      real tintp(ln2)
      real rg(ln2)
      real blg(ln2)
      real z4(ln2)
      real wg(ln2)
      real wg2(ln2) 
      real snowd(ln2)
      real rcondx(ln2)
      real pg(ln2)
      real ttg(ln2,nl)
      real als(ln2)
      real qtg(ln2,nl)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real cls(ln2)
      real eg(ln2)
      real fg(ln2)
      real tg(ln2)
      real tb2(ln2)
      real tb3(ln2)
      real wgp(ln2)
      real tddd(ln2)
      real totpev(ln2)
      real scalev(ln2)
      real hfrm(ln2)
      real he(ln2)
      real mcmax(ln2)
      real runoff(ln2)
      real perc(ln2)
      real vmod(ln2)  
c temporary variables to hold the SURF1.f common block for output from CABLE
      real tgg_temp   (ln2) 
      real tgf_temp   (ln2)
      real mc_temp    (ln2)
      real osnowd_temp(ln2)
      real snage_temp (ln2)
      real ssdnn_temp (ln2)
      real gflux_temp (ln2)
      real sgflux_temp(ln2)
      real wb_temp    (ln2,ms)
      real wbice_temp (ln2,ms)
      real tggsl_temp (ln2,ms)
      real tggsn_temp (ln2,3)
      real smass_temp (ln2,3)
      real ssdn3_temp (ln2,3)
      integer isflag_temp(ln2)
c end temporary variables for SURF1.f common block
      real rpreci(ln2) 
      real zmin(ln2)
      real cabTscrn(ln2)
C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'LOGIMSL.f'
      include 'RELATE.f'

C Global data blocks
      include 'FEWFLAGS.f'
      include 'TIMEX.f'  ! made available to retrieve time details for cable
      include 'RADISW.f' ! made avaialble to retreive CO2 concentration
      include "LSMI.f"   ! the array imsl is made available
      include "MASIV3.f" ! rgsav is made available
      include "SURF1.f"  ! some output variables (tgg,tgf..) are made available
      include "GIANT5.f" ! phnt made available for za calculations in CABLE
      include "CNSTA.f"  ! algf made available to be passed to CABLE for rough%za calcuations

C Local work arrays and variables
      logical flag1

      integer mg
      integer ns

C Local data, functions etc

C Start code : ----------------------------------------------------------

C****
C**** ROUTINE TO DO THE FIRST PART OF THE UPDATE OF THE SURFACE VALUES
C**** (Note: final updating done in surfupb.f)
C****

C****
C**   UPDATE LAND POINTS.
C****
!      IF (lg .eq. 23) THEN
!        WRITE(45,*) 'surfupa  = ',wb(4,1,lg),wbice(4,1,lg)
!        WRITE(45,
!     &     '(3f10.2,8i10,2f10.2,2x,a8,f10.2,i10,24f10.2,i10,6f10.2)')
!     &  sg(4),rcondx(4),pg(4),iyear,month,kdays,mins,
!     &  ln2,lat,nl,nstepsa,vmod(4),zmin(4),qcloud,rpreci(4),mstep,
!     &  rg(4),scalev(4),tg(4),fg(4),eg(4),tb2(4),tb3(4),totpev(4),
!     &  mcmax(4),tddd(4),tgg(4,lg),tgf(4,lg),mc(4,lg),osnowd(4,lg),
!     &  snage(4,lg),ssdnn(4,lg),gflux(4,lg),sgflux(4,lg),tggsl(4,1,lg),
!     &  tggsn(4,1,lg),wb(4,1,lg),wbice(4,1,lg),smass(4,1,lg),
!     &  ssdn3(4,1,lg),
!     &  isflag(4,lg),als(4),snowd(4),runoff(4),wg(4),wg2(4),cabTscrn(4)
!      ENDIF

      select case (lsm_type) 
      case ("nsib ")
        call surfa(lg,sg,rg,
     &           snowd,rcondx,pg,ttg,qtg,qlg,qfg,
     &           cls,eg,fg,tg,tb2,
     &           tb3,totpev,scalev,wg,wg2,tddd,he,mcmax)

      case ("cable") ! the CABLE execution step part is called
        tgg_temp   (:) = tgg   (:,lg)
        tgf_temp   (:) = tgf   (:,lg)
        mc_temp    (:) = mc    (:,lg)
        osnowd_temp(:) = osnowd(:,lg)
        snage_temp (:) = snage (:,lg)
        ssdnn_temp (:) = ssdnn (:,lg)
        gflux_temp (:) = gflux (:,lg)
        sgflux_temp(:) = sgflux(:,lg)
        wb_temp   (:,:) = wb   (:,:,lg)
        wbice_temp(:,:) = wbice(:,:,lg)
        tggsl_temp(:,:) = tggsl(:,:,lg)
        tggsn_temp(:,:) = tggsn(:,:,lg)
        smass_temp(:,:) = smass(:,:,lg)
        ssdn3_temp(:,:) = ssdn3(:,:,lg)
        isflag_temp(:) = isflag(:,lg)
        call CABLE_run(lg,sg,rgsav,rcondx,pg,ttg,iyear,month,kdays,mins,!INPUTS
     &       ln2,lat,nl,imsl,nstepsa,vmod,zmin,qcloud,qtg,rpreci,mstep, !INPUTS
     &              rg,scalev,tg,fg,eg,tb2,tb3,totpev,mcmax,tddd,       !OUTPUTS
     &              tgg_temp,tgf_temp,mc_temp,osnowd_temp,snage_temp,   !OUTPUTS
     &              ssdnn_temp,gflux_temp,sgflux_temp,tggsl_temp,       !OUTPUTS
     &              tggsn_temp,wb_temp,wbice_temp,smass_temp,ssdn3_temp,!OUTPUTS
     &              isflag_temp,als,snowd,runoff,perc,wg,wg2,cabTscrn)  !OUTPUTS
        write(38,'(a21,i3,128f9.3)') 'lg, cable_snowd: ', lg, snowd
        ! mapping the one dimensional output to the global arrays
        tgg   (:,lg) = tgg_temp   (:)
        tgf   (:,lg) = tgf_temp   (:)
        mc    (:,lg) = mc_temp    (:)
        osnowd(:,lg) = osnowd_temp(:)
        snage (:,lg) = snage_temp (:)
        ssdnn (:,lg) = ssdnn_temp (:)
        gflux (:,lg) = gflux_temp (:)
        sgflux(:,lg) = sgflux_temp(:)
        wb   (:,:,lg) = wb_temp   (:,:)
        wbice(:,:,lg) = wbice_temp(:,:)
        tggsl(:,:,lg) = tggsl_temp(:,:)
        tggsn(:,:,lg) = tggsn_temp(:,:)
        smass(:,:,lg) = smass_temp(:,:)
        ssdn3(:,:,lg) = ssdn3_temp(:,:)
        isflag(:,lg) = isflag_temp(:)

      case default ! old land surface scheme removed
        print *,'***  Invalid Input: Program Exiting'
        stop
      end select
      write(54,*) 'After CABLE_run, lg = ', lg

C****
C**    UPDATE MLO POINTS
C****
      zfrm(1)=0.0
      zfrm(2)=0.0

      if(.not.lcouple)then

        call checkl(mlo,ln2,flag1)
        if(flag1)then
          call updmlo1(lg,sg,blg,degdt,dfgdt,z4,tintp, ! Input
     &                 tg,rg,eg,fg,             ! Input/Output
     &                 hfrm)                    ! Output
        endif

      endif

C****
C**   UPDATE THERMODYNAMIC ICE POINTS.
C**   (If semice and leads, then updates done via seaice.f and surfupl.f)
C****
      if(.not.semice)then

        call checkl(cice,ln2,flag1)
        if(flag1)then
          call updice1(lg,sg,blg,degdt,dfgdt, ! Input
     &                 cie,gamms,cls,         ! Input
     &                 tg,rg,eg,fg)           ! Input/Output
        endif

      endif
C****
C**   END OF FIRST PART OF SURFACE UPDATING
C****

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,i1)')'After surfupa. IPASS = 1'
          write(25,1)'eg ',eg(mg),' fg ',fg(mg),' tg ',tg(mg)
          write(25,1)'tb2 ',tb2(mg),' tb3 ',tb3(mg)
          write(25,*)
        endif
      endif
 1    format(3(a,f8.2))

      return
      end
C---------------------------------------------------------------------
      subroutine updmlo1(lg,sg,blg,degdt,dfgdt,z4,tintp,
     &                   tg,rg,eg,fg,
     &                   hfrm)

      implicit none

!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /RELATE/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg
      real sg(ln2)
      real blg(ln2)
      real degdt(ln2)
      real dfgdt(ln2)
      real z4(ln2)
      real tintp(ln2)
      real tg(ln2)
      real rg(ln2)
      real eg(ln2)
      real fg(ln2)
      real hfrm(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'LOGIMSL.f'
      include 'RELATE.f'

C Global data blocks
      include 'FEWFLAGS.f'
      include 'QFLXDAT.f'
      include 'TIMEX.f'

      real hmo,hmo2
      common/hm/hmo(ln2,lat),hmo2(lon,lat2,2)

C Local work arrays and variables

      integer ma
      integer mg
      integer ns

      real b1
      real deltat
      real dirad
      real frmix
      real gtop
      real gbot
      real hcap1
      real sheat
      real tratlm
      real wt5day
      real zfrmx

C Local data, functions etc
      real hcap
      data hcap/2.095e8/

C Start code : ----------------------------------------------------------

C****
C**   UPDATE MLO POINTS.
C****
      IF ( qflux ) THEN

c     QFLUX case
       wt5day=1.0/(5.0*(1440.0/mstep))
       newqfluxrun=.false.
Cfix  newqfluxrun=mins.eq.?????? ! change this for new qflux run
       if(newqfluxrun)then !Use weight of zero at start of new qflux run
         wt5day=1.0
       endif

       do mg=1,ln2
        if ( mlo(mg) ) then
C
C       MLO : COMPUTE SURFACE TEMPERATURE
C       sheat=SURFACE HEATING IN WATTS/M**2
C
         gtop=sg(mg)-rg(mg)-eg(mg)-fg(mg)
         dirad=4.0*blg(mg)/tg(mg)
         b1=dirad+degdt(mg)+dfgdt(mg)
         ochfa(mg,lg)=(1.0-wt5day)*ochfa(mg,lg) + wt5day*gtop
         tratlm=1.0-ratlm
         sheat=ochfa(mg,lg)+tratlm**2 * gm1cur(mg,lg)
     &    +tratlm*gm2cur(mg,lg)+gm3cur(mg,lg) !IGW's variable MLO depth
         hcap1=hcap/50.*hmo(mg,lg)
         gbot=(hcap1/dt)+b1*wt5day
C**     REMOVE Z4 FACTOR FROM MLO T* (SURFACE FLUX COMPUTATIONS
C**     ARE NOW COMPLETED)
         tg(mg)=(tg(mg)+z4(mg)*0.0065)  + sheat/gbot
	end if
       end do

      ELSE

c     Interpolated SSTs as a background value
       frmix=0.5e-06

       do 383 mg=1,ln2
        if ( mlo(mg) ) then
C
C     MLO : COMPUTE SURFACE TEMPERATURE,CHECK FOR STATUS CHANGE.
C       g=SHEAT=SURFACE HEATING IN WATTS/M**2
C       NOW INCLUDES MIXING/DIFFUSION TO LIMIT TEMP CHANGE:
C       1/FRMIX GIVES E-FOLDING TIME(SECS) FOR MLO TEMP DIFF.
C
c     remove elevation correction set in surfset
CZZZZ Change to 50 m for MLO Feb 94, Ice model Tmix will have 150m
c        hcap1=hcap/50.*hmo(mg,lg)
CZZZZ Change to 100 m for MLO
         hcap1=hcap*2.

         tg(mg)=tg(mg)+z4(mg)*0.0065
         hfrm(mg)=frmix*hcap1*(tintp(mg)-tg(mg))
         gtop=sg(mg)-rg(mg)-eg(mg)-fg(mg)+hfrm(mg)
         dirad=4.0*blg(mg)/tg(mg)
         b1=dirad+degdt(mg)+dfgdt(mg)
         gbot=(hcap1/dt)+b1+frmix*hcap1
         deltat=gtop/gbot
         tg(mg)=tg(mg)+deltat
         fg(mg)=fg(mg)+deltat*dfgdt(mg)
         eg(mg)=eg(mg)+deltat*degdt(mg)
         rg(mg)=rg(mg)+deltat*dirad
         hfrm(mg)=frmix*hcap1*(tintp(mg)-tg(mg))
        endif
  383  continue

       do 3832 ns=1,2
        ma=(ns-1)*lon
        zfrmx=0.0
        do 3831 mg=1+ma,lon+ma
 3831    if(mlo(mg))zfrmx=zfrmx+hfrm(mg)
 3832  zfrm(ns)=zfrmx

      END IF

      return
      end
C---------------------------------------------------------------------
      subroutine updice1(lg,sg,blg,degdt,dfgdt,
     &                   cie,gamms,cls,
     &                   tg,rg,eg,fg)

      implicit none

!$OMP THREADPRIVATE ( /LOGIMSL/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg
      real sg(ln2)
      real blg(ln2)
      real degdt(ln2)
      real dfgdt(ln2)
      real cie(ln2)
      real gamms(ln2)
      real cls(ln2)
      real tg(ln2)
      real rg(ln2)
      real eg(ln2)
      real fg(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'LOGIMSL.f'

C Global data blocks
      include 'FREEZE.f'
      include 'TIMEX.f'

C Local work arrays and variables

      real b1
      real deltat
      real dirad
      real gtop
      real gbot
      real hicf

      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

C****
C**   UPDATE SEA-ICE POINTS.
C****
      do 385 mg=1,ln2
      if(cice(mg))then
        hicf=cie(mg)*(tg(mg)-tfi)
        gtop=sg(mg)-rg(mg)-cls(mg)*eg(mg) -fg(mg)-hicf
        dirad=4.0*blg(mg)/tg(mg)
        b1=dirad+cls(mg)*degdt(mg)+dfgdt(mg)+cie(mg)
        gbot=(gamms(mg)/dt)+b1
        deltat=gtop/gbot
        tg(mg)=tg(mg)+deltat
        fg(mg)=fg(mg)+deltat*dfgdt(mg)
        eg(mg)=eg(mg)+deltat*degdt(mg)
        rg(mg)=rg(mg)+deltat*dirad
      endif
  385 continue

      return
      end
