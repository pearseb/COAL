c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Removed DEFAULT(NONE) clauses from the OpenMP directives, in order to fix
c problems on both the SC and blizzard.
c SJP 2004/03/16
c
c Modified to make use of DYNAMIC scheduling, which leads to a considerable
c improvement in performance
c SJP 2003/03/29
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13

      subroutine physgm

c Calculate global means from physics zonal averages 

      implicit none
C Global parameters
      include 'PHYSPARAMS.f'
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'
      include 'GAUSL.f'
      include 'GLMEAN.f'   !contains common block glmean used for global
      include 'TIMEX.f'
      include 'WORKA.f'
      include 'ZMEAN.f'

C Local work arrays and variables
      real edbar(nl)
      logical radstep

      integer k
      integer lg

      real rofbar
      real x3

C Local data, functions etc

C Start code : ----------------------------------------------------------

      if(nrad.eq.0)then
        radstep=.false.
      else
        radstep = mod(mins/int(mstep, 8),int(nrad, 8)).eq.0_8
      endif

      do 525 k=1,nl
      edbar(k)=0.0
      emean(k)=0.0
  525 tbar(k)=0.0
      tempbar=0.0
      rofbar=0.0
      gsbal=0.0
      sume=0.
      tsbar=0.0
      rgbar=0.0
      sgbar=0.0
      evbar=0.0
      rnbar=0.0
      hfbar=0.0
      sndbar=0.0
      sidbar=0.0
c     seabar=0.0
c     Only calculate radiative diagnostics on a radiation step
      if (radstep) then
        rtbar=0.0
        sinbar=0.0
        soubar=0.0
        clforbar=0.0
        albbar=0.0
        clhbar=0.0
        clmbar=0.0
        cllbar=0.0
        clbar=0.0
      end if

      global_tmin=999.
      global_tmax=-999.

C**** FORM GLOBAL MEAN TEMPS ETC FROM THE LATITUDINAL PARTS

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (k ,lg)
!$OMP& SHARED  (edbar, edifzk, emean, ezonk, tbar, tzonk, w)

      do 82 k=1,nl
          do 83 lg=1,lat
            edbar(k)=edbar(k)+(edifzk(lg,1,k)+edifzk(lg,2,k))*w(lg)
            emean(k)=emean(k)+(ezonk(lg,1,k)+ezonk(lg,2,k))*w(lg)
            tbar(k)=tbar(k)+(tzonk(lg,1,k)+tzonk(lg,2,k))*w(lg)
 83       continue
        edbar(k)=edbar(k)*0.5/lon
        emean(k)=(emean(k)*0.5/lon)/eradsq
        tbar(k)=tbar(k)*0.5/lon
 82   continue

!$OMP END PARALLEL DO

*PDIR SERIAL
      do lg = 1, lat
        x3=0.5*w(lg)/lon
        global_tmin = min(global_tmin,z_tmink(lg,1),z_tmink(lg,2)) !tmin
        global_tmax = max(global_tmax,z_tmaxk(lg,1),z_tmaxk(lg,2)) !tmax
        rofbar = rofbar + (zrunofk(lg,1)+zrunofk(lg,2))*x3 !mean runoff
        gsbal = gsbal + (zsbalk(lg,1)+zsbalk(lg,2))*x3     !mean (Sg-Rg-Ev-Hf)
        tsbar = tsbar + (ztsk(lg,1)+ztsk(lg,2))*x3         !mean surf temp
        rgbar = rgbar + (zrgk(lg,1)+zrgk(lg,2))*x3         !mean lw at ground
        sgbar = sgbar + (zsgk(lg,1)+zsgk(lg,2))*x3         !mean sw at ground
        evbar = evbar + (zevk(lg,1)+zevk(lg,2))*x3         !mean evap
        rnbar = rnbar + (zrnk(lg,1)+zrnk(lg,2))*x3         !mean rain
        hfbar = hfbar + (zhfk(lg,1)+zhfk(lg,2))*x3         !mean hflux
        sndbar = sndbar + (zsndk(lg,1)+zsndk(lg,2))*x3     !mean snow depth
        sidbar = sidbar + (zsidk(lg,1)+zsidk(lg,2))*x3     !mean sea ice depth
c       seabar = seabar + (zfrack(lg,1)+zfrack(lg,2))*x3   !mean sea fraction

c Accumulate the following on a radiation step
        if (radstep) then
          rtbar = rtbar + (zrtk(lg,1)+zrtk(lg,2))*x3      !mean lw at toa
          sinbar = sinbar + (zsink(lg,1)+zsink(lg,2))*x3  !mean solar in at toa
          soubar = soubar + (zsouk(lg,1)+zsouk(lg,2))*x3 !mean solar out at toa
                                              !mean net cloud radiative forcing
          clforbar = clforbar + (zclfork(lg,1)+zclfork(lg,2))*x3
                                            !mean albedo (yes, *lon is needed!)
          albbar = albbar + (zalbk(lg,1)+zalbk(lg,2))*x3*lon
          clhbar = clhbar + (zclhk(lg,1)+zclhk(lg,2))*x3  !mean high cloud
          clmbar = clmbar + (zclmk(lg,1)+zclmk(lg,2))*x3  !mean mid cloud
          cllbar = cllbar + (zcllk(lg,1)+zcllk(lg,2))*x3  !mean low cloud
          clbar = clbar + (zclk(lg,1)+zclk(lg,2))*x3      !mean total cloud
	endif

      enddo

      do 230 k=1,nl
        sume=sume+edbar(k)*dsk(k)
        tempbar=tempbar+tbar(k)*dsk(k) !mean temp (all levels)
 230  continue 

      grunof=grunof+rofbar

*PDIR ENDSERIAL

!CMIC$ GUARD
!*PDIR CRITICAL
!      print *,'sume=',sume
!      print *,'gsbal=',gsbal
!      print *,'tsbar=',tsbar
!      print *,'rgbar=',rgbar
!      print *,'sgbar=',sgbar
!      print *,'evbar=',evbar
!      print *,'rnbar=',rnbar
!      print *,'hfbar=',hfbar
!      print *,'sndbar=',sndbar
!      print *,'sidbar=',sidbar
!      print *,'seabar=',seabar
!      print *,'tempbar=',tempbar
!      if(lg*ns.ne.77)stop
!CMIC$ ENDGUARD
!*PDIR ENDCRITICAL

      return
      end
