c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c the old land surface scheme removed updlnd2
c nsib comparison is modified to lsm_type = "nsib "
c the if-else statement for lsm_type is modified to case statement
c AJA 2009/01/22
c
c (1) Modified to enable ATHFA, the heat flux into the ocean with zero sub-ice
c     heat input, to be saved.
c (2) HFACA removed from /FICECON/, as this array is never used.
c SJP 2003/06/20
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: surfupb.f,v $
c Revision 1.65  2001/10/12 02:06:58  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.64  2001/02/22 05:34:44  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.63  2001/02/12 05:39:55  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.62  2000/06/20 02:08:34  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.61  1999/05/20 06:23:56  rot032
c HBG changes to V5-2
c
c Revision 1.60  1998/12/10  00:55:53  ldr
c HBG changes to V5-1-21
c
c Revision 1.59  1997/12/23  00:23:38  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.58  1997/12/17  23:22:58  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.57  1997/12/03  05:35:42  ldr
c Put IGW's qflux code with variable MLO depth into main version of model.
c
c Revision 1.56.2.1  1997/12/19  02:03:18  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.56  1996/06/13  02:08:34  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.55  1996/03/21  03:19:09  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.54  1995/11/16  03:54:16  ldr
c New 9 soil types from EAK.
c
c Revision 1.53  1995/08/31  04:30:49  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.52  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.45.2.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.51.1.1  1995/11/10  05:23:03  ldr
c EAK's changes to soil types (increase to 9 soil types).
c
c Revision 1.51  1995/05/02  07:49:15  ldr
c Merge of HBG stuff in V4-5-43h3 into latest version.
c
c Revision 1.50  1995/03/14  01:21:44  ldr
c Changes to qcloud: incorporate hlfusion in snow that hits ground; use
c preset vice in slow icefall calc. - this is run E85, still with small
c fluxi array in icefall.f.
c
c Revision 1.49  1994/09/12  16:33:15  ldr
c Comment out an unnecessary diagnostic.
c
c Revision 1.48  94/08/08  17:22:47  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.47  94/08/08  13:16:30  ldr
c Push the debug prints down into individual routines.
c 
c Revision 1.46  94/08/04  16:56:42  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.46.1.1  1995/05/02  07:31:32  ldr
c Extra changes from HBG to files which already went in as V4-5-43h.
c
c Revision 1.45  94/03/30  12:35:21  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.44  93/12/23  15:31:27  ldr
c Changes to V4-4-54l from HBG for coupled model of V4-5
c 
c Revision 1.43  93/12/20  16:22:07  ldr
c Minor changes to V4-4-53l from HBG for coupled model
c 
c Revision 1.42  93/12/17  15:33:59  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.41  93/12/17  12:06:00  ldr
c Merge of HBG and LDR changes.
c 
c Revision 1.40  93/12/06  17:32:35  ldr
c Get rid of some irritating warnings on the VP.
c 
c Revision 1.39.1.1  93/12/17  11:51:50  ldr
c Changes to V4-4-45l from HBG for coupled model
c 
      subroutine surfupb(lg,tintp,ttg,degdt,dfgdt,eg,fg,ci, !Inputs
     &                   gamms,sg,rg,blg,condx,rcondx,preci,rpreci,he, !Inputs
     &                   snowd,siced,tg,cie,bs,wgp,                !I & O
     &                   tstar,sublmi,snowmi,il,perc,ifroz,wg,wg2, !Outputs
     &                   runoff)                                   !Outputs

      implicit none

!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /RELATE/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg
      real tintp(ln2)
      real ttg(ln2,nl)
      real degdt(ln2)
      real dfgdt(ln2)
      real eg(ln2)
      real fg(ln2)
      real ci(ln2)
      real gamms(ln2)
      real sg(ln2)
      real rg(ln2)
      real blg(ln2)
      real condx(ln2)
      real rcondx(ln2)
      real preci(ln2),rpreci(ln2) !Tau and tau-1 ice precip if (qcloud)
      real he(ln2)
      real snowd(ln2)
      real siced(ln2)
      real tg(ln2)
      real cie(ln2)
      real bs(ln2)
      real wgp(ln2)
      real tstar(ln2)
      real sublmi(ln2)
      real snowmi(ln2)
      integer il(ln2)
      real perc(ln2)
      integer ifroz(2)
      real wg(ln2)
      real wg2(ln2)
      real runoff(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'LOGIMSL.f'
      include 'RELATE.f'

C Global data blocks
      include 'CHMAP.f'
      include 'FEWFLAGS.f'
      include 'FREEZE.f'
      include 'LSMI.f'
      include 'QFLXDAT.f'
      include 'SURF1.f'
      include 'TIMEX.f'

      real tsst1,tsst2,tdtm1,tdtm2
      common/compsst/tsst1(ln2,lat),tsst2(ln2,lat)
     &              ,tdtm1(ln2,lat),tdtm2(ln2,lat)

C Local work arrays and variables
      logical tester1(ln2),flag1
      real tintpx(ln2)

      integer idiff
      integer mg
      integer ns

      real tadjum

C Local data, functions etc
      character*1 c10d(21)
      data c10d/'w','j','i','h','g','e','d','c','b','a','-','1',
     & '2','3','4','5','6','7','8','9','t'/
      real tfrz
      data tfrz/273.1/

C Start code : ----------------------------------------------------------

C****
C**** ROUTINE TO DO THE LAST PART OF THE UPDATE OF THE SURFACE VALUES
C**** (Note: basic updating of tg for land, cice, mlo done in surfupa.f)
C****

c  Set il=0 for all points - later set to 1 for cice and newly formed cice

      do mg=1,ln2
        il(mg)=0
      enddo

C****
C**    SEA POINTS
C****

      do 382 mg=1,ln2
c     remove elevation correction
        if(sea(mg))tg(mg)=tintp(mg)
  382 continue

C****
C**    CHECK MLO POINTS FOR FREEZING (CHANGE TO SEA-ICE)
C****

      do mg=1,ln2
        tester1(mg)=
     &    (mlo(mg).or.(lcouple.and.sea(mg))).and.(tg(mg).lt.tfi)
      enddo
      call checkl(tester1,ln2,flag1)
      if(flag1)then
        call updmlo2(lg,tester1,               ! Input
     &               tg,                       ! Input/Output
     &               siced,gamms,cie,il,ifroz) ! Output
      endif

C****
C**   UPDATE LAND POINTS.
C****
      select case (lsm_type)
      case("nsib ")
        call surfb(lg,ttg,snowd,tg,rcondx,eg,perc,rpreci,he,
     &             wg, wg2,runoff)

      end select

C****
C**   UPDATE ICE POINTS.
C****

      zflxi(1)=0.0
      zflxi(2)=0.0

C**** CHECK IF ANY ICE AT THIS LAT
      do 3850 mg=1,ln2
 3850   tester1(mg)=imsl(mg,lg).eq.1

      call checkl(tester1,ln2,flag1)

      IF(flag1)THEN !Do ice check

       if(semice)then
        call seaice(lg,sg,rg,blg,degdt,dfgdt,eg,fg,
     &              condx,gamms,preci,
     &              tg,siced,snowd,sublmi,snowmi,
     &              il)

       else  !  use old sea-ice scheme

        call updice2(lg,blg,degdt,dfgdt,gamms, ! Input
     &               tintp,eg,ci,condx,preci,  ! Input
     &               tg,cie,bs,snowd,siced,    ! Input/Output
     &               sublmi,snowmi)            ! Output
       endif

      END IF ! END OF ICE CHECK (flag1=true)

C****
C**   UPDATE SURFACE TEMPERATURE.
C****
      do 286 mg=1,ln2
        tstar(mg)=tg(mg)
  286   snowd(mg)=max(0.0,snowd(mg))

C****
C**   END OF SURFACE UPDATING
C****

C****
C**   Optional diagnostics and mapping next
C****

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,3i3)')'IPASS=1, after surfupb. INS,LG,MG= ',
     &         ns,lg,mgdebug
          write(25,1)'wb1',wb(mg,1,lg),' wb ',wb(mg,ms,lg),
     &               ' tg ',tg(mg)
          write(25,1)' snowd ',snowd(mg),' siced ',siced(mg)
          write(25,*)
        endif
      endif
 1    format(3(a,f7.2))

C****
C**   UPDATE TMLO DIFFERENCE MAP AT END OF DAY ONLY
C****
      if(mod(mins+int(mstep, 8),1440_8).ne.0_8)return

      do 3831 mg=1,ln2
 3831 itdiff(mg,lg)=9999

      IF (lcouple) THEN

c.... all points are sea points (no mlo points)
c.... compare temperature supplied to atmospheric model by the 
c.... ocean model (which has Delta-T corrections and dtm correction
c.... terms added) to the apparent SSTs used to drive the atmos 
c.... model i.e. the interpolated climatological SSTs plus dtm
      do 3833 mg=1,ln2
      if(sea(mg).or.cice(mg))then
       if(fluxadj)then ! Delta-T correction
        tadjum=tdtm1(mg,lg)*(1.0-ratlm)+tdtm2(mg,lg)*ratlm
       else
        tadjum=0.0
       endif
       tintpx(mg)=tsst1(mg,lg)*(ratlm-1.0)-tsst2(mg,lg)*ratlm
     &   +tadjum
       itdiff(mg,lg)=nint(tintp(mg)-tintpx(mg))
      end if
 3833 continue

      if(mlomap)then
      do 3835 mg=1,ln2
      if(sea(mg).or.cice(mg))then
       idiff=max(-10,min(10,int((tintp(mg)-tintpx(mg))*10.0)))
       i10d(mg,lg)=c10d(idiff+11)
      end if
 3835 continue
      endif

      ELSE

c.... if not coupled model, check mlo points only
      do 3834 mg=1,ln2
      if(mlo(mg))then
       itdiff(mg,lg)=nint(tg(mg)-tintp(mg))
      elseif(cice(mg))then
       itdiff(mg,lg)=150
      end if
 3834 continue

      if(mlomap)then
      do 3836 mg=1,ln2
      if(mlo(mg))then
       idiff=max(-10,min(10,int((tg(mg)-tintp(mg))*10.0)))
       i10d(mg,lg)=c10d(idiff+11)
      end if
 3836 continue
      endif

      END IF

      if(mlomap)then
      do 3832 mg=1,ln2
      if(cice(mg))then
       i10d(mg,lg)='m'
       if(tg(mg).lt.tfrz)i10d(mg,lg)='f'
      end if
 3832 continue
      endif

      return
      end
C---------------------------------------------------------------------
      subroutine updmlo2(lg,tester1,
     &                   tg,
     &                   siced,gamms,cie,il,ifroz)

      implicit none

!$OMP THREADPRIVATE ( /FLOES/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg
      logical tester1(ln2)
      real siced(ln2)
      real tg(ln2)
      real gamms(ln2)
      real cie(ln2)
      integer il(ln2)
      integer ifroz(2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      real dsn,dic,tn0,tn1,tn2,tni,sto,pl,tmix
     & ,plsum,vol,tflux,sflux,fflux
      common/floes/dsn(ln2),dic(ln2),tn0(ln2),
     &tn1(ln2),tn2(ln2),tni(ln2),sto(ln2),pl(ln2),tmix(ln2),
     &plsum(ln2),vol(ln2),tflux(ln2),sflux(ln2),fflux(ln2)
      include 'LOGIMSL.f'

C Global data blocks
      include 'FEWFLAGS.f'
      include 'FREEZE.f'
      include 'LSMI.f'
      include 'TIMEX.f'
c.... From ocean model - need depth of 1st layer
      include 'TTFC.f'

      real drunof,dsis,dsisb,dsfw,athf,athfa,atsi
      common/aforce/drunof(ln2,lat),dsis(ln2,lat),dsisb(ln2,lat)
     & ,dsfw(ln2,lat),athf(ln2,lat),athfa(ln2,lat),atsi(ln2,lat)

      real hmo,hmo2
      common/hm/hmo(ln2,lat),hmo2(lon,lat2,2)

C Local work arrays and variables
      real hcap1(ln2)

      real dicx
      real dz1
      real plx
      real tfix

      integer mg

C Local data, functions etc
      real qice,hcap
      data qice,hcap/3.35e8,2.095e8/

C Start code : ----------------------------------------------------------

C**** Mixed Layer Ocean (MLO) points (or sea points when coupled)
C****  are checked if they are cold enough to freeze.
C****  If so, convert status to that of an ice point.

ch... Set the mixed layer depth at 100m or depth of 1st ocean level
ch... if coupled model.
      if(qflux)then
        do mg=1,ln2
          dz1=hmo(mg,lg)
          hcap1(mg)=hcap/50.*dz1
        enddo
      else
        if(lcouple)then
          dz1=dzlev1
        else
          dz1=100.0
        endif
        do mg=1,ln2
          hcap1(mg)=hcap/50.*dz1
        enddo
      endif

      do 383 mg=1,ln2
c sea is added in for coupled points and it is written this way to
c prevent land points being included
      if(tester1(mg))then
c--     tester1(mg)=
c--  &    (mlo(mg).or.(lcouple.and.sea(mg))).and.(tg(mg).lt.tfi)

C
C        START SEA ICE FORMATION,100% OF GRID AREA BECOMES ICE COVERED
C        SEA ICE THICKNESS IN M
C

       siced(mg)=(tfi-tg(mg))*hcap1(mg)/qice
       tfix=tfi
       if(siced(mg).gt.0.15)then
         tfix=tg(mg)+(0.15*qice/hcap1(mg))
         siced(mg)=0.15
       endif
       if(semice)then
c..... start with a minimum cover of 4% ice if possible
         plx=0.96
         dicx=siced(mg)/(1.0-plx)
         if(dicx.gt.0.1)then
           dicx=0.1
           plx=1.0-siced(mg)/dicx
         endif
         dic(mg)=dicx
         pl(mg)=plx
         il(mg)=1
         dsn(mg)=0.0
c set inital internal snow temperature to melt value for snow used in seaice.f
         tn0(mg)=273.00
         tni(mg)=tfix
         tn1(mg)=tfix
         tn2(mg)=tfix
         sto(mg)=0.0
       endif
         athf(mg,lg)=athf(mg,lg)+(tfix-tg(mg))*hcap1(mg)/dt
         athfa(mg,lg)=athfa(mg,lg)+(tfix-tg(mg))*hcap1(mg)/dt
         tg(mg)=tfix
         tmix(mg)=tfix
         mlo(mg)=.false.
         sea(mg)=.false.
         cice(mg)=.true.
         imsl(mg,lg)=1
         just_frozen(mg,lg)=.true.
         ifroz((mg-1+lon)/lon)=1
         gamms(mg)=3.471e+05
         cie(mg)=2.04/(siced(mg)+0.001)
       endif
  383 continue

      return
      end
C---------------------------------------------------------------------
      subroutine updice2(lg,blg,degdt,dfgdt,gamms,
     &                   tintp,eg,ci,condx,preci,
     &                   tg,cie,bs,snowd,siced,
     &                   sublmi,snowmi)

      implicit none

!$OMP THREADPRIVATE ( /FLOES/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /RELATE/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg
      real blg(ln2)
      real degdt(ln2)
      real dfgdt(ln2)
      real gamms(ln2)
      real tintp(ln2)
      real eg(ln2)
      real ci(ln2)
      real condx(ln2)
      real preci(ln2) ! ice precip if (qcloud)
      real tg(ln2)
      real cie(ln2)
      real bs(ln2)
      real snowd(ln2)
      real siced(ln2)
      real sublmi(ln2)
      real snowmi(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      real dsn,dic,tn0,tn1,tn2,tni,sto,pl,tmix
     & ,plsum,vol,tflux,sflux,fflux
      common/floes/dsn(ln2),dic(ln2),tn0(ln2),
     &tn1(ln2),tn2(ln2),tni(ln2),sto(ln2),pl(ln2),tmix(ln2),
     &plsum(ln2),vol(ln2),tflux(ln2),sflux(ln2),fflux(ln2)
      include 'LOGIMSL.f'
      include 'RELATE.f'

C Global data blocks
      include 'DATICE1.f'
      include 'FEWFLAGS.f'
      include 'FREEZE.f'
      include 'LSMI.f'
      include 'QFLXDAT.f'
      include 'SURF1.f'
      include 'TIMEX.f'

      real flxia
      common/ficecon/flxia(lat,2)

C Local work arrays and variables
      real hbotx(ln2)
      real flx(ln2)

      integer ma
      integer mg
      integer ns

      real b1
      real dirad
      real dsic
      real dsiced
      real rainf
      real subl
      real tboc
      real tice
      real tgtf
      real smelt
      real zflxix

C Local data, functions etc
      real qice,tfrz
      data qice,tfrz/3.35e8,273.1/

C Start code : ----------------------------------------------------------

C**** Code for the simple thermodynamic sea-ice model (Not semice or leads)

      do 3851 mg=1,ln2
      if(cice(mg))then
        dirad=4.0*blg(mg)/tg(mg)
        b1=dirad+1.13344*degdt(mg)+dfgdt(mg)+cie(mg)
        hbotx(mg)=(gamms(mg)/dt)+b1
        snowd(mg)=snowd(mg)+preci(mg)
        if(tg(mg).le.273.1)then
          rainf=condx(mg)-preci(mg)
          snowd(mg)=snowd(mg) + rainf
          tg(mg)=tg(mg)+rainf*hlf/(dt*hbotx(mg))
          if(snowd(mg).gt.0.0)then
            bs(mg)=0.31e+03/(snowd(mg)+1.0)
            cie(mg)=bs(mg)*ci(mg)/(bs(mg)+ci(mg))
          endif
        endif
      endif
 3851 continue

c**** Sub-ice flux of heat from ocean 
c**** The flxia values are preset in icecon.f
c**** and hoice is usually zero if not semice model.
      do mg=1,lon
        flx(mg)=flxia(lg,1)+hoice(mg,lg)
        flx(mg+lon)=flxia(lg,2)+hoice(mg+lon,lg)
      enddo

      do 385 mg=1,ln2

      if(cice(mg))then

      tice=tg(mg)
chal  old ice model method had temp under ice always at freezing (tfi)
chal  Changed 10/11/93 - untested. 
chal  tboc=tfi
      tboc=tintp(mg)
      if(snowd(mg).gt.0.0)
     &  tice=(bs(mg)*tg(mg)+ci(mg)*tboc)/(bs(mg)+ci(mg))
      subl=dt*eg(mg)/hl
C     NOTE : SUBL IN MMS OF WATER, SNOWD IN CMS OF SNOW=MMS OF WATER
C          : SICED IN M OF ICE (FACTOR 1000 FOR MMS OF WATER)
C         SUBL CHANGES SNOWD IF SNOW COVER, ELSE CHANGES SICED
      if(tg(mg).gt.tfrz)go to 471
C
C     FREEZING CONDITIONS
C
      if(snowd(mg).le.0.0)go to 47
      snowd(mg)=snowd(mg)-subl
      subl=0.0
      go to 47
  471 continue
C
C     MELTING CONDITIONS
C
      tgtf=tg(mg)-tfrz
      tg(mg)=tfrz
      if(snowd(mg).gt.0.0) go to 46
C
C     SNOW FREE ICE,COMPUTE SURFACE MELTING OF ICE.
C
      siced(mg)=siced(mg)-(dt*tgtf*hbotx(mg)/qice)
      tice=tfrz
      if(siced(mg).gt.0.0) go to 47
      go to 48
C
C     SNOW COVERED ICE,COMPUTE SNOW MELT
c     (computed as kg water/snow, = mm of water, approx = cm of snow)
   46 smelt=dt*tgtf*hbotx(mg)/hlf
      tice=(bs(mg)*tfrz+ci(mg)*tboc)/(bs(mg)+ci(mg))
      snowd(mg)=snowd(mg)-smelt
      snowmi(mg)=smelt*0.001
      if(snowd(mg).lt.0.0)go to 47
      snowd(mg)=snowd(mg)-subl
      subl=0.0
C
C     NOW COMPUTE CHANGE IN ICE THICKNESS UNDER ICE
C
   47 continue 
      dsiced=min(0.01,dt*(ci(mg)*(tboc-tice) - flx(mg))/qice)
      sublmi(mg)=subl*0.001
      siced(mg)=siced(mg)+dsiced-sublmi(mg)
C---- SNOW TO ICE CONVERSION AND LIMITATION
      if(snowd(mg).gt.200.0)then
        siced(mg)=siced(mg)+(snowd(mg)-200.0)*0.001
        snowd(mg)=200.0
      end if
c       siced(mg)=min( siced(mg),6.0 )
      if(siced(mg).gt.6.0)then
        dsic=(siced(mg)-6.0)*1000.0*(1-pl(mg))
c----  change local tg to account for energy - clearly not best method
        tg(mg)=tg(mg)-dsic*hlf/(dt*hbotx(mg))
        siced(mg)=6.0
      end if
C----
      if(siced(mg).gt.0.0)go to 385
C
C     SEA ICE HAS MELTED,POINT CHANGES STATUS TO MIXED LAYER OCEAN
C
   48 siced(mg)=0.0
      snowd(mg)=0.0
      tg(mg)=tintp(mg)
      cice(mg)=.false.
      mlo(mg)=.true.
C**   CHANGE STATUS OF POINT AT NEXT LATITUDE IF MLO
      imsl(mg,lg)=2
      just_melted(mg,lg)=.true.
c      if( .not.qflux .and. imsl(mg,lg+1).eq.2)imsl(mg,lg+1)=3

      endif ! cice

  385 continue

      do 387 ns=1,2
       ma=(ns-1)*lon
       zflxix=0.0
       do 386 mg=1+ma,lon+ma
  386   if(cice(mg))zflxix=zflxix+flx(mg)
  387  zflxi(ns)=zflxix

      return
      end
