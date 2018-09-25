c Add a loop after the call to SMOISTURE to detect and replace nonsensical
c values for the soil moisture. See the comments within the source code for a
c more detailed explanation.
c SJP 2009/04/30
c
c (1) Fixing syntax errors identified by the g95 Fortran compiler.
c (2) Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: surfb.f,v $
c Revision 1.40  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.39  2001/02/28 04:36:38  rot032
c Further tidy ups from HBG
c
c Revision 1.38  2001/02/22 05:34:40  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.37  2001/02/12 05:39:49  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.36  2000/11/14 03:11:36  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.35  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.34  1998/12/10  00:55:37  ldr
c HBG changes to V5-1-21
c
c Revision 1.33  1998/05/26  05:29:35  ldr
c Final changes (inc. ACH salinity stuff) for 10k run.
c
c Revision 1.32  1997/12/23  06:43:54  ldr
c Changes to allow array bounds checking.
c
c Revision 1.31  1997/12/23  01:32:32  ldr
c Add a couple of TASKLOCAL directives for NEC.
c
c Revision 1.30  1997/12/23  00:23:35  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.29  1997/12/17  23:22:48  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.28.2.2  1997/12/23  00:01:07  ldr
c Hal's restart fixes.
c
c Revision 1.28.2.1  1997/12/19  02:03:13  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.28  1997/01/29  23:18:12  ldr
c Reinstate limiting of snowd to 400cm.
c
c Revision 1.27  1996/12/23  23:34:43  ldr
c Comment out line restricting snow depth to 400 (requested by EAK).
c
c Revision 1.26  1996/12/19  23:40:50  ldr
c EAK's snow albedo changes merged and debugged.
c
c Revision 1.25  1996/10/24  01:03:17  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.24.1.1  1996/12/05  03:42:33  ldr
c Fix to snow albedo from EAK.
c
c Revision 1.24  1996/06/13  02:08:25  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.23  1996/03/21  03:19:07  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.22  1995/11/16  03:54:16  ldr
c New 9 soil types from EAK.
c
c Revision 1.21  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.18.2.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.20  1995/03/14  01:21:44  ldr
c Changes to qcloud: incorporate hlfusion in snow that hits ground; use
c preset vice in slow icefall calc. - this is run E85, still with small
c fluxi array in icefall.f.
c
c Revision 1.19  1995/02/08  23:18:11  mrd
c Improvements to soil moisture scheme (from mrd and eak).
c
c Revision 1.18  1994/08/08  17:22:43  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.20.1.1  1995/11/10  05:23:03  ldr
c EAK's changes to soil types (increase to 9 soil types).
c
c Revision 1.17  94/03/30  12:35:14  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.16  93/12/17  15:33:55  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.15  93/09/13  14:33:59  mrd
c Addded deep soil moisture percolation diagnostic
c 
c Revision 1.14  93/08/31  17:22:27  ldr
c Put /surf1 in include file to avoid problems with declaration of mc.
c 
c Revision 1.13  93/08/19  15:10:55  ldr
c Minor cosmetic changes.
c 
c Revision 1.12  93/07/06  16:38:47  ldr
c Only calculate maps of wg, wb at end of day (HBG).
c 
c Revision 1.11  93/06/28  17:14:56  ldr
c Added code to set up HBG's instantaneous wg, wb maps scaled correctly
c for NSIB scheme.
c 
c Revision 1.10  93/06/18  13:07:17  ldr
c Move /relate for include file.
c
c     INPUT/OUTPUT
c     Input:   from common/fewflags in FEWFLAGS.f
c                  qcloud - T if using prognostic cloud scheme
c
c              from common/logimsl in LOGIMSL.f
c                  land - logical variable for surface type
c
c              from common/surf1 in SURF1.f
c                  tgf - canopy temperature
c
c              from common/timex in TIMEX.f
c                  dt - time step in seconds
c                  mins - current model time in mins
c                  mstep - timestep in minutes
c
c              from common/printt in PRINTT.f
c                  gwicm - if true, prints maps of snow depth and ice depth
c
c     In/Out:  from common/relate in RELATE.f
c                  zgw  - zonal dip moisture
c                  ztsl - zonal surface temperature
c                  zrunof - zonal mean runoff
c 
c This routine is used for NSIB land surface scheme. The routine computes: 
c - the soil moisture in top and bottom layers
c - snow depth, evaporation, melting and the appropriate temperature adjustment
c
c INPUTS:
c ttg   - air temperature
c snowd - snow depth at the previous time step
c tg    - combined soil-cannopy temperature of the surface layer at the
c         current time step 
c rcondx - precipitation 
c eg    - evapotranspiration at the current time step
c wb    - soil moisture (ms layers)
c 
c OUTPUTS:
c tg    - adjusted combined soil-canopy temperature of the surface layer 
c wb    - soil moisture
c wg    - soil moisture available for evapotranspiration at the top layer
c wg2   - soil moisture available for evapotranspiration at the lowest layer
c snowd - snow depth
c runoff- groud runoff
c 
c
      subroutine surfb(lg,ttg,snowd,tg,rcondx,eg,
     &                 perc,rpreci,he,wg,wg2,runoff)

      implicit none

!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /RELATE/ )
!$OMP THREADPRIVATE ( /SOILPR/ )
!$OMP THREADPRIVATE ( /SURFBC/ )
!$OMP THREADPRIVATE ( /VEGDAT/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      real snmin,d1land,rhosnow
      parameter (snmin=0.15)
      parameter (d1land=.03)
      parameter (rhosnow=200.)

C Argument list
      integer lg
      real ttg(ln2,nl)
      real snowd(ln2)
      real tg(ln2)
      real rcondx(ln2)
      real eg(ln2)
      real perc(ln2)
      real rpreci(ln2)!Tau-1 frozen precipitation in mm water/tdt if (qcloud)
      real he(ln2)
      real wg(ln2)
      real wg2(ln2)
      real runoff(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'LOGIMSL.f'
      include 'RELATE.f'

      logical lakesg,lakeind
      common/glakes/lakesg(ln2,lat),lakeind(lat)

      real tgsoil,ssatcur,gamm,coef,scond1s,cap1s,sgamm,sdepth
     & ,hs1tm,hs1tp
      common/soilpr/tgsoil(ln2,ms),ssatcur(ln2,ms+1),gamm(ln2,ms)
     &,coef(ln2,ms+4),scond1s(ln2,3),cap1s(ln2,3)
     &,sgamm(ln2,3),sdepth(ln2,3),hs1tm(ln2,9),hs1tp(ln2,9)

      real egg,condxg,tsigmf,evapxf,ewwp,pmcmax
      common/surfbc/egg(ln2),condxg(ln2),tsigmf(ln2),evapxf(ln2)
     &             ,ewwp(ln2),pmcmax(ln2)

      real rsmin,z0m,sigmf,vegt
      common/vegdat/ rsmin(ln2),z0m(ln2),sigmf(ln2),vegt(ln2)

C Global data blocks
      include 'CHMAP.f'
      include 'FEWFLAGS.f'
      include 'PRINTT.f'
      include 'SURF1.f'
      include 'TIMEX.f'

      real prhos,pcss,pcnsd,pbch,pssat,pswilt,phsbh,psfc,psucs
     & ,phyds,pcdr3,prlaim,prlais,pslveg,pscveg
      integer ip2bp3,ipbp2
      common/psoilvegpr/prhos(ln2,lat),pcss(ln2,lat),
     &pcnsd(ln2,lat),pbch(ln2,lat),pssat(ln2,lat),
     &pswilt(ln2,lat),phsbh(ln2,lat),ip2bp3(ln2,lat),
     &ipbp2(ln2,lat),psfc(ln2,lat),psucs(ln2,lat),
     &phyds(ln2,lat),pcdr3(ln2,lat),prlaim(ln2,lat),
     &prlais(ln2,lat),pslveg(ln2,lat),pscveg(ln2,lat)

      real zs,zsh,zsdzs,zsfrac
      common/soilzs/zs(ms),zsh(ms+1),zsdzs(ms),zsfrac(ms)

      real aftfh,aftfhg,ocondx,opreci,pmc
      integer isoilm
      common/surf/aftfh(ln2,lat),aftfhg(ln2,lat),ocondx(ln2,lat)
     &          ,isoilm(ln2,lat),opreci(ln2,lat)
     &          ,pmc(ln2,lat)

C Local work arrays and variables
      real evapsn(ln2)
      real rnof1(ln2),rnof2(ln2)
      real smelt(ln2)
      real twf(ln2)
      real wblf(ln2,ms)
      logical ism(ln2)
      logical tester1(ln2),tester2(ln2),flag1,flag2
      logical tester3(ln2)

      integer k
      integer kx
c     integer lgns
      integer ma
      integer mg
      integer ns

      real dwb
      real fracs
c     real pmcold
      real rnof5
      real segg
      real sfall
      real sinfil
      real smelt1
      real snowflx
      real tfrzx
      real totwet
      real totwetk
      real wbknew
      real xsinfl
      real zgwx
      real zrunofx
      real ztslx
      real wblfint

C Local data, functions etc
      real tfrz
      data tfrz/273.15/   ! BP aug2010
      save tfrz
      character*1 cint(12)
      data cint/'0','1','2','3','4','5','6','7','8','9','t','t'/

C Start code : ----------------------------------------------------------

      do 302 mg=1,ln2
       ism(mg)=.true.
       tester1(mg)=.false.
       tester2(mg)=.false.
       runoff(mg)=0.0
       perc(mg)=0.
      if(land(mg))then

       evapsn(mg)=0.0
       smelt(mg)=0.0
       osnowd(mg,lg)=snowd(mg)
      if(qcloud.and.rcondx(mg).gt.0.)then  ! Snowfall calculation
        fracs=rpreci(mg)/rcondx(mg) !Fraction of precip falling as snow
        sfall=fracs*((1.-tsigmf(mg))*rcondx(mg)+tsigmf(mg)*condxg(mg))
        snowd(mg) = max (0.0, snowd(mg) + sfall)
        condxg(mg)=condxg(mg)*(1-fracs)
        rcondx(mg)=rcondx(mg)*(1-fracs)
      endif
       if(isflag(mg,lg).eq.0)tester1(mg)=.true.
       if(isflag(mg,lg).ne.0)tester2(mg)=.true.

      endif
  302 continue
c                                          snow accumulation

      call checkl(tester1,ln2,flag1)
      if(flag1)then

      do 306 mg=1,ln2
      if(tester1(mg))then
        if(ttg(mg,2).lt.tfrz.and.tggsl(mg,1,lg).lt.tfrz.and.
     &                                             rcondx(mg).gt.0.)then
          snowflx=(1.-tsigmf(mg))*rcondx(mg)+tsigmf(mg)*condxg(mg)
          snowd(mg)=max(snowd(mg) + snowflx, 0.0)
c         update surface temperatures explicitly for extra snow g fluxes
c         to allow for freezing of rain at surface
          snowflx=snowflx*hlf/dt
          tggsl(mg,1,lg)=tggsl(mg,1,lg)+snowflx*dt/gamm(mg,1)
          rcondx(mg)=0.0
          condxg(mg)=0.0
        endif
c
        if(snowd(mg).le.0.0) go to 44
c
c       snow covered land
        if( snowd(mg).gt.0.1) 
     &  evapsn(mg)=min(snowd(mg),dt*(eg(mg)-tsigmf(mg)*evapxf(mg))/
     &       hl)
        snowd(mg)=snowd(mg)-evapsn(mg)

c---- To prevent small amounts of snow from remaining present
c---- at the temperature (tfrz) for soil-ice being melted (in surfa)
c---- then set the condition for melting snow to be slightly less 
c---- tfrz :
        tfrzx=tfrz-0.1
        if(tggsl(mg,1,lg).ge.tfrzx)then
c**       land,snow,melting
          snowflx=(tggsl(mg,1,lg)-tfrzx)*gamm(mg,1)
c         prevent snow depth going negative
          smelt(mg)= min(snowflx/hlf ,snowd(mg)) 
          snowd(mg)=snowd(mg)-smelt(mg)
          tggsl(mg,1,lg)= tggsl(mg,1,lg)-smelt(mg)*hlf
     &                                            /gamm(mg,1)
        endif
   44 continue
      endif
  306 continue

      endif ! if(flag1)

      call checkl(tester2,ln2,flag2)
      if(flag2)then

      do 308 mg=1,ln2
      if(tester2(mg))then
        if(ttg(mg,2).lt.tfrz.and.rcondx(mg).gt.0.)then
          snowflx=(1.-tsigmf(mg))*rcondx(mg)+tsigmf(mg)*condxg(mg)
          snowd(mg)=max(snowd(mg) + snowflx, 0.0)
c         update surface temperatures explicitly for extra snow g fluxes
c         to allow for freezing of rain at surface
          snowflx=snowflx*hlf/snowd(mg)
          do k=1,3
          tggsn(mg,k,lg)=tggsn(mg,k,lg)+snowflx/sgamm(mg,k)*
     &                       smass(mg,k,lg)
          enddo
          rcondx(mg)=0.0
          condxg(mg)=0.0
        endif
        evapsn(mg)=min(snowd(mg),dt*(eg(mg)-tsigmf(mg)*evapxf(mg))/
     &       hl)
        snowd(mg)=snowd(mg)-evapsn(mg)
        smelt(mg)=0.0
      endif
  308 continue

      do k=1,3
      do 310 mg=1,ln2
      if(tester2(mg))then
c**     land,snow,melting
        if( tggsn(mg,k,lg).gt.tfrz ) then
          snowflx=(tggsn(mg,k,lg)-tfrz)*sgamm(mg,k)
          smelt1= min(snowflx/hlf ,0.9*smass(mg,k,lg)) 
          smass(mg,k,lg) = smass(mg,k,lg) - smelt1
          snowd(mg)=snowd(mg)-smelt1
          tggsn(mg,k,lg)= min(tggsn(mg,k,lg)-smelt1*hlf/
     &                      sgamm(mg,k), tfrz ) 
          smelt(mg)=smelt(mg)+smelt1
        endif
      endif
  310 continue
      enddo

      endif ! if(flag2)

      do mg=1,ln2
        tester3(mg)=land(mg).and.(.not.lakesg(mg,lg))
      enddo

      do k=1,ms
      do 312 mg=1,ln2
      if(tester3(mg))then
        if((tggsl(mg,k,lg).gt.(tfrz-2.0)).or.(snowd(mg).lt.1.0))
     &                                    ism(mg)=.false.
c       ssatcur(mg,k)=pssat(mg,lg)*(1.-wbice(mg,k,lg)/
c    &  pssat(mg,lg)) ! Rewritten : see next line
        ssatcur(mg,k)=pssat(mg,lg)-wbice(mg,k,lg)
        wblf(mg,k)=(wb(mg,k,lg)-wbice(mg,k,lg))/ssatcur(mg,k)
      endif
  312 continue
      enddo

      
      do 314 mg=1,ln2
      if(tester3(mg))then
       totwet=(1.-tsigmf(mg))*rcondx(mg)+smelt(mg)+tsigmf(mg)*condxg(mg)
       pmc(mg,lg)=max(pmc(mg,lg)+totwet-ewwp(mg)*dt/hl, 0.) !-puddle evap.
       totwet=pmc(mg,lg)
       xsinfl=0.00001 + 1.0e-04*min(he(mg),2000.0)
       sinfil=(1.0-xsinfl)*
     &  min((pssat(mg,lg)-wb(mg,1,lg))*zs(1)*1000.,totwet)
       segg=egg(mg)/hl*dt*(1.-tsigmf(mg)) ! Evap from surface (mms/step)
       if(evapsn(mg).ne.0.0) segg=0.0     ! None if evap is from snow 
       twf(mg)=(sinfil-segg)/dt           ! mms/sec going into ground
       pmc(mg,lg)=max(pmc(mg,lg)-sinfil,0.0)
       rnof1(mg)=max(0.,pmc(mg,lg)-pmcmax(mg))
       pmc(mg,lg)=min(pmc(mg,lg),pmcmax(mg))
       if(ism(mg)) rnof1(mg)=rnof1(mg)+(sinfil-segg)
      endif
  314 continue

      if(newriver.and.lakeind(lg))then

       do mg=1,ln2
       if(lakesg(mg,lg))then
        ism(mg)=.true.
        totwetk=rcondx(mg)+smelt(mg)
        totwet=totwetk-evapxf(mg)*dt/hl !-(ground+puddle) evap.
c       pmcold=pmc(mg,lg)
        pmc(mg,lg)=max(pmc(mg,lg)+totwet, 0.)
        twf(mg)=0.0
        rnof1(mg)=max(0.,pmc(mg,lg)-pmcmax(mg))
        pmc(mg,lg)=min(pmc(mg,lg),pmcmax(mg))
c        write(6,6767)mg,lg,pmcold,pmc(mg,lg),totwetk,evapxf(mg)*dt/hl
c6767 format(2i4,7(1x,e13.6))
       endif
       enddo

      endif

      call smoisture(lg ,wblf ,twf ,dt , ism)

c...  If negative values of WBLF have been returned at any gridpoint, replace
c...  the values at that gridpoint with the volume-weighted vertical mean. This
c...  eliminates nonsensical values of WBLF whilst also conserving freshwater.
c...
c...  Diagnostic evaluation of the model has indicated that numerical
c...  instabilities only occur at the rate of O(1) event per model year.
c...  However, the treatment of nonsensical values of WBLF by the subsequent
c...  section of code is non-conservative and can give rise to large pulses of
c...  runoff into the ocean. By eliminating nonsensical values here, we
c...  eliminate these anomalous "pulses" without having any significant impact
c...  upon the climatology of the model.
c...
c...  SJP 2009/04/30

      do mg = 1, ln2
        if (minval(wblf(mg, :)) .lt. 0.0) then
          wblfint = 0.0
          do k = 1, ms
            wblfint = wblfint + wblf(mg, k) * zs(k)
          end do
          wblfint = wblfint / sum(zs)
          do k = 1, ms
            wblf(mg, k) = wblfint
          end do
        end if
      end do

      do k=1,ms
      do 316 mg=1,ln2
      if(tester3(mg))then
       wbknew=max(wblf(mg,k)*ssatcur(mg,k)+wbice(mg,k,lg),
     &                 0.00001)
       if(wbknew.lt.wbice(mg,k,lg))
     &      wbice(mg,k,lg)=max(0.00001,wbknew-0.01)
       rnof1(mg)=rnof1(mg)+max(wbknew-pssat(mg,lg),0.)
     & *1000*zs(k)
       wb(mg,k,lg)=min(wbknew,pssat(mg,lg))
      endif
  316 continue
      enddo

      do 385 mg=1,ln2
      if(land(mg))then
       dwb=max( min(wb(mg,ms,lg)-psfc(mg,lg),.9*wb(mg,ms,lg)-
     &     wbice(mg,ms,lg))*pcdr3(mg,lg)*
     &           zs(ms)*1000., 0.0)/1000.0/zs(ms)/86400.0
c       dwb=max( (wb(mg,ms,lg)-(psfc(mg,lg)))*pcdr3(mg,lg)*
c     &           zs(ms)*1000., 0.0)/1000.0/zs(ms)/86400.0
       dwb=dwb*0.5
          if(lakesg(mg,lg))dwb=0.0
       rnof2(mg)=zs(ms)*1000.0*dwb*dt
       wb(mg,ms,lg)=wb(mg,ms,lg)-dwb*dt


       runoff(mg)=rnof1(mg)+rnof2(mg)
       perc(mg)=rnof2(mg)

      if( isflag(mg,lg).eq.0) then 
      tg(mg)=tsigmf(mg)*tgf(mg,lg)+
     &    (1.-tsigmf(mg))*tggsl(mg,1,lg)
      else
      tg(mg)=tsigmf(mg)*tgf(mg,lg)+(1.-tsigmf(mg))*
     &       tggsn(mg,1,lg)
      endif

C---- GLACIER FORMATION
c     SNOWD(MG)=MIN( SNOWD(MG), 400.0)
      if(snowd(mg).gt.400.0)then
        rnof5=snowd(mg)-400.0
        runoff(mg)=runoff(mg)+rnof5
c----  change local tg to account for energy - clearly not best method
        tggsn(mg,1,lg)=tggsn(mg,1,lg)-rnof5*hlf/sgamm(mg,1)
        snowd(mg)=400.0
      end if

       wg(mg)=max(0.,wb(mg,1,lg)-pswilt(mg,lg))
       wg2(mg)=0.0
      endif
 385  continue

c     if(newriver)then

c      do ns=1,2
c        lgns=(ns-1)*lg+(lat2+1-lg)*(2-ns)
c        do mg=1,lon
c         ma=mg+(ns-1)*lon
c         if(lakesg(mg,ns,lg))then
c           write(6,6666)mg,lgns,rnof1(mg),rnof2(mg)
c    & ,wbice(ma,ms,lg),twf(mg)*dt
c           write(6,6666)mg,lgns,(wb(ma,k,lg),k=1,ms)
c6666 format(2i4,7(2x,e13.6))
c         endif
c        enddo
c      enddo

c     endif

      do k=1,ms
      do mg=1,ln2
        if(land(mg))then
          wg2(mg)=wg2(mg)+zsfrac(k)
     &       *max(0.,wb(mg,k,lg)-pswilt(mg,lg))
        endif
      enddo
      enddo

      do 3852 ns=1,2
      ma=(ns-1)*lon
      zgwx=0.0
      ztslx=0.0
      zrunofx=0.0
      do 3851 mg=1+ma,lon+ma
      if(land(mg))then
      zgwx=zgwx+wb(mg,1,lg)/pssat(mg,lg)
      ztslx=ztslx+tg(mg)
      zrunofx=zrunofx+runoff(mg)
      endif
 3851 continue
      zgw(ns)=zgwx
      ztsl(ns)=ztslx
 3852 zrunof(ns)=zrunofx

c Set up HBG's instantaneous soil moisture maps (scaled correctly for NSIB)
c This overwrites the values set up in surfset for old soil scheme

      if ( (mod(mins+int(mstep, 8),1440_8).eq.0_8) .and. gwicm ) then
        do 388 mg=1,ln2
          if(land(mg))then
            kx=1+nint(10*wg(mg)/pssat(mg,lg))
            chwg(mg,lg)=cint(kx)
            kx=1+nint(10*wg2(mg)/pssat(mg,lg))
            chwb(mg,lg)=cint(kx)
          endif
 388    continue
      endif

      return
      end
c
      subroutine smoisture(lg, wb, fwtop, dt, ism)

      implicit none

!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /SOILPR/ )
!$OMP THREADPRIVATE ( /SURFBC/ )

C Global parameters
      include 'PARAMS.f'
      integer nmeth
      parameter (nmeth=1)  ! 1 for full implicit,
c     parameter (nmeth=2)    2 for simple implicit D, implicit K
c      NOTE : nmeth=2 has some instability : not finalized : do not use

C Argument list
      integer lg
      real wb(ln2,ms)
      real fwtop(ln2)
      real dt
      logical ism(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'LOGIMSL.f'

      real tgsoil,ssatcur,gamm,coef,scond1s,cap1s,sgamm,sdepth
     & ,hs1tm,hs1tp
      common/soilpr/tgsoil(ln2,ms),ssatcur(ln2,ms+1),gamm(ln2,ms)
     &,coef(ln2,ms+4),scond1s(ln2,3),cap1s(ln2,3)
     &,sgamm(ln2,3),sdepth(ln2,3),hs1tm(ln2,9),hs1tp(ln2,9)

      real egg,condxg,tsigmf,evapxf,ewwp,pmcmax
      common/surfbc/egg(ln2),condxg(ln2),tsigmf(ln2),evapxf(ln2)
     &             ,ewwp(ln2),pmcmax(ln2)

C Global data blocks
      real prhos,pcss,pcnsd,pbch,pssat,pswilt,phsbh,psfc,psucs
     & ,phyds,pcdr3,prlaim,prlais,pslveg,pscveg
      integer ip2bp3,ipbp2
      common/psoilvegpr/prhos(ln2,lat),pcss(ln2,lat),
     &pcnsd(ln2,lat),pbch(ln2,lat),pssat(ln2,lat),
     &pswilt(ln2,lat),phsbh(ln2,lat),ip2bp3(ln2,lat),
     &ipbp2(ln2,lat),psfc(ln2,lat),psucs(ln2,lat),
     &phyds(ln2,lat),pcdr3(ln2,lat),prlaim(ln2,lat),
     &prlais(ln2,lat),pslveg(ln2,lat),pscveg(ln2,lat)

      real zs,zsh,zsdzs,zsfrac
      common/soilzs/zs(ms),zsh(ms+1),zsdzs(ms),zsfrac(ms)

C Local work arrays and variables
      real dtt(ln2,ms)
      real swb(ln2,ms)
      real  c1(ln2,ms), c2(ln2,ms), c3(ln2,ms), rhs(ln2,ms)
      real wc1(ln2,ms),wc2(ln2,ms),wc3(ln2,ms),wrhs(ln2,ms)
      real z1(ln2,ms+1),z2(ln2,ms+1),z3(ln2,ms+1)
      logical tester1(ln2),flag1
      INTEGER INDEX1(ln2) ! INDEX FOR COMPRESS AND EXPAND
      real pwb_min(ln2),z1mult(ln2,ms+1)

      integer i
      integer k
      integer mg
      integer ml
      integer ms0
      integer npts

      real fact
      real fact2
      real hydss
      real pwb
      real pwbb
      real wbh

C Local data, functions etc
      real rhowat
      data rhowat /1000.0/
      real eps
      data eps /1.0e-9/
      save eps

C Start code : ----------------------------------------------------------
c
c        solves implicit soil moisture equation
c
c        fwtop  - water flux into the surface (precip/evap)
c        dt   - time step 
c        ism    - true if all soil layers frozen, otherwise false
c

      do 40 mg=1,ln2
      tester1(mg)=.false.
      if(land(mg))then
        if(.not.ism(mg)) tester1(mg)=.true.
      endif
   40 continue

      call checkl(tester1,ln2,flag1)
      if(.not.flag1)return

      do 52 mg=1,ln2
      if(tester1(mg))then
        dtt(mg,1)=1./(zs(1)*ssatcur(mg,1))*dt
        z1(mg,1)=0.0
        z2(mg,1)=0.0
        z3(mg,1)=0.0
        z1(mg,ms+1)=0.0
        z2(mg,ms+1)=0.0
        z3(mg,ms+1)=0.0
c       wbh(mg,1) = min(1.,wb(mg,1))
c       wbh(mg,ms+1)=min(1.,wb(mg,ms))
      endif
   52 continue
c
      if(nmeth.eq.1)then  ! full implicit method

       do k=2,ms
       c1(1,k)=0.0 ! Dummy to force NEC vectorizing
      do 54 mg=1,ln2
      if(tester1(mg))then
         dtt(mg,k)=1./(zs(k)*ssatcur(mg,k))*dt

       wbh=min(1.,zsdzs(k)*wb(mg,k-1)+(1.-zsdzs(k))*wb(mg,k))
       pwb = max(phsbh(mg,lg)*(wbh**(ipbp2(mg,lg)-1)),eps)
       hydss=phyds(mg,lg)
       pwbb  = pwb*wbh
       z1(mg,k)=(-1.+ip2bp3(mg,lg))*hydss*(wbh**ip2bp3(mg,lg))
     &   -ipbp2(mg,lg)*pwb*wbh*(wb(mg,k)-wb(mg,k-1))/zsh(k)
       z2(mg,k)=-ip2bp3(mg,lg)*hydss*(wbh**(ip2bp3(mg,lg)-1))
     &   + ipbp2(mg,lg)*pwb*(wb(mg,k)-wb(mg,k-1))/zsh(k)
       z3(mg,k)=pwbb/zsh(k)
      endif
   54 continue
       enddo      
c
c For subsequent use of routine trim3x : set values 
c
      do  k=1,ms
      do 56 mg=1,ln2
      if(tester1(mg))then
        c1(mg,k) = dtt(mg,k)*(z2(mg,k)*0.5*zs(k)/zsh(k)-z3(mg,k) )
        rhs(mg,k) = wb(mg,k)+dtt(mg,k)*( z1(mg,k+1) - z1(mg,k) )
      endif
   56 continue
      enddo

      do 58 mg=1,ln2
      if(tester1(mg))then
        rhs(mg,1)=rhs(mg,1)+dtt(mg,1)*(fwtop(mg))/rhowat
        k=ms
        c2(mg,k)=1.+dtt(mg,k)*( z2(mg,k)*.5*zs(k-1)/zsh(k)+z3(mg,k) )
      endif
   58 continue
c
      do k=1,ms-1
        ml = max (k-1,1)
      do 60 mg=1,ln2
      if(tester1(mg))then
        c3(mg,k)=dtt(mg,k)*(-z2(mg,k+1)*.5*zs(k)/zsh(k+1)-z3(mg,k+1))
        c2(mg,k)=1.+dtt(mg,k)*(-z2(mg,k+1)*.5*zs(k+1)/zsh(k+1)
     *     + z2(mg,k)*.5*zs(ml)/zsh(k)+z3(mg,k+1)+z3(mg,k) )
      endif
   60 continue
      enddo

      else                         ! (nmeth.eq.2)  
c                                     simple implicit D, implicit K 

      do mg=1,ln2
      if(tester1(mg))then
        do k=1,ms
          dtt(mg,k)=1./(zs(k)*ssatcur(mg,k))*dt
          c1(mg,k)=0.
          c2(mg,k)=1.
          c3(mg,k)=0.
        enddo
        pwb_min(mg)=(pswilt(mg,lg)/pssat(mg,lg))**ipbp2(mg,lg)
        z1mult(mg,1)=0.      ! corresponds to 2b+3
        z1mult(mg,ms+1)=0.   ! corresponds to 2b+3
      endif
      enddo

      do  mg=1,ln2
      if(tester1(mg))then

       do k=2,ms
        z1mult(mg,k)=ip2bp3(mg,lg)   ! corresponds to 2b+3
        wbh=min(1.,zsdzs(k)*wb(mg,k-1)+(1.-zsdzs(k))*wb(mg,k))
        fact=wbh**(ipbp2(mg,lg)-1)
        pwb = phsbh(mg,lg)*max(pwb_min(mg),wbh*fact)
        fact2=fact*fact
        z1(mg,k)=phyds(mg,lg)*fact2
        z3(mg,k)=pwb/zsh(k)
        c1(mg,k)=-dtt(mg,k)*z3(mg,k)
        c3(mg,k-1)=-dtt(mg,k-1)*z3(mg,k)
       enddo   !ms

       do k=1,ms
        c2(mg,k)=1.-c1(mg,k)-c3(mg,k)
       enddo  !ms

       do k=2,ms
        c1(mg,k)=c1(mg,k)
     &          -dtt(mg,k)*z1mult(mg,k)*z1(mg,k)*zs(k)/(zs(k)+zs(k-1))
        c3(mg,k-1)=c3(mg,k-1)
     &      -dtt(mg,k-1)*z1mult(mg,k)*z1(mg,k)*zs(k-1)/(zs(k)+zs(k-1))
       enddo   !ms

       do k=1,ms
        c2(mg,k)= c2(mg,k)
     &  -dtt(mg,k)*z1mult(mg,k)*z1(mg,k)*zs(k-1)/(zs(k)+zs(k-1))
     &  +dtt(mg,k)*z1mult(mg,k+1)*z1(mg,k+1)*zs(k+1)/(zs(k)+zs(k+1))
       enddo   !ms

       do k=1,ms
        wbh=min(1.,zsdzs(k)*wb(mg,k-1)+(1.-zsdzs(k))*wb(mg,k))
        z1(mg,k)= z1(mg,k)*wbh
       enddo  !ms

!          the following top & bottom b.c.'s will preserve a uniform column
       z1(mg,1) =min(z1(mg,2),z1(mg,ms))   !  N.B. z1 are here +ve
       z1(mg,ms+1)=z1(mg,1)
       do k=1,ms
        rhs(mg,k) = wb(mg,k)+dtt(mg,k)*( (z1mult(mg,k+1)-1.)*z1(mg,k+1)
     &              - (z1mult(mg,k)-1.)*z1(mg,k) )
       enddo  !ms

       rhs(mg,1)=rhs(mg,1)+dtt(mg,1)*fwtop(mg)/rhowat

       endif !tester
       enddo ! mg

      endif ! nmeth
       
      npts = 0
      do mg=1,ln2
        if(tester1(mg))then
          npts = npts + 1
          INDEX1(npts) = mg
        endif
      enddo
c
      ms0 = ms
      do k=1,ms0
      do i=1,npts
        wc1(i,k)  = c1(INDEX1(i),k)
        wc2(i,k)  = c2(INDEX1(i),k)
        wc3(i,k)  = c3(INDEX1(i),k)
        wrhs(i,k)  = rhs(INDEX1(i),k)
      enddo
      enddo
 
      call trim3x(wc1,wc2,wc3,wrhs,swb,ms0,npts)

      do k=1,ms
      do 10 i=1,npts
        wb(INDEX1(i),k)=swb(i,k)
10    continue
      enddo

      return
      end
