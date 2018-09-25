c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified by HBG to prevent negative snowmelt when sublimation removes all
c available snow.
c SJP 2005/02/24
c
c Further changes by Hal Gordon to improve freshwater conservation.
c SJP 2003/08/16
c
c Conservation error fixed: when white ice forms, the model currently converts
c part of the snow into ice AND melts it. The former of these is correct, the
c latter is not. The line in which SNOWMI is incremented is therfore
c commented out.
c SJP 2003/06/28
c
c HFACA removed from /FICECON/, as this array is never used.
c SJP 2003/06/20
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: seaice.f,v $
c Revision 1.46  2001/10/12 02:06:58  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.45  2001/02/22 05:34:42  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.44  2001/02/12 05:39:53  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.43  2000/11/14 03:11:39  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.42  1999/05/20 06:23:55  rot032
c HBG changes to V5-2
c
c Revision 1.41  1998/12/10  00:55:49  ldr
c HBG changes to V5-1-21
c
c Revision 1.40  1997/12/23  00:23:37  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.39  1997/12/17  23:22:54  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.38.1.1  1997/12/19  02:03:16  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.38  1996/10/24  01:03:13  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.37  1996/06/13  02:08:09  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.36  1996/03/21  03:19:03  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.35  1995/08/14  05:22:36  ldr
c Inserted definition of jj so that groice diagnostic works correctly.
c
c Revision 1.34  1995/07/31  05:31:20  ldr
c Fix up the "Growth of ice" diagnostic.
c
c Revision 1.33  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.29.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.32  1995/03/14  01:21:44  ldr
c Changes to qcloud: incorporate hlfusion in snow that hits ground; use
c preset vice in slow icefall calc. - this is run E85, still with small
c fluxi array in icefall.f.
c
c Revision 1.31  1994/08/08  17:22:30  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.30  94/08/04  16:56:33  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.29  94/03/30  12:34:59  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.28  93/12/17  15:33:42  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.27  93/12/17  11:51:44  ldr
c Changes to V4-4-45l from HBG for coupled model
c 
c Revision 1.26  93/11/29  11:38:37  ldr
c Changes to V4-4-32l from HBG for coupled model
c 
c Revision 1.25  93/11/03  11:44:30  ldr
c Changes to V4-4-24l from HBG for coupled model
c 
c Revision 1.24  93/08/19  15:10:41  ldr
c Minor cosmetic changes.
c 
c Revision 1.23  93/08/13  14:43:39  ldr
c Minor changes to get coupled model to run on SGI.
c 
c Revision 1.22  93/08/10  15:27:58  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c
c     INPUT/OUTPUT
c     Input:   from common/datice in DATICE1.f
c                  hoice - heat flux from ocean to ice
c
c              from common/ficecon in this subroutine
c                  flxia - subice heat input
c
c              from common/floes in this subroutine
c                  tmix - temp at which adjustment of computed ice melt occurs
c
c              from common/timex in TIMEX.f
c                  dt - time step in seconds
c
c              from arguments
c                  blg - Stefan Boltzmann const*(tg)**4 
c                  condx - total precipiation 
c                  gamms - density*specific heat*depth (J/m**2/K)
c                  degdt - derivative of eg with respect to temp
c                  dfgdt - derivative of fg with respect to temp
c                  preci - ice precipitation
c                  sg - solar absorbed in ground 
c                  lg - latitude index
c
c     Output:  from common/logimsl in LOGIMSL.f
c                  mlo - mixed layer ocean
c
c              from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c                  just_melted - if T, whole grid point just melted
c
c              from arguments
c                  il - longitude grid point
c                  snowmi - snow melt  sublmi - sublimation
c
c     In/Out:  from common/floes in this subroutine
c                  dic - ice thickness  dsn - snow thickness
c                  sto - stored internal heat
c                  tn0,tn1,tn2,tni - internal ice temps at varying levels
c                  pl - proportion of leads
c
c              from common/logimsl in LOGIMSL.f
c 
      subroutine seaice(lg,sg,rg,blg,degdt,dfgdt,eg,fg,
     &                  condx,gamms,preci,
     &                  tg,siced,snowd,sublmi,snowmi,
     &                  il)

      implicit none

!$OMP THREADPRIVATE ( /FLOES/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      integer nkmax
      parameter (nkmax=2)

C Argument list
      integer lg
      real sg(ln2)
      real rg(ln2)
      real blg(ln2)
      real degdt(ln2)
      real dfgdt(ln2)
      real eg(ln2)
      real fg(ln2)
      real condx(ln2)
      real gamms(ln2)
      real preci(ln2) !Frozen precipitation in mm water/tdt if (qcloud)
      real tg(ln2)
      real siced(ln2)
      real snowd(ln2)
      real sublmi(ln2)
      real snowmi(ln2)
      integer il(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      real dsn,dic,tn0,tn1,tn2,tni,sto,pl,tmix
     & ,plsum,vol,tflux,sflux,fflux
      common/floes/dsn(ln2),dic(ln2),tn0(ln2),
     &tn1(ln2),tn2(ln2),tni(ln2),sto(ln2),pl(ln2),tmix(ln2),
     &plsum(ln2),vol(ln2),tflux(ln2),sflux(ln2),fflux(ln2)
      include 'LOGIMSL.f'

C Global data blocks
      include 'DATICE1.f'
      include 'LSMI.f'
      include 'SURF1.f'
      include 'TIMEX.f'
      real flxia
      common/ficecon/flxia(lat,2)
      real redice,groice,hav
      common/red/redice(ln2,lat),groice(ln2,lat),hav(ln2,lat)

C Local work arrays and variables
      real flx(ln2),dsict(ln2)
      real whitei(ln2)
      integer nkx(ln2)
      real t(2),tp(0:3),fl(0:2),z(0:2),zp(0:2)
      logical flagt(2)

      integer k
      integer kice
      integer ma
      integer mg
      integer nk
      integer nkm1
      integer ns

      real bot
      real cls
      real con
      real con0
      real con1
      real con2
      real deltat
      real dfcond
      real dfdt
      real dhb
      real dhi
      real dhs
      real dirad
      real en
      real excess
      real eye
      real eyesav
      real f
      real fa
      real fact
      real fb
      real fs
      real f0
      real hb
      real hi
      real hinew
      real hn
      real hs
      real hsold
      real htdown
      real htup
      real imelt
      real qmax
      real qneed
      real qstor
      real rhin
      real rhi2n
      real shir
      real smelt
      real snowfall
      real subl
      real sumt
      real tb
      real ti
      real titm
      real tpast
      real ts
      real tstz
      real t0
      real xxx

cx    real dsnold(ln2),sublsn,xxc1,xxc2 ! for checking snow conservation

C Local data, functions etc
c     qs,qi,qb are latent heats of fusion in J/m**3
c     (hlf=3.35e5 in other code is in J/Kgm)
      real sk,yk,qs,qi,qb,hmin,tz,tm
      data sk/0.30976/,yk/2.03439/,qs/3.35e8/,
     & qi/3.35e8/,qb/3.35e8/,hmin/0.10/,
     & tz/0.0/,tm/-0.15/
      real shs,shi
      data  shs/6.9069e5/,shi/1.8837e6/
      real pw,rhowi,ps
      data  pw/1025./,rhowi/900./,ps/330./
      real tfrz
      data  tfrz/273.15/

C Start code : ----------------------------------------------------------

c....
c.... model recoded in degrees C instead of degrees K
c....
c   ***********************************************
c   dictionary of names
c   ***********************************************
c   stefbo   stefan-boltzmann constant
c   sk       conductivity of snow
c   yk       conductivity of ice
c   qs       heat of fusion of snow (Joules/m**3)
c   qi       heat of fusion of surface ice
c   qb       heat of fusion of bottom ice
c   hl       heat of evaporation (Joules/Kgm)
c   hi       ice thickness
c   hs       snow thickness
c   ts       surface temperature
c   ti       temperature at snow-ice interface
c   tb       temperature at bottom of ice
c   hmin     minimum allowable layer thickness
c   eye      fraction of penetrating solar radiation
c   fb       oceanic heat flux upward
c   shs,shi  in Joules/m**3/K (rho*Cp)
c   cls      =1+(hlf/hl) - factor for latent heat of fusion
c   ***********************************************

c Set nominal eye value (eyesav)
      eyesav=0.35

c  For evaporation to occur the snow must have been melted. The latent
c  heat of fusion has not yet been included in the surface temperature.
c  Note that evap is calculated in mms per step using the
c  latent heat of evaporation (hl=2.5104E+06 Joules/Kgm).
      cls=1.0+(hlf/hl)
c---- precompute some values
         do 11 mg=1,ln2
c        tf(mg)=0.
c        sf(mg)=0.0
c        ff(mg)=0.0
         snowmi(mg)=0.0
         sublmi(mg)=0.0
         whitei(mg)=0.0
cx       dsnold(mg)=dsn(mg)
  11     dsict(mg)=0.0

c**** Sub-ice flux of heat from ocean 
c**** The flxia values are preset in ICECON
      do 14 mg=1,lon
      flx(mg)=flxia(lg,1)+hoice(mg,lg)
   14 flx(mg+lon)=flxia(lg,2)+hoice(mg+lon,lg)

c.... set up nk values
         do 16 mg=1,ln2
         if(cice(mg))then
           il(mg)=1
           hi=dic(mg)
c          nk=2
c          if (hi.lt.2.*hmin) nk=1
c          if (hi.lt.hmin) nk=0
c          nk=min(nk,nkmax)
           nk=min(int(hi/hmin),nkmax)
           if (hi.eq.0.0) nk=-1
           nkx(mg)=nk
         end if
  16     continue

      shir=1./shi
c---- end of pre-computes

c---- Main ice updating loop
      flagt(1)=.false.
      flagt(2)=.false.
      do 5000 mg=1,ln2
      if(.not.cice(mg))go to  5000

      tb=tmix(mg)-tfrz
      tpast=tg(mg)-tfrz
      dirad=0.0
cx    sublsn=0.0
      snowfall=0.0

      if(pl(mg).ge.1.00)goto 7090 ! No ice left

      fb=flx(mg)
      zp(0)=0.0
      z(0)=0.0
      z(1)=0.0
      z(2)=0.0
      smelt=0.0  
c---- condx(mg) = precipitation in mms of water per step
c     Note that even if running qcloud, condx is total precip (LDR)
      snowfall=condx(mg)*0.01
c---- evap = evaporation in mms of water per step
c---- evap =dt*eg(mg)/hl
      subl=(dt*eg(mg)/hl)*0.01
      t(1)=tn1(mg)-tfrz
      t(2)=tn2(mg)-tfrz
      osnowd(mg,lg)=snowd(mg) ! for snow age/albedo
      t0=tn0(mg)-tfrz
      ti=tni(mg)-tfrz
      qstor=sto(mg)

      hi=dic(mg)
      hs=dsn(mg)
c Add any precipitation as snowfall. (Any rain is assumed to be
c snow, and the latent heat release is added later)
      if(snowfall.gt.0.0)then
        if(hs.eq.0.0)t0=ti
        hs=hs+snowfall
      endif
      nk=nkx(mg)
      en=float(nk)
      nkm1=nk-1

      if(nk.le.0) go to 1355 ! Go to zero layer & no ice code

      rhin=en/hi
      rhi2n=2.*rhin
      qmax=qi*0.5*(hi-hmin)
      if(nk.eq.1)then
        if(qstor.gt.0)then
c     adjust ice thickness by the remaining reservoir in qstor
          hi=max(0.0,hi-qstor/qi)
        endif
        qstor=0.0
        if(hi.lt.hmin+0.001)t(1)=0.5*(ti+tb)
        t(2)=tb
      endif

c   jump to 200 if there is no snow on the ice
      if(hs.eq.0.00) go to 200

      con=2.*sk/hs
      if(hs.lt.0.05) then
        con=sk/(hs+(hi*sk)/(yk*en*2.))
        t0=t(1)
      endif
c
c   Solve for snow surface temperature (implicit time scheme)
c
c     Solar Sg(in), LWR Rg(out), Heat Flux Fg(out), Evap Eg(out)
c
c     ####################### Surface energy flux f ##### and fa(up)##
c     ...ts = Snow temp for thin layer near surface................. S
c     -------------------------------------------------------------- N
c                                                             fs(up) O
c     ...t0 = Mid point temp of snow layer.......................... W
c
c     ## ti = Temp at snow/ice interface #############################
c                                                             f0(up) I
c     ...t(1) = Mid point temp of first ice layer................... C
c                                                                    E
c     ########################################################fl(1)###
c                                                                    I
c     ...t(2) = Mid point temp of second ice layer.................. C
c                                                                    E
c     ################################################################
c
c Flux of energy into surface
      f=sg(mg)-rg(mg)-fg(mg)-cls*eg(mg)
c If(qcloud), condx-preci gives the liquid precip.
c If(.not.qcloud), preci is zero.
c Add flux of heat due to converting any rain to snowfall
      dfcond=hlf*(condx(mg)-preci(mg))/dt ! rain (mms/step) to W/m**2
      f=f+dfcond
c Get total flux into thin surface layer by adding flux from below
      f=f+con*(t0-tpast)
c Implicit correction terms
      dirad=4.0*blg(mg)/tg(mg)
c     dfdt=dirad+cls*degdt(mg)+dfgdt(mg)+con
      dfdt=dirad+con
      bot=(gamms(mg)/dt)+dfdt
      ts=tpast+ f/bot  
      tstz=max(0.0,ts-tz) ! >0 only when melting
      ts=min(ts,tz)       ! melting condition ts=tz

c   find temperature at snow-ice interface
      ti=(hi*sk*t0+en*hs*yk*t(1))/(hi*sk+en*hs*yk)
c   Compute flux of energy leaving surface using new ts + implicit terms
c   (but don't include the flux from below as was done for f)
      fa=rg(mg)+fg(mg)+cls*eg(mg)-sg(mg)-dfcond+(dfdt-con)*(ts-tpast)
c   Upward fluxes between various levels below surface
      fs=con*(t0-ts)            ! Middle of snow layer to surface
      f0=rhi2n*yk*(t(1)-ti)     ! Middle of first ice layer to snow interface
      fl(1)=yk*(t(2)-t(1))*rhin ! Between ice layers 2 and 1

      if(hs.lt.0.05)then
c   Ammendments if surface layer < 5cms
        hs=hs+0.001                     ! Add some snow
        snowmi(mg)=snowmi(mg)-0.001*0.1 ! Accounts for added snow
        ti=(hi*sk*ts+2.*en*hs*yk*t(1))/(hi*sk+2.*en*hs*yk)
        fs=sk*(ti-ts)/hs
        f0=fs
        t0=(ti+ts)*0.5
      endif

c   Surface evap/sublimation (can be >0 or <0)
c   Note : dt*eg/hl in Kgm/m**2 => mms of water
c    0.001*dt*eg/hl => m of water,
c     ==> 10*0.001*dt*eg/hl gives m of "snow"
c     subl=0.01*dt*(eg(mg)+(ts-tpast)*degdt(mg))/hl
      subl=0.01*dt*eg(mg)/hl          ! m of "snow"
cx    sublsn=sublsn+subl
c   Snow melt (smelt >= 0)
      smelt=10.0*tstz*gamms(mg)/qs    ! m of "snow"
c   Change the snow thickness
      dhs=smelt+subl
      hsold=hs

      if(dhs.gt.hs)then
        hi=hi-0.1*(dhs-hs)
        smelt=max((hs-subl),0.0)
        hs=0.0
      else
        hs=hs-dhs
      endif

      snowmi(mg)=snowmi(mg)+smelt*0.1 ! m of water

      subl=0.0 
      smelt=0.0
      dhs=0.0
      dhi=0.00

c   Update the mid-level snow temperature
      if(hsold.ge.0.05) t0=t0+dt*(f0-fs)/(hs*shs)
c   Remove very thin snow
      if(hs.lt.0.001)then
        snowmi(mg)=snowmi(mg)+hs*0.1
        hs=0.0         ! Need heat release
      endif
      if(hs.eq.0.0) t0=273.05-tfrz

      t0=min(t0,tz)
c  (Should melt snow if t0>tz, but apparently rare)

      fl(0)=f0
      if(nk.eq.1) go to 510

c   update temperature in top layer of ice
      tp(1)=min((t(1)+dt*(fl(1)-f0)*rhin*shir),tm)
c  (Should melt ice if tp(1)>tm, but apparently rare)
      if(qstor.eq.0..or.tp(1).eq.tm) go to 400
c   use stored heat in brine pockets to keep temperature at -0.1 until
c   heat is used up
      qneed=(tm-tp(1))*shi*(hi/en) ! J/m**2
      if(qstor.le.qneed)then
        fact=qstor/qneed
        tp(1)=tm*fact+tp(1)*(1.-fact)
        qstor=0.
      else
        qstor=qstor-qneed
        tp(1)=tm
      endif
      go to 400

c   come here when ice is snow free
 200  continue

c   allow penetrating radiation to be stored up to some maximum value
      eye=eyesav
      if(qstor.gt.qmax) eye=0.
      qstor=qstor+dt*eye*sg(mg)
      con=yk*rhi2n
c
c   Solve for surface ice temperature (implicit time scheme)
c
c     Solar Sg(in), LWR Rg(out), Heat Flux Fg(out), Evap Eg(out)
c
c     ####################### Surface energy flux f ##### and fa(up)##
c     ...ti = Ice temp for thin layer near surface..................
c     -------------------------------------------------------------- I
c                                                             f0(up) C
c     ...t(1) = Mid point temp of first ice layer................... E
c                                                                    
c     ########################################################fl(1)###
c                                                                    I
c     ...t(2) = Mid point temp of second ice layer.................. C
c                                                                    E
c     ################################################################
c
c Flux of energy into surface
      f=(1-eye)*sg(mg)-rg(mg)-cls*eg(mg)-fg(mg)
c Get total flux into thin surface ice layer by adding flux from below
      f=f+con*(t(1)-tpast)
c Implicit correction terms
      dirad=4.0*blg(mg)/tg(mg)
c     dfdt=dirad+cls*degdt(mg)+dfgdt(mg)+con
      dfdt=dirad+con
      bot=(gamms(mg)/dt)+dfdt
      ti=tpast + f/bot  

      titm=max(0.0,ti-tm) ! >0 only when melting
      ti=min(ti,tm)       ! melting condition ti=tm
      ts=ti
c   Compute flux of energy leaving surface using new ti + implicit terms
c   (but don't include the flux from below as was done for f)
      fa=rg(mg)+cls*eg(mg)+fg(mg)-(1.-eye)*sg(mg)
     & +((dfdt-con)*(ti-tpast))
c   Upward fluxes between various levels below surface
      f0=con*(t(1)-ti)          ! Middle of first ice layer to thin surface layer
      fl(0)=fa ! ?? f0 ??
      fl(1)=yk*(t(2)-t(1))*rhin ! Between ice layers 2 and 1

c   Surface evap/sublimation (can be >0 or <0)
c     subl=0.01*dt*(eg(mg)+(ti-tpast)*degdt(mg))/hl
      subl=0.01*dt*eg(mg)/hl
      subl=subl*0.1                ! m of ice
      sublmi(mg)=subl
c   Ice melt (imelt >= 0)
      imelt=titm*gamms(mg)/qi      ! m of ice

      dhi=-(subl+imelt)
      if(qstor.gt.qmax)then
        dhi=dhi-(qstor-qmax)/qi
        qstor=qmax
      endif
      subl=0.0
      imelt=0.0

      if(dhi.lt.0.0)go to 350 ! melting

      if(nk.eq.1) go to 510

c   go to update upper ice layer temperature (for nk>=2)
      go to 390

 350  continue
      if(nk.ge.2) go to 390
c   if ice is melting, update thickness when single level
      hi=hi+dhi ! ? Check dhi does not exceed hi
      dhi=0.
      rhin=1./hi
      rhi2n=2.*rhin
      go to 510

 390  continue
      tp(1)=t(1)+dt*shir*(fl(1)-f0)/(hi/en+dhi)

 400  continue
      if(qstor.eq.0..or.tp(1).eq.tm) go to 490
      qneed=(tm-tp(1))*shi*(hi/en)
      if(qstor.le.qneed)then
        fact=qstor/qneed
        tp(1)=tm*fact+tp(1)*(1.-fact)
        qstor=0.
      else
        qstor=qstor-qneed
        tp(1)=tm
      endif

 490  continue
      tp(1)=min(tp(1),tm)
c  (Should melt ice if tp(1)>tm, but apparently rare)

      if(nk.gt.2)then
c   update temperatures of internal ice layers
        do k=2,nkm1
          fl(k)=yk*(t(k+1)-t(k))*rhin
        enddo
        do k=2,nkm1
          tp(k)=t(k)+dt*(fl(k)-fl(k-1))*rhin*shir
        enddo
      endif

 510  continue
c   determine amount of bottom ablation or accretion
      fl(nk)=yk*(tb     -t(nk))*rhi2n
ch    limit the ice growth for cases of poor initial conditions
      dhb=min(0.01,dt*(fl(nk)-fb)/qb)

      hinew=hi+dhi+dhb
c   determine depths below surface of interim and final interfaces
      con0=hinew/en
      do k=1,nk
        z(k)=z(k-1)-con0
      enddo
      con0=hi/en
      zp(1)=-con0-dhi
      if(nk.gt.2)then
        do k=2,nkm1
          zp(k)=zp(k-1)-con0
        enddo
      endif
      if(dhb.le.0.) go to 600
c   construct final values of temperature for case of bottom accretion
      zp(nk)=zp(nkm1)-con0
      tp(nk)=t(nk)+dt*(fl(nk)-fl(nkm1))*rhin*shir
      if(tp(nk).gt.tm)tp(nk)=tm
      tp(nk+1)=tb
      do k=1,nk
        t(k)=((z(k-1)-zp(k))*tp(k)+(zp(k)-z(k))*tp(k+1))*en/hinew
      enddo
      go to 705
 600  continue
c   construct final values of temperature for case of bottom ablation
      zp(nk)=z(nk)
      tp(nk)=t(nk)+dt*shir*(fl(nk)-fl(nkm1))/(hi/en+dhb)
      if(tp(nk).gt.tm)tp(nk)=tm
      tp(0)=ti
      tp(nk+1)=tb
      do k=1,nk
        if(z(k).lt.zp(k).and.zp(k).lt.z(k-1))then
          t(k)=((z(k-1)-zp(k))*tp(k)+(zp(k)-z(k))*tp(k+1))*en/hinew
        elseif(z(k).lt.zp(k-1).and.zp(k-1).lt.z(k-1))then
          t(k)=((z(k-1)-zp(k-1))*tp(k-1)+(zp(k-1)-z(k))*tp(k))*en/hinew
        else
          t(k)=tp(k)
        endif
      enddo
 705  hi=hinew
c   test whether to change number of layers
      htup=float(nk+1)*hmin
      if(hi.ge.htup) go to 1100
      htdown=float(nk)*hmin
      if(hi.le.htdown) go to 1300
      go to 1500
 1100 if(nk.eq.nkmax) go to 1500
c   code to increase number of layers
      nk=nk+1
      nkm1=nk-1
      en=float(nk)
      tp(1)=t(1)
      tp(nk)=t(nkm1)
      if(nk.eq.2) go to 1200
      con1=hi/float(nkm1)
      con2=hi/en
      do k=1,nk
        z(k)=z(k-1)-con1
        zp(k)=zp(k-1)-con2
      enddo
      do k=2,nkm1
        tp(k)=((zp(k-1)-z(k-1))*t(k-1)+(z(k-1)-zp(k))*t(k))/con2
      enddo
 1200 continue
      do k=1,nk
        t(k)=tp(k)
      enddo
      go to 1500

 1300 continue
c   code to decrease number of layers
      nk=nk-1
      nkm1=nk-1
      en=float(nk)
      if(nk.eq.0) go to 1500
      con1=hi/float(nk+1)
      con2=hi/en
      do k=1,nk
        z(k)=z(k-1)-con1
        zp(k)=zp(k-1)-con2
      enddo
      do k=1,nk
        tp(k)=((zp(k-1)-z(k))*t(k)+(z(k)-zp(k))*t(k+1))/con2
      enddo
      do k=1,nk
        t(k)=tp(k)
      enddo

 1500 continue
      go to 1501

 1355 continue

c   when no ice is present, go apply fluxes to upper ocean layer
      if(nk.eq.-1)go to 7090 

c   jump to 2001 if there is no snow on the ice
      if(hs.eq.0.00) go to 2001

      con=sk/(hs+(hi*sk/yk))
c
c   Solve for snow surface temperature (implicit time scheme)
c
c     Solar Sg(in), LWR Rg(out), Heat Flux Fg(out), Evap Eg(out)
c
c     ####################### Surface energy flux f ##### and fa(up)##
c     ...ts = Snow temp for thin layer near surface................. S
c     -------------------------------------------------------------- N
c                                                             fs(up) O
c     ...t0 = Mid point temp of snow layer.......................... W
c
c     ## ti = Temp at snow/ice interface #############################
c                                                             f0(up) I
c     ...t(1) = Mid point temp of ice layer......................... C
c                                                                    E
c     ## tb ##################################################fl(1)###
c
c Flux of energy into surface
      f=sg(mg)-rg(mg)-cls*eg(mg)-fg(mg)
c If(qcloud), condx-preci gives the liquid precip.
c If(.not.qcloud), preci is zero.
c Add flux of heat due to converting any rain to snowfall
      dfcond=hlf*(condx(mg)-preci(mg))/dt ! rain (mms/step) to W/m**2
      f=f+dfcond
c Get total flux into thin surface layer by adding flux from below
      f=f+con*(tb-tpast)
c Implicit correction terms
      dirad=4.0*blg(mg)/tg(mg)
c     dfdt=dirad+cls*degdt(mg)+dfgdt(mg)+con
      dfdt=dirad+con
      bot=(gamms(mg)/dt)+dfdt
      ts=tpast+f/bot  
      tstz=max(0.0,ts-tz) ! >0 only when melting
      ts=min(ts,tz)       ! melting condition ts=tz

c   find temperature at snow-ice interface
      ti=(hi*sk*ts +  hs*yk*tb  )/(hi*sk+   hs*yk)
c   Compute flux of energy leaving surface using new ts + implicit terms
c   (but don't include the flux from below as was done for f)
      fa=rg(mg)+cls*eg(mg)+fg(mg)-sg(mg)-dfcond+(dfdt-con)*(ts-tpast)
c   Upward fluxes between various levels below surface
      fs=sk*(ti-ts)/hs
      f0=fs

c   Surface evap/sublimation (can be >0 or <0)
c     subl=0.01*dt*(eg(mg)+(ts-tpast)*degdt(mg))/hl
      subl=0.01*dt*eg(mg)/hl          ! m of "snow"
cx    sublsn=sublsn+subl
c   Snow melt (smelt >= 0)
      smelt=10.0*tstz*gamms(mg)/qs    ! m of "snow"
c   Change the snow thickness
      dhs=smelt+subl
      hsold=hs

      if(dhs.gt.hs)then
        hi=hi-0.1*(dhs-hs)
        smelt=max((hs-subl),0.0)
        hs=0.0
      else
        hs=hs-dhs
      endif

      snowmi(mg)=snowmi(mg)+smelt*0.1 ! m of water

      subl=0.0
      smelt=0.0
      dhs=0.0
      dhi=0.0

c   Remove very thin snow
      if(hs.lt.0.001)then
        snowmi(mg)=snowmi(mg)+hs*0.1
        hs=0.0         ! Need heat release
      endif
      go to 5100

c   come here when ice is snow free
 2001 continue

c   solve for surface ice temperature
      con=yk/hi
      eye=eyesav
      if(nkmax.ne.0) eye=0.
      if(nkmax.eq.0) eye=eye*0.4
c
c   Solve for surface ice temperature (implicit time scheme)
c
c     Solar Sg(in), LWR Rg(out), Heat Flux Fg(out), Evap Eg(out)
c
c     ####################### Surface energy flux f ##### and fa(up)##
c     ...ti = Ice temp for thin layer near surface..................
c     -------------------------------------------------------------- I
c                                                             f0(up) C
c     ...t(1) = Mid point temp of ice layer......................... E
c
c     ## tb ##################################################fl(1)###
c
c Flux of energy into surface
      f=(1-eye)*sg(mg)-rg(mg)-cls*eg(mg)-fg(mg)
c Get total flux into thin surface ice layer by adding flux from below
      f=f+con*(tb-tpast)
c Implicit correction terms
      dirad=4.0*blg(mg)/tg(mg)
c     dfdt=dirad+cls*degdt(mg)+dfgdt(mg)+con
      dfdt=dirad+con
      bot=(gamms(mg)/dt)+dfdt
      ti=tpast+ f/bot  

      titm=max(0.0,ti-tm) ! >0 only when melting
      ti=min(ti,tm)       ! melting condition ti=tm
      ts=ti
c   Compute flux of energy leaving surface using new ti + implicit terms
c   (but don't include the flux from below as was done for f)
      fa=rg(mg)+cls*eg(mg)+fg(mg)-(1.-eye)*sg(mg)
     & +(dfdt-con)*(ti-tpast)
c   Upward flux between various levels below surface
      f0=yk*(tb-ti)/hi

c   Surface evap/sublimation (can be >0 or <0)
c     subl=0.01*dt*(eg(mg)+(ti-tpast)*degdt(mg))/hl
      subl=0.01*dt*eg(mg)/hl
      subl=subl*0.1                ! m of ice
      sublmi(mg)=subl
c   Ice melt (imelt >= 0)
      imelt=titm*gamms(mg)/qi      ! m of ice

      dhi=-(subl+imelt)
      subl=0.0
      imelt=0.0

 5100 continue
c   determine amount of bottom ablation or accretion
ch    limit the ice growth for cases of poor initial conditions
      dhb=min(0.01,dt*(f0-fb)/qb)
c   update ice thickness
      hi=hi+dhi+dhb
      if(hi.lt.0.001)hi=0.0
      if(hi.lt.hmin.or.nkmax.eq.0) go to 7070
      nk=1
      nkm1=0
      en=1.0
      t0=0.5*(ts+ti)
      t(1)=0.5*(ti+tb)
      t(2)=tb
      qstor=0.0
      go to 1501

 7070 continue
c  Remove very thin snow
      if(hs.lt.0.001)then
        snowmi(mg)=snowmi(mg)+hs*0.1
        hs=0.0         ! Need heat release
      endif
      if(hi.le.0.0)go to 7090
      go to 1501

 7090 continue
c   when ice disappears, predict ocean temperature until freezing again
      snowmi(mg)=snowmi(mg)+hs*0.1
      hs=0.00
      hi=0.00
      deltat=tb-tpast
      tg(mg)=tb+tfrz
c     eg(mg)=eg(mg)+deltat*degdt(mg)
c     fg(mg)=fg(mg)+deltat*dfgdt(mg)
      rg(mg)=rg(mg)+deltat*dirad
      pl(mg)=0.0
      il(mg)=0
      cice(mg)=.false.
      mlo(mg)=.true.
      imsl(mg,lg)=2
      just_melted(mg,lg)=.true.
      go to 4980

 1501 continue
      tn1(mg)=t(1)+tfrz
      tn2(mg)=t(2)+tfrz
      tn0(mg)=t0+tfrz
      tni(mg)=ti+tfrz
      sto(mg)=qstor

c    White ice formation
c  Compute buoyancy - if weight of snow is heavy enough, the
c  snow+ice will sink, and part of the snow will become ice.
c  When part of the snow is converted to ice, the snow will
c   need salinity (at value for ice = 0.01 = Sice). This is
c   equivalent to the reverse of sublimation at the surface
c   which drops the salt content of the sublimated ice into
c   the ocean below. Thus use whitei() to account for the salt
c   needed from the ocean to change the salt-free snow into
c   salty ice. Must therefore remove salt from ocean, i.e.
c   make it fresher (put amount into whitei() ).
c This should conserve ocean salinity
      xxx=hi+hs-(ps*hs+rhowi*hi)/pw
      if(xxx.lt.hs) then
        excess=(hs-xxx)*0.1          ! m of water
        whitei(mg)=whitei(mg)+excess ! m of water
        hn=hi+excess
        hb=(0.5*hi)+excess ! Assume 2 levels of ice, hence 0.5*hi
        tn1(mg)=((0.5*hi)*t(1)+excess*t0)/hb + tfrz
        hs=xxx
        hi=hn
      endif 

      dsn(mg)=hs
      groice(mg,lg)=groice(mg,lg)+(hi-dic(mg))
      dic(mg)=hi
      tg(mg)=ts+tfrz
      deltat=ts-tpast
c     eg(mg)=eg(mg)+deltat*degdt(mg)
c     fg(mg)=fg(mg)+deltat*dfgdt(mg)
      rg(mg)=rg(mg)+deltat*dirad
      
      siced(mg)=hi
      snowd(mg)=hs*100.0

c---- Snow depth limitation and conversion to ice
      if(snowd(mg).gt.200.0)then
c
c   snowd in mms of water (equivalent to cms of "snow")
c Excess snow over 200mms (of water = 2m of snow) removed, 
c        and replaced by ice
c   See comments on "White ice" formation above:
c   Use whitei() to account for the salt
c   needed from the ocean to change the salt-free snow into
c   salty ice. Must therefore remove salt from ocean, i.e.
c   make it fresher (put amount into whitei() ).
c This should conserve ocean salinity
c
        excess=(snowd(mg)-200.0)*0.001  ! m water
c       print *,'Removing excess snow(m)=',excess
        siced(mg)=siced(mg)+excess
        whitei(mg)=whitei(mg)+excess    ! m of water
        snowd(mg)=200.0
        dsn(mg)=2.0
        dic(mg)=siced(mg)
      end if

c---- Ice depth limitation
      if(siced(mg).gt.6.0)then
c.... dsict is change in temp at point due to excess ice removal
        dsict(mg)=(siced(mg)-6.0)*1000.0*(1-pl(mg))*hlf/(dt*bot)
c.... limit this change (due to some initial values being >> 6.0m)
        dsict(mg)=min(dsict(mg),1.0)
        ns=(mg-1+lon)/lon
        flagt(ns)=.true.
        siced(mg)=6.0
        dic(mg)=siced(mg)
      endif
c----
 4980 continue

cx Check snow conservation
cx     xxc1=(dsnold(mg)+snowfall-sublsn-dsn(mg))*0.1
cx     xxc2=snowmi(mg)+whitei(mg)
cx     if(abs(xxc1-xxc2).gt.1.0e-10)then
cx      write(60,888)mg,lg,xxc1,xxc2,abs(xxc1-xxc2)
cx888 format(2i4,4(2x,e13.6))
cx     endif

c  Sublimation(+ive) from bare ice (sublmi): Drops ice-salt into ocean
c   making the ocean saltier.
C  White ice formation (whitei): This changes snow (zero salt) to
c   ice, and thus the appropriate amount of ice-salt is to be removed
c   from the ocean. This will make it less salty.
c  For ease of coding, combine these 2 components into sublmi() :
      sublmi(mg)=sublmi(mg)-whitei(mg)

 5000 continue

c---- check if excess ice removed : change Tg for all ice points at
c---- this latitude to account for energy (better than local change?).
      do ns=1,2
      If(flagt(ns))then
        ma=(ns-1)*lon
        sumt=0.0
        kice=0
        do 6000 mg=1+ma,lon+ma
        if(cice(mg))kice=kice+1
 6000   sumt=sumt+dsict(mg)
        sumt=sumt/max(1,kice)
c       print *,lg,ns,kice,sumt
        do 6010 mg=1+ma,lon+ma
 6010   if(cice(mg))tg(mg)=tg(mg)-sumt
      End if
      enddo
c----
      return
      end
