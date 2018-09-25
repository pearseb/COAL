c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Loop 225 replaced with modified version supplied by Leon Rotstayn. This
c removes a CFL problem that prevents the use of a 30-minute timestep.
c SJP 2006/06/02
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: surfa.f,v $
c Revision 1.51  2001/10/12 02:06:57  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.50  2001/02/28 04:36:38  rot032
c Further tidy ups from HBG
c
c Revision 1.49  2001/02/22 05:34:40  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.48  2001/02/12 05:39:50  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.47  2000/11/14 03:42:55  rot032
c Make sure kaos bugfix gets into HBG version!
c
c Revision 1.46  2000/11/14 03:11:37  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.45  2000/10/04 23:25:32  rot032
c Really fix the kaos bug this time!
c
c Revision 1.44  2000/08/22 02:38:28  rot032
c Hard code the call to fle for better vectorization.
c
c Revision 1.43  2000/08/21 07:39:25  rot032
c Little bugfix that was causing inconsistent results on kaos.
c
c Revision 1.42  2000/06/20 02:54:07  rot032
c Merge of HBG and LDR changes
c
c Revision 1.41  2000/06/20 02:08:33  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.40.1.1  2000/06/20 02:27:10  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.40  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.39  1998/12/10  00:55:38  ldr
c HBG changes to V5-1-21
c
c Revision 1.38  1998/05/26  05:29:36  ldr
c Final changes (inc. ACH salinity stuff) for 10k run.
c
c Revision 1.37  1997/12/23  01:40:46  ldr
c Add a couple of TASKLOCAL directives for NEC.
c
c Revision 1.35  1997/12/23  00:23:36  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.34  1997/12/17  23:22:48  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.33  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.32  1996/12/20  01:36:42  ldr
c Tidy-ups to EAK snow stuff.
c
c Revision 1.31  1996/12/19  23:40:50  ldr
c EAK's snow albedo changes merged and debugged.
c
c Revision 1.30  1996/10/24  01:03:16  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.29.1.1  1996/12/05  03:42:33  ldr
c Fix to snow albedo from EAK.
c
c     INPUT/OUTPUT
c     Input:   from common/hybrpr in HYBRPR.f
c                  pdpsk - (pressure/surface pressure)**cappa
c
c              from common/logimsl in LOGIMSL.f
c                  land - logical variable for surface type
c
c              from common/rshead in this subroutine
c                  header2 - second header on restart file
c
c              from common/surf in this subroutine
c                  aftfh  - reciprocal of air resistance for ground
c                  aftfhg - reciprocal of air resistance for bare ground
c                  isoilm - soil types
c
c              from common/timex in TIMEX.f
c                  dt    - time step in seconds
c                  mstep - timestep in minutes
c
c              from arguments
c                  cls  - factor for latent heat of fusion
c                  rcondx- precipitation in mms of water per step
c                  1  - hemisphere index       lg - latitude index
c                  pg   - surface pressure(mbs)  qtg - mixing ratio
c
c     Output:  from arguments
c                  rg     - net long wave heating at ground
c                  scalev - scaling evaporation
c                  tddd   - total wet evaporation
c
c     In/Out:  from common/surf1 in SURF1.f
c                  mc   - moisture depth on canopy
c
c     In:      from common/surf1 in SURF1.f
c                  ssdn   - snow density
c
c
c This routine is used for NSIB land surface scheme. The routine computes
c the soil and canopy temperatures and sensble and latent heat fluxes.
c The details of the scheme are described in the DAR Technical Paper No.23.
c
c INPUTS:
c sg    - incoming short-wave
c rgsav - incoming long-wave
c wb    - soil moisture at the previous time step (ms layers)
c tgg   - ground temperature at the previous time step (ms layers)
c tgf   - vegetation covered ground temperature at the previous time step
c ttg   - air temperature
c snowd - snow depth
cc
c OUTPUTS:
c tg    - combined soil-canopy temperature of the surface layer
c fg,eg - sensible and latent heat fluxes for the whole grid
c tgf   - vegetation covered ground temperature
c fgf,evapxf - sensible and latent heat fluxes for the canopy surface
c tgg   - ground temperature
c fgg,egg - sensible and latent heat fluxes for the bare ground surface
c tb2   - temperature of the second layer
c tb3   - temperature of the lowest layer
c totpev - potential evaporation in mm/timestep.
c scalev - scaling evaporation in mm/timestep.
c tddd - wet evaporation
c
      subroutine surfa(lg,sg,rg,
     &                 snowd,rcondx,pg,ttg,qtg,qlg,qfg,
     &                 cls,eg,fg,tg,tb2,tb3,
     &                 totpev,scalev,wg,wg2,tddd,he,mcmax)

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /SOILPR/ )
!$OMP THREADPRIVATE ( /SURFBC/ )
!$OMP THREADPRIVATE ( /VEGDAT/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      real snmin,cgsnow
      parameter (snmin=0.15)
      parameter (cgsnow=2105.)
      integer nalph
c     parameter (nalph=1) ! Mk3 model
      parameter (nalph=2) ! Mk3.1 model

C Argument list
      integer lg
      real sg(ln2)
      real rg(ln2)
      real snowd(ln2)
      real rcondx(ln2)
      real pg(ln2) 
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real cls(ln2)
      real eg(ln2)
      real fg(ln2)
      real tg(ln2)
      real tb2(ln2)
      real tb3(ln2)
      real totpev(ln2)
      real scalev(ln2)
      real wg(ln2)
      real wg2(ln2)
      real tddd(ln2)
      real he(ln2)
      real mcmax(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'
      include 'LOGIMSL.f'

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
      include 'FEWFLAGS.f'
      include 'MASIV3.f'
      include 'SURF1.f'
      include 'TIMEX.f'


      logical lakesg,lakeind
      common/glakes/lakesg(ln2,lat),lakeind(lat)

      real prhos,pcss,pcnsd,pbch,pssat,pswilt,phsbh,psfc,psucs
     & ,phyds,pcdr3,prlaim,prlais,pslveg,pscveg
      integer ip2bp3,ipbp2
      common/psoilvegpr/prhos(ln2,lat),pcss(ln2,lat),
     &pcnsd(ln2,lat),pbch(ln2,lat),pssat(ln2,lat),
     &pswilt(ln2,lat),phsbh(ln2,lat),ip2bp3(ln2,lat),
     &ipbp2(ln2,lat),psfc(ln2,lat),psucs(ln2,lat),
     &phyds(ln2,lat),pcdr3(ln2,lat),prlaim(ln2,lat),
     &prlais(ln2,lat),pslveg(ln2,lat),pscveg(ln2,lat)

      real resvr1,resrem,srunoff,sresvr1
      common/rivers/resvr1(0:lon+1,0:lat2),resrem(lon,lat2)
     & ,srunoff(lon,lat2),sresvr1(lon,lat2)

      character*50 header2
      common/rshead/header2

      real zs,zsh,zsdzs,zsfrac
      common/soilzs/zs(ms),zsh(ms+1),zsdzs(ms),zsfrac(ms)

      real aftfh,aftfhg,ocondx,opreci,pmc
      integer isoilm
      common/surf/aftfh(ln2,lat),aftfhg(ln2,lat),ocondx(ln2,lat)
     &          ,isoilm(ln2,lat),opreci(ln2,lat)
     &          ,pmc(ln2,lat)

C Local work arrays and variables
      real airr(ln2)
      real bbb(ln2)
      real cc(ln2)
c     real c1(ln2)
      real degdt(ln2),dfgdt(ln2),dirad(ln2),dgdtg(ln2)
      real esatta(ln2)
     & ,evapfb1(ln2),evapfb2(ln2)
     & ,evapfb3(ln2),evapfb4(ln2),evapfb5(ln2),evapxw(ln2)
     & ,eggtot(ln2)
c     real evap(ln2)
      real f124(ln2),fgf(ln2),fgg(ln2),ftgsoill(ln2)
      real ga(ln2)
c     real hstm(ln2),hstp(ln2)
      real omc(ln2),otggsn(ln2),otgsoil(ln2)
      real pevap(ln2)
      real qt1(ln2)
      real rdg(ln2),rgg(ln2),rlai(ln2),extin(ln2)
c    & ,residg(ln2)
      real theta(ln2),trho(ln2)
c     real tpevap(ln2)
      real wba(ln2),wbfice(ln2,ms),wblf(ln2,ms)
      real zsmod(ln2,ms) ! =zs (land) , or zs modified for lakes
      logical tester1(ln2),tester2(ln2),flag1,flag2

      integer icount
      integer lg60s
      integer k
      integer kk
      integer mg

      real bbbveg
      real beta
      real betetrdt
      real calgamm0
      real ccc
      real conw_fh
      real ccoef
      real deg
      real deltat
      real den
      real devf
      real dirad1
      real dqg
      real dqsttg
      real easss
      real eggx
      real eg1
      real eg2
      real eg2ep
      real etr
      real evapfb
      real ew
      real ewi
      real f
      real f1
      real f2
      real f4
      real f3
      real ftgsoillx
      real hlrvap
      real pm
      real prz
      real qsatgf
c      real qstta
      real qsttg
      real qtgair
      real qtgairep
      real qtgnet
      real res
      real residf
      real residp
      real rsi
      real sfl
      real sicefreeze
      real sicemelt
      real snadd
      real snowdd
      real srlai
      real sstar
      real stgf
      real tempmc
      real tgfnew
      real tggd
      real tgsoill
      real tgss
      real tsigmfx
      real tr1
      real tr11
      real wbav
      real wbicep
      real wbliqp
      real wbtest
      real wetfct

C Local data, functions etc
      real root(5)
c     data root/0.20,0.45,0.20,0.10,0.05/ 
      data root/0.05,0.10,0.35,0.40,0.10/ 
      save root
      real csice,cswat,rhowat,tfrz
      data csice /2.100e3/, cswat /4.218e3/, rhowat /1000.0/
      data tfrz/273.15/       ! (BP aug2010)
      include 'ESTABL.f'

C Start code : ----------------------------------------------------------

      hlrvap= hl/rvap
c**
c**   Check for any land points
c** 
c     call checkl(land,ln2,flag1)
c     if(.not.flag1)return

c**
c**   update land points.
c**
      if(header2(50:50).ne.'M')then
c.... Initialize Mk3 nsib data from Mk2 nsib data
      do 3800 mg=1,ln2
      if( .not.land(mg)) go to 3800
      tggsl(mg,2,lg)=tb2(mg)
      tggsl(mg,3,lg)=tb2(mg)
      isflag(mg,lg)=0
      snadd=0.0
      if(snowd(mg).gt.100.)snadd=100.0
      ssdnn(mg,lg) =140.+snadd
      ssdn3(mg,1,lg)=140.+snadd
      ssdn3(mg,2,lg)=140.+snadd
      ssdn3(mg,3,lg)=140.+snadd
      osnowd(mg,lg)=snowd(mg)
      wb(mg,1,lg)=wg(mg)
      wb(mg,2,lg)=wg2(mg)
      wb(mg,3,lg)=wg2(mg)
      gflux(mg,lg)=0.0
      sgflux(mg,lg)=0.0
      do kk=4,ms
      tggsl(mg,kk,lg)=tb3(mg)
      wb(mg,kk,lg)=wg2(mg)
      enddo
      do kk=1,ms
      if(wb(mg,kk,lg).lt.pswilt(mg,lg)) wb(mg,kk,lg)=
     &                                          pswilt(mg,lg)
      wbice(mg,kk,lg)=0.
      if(tggsl(mg,kk,lg).lt. tfrz-5.0) 
     &     wbice(mg,kk,lg)=0.85*wb(mg,kk,lg)
      enddo
 3800 continue
      lg60S=int(0.31*lat) ! For latitude row near 60S
      if(lg.le.lg60S)then
      do 3801 mg=1,lon
       if(land(mg+lon))tggsl(mg+lon,3,lg)=tb3(mg+lon)
 3801 continue
      endif
      endif

      do 3840 mg=1,ln2
        tddd(mg) = 0.0
        totpev(mg)=0.0
        scalev(mg)=0.0
        if(qcloud)then
          qt1(mg)=qtg(mg,1)+qfg(mg,1)+qlg(mg,1)
        else
          qt1(mg)=qtg(mg,1)
        endif
 3840 continue
    
      do 3841 mg=1,ln2
      if(land(mg))then
        tgsoill=0.5*(0.3333*tggsl(mg,2,lg)+0.6667*tggsl(mg,3,lg)
     & +0.95*tggsl(mg,4,lg) +  0.05*tggsl(mg,5,lg))
        ftgsoillx=max(0.0,1.0-0.0016*(298.0-tgsoill)**2)
        if( tgsoill .ge. 298.0 ) ftgsoillx=1.0
        tsigmfx=max(0.0,sigmf(mg)-pscveg(mg,lg)*(1.-ftgsoillx))
        ftgsoill(mg)=ftgsoillx
        ccc=min(10.,snowd(mg)/max(2.5*100.*z0m(mg),1.))
        pp=exp(ccc)
        pm=exp(-ccc)
c                                                tanh=(pp-pm)/(pp+pm)
        tsigmf(mg)=(1.-(pp-pm)/(pp+pm))*tsigmfx
        condxg(mg) = 0.0
        trho(mg) = pg(mg)/ttg(mg,1)*100.0/rdry
        rlai(mg)= max(0.1,prlaim(mg,lg)-pslveg(mg,lg)*(1.0-
     &  ftgsoill(mg)))
        extin(mg)=exp(-0.6*max(2.7,rlai(mg)))
      endif
 3841 continue

      do 3842 mg=1,ln2
      if (land(mg))then
       do k=1,ms
         tgsoil(mg,k)=tggsl(mg,k,lg)
         wblf(mg,k)=max(0.001,(wb(mg,k,lg)-wbice(mg,k,lg))
     &             /pssat(mg,lg))
         wbfice(mg,k)=max(0.0,wbice(mg,k,lg)/pssat(mg,lg))
       enddo
      endif
 3842 continue

      do mg=1,ln2
        pmcmax(mg)=(1.0-4.0e-04*min(he(mg),2000.))*20. !max puddle 20mm?
      enddo
      do k=1,ms
      do mg=1,ln2
        zsmod(mg,k)=zs(k) ! soil depth (modified next if lakes)
      enddo
      enddo

c Set lakes properties. Fill up puddle water from reservoirs
      if(newriver.and.lakeind(lg))call surfa_riv(1,lg,zsmod)

      do 3843 mg=1,ln2
      if(land(mg))then
c
c     bare ground
c
        tgss=tggsl(mg,1,lg)
        if( isflag(mg,lg).eq.1 ) tgss=tggsn(mg,1,lg)
        esatta(mg)=establ(ttg(mg,1))
        qsttg=qsat(100*pg(mg), tgss)
        hlrvap=hl/rvap
        dqsttg=qsttg*hlrvap/tgss**2
c        qstta=qsat(100*pg(mg), ttg(mg,1))

c        sss=abs(qsttg-qstta)/max(abs(tgss-ttg(mg,1)),
c     &      0.000001)
c        c1(mg)=sss/(sss+cp/hl)
c        c2=(cp/hl)/(sss+cp/hl)
        bbb(mg) =  stefbo*tgss**4
c        rgg(mg) = rgsav(mg,lg) + bbb(mg)
        bbbveg =  stefbo*tgf(mg,lg)**4
        rgg(mg) = rgsav(mg,lg)*(1.- tsigmf(mg)) + bbb(mg)
     &             - tsigmf(mg)*bbbveg
        dirad(mg)=4.0*bbb(mg)  /tgss
        theta(mg) = ttg(mg,1)/pdpsk(mg,1)

c        ssgfl=gflux(mg,lg)
c        if( isflag(mg,lg).eq.1 ) ssgfl=sgflux(mg,lg)
c        term1=c1(mg)*(sg(mg)-rgg(mg)-ssgfl)
c        vec1=c2*trho(mg)*hl*(qstta-qt1(mg)+sss*(ttg(mg,1)-theta(mg)))
c        term2=aftfhg(mg,lg)*vec1
c        stc=wb(mg,1,lg)/pssat(mg,lg)
c        pevap(mg)=aftfhg(mg,lg)*trho(mg)*hl*(qsttg-qt1(mg))
c        tpevap(mg)=term1+term2
        conw_fh = aftfhg(mg,lg)*trho(mg)*hl
        wba(mg)=max(0.0,wb(mg,1,lg)-wbice(mg,1,lg))
c        wetfct=fle(mg,lg,wba(mg),1.,tsigmf(mg)) !Use next line instead
        wetfct=max(0.,min(1.0,         !Hard coding the call to fle vectorizes...
     &       (wba(mg)-tsigmf(mg)*pswilt(mg,lg))/
     &       (psfc(mg,lg)-tsigmf(mg)*pswilt(mg,lg))))
        if(lakesg(mg,lg))wetfct=1.0 ! Max evap for lakes (== bare ground)

        if(nalph.eq.1) then ! alpha scheme
          pevap(mg)=conw_fh*(qsttg-qt1(mg))
          eggx=wetfct*pevap(mg)
          degdt(mg)=wetfct*conw_fh*dqsttg
        else
c                             following is beta scheme
         qtgnet=qsttg*wetfct -qt1(mg)
         qtgair=qsttg*wetfct-max(qtgnet,.1*qtgnet)
         qtgairep=qsttg -max(qtgnet,.1*qtgnet)
         eg2=-conw_fh*qtgair
         eg2ep = -conw_fh*qtgairep
         eg1=conw_fh*qsttg
c                                 evaporation from the bare ground
         eggx=eg1*wetfct + eg2
         pevap(mg)=eg1 + eg2ep
         deg=wetfct*conw_fh*dqsttg
c        following reduces degdt by factor of 10 for dew
         degdt(mg)=.55*deg+sign(.45*deg,qtgnet)
       endif  ! (nalph.eq.1)

c        eggx=fle(mg,lg,wba(mg),tpevap(mg))
    
        eggx=min(eggx,wb(mg,1,lg)*1000.*zsmod(mg,1)*hl/dt)
        if(qsttg .lt. qt1(mg)) eggx = pevap(mg)
        if(snowd(mg).gt.1.) eggx=min(snowd(mg)/dt*hl,pevap(mg))
        if((pevap(mg).lt.0.0).and.(qsttg.ge.qt1(mg))) eggx=0.0
c
        beta=min(pmc(mg,lg)/pmcmax(mg),0.1)
        ewwp(mg)=min(pmc(mg,lg)/dt*hl,max(beta*pevap(mg),0.)) ! not for dew
        if(qsttg .lt. qt1(mg))beta=0.0
        egg(mg)=eggx*(1-beta)
        eggtot(mg)=egg(mg)+ewwp(mg)

        dfgdt(mg)=aftfhg(mg,lg)*trho(mg)*cp
        fgg(mg)=dfgdt(mg)*(tgss-theta(mg))
        
c        ga(mg)=sg(mg)-rgg(mg)-fgg(mg)-eggx*cls(mg)
c        ga(mg)=(sg(mg)-rgg(mg)-fgg(mg)-eggx*cls(mg))*(1.-tsigmf(mg))
        ga(mg)=(sg(mg)-fgg(mg)-eggtot(mg)*cls(mg))*
     &        (1.-tsigmf(mg)) - rgg(mg) + tsigmf(mg)*extin(mg)*sg(mg) 

c        dgdtg(mg)=-(dirad(mg)+dfgdt(mg)+cls(mg)*degdt(mg))*
        dgdtg(mg)=-(dfgdt(mg)+cls(mg)*degdt(mg))*
     &             (1.-tsigmf(mg))-dirad(mg)       

        otgsoil(mg)=tgsoil(mg,1)

      endif
 3843 continue

      do 3702 mg=1,ln2
 3702 tester1(mg)=land(mg)

      do 3844 mg=1,ln2
      if(tester1(mg).and.(snowd(mg).le.0.)) then
       isflag(mg,lg)=0
       ssdn3(mg,1,lg)=140.
       ssdnn(mg,lg)=140.
       tester1(mg)=.false.
      endif
 3844 continue

      call checkl(tester1,ln2,flag1)
      if(flag1)then

      do 3846 mg=1,ln2
      scond1s(mg,1)=0.
      if(tester1(mg).and.(snowd(mg).gt.0.))then

       if(snowd(mg)/ssdnn(mg,lg).lt. snmin) then
       ccoef=0.0
       if(ssdn3(mg,1,lg) .ge.150.0) ccoef=4.6e-2
       tggd=min(tfrz,tggsl(mg,1,lg))
       if(isflag(mg,lg).eq.1)tggd=min(tfrz,tggsn(mg,1,lg))
       if(isflag(mg,lg).eq.1) ssdn3(mg,1,lg)=ssdnn(mg,lg) 
       ssdn3(mg,1,lg) = max(140.,ssdn3(mg,1,lg)+
     &   dt*ssdn3(mg,1,lg)*2.8e-6* exp(-3.e-2*
     &    (273.1-tggd)-ccoef*(ssdn3(mg,1,lg)-150.0)))
       tr1 =  max(0.,snowd(mg)-osnowd(mg,lg))
       tr11 = tr1/snowd(mg)
       ssdn3(mg,1,lg) = tr11*140.+ (1.-tr11)*ssdn3(mg,1,lg) 
       scond1s(mg,1) = max(0.2,min(2.576e-6*ssdn3(mg,1,lg)*
     &                 ssdn3(mg,1,lg)+0.074,0.8))
c       cap1s(mg,1) = ssdn3(mg,1,lg) * (92.88+7.364*tggd)
       cap1s(mg,1)=ssdn3(mg,1,lg)*2105.
       ssdnn(mg,lg)=ssdn3(mg,1,lg)
       isflag(mg,lg)=0
       sdepth(mg,1)=snowd(mg)/ssdn3(mg,1,lg)
       tester1(mg)=.false.
       endif

      endif
 3846 continue

      call checkl(tester1,ln2,flag2)
      if(flag2)call snowpr(lg,otggsn,dt,snowd,tester1)

      endif ! if flag1

      call stemp (lg,wblf,wbfice,ga,dgdtg,dt,snowd,zsmod)

      do k=1,ms
      do 3851 mg=1,ln2
      if(land(mg)) tggsl(mg,k,lg)=tgsoil(mg,k)
 3851 continue
      enddo

      do 3908 k=1,ms

      do 3902 mg=1,ln2
 3902 tester1(mg)=land(mg)

      do 3904 mg=1,ln2
      if(tester1(mg))then

       wbtest=0.9*wb(mg,k,lg)-wbice(mg,k,lg)
       if((tggsl(mg,k,lg).lt.tfrz).and.(wbtest.gt.0.001)) then

         sfl=(tfrz-tggsl(mg,k,lg))*gamm(mg,k)
         sicefreeze=min(max(0.,wbtest)*zsmod(mg,k)*1000. ,sfl/hlf)
         wbicep=min(wbice(mg,k,lg)+sicefreeze/
     &               (zsmod(mg,k)*1000.0),0.9*wb(mg,k,lg))
         wbicep=max(wbicep,0.)
         wbice(mg,k,lg)=wbicep
     
         wbliqp=wb(mg,k,lg)-wbicep
         snowdd=snowd(mg)
calg     gamm(mg,k)=calgamm(mg,lg,k,wbliqp,wbicep,snowdd,zsmod(mg,k))
      calgamm0 = max( (1.-pssat(mg,lg))*pcss(mg,lg)*
     & prhos(mg,lg) + wbliqp*cswat*rhowat + wbicep*csice*rhowat*0.9,
     & pcss(mg,lg)*prhos(mg,lg) ) * zsmod(mg,k)
      if(k.eq.1.and.isflag(mg,lg).eq.0) calgamm0 = calgamm0+
     &   cgsnow*snowdd
         gamm(mg,k)= calgamm0
         tggsl(mg,k,lg)=tggsl(mg,k,lg)+sicefreeze*hlf/calgamm0
         if((0.9*wb(mg,k,lg)-wbicep).gt.0.001) 
     &                tggsl(mg,k,lg)= max(tfrz,tggsl(mg,k,lg))

        tester1(mg)=.false.
       endif
      endif
 3904 continue

      call checkl(tester1,ln2,flag1)
      if(flag1)then
      
      do 3906 mg=1,ln2
      if(tester1(mg))then
       if( tggsl(mg,k,lg).gt.tfrz.and.wbice(mg,k,lg)
     &        .gt.0.0) then

         sfl=(tggsl(mg,k,lg)-tfrz)*gamm(mg,k)
         sicemelt=min(wbice(mg,k,lg)*zsmod(mg,k)*1000.,sfl/hlf)
         wbicep=wbice(mg,k,lg)-sicemelt/(zsmod(mg,k)*1000.)
         if(wbicep.lt.0.001) wbicep=0.
         wbice(mg,k,lg)=wbicep

         wbliqp=wb(mg,k,lg)-wbicep
         snowdd=snowd(mg)
calg     gamm(mg,k)=calgamm(mg,lg,k,wbliqp,wbicep,snowdd,zsmod(mg,k))
      calgamm0 = max( (1.-pssat(mg,lg))*pcss(mg,lg)*
     & prhos(mg,lg) + wbliqp*cswat*rhowat + wbicep*csice*rhowat*0.9,
     & pcss(mg,lg)*prhos(mg,lg) ) * zsmod(mg,k)
      if(k.eq.1.and.isflag(mg,lg).eq.0) calgamm0 = calgamm0+
     &   cgsnow*snowdd
         gamm(mg,k)= calgamm0
         tggsl(mg,k,lg)=tggsl(mg,k,lg)-sicemelt*hlf/gamm(mg,k)
         if(wbicep.gt.0.) tggsl(mg,k,lg)=
     &                               min(tfrz,tggsl(mg,k,lg))

       endif
      endif
 3906 continue

      endif ! if flag1

 3908 continue

      if(newriver.and.lakeind(lg))call surfa_riv(2,lg,zsmod)
c     if(lg.eq.lat)stop

      do 385 mg=1,ln2
      if(land(mg)) then
        if(isflag(mg,lg).eq.0) then
         deltat=tggsl(mg,1,lg)-otgsoil(mg)
         egg(mg)=egg(mg)+deltat*degdt(mg)
         egg(mg)=min(egg(mg),wb(mg,1,lg)*1000.*zsmod(mg,1)*hl/dt)
        else
         deltat=tggsn(mg,1,lg)-otggsn(mg)
         egg(mg)=egg(mg)+deltat*degdt(mg)
chbg     egg(mg)=min(egg(mg),snowd(mg)/dt*hl) ! not worth checking
        endif
        eggtot(mg)=egg(mg)+ewwp(mg)
        fgg(mg)=fgg(mg)+deltat*dfgdt(mg)
        rgg(mg)=rgg(mg)+deltat*dirad(mg)
      endif
385   continue

c     vegetation

      do 3854 mg=1,ln2
      if(land(mg).and.(tsigmf(mg).le.0.01)) then
        fgf(mg)  = fgg(mg)
        rdg(mg)=rgg(mg)
        tgf(mg,lg) = tggsl(mg,1,lg)
        evapxf(mg) = eggtot(mg)
        evapxw(mg)=0.0
      endif
 3854 continue

      do 392 mg=1,ln2
        tester1(mg)=land(mg).and.(tsigmf(mg).gt.0.01)
        tester2(mg)=tester1(mg)
  392 continue

      call checkl(tester1,ln2,flag1)
      if(flag1)then

      do 391 mg=1,ln2
      if(tester1(mg)) then 
c      iveg=vegt(mg)
c      rlai(mg)= max(0.1,prlaim(mg,lg)-pslveg(mg,lg)*(1.0-
c     &ftgsoill(mg)))
      srlai=rlai(mg)+prlais(mg,lg)
      sstar = 30.0
      if( z0m(mg) .lt. 0.5 ) sstar = 150.0
      f= 1.1*sg(mg)/(rlai(mg)*sstar)
c     rsi = rsmin(mg) * rlai(mg)
      rsi = rsmin(mg)
      f1= (1+f)/(f+rsi/5000.0)
      den=psfc(mg,lg)-pswilt(mg,lg)
      wbav=max(0.,root(1)*(wb(mg,1,lg)-pswilt(mg,lg))/den)+
     &     max(0.,root(2)*(wb(mg,2,lg)-pswilt(mg,lg))/den)+
     &     max(0.,root(3)*(wb(mg,3,lg)-pswilt(mg,lg))/den)+
     &     max(0.,root(4)*(wb(mg,4,lg)-pswilt(mg,lg))/den)+
     &     max(0.,root(5)*(wb(mg,5,lg)-pswilt(mg,lg))/den)
      f2=max(1.0,0.5/ max( wbav,0.0000001))
      f4=max(1.-0.0016*(298.0-ttg(mg,1))**2,0.20)
      f124(mg)=f1*f2/f4
      airr(mg) = 1.0/aftfh(mg,lg)
      cc(mg) = rcondx(mg) ! mms/timestep
chbg  if(rcondx(mg)*(1440.0/mstep).gt.4.0) cc(mg) = 0.08333
      if(rcondx(mg)*(1440.0/mstep).gt.4.0) cc(mg) = 4.0*mstep/1440.0
c     mcmax(mg) = max(srlai,0.5)*0.1
      mcmax(mg) = max(srlai,0.5)*0.2
      omc(mg) = mc(mg,lg)
      endif
 391  continue

      icount=0
  100 icount=icount+1

      do 3921 mg=1,ln2
      if(tester1(mg)) then 

      tempmc = omc(mg)
      qsatgf=qsat(100.0*pg(mg), tgf(mg,lg))
      easss=qt1(mg)*100.0*pg(mg)/0.622
      f3=max(1.-0.00025*(esatta(mg)-easss),0.5)
c
c     res= airr(mg) + max((rsmin(mg)*f124(mg)/f3),20.)
      res= airr(mg) + max(((rsmin(mg)/rlai(mg))*f124(mg)/f3),20.)

      ew = trho(mg)/airr(mg) *(qsatgf-qt1(mg)) 
      ewi=ew
      if( qsatgf .ge. qt1(mg) ) 
     &  ewi  = min(tempmc/dt,tempmc/mcmax(mg) * ew)
      tempmc = omc(mg)+cc(mg)-ewi*dt

      condxg(mg)=max( rcondx(mg)-cc(mg)+
     &   max(0.,tempmc-mcmax(mg)),0.0)
      tempmc = min( max(tempmc,0.0), mcmax(mg))
      beta =  min(tempmc/mcmax(mg),1.0)
      if(tempmc.lt.1.0e-10)tempmc=0.0
      mc(mg,lg)=tempmc
      if( qsatgf .lt. qt1(mg) ) beta = 1 
      etr=max(0.,trho(mg)/res*(qsatgf-qt1(mg)))
     
      prz = trho(mg)/airr(mg)*cp
      fgf(mg) = prz*(tgf(mg,lg)-theta(mg)) 
      betetrdt =(1.-beta)*etr*dt 
      evapfb1(mg)=min(betetrdt*root(1),max(0.,
     &            (wb(mg,1,lg)-pswilt(mg,lg))*zsmod(mg,1)*1000.))
      evapfb2(mg)=min(betetrdt*root(2),max(0.,
     &            (wb(mg,2,lg)-pswilt(mg,lg))*zsmod(mg,2)*1000.))
      evapfb3(mg)=min(betetrdt*root(3),max(0.,
     &            (wb(mg,3,lg)-pswilt(mg,lg))*zsmod(mg,3)*1000.))
      evapfb4(mg)=min(betetrdt*root(4),max(0.,
     &            (wb(mg,4,lg)-pswilt(mg,lg))*zsmod(mg,4)*1000.))
      evapfb5(mg)=min(betetrdt*root(5),max(0.,
     &            (wb(mg,5,lg)-pswilt(mg,lg))*zsmod(mg,5)*1000.))

      evapfb=evapfb1(mg)+evapfb2(mg)+evapfb3(mg)+evapfb4(mg)+
     &           evapfb5(mg)

      evapxf(mg) = (evapfb/dt + ewi)*hl
      evapxw(mg) =  ewi*hl
      tddd(mg) = max(ewi*dt,0.0)*tsigmf(mg)   

      bbbveg     =  stefbo*tgf(mg,lg)**4
c      rdg(mg) = rgsav(mg,lg) + bbbveg
      rdg(mg) = rgsav(mg,lg) + 2.*bbbveg - bbb(mg)
      residf   = sg(mg)- rdg(mg)-fgf(mg)-evapxf(mg)*cls(mg)
     &           -extin(mg)*sg(mg)

      if( abs(residf    ) .lt. 3.0 ) tester1(mg)=.false.
      if(tester1(mg))then
       dirad1 = 4.0*bbbveg    /tgf(mg,lg)
c                                               calculate de/dtg
       dqg=qsatgf*hlrvap/tgf(mg,lg)**2
       devf= hl*(((1.0-beta)*trho(mg)*1.0/res
     +    *dqg)+beta*trho(mg)/airr(mg)*dqg )
       residp = -(2*dirad1 + devf + prz)
       stgf = tgf(mg,lg)
       tgfnew = stgf + residf    /(-residp)

       tgfnew=min(tgfnew,stgf + 1.4)
       tgfnew=max(tgfnew,stgf - 1.4)
       tgf(mg,lg) = tgfnew
      endif

      endif
 3921 continue

c Check if any land points still need iterating
      call checkl(tester1,ln2,flag2)

      if(flag2.and.(icount.lt.5)) go to 100

      do 3923 mg=1,ln2
      if(tester2(mg)) then 

       wb(mg,1,lg)=wb(mg,1,lg)-evapfb1(mg)
     &    /(zsmod(mg,1)*1000.)*tsigmf(mg)
       wb(mg,2,lg)=wb(mg,2,lg)-evapfb2(mg)
     &    /(zsmod(mg,2)*1000.)*tsigmf(mg)
       wb(mg,3,lg)=wb(mg,3,lg)-evapfb3(mg)
     &    /(zsmod(mg,3)*1000.)*tsigmf(mg)
       wb(mg,4,lg)=wb(mg,4,lg)-evapfb4(mg)
     &    /(zsmod(mg,4)*1000.)*tsigmf(mg)
       wb(mg,5,lg)=wb(mg,5,lg)-evapfb5(mg)
     &    /(zsmod(mg,5)*1000.)*tsigmf(mg)

       endif
 3923 continue

      endif ! if flag1

      do 3940 mg=1,ln2
c       hstm(mg)=0.
c       hstp(mg)=0.
        tester1(mg)=land(mg).and.(isflag(mg,lg).eq.0)
        tester2(mg)=land(mg).and.(isflag(mg,lg).ne.0)
 3940 continue

       do 3942 mg=1,ln2
       if(tester1(mg)) then
         if( tsigmf(mg) .le. 0.01 ) then
           tg(mg) = tggsl(mg,1,lg)
           tgf(mg,lg) = tggsl(mg,1,lg)
         else
           tg(mg)=tsigmf(mg)*tgf(mg,lg)+(1.0-tsigmf(mg))*
     &            tggsl(mg,1,lg)
         endif
       endif
 3942 continue
c      do 3943 mg=1,ln2
c      if(tester1(mg)) then
c        do k=1,6
c         hs1tp(mg,k)=gamm(mg,k)*tgsoil(mg,k)
c        enddo 
c        do k=1,6
c         hstm(mg)=hstm(mg)+hs1tm(mg,k)
c         hstp(mg)=hstp(mg)+hs1tp(mg,k)
c        enddo
c      endif
c3943 continue

       do 3944 mg=1,ln2
       if(tester2(mg)) then
         if( tsigmf(mg) .le. 0.01 ) then
           tg(mg) = tggsn(mg,1,lg)
           tgf(mg,lg) = tggsn(mg,1,lg)
         else
           tg(mg)=tsigmf(mg)*tgf(mg,lg)+(1.0-tsigmf(mg))*
     &            tggsn(mg,1,lg)
         endif
       endif
 3944 continue
c      do 3945 mg=1,ln2
c      if(tester2(mg)) then
c        do k=1,3
c         hs1tp(mg,k)=sgamm(mg,k)*tggsn(mg,k,lg)
c        enddo
c        do k=1,6
c         hs1tp(mg,k+3)=gamm(mg,k)*tgsoil(mg,k)
c        enddo 
c        do k=1,9
c         hstm(mg)=hstm(mg)+hs1tm(mg,k)
c         hstp(mg)=hstp(mg)+hs1tp(mg,k)
c        enddo
c      endif
c3945 continue

      do 3970 mg=1,ln2
      if(land(mg))then
c---- compute change & ensure step conservation
c     dxxx=(hstp(mg)-hstm(mg))/dt
c    &    -(sg(mg)-rgg(mg)-cls(mg)*egg(mg)-fgg(mg))
c     arg=abs(rgg(mg))
c     aeg=abs(egg(mg))
c     afg=abs(fgg(mg))
c     sumab=arg+cls(mg)*aeg+afg+0.1e-10
c      eg(mg) = tsigmf(mg)*evapxf(mg) + (1.0 - tsigmf(mg))*egg(mg)
      eg(mg) = tsigmf(mg)*evapxf(mg) + (1.0 - tsigmf(mg))*eggtot(mg)
      if(snowd(mg).gt.1.0) eg(mg)=tsigmf(mg)*evapxf(mg) + eggtot(mg)
      fg(mg) = tsigmf(mg)*fgf(mg)+(1.0 - tsigmf(mg))*fgg(mg)
      rg(mg) = tsigmf(mg)*rdg(mg)+rgg(mg)
c      rg(mg) = tsigmf(mg)*rdg(mg)+(1.0 - tsigmf(mg))*rgg(mg)
      totpev(mg)=(tsigmf(mg)*evapxw(mg)+(1.-tsigmf(mg))*pevap(mg))/
     &    hl*dt
      scalev(mg)=totpev(mg)
c
c     residg(mg)=sg(mg)-rgg(mg)-fgg(mg)-egg(mg)-gflux(mg,lg)
c     if(isflag(mg,lg).eq.1) 
c    &residg(mg)=sg(mg)-rgg(mg)-fgg(mg)-egg(mg)-sgflux(mg,lg)
c
c**   evap=evaporation/timestepin kgm/m**2(=mms)
c     evap(mg)=dt*eg(mg)/hl
c     temporary fix
      tb2(mg)=tggsl(mg,2,lg)
      tb3(mg)=tggsl(mg,ms,lg)
      endif
 3970 continue

      return
      end
c
      function fle(mg,lg,g,ep,frac)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      real fle
      integer mg
      integer lg
      real g
      real ep
      real frac

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

C Local work arrays and variables

C Local data, functions etc
   
C Start code : ----------------------------------------------------------

       fle=min(1.0,(g-frac*pswilt(mg,lg))/(psfc(mg,lg)-
     &                                 frac*pswilt(mg,lg)))
       fle = max( fle,0.0)*ep
      return
      end
c
      subroutine stemp(lg,owb,owib,tbhfl,dgdtg,dt,snowd,zsmod)

      implicit none

!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /SOILPR/ )
!$OMP THREADPRIVATE ( /SURFBC/ )

C Global parameters
      include 'PARAMS.f'
      real cgsnow
      parameter (cgsnow=2105.)

C Argument list
      integer lg
      real owb(ln2,ms)
      real owib(ln2,ms)
      real tbhfl(ln2)
      real dgdtg(ln2)
      real dt
      real snowd(ln2)
      real zsmod(ln2,ms)

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
      include 'SURF1.f'

      logical lakesg,lakeind
      common/glakes/lakesg(ln2,lat),lakeind(lat)

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
      real cnsw(ln2,ms)
      real ccnsw(ln2,ms+3),tzs(ln2,ms+3),ttgsoil(ln2,ms+3)
      real  c1(ln2,ms+3), c2(ln2,ms+3), c3(ln2,ms+3), rhs(ln2,ms+3)
      real wc1(ln2,ms+3),wc2(ln2,ms+3),wc3(ln2,ms+3),wrhs(ln2,ms+3)
      logical tester1(ln2),tester2(ln2)
      logical flag1,flag2
      INTEGER INDEX1(ln2) ! INDEX FOR COMPRESS AND EXPAND

      integer i
      integer k
      integer mg
      integer ms0
      integer ms3
      integer npts

      real calgamm0
      real ccf
      real cfr1
      real dtg
      real ei
      real ew
      real eww
      real snowdd
      real wbliqp
      real wbicep
      real xx
      real xy


C Local data, functions etc
      real csice,cswat,rhowat
      data csice /2.100e3/, cswat /4.218e3/, rhowat /1000.0/

C Start code : ----------------------------------------------------------

c        calculates temperatures of the soil 
c
c        tgsoil - new soil/ice temperature
c        tbhfl - heat flux from the atmosphere
c        cnsw - soil conductivity 
c        dt  - time step 

      do 200 mg=1,ln2
      tester1(mg)=.false.
      tester2(mg)=.false.
      if(land(mg))then
        tester1(mg)=isflag(mg,lg).eq. 0
        tester2(mg)=isflag(mg,lg).ne. 0
      endif
  200 continue

      call checkl(tester1,ln2,flag1)
      call checkl(tester2,ln2,flag2)

      if(flag1)then

       do k=1,ms
      do 210 mg=1,ln2
      if(tester1(mg))then
        ew=max(0.05,owb(mg,k)*pssat(mg,lg))
        eww=min(ew,pssat(mg,lg)/2.)
        ccf=max(1.,sqrt(min(2.0,pssat(mg,lg)/(2.*eww))))
        ei=owib(mg,k)*pssat(mg,lg)
        cnsw(mg,k)= min(pcnsd(mg,lg)*(60.)**ew*(250.)**ei,2.2)*ccf
        if(pcnsd(mg,lg).gt.2.51) cnsw(mg,k)=2.51
          if(lakesg(mg,lg))cnsw(mg,k)=6.0*2.51
        tzs(mg,k)=zsmod(mg,k)
      endif
  210 continue
       enddo

      do 211 mg=1,ln2
      if(tester1(mg))then
       xx=snowd(mg)/ssdnn(mg,lg)
       xy=zsmod(mg,1)/(zsmod(mg,1)+xx)
       if(snowd(mg).eq.0.)xy=1. !Needed for Crays
       cnsw(mg,1)=cnsw(mg,1)*xy + (1.-xy)*scond1s(mg,1)
       tzs(mg,1)=zsmod(mg,1)+snowd(mg)/ssdnn(mg,lg)
       coef(mg,1)=0.
       coef(mg,ms+1)=0.
      endif
  211 continue

       do k=2,ms
      do 212 mg=1,ln2
      if(tester1(mg))then
         coef(mg,k)=1./
     &      (.5*(tzs(mg,k-1)/cnsw(mg,k-1)+tzs(mg,k)/cnsw(mg,k)))
      endif
  212 continue
       enddo
c
c For subsequent use of routine trim3x : set values 
c
      do 213 mg=1,ln2
       do k=1,ms
      if(tester1(mg))then
          wbliqp=owb(mg,k)*pssat(mg,lg)
          wbicep=owib(mg,k)*pssat(mg,lg)
          snowdd=snowd(mg)
calg      gamm(mg,k)=calgamm(mg,lg,k,wbliqp,wbicep,snowdd,zsmod(mg,k))
      calgamm0 = max( (1.-pssat(mg,lg))*pcss(mg,lg)*
     & prhos(mg,lg) + wbliqp*cswat*rhowat + wbicep*csice*rhowat*0.9,
     & pcss(mg,lg)*prhos(mg,lg) ) * zsmod(mg,k)
      if(k.eq.1.and.isflag(mg,lg).eq.0) calgamm0 = calgamm0+
     &   cgsnow*snowdd
c         if(lakesg(mg,lg))then
c          write(6,6666)mg,lg,k,pssat(mg,lg),wbliqp,wbicep
c    &    ,wbliqp*cswat*rhowat,wbicep*csice*rhowat*0.9
c    &    ,pcss(mg,lg)*prhos(mg,lg),zsmod(mg,k)
c6666 format(i3,1x,i2,1x,i1,3(2x,f6.3),3(2x,e10.3),2x,f7.3)
c         endif
          gamm(mg,k)= calgamm0
          hs1tm(mg,k)=gamm(mg,k)*tgsoil(mg,k)
          dtg=dt/gamm(mg,k)
          rhs(mg,k) = tgsoil(mg,k)
          c1(mg,k)= -dtg*coef(mg,k)
          c3(mg,k) = -dtg*coef(mg,k+1)
          c2(mg,k)= 1.-c1(mg,k)-c3(mg,k)
      endif
       enddo
  213 continue
c     if(lg.eq.lat)stop

      do 214 mg=1,ln2
      if(tester1(mg))then
c        rhs(mg,1) = rhs(mg,1)+tbhfl(mg)*dt/gamm(mg,1)
        rhs(mg,1) = rhs(mg,1)+dt/gamm(mg,1)*(tbhfl(mg)-tgsoil(mg,1)*
     &               dgdtg(mg))
        c2(mg,1)=c2(mg,1)-dt/gamm(mg,1)*dgdtg(mg)
      endif
  214 continue

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

      call trim3x(wc1,wc2,wc3,wrhs,ttgsoil,ms0,npts)

       do k=1,ms
      do 215 i=1,npts
        tgsoil(INDEX1(i),k)=ttgsoil(i,k)
  215 continue
       enddo

      do 216 mg=1,ln2
      if(tester1(mg))then
       sgflux(mg,lg)=0.
       gflux(mg,lg)=coef(mg,2)*(tgsoil(mg,1)-tgsoil(mg,2))
      endif
  216 continue

      endif ! if flag1
c------------------------------------------------------------------------

      if(flag2)then

      do 220 mg=1,ln2
      if(tester2(mg))then
        tzs(mg,1)=sdepth(mg,1)
        tzs(mg,2)=sdepth(mg,2)
        tzs(mg,3)=sdepth(mg,3)
        ccnsw(mg,1)=scond1s(mg,1)
        ccnsw(mg,2)=scond1s(mg,2)
        ccnsw(mg,3)=scond1s(mg,3)
        coef(mg,1)=0.
        coef(mg,ms+4)=0.
      endif 
  220 continue

       do k=1,ms
      do 221 mg=1,ln2
      if(tester2(mg))then
         tzs(mg,k+3)=zsmod(mg,k)
         ew=owb(mg,k)*pssat(mg,lg)
         eww=min(ew,pssat(mg,lg)/2.)
         ccf=max(1.,sqrt(min(2.0,pssat(mg,lg)/(2.*eww))))
         ei=owib(mg,k)*pssat(mg,lg)
         ccnsw(mg,k+3)= 
     &          min(pcnsd(mg,lg)*(60.)**ew*(250.)**ei,2.2)*ccf
         if(pcnsd(mg,lg).gt.2.51) ccnsw(mg,k+3)=2.51
          if(lakesg(mg,lg))ccnsw(mg,k+3)=6.0*2.51
      endif 
  221 continue
       enddo

       do k=2,ms+3
      do 222 mg=1,ln2
      if(tester2(mg))then
          coef(mg,k)=1./
     &       (.5*(tzs(mg,k-1)/ccnsw(mg,k-1)+tzs(mg,k)/ccnsw(mg,k)))
      endif 
  222 continue
       enddo
c
c For subsequent use of routine trim3x : set values 
c
       do  k=1,3
      do 223 mg=1,ln2
      if(tester2(mg))then
          dtg=dt/sgamm(mg,k)
          rhs(mg,k) = tggsn(mg,k,lg)
          c1(mg,k)= -dtg*coef(mg,k)
          c3(mg,k) = -dtg*coef(mg,k+1)
          c2(mg,k)= 1.-c1(mg,k)-c3(mg,k)
          hs1tm(mg,k)=sgamm(mg,k)*tggsn(mg,k,lg)
      endif 
  223 continue
       enddo

       do k=1,ms
      do 224 mg=1,ln2
      if(tester2(mg))then
          wbliqp=owb(mg,k)*pssat(mg,lg)
          wbicep=owib(mg,k)*pssat(mg,lg)
          cfr1=0.
          snowdd=cfr1
calg      gamm(mg,k)=calgamm(mg,lg,k,wbliqp,wbicep,snowdd,zsmod(mg,k))
      calgamm0 = max( (1.-pssat(mg,lg))*pcss(mg,lg)*
     & prhos(mg,lg) + wbliqp*cswat*rhowat + wbicep*csice*rhowat*0.9,
     & pcss(mg,lg)*prhos(mg,lg) ) * zsmod(mg,k)
      if(k.eq.1.and.isflag(mg,lg).eq.0) calgamm0 = calgamm0+
     &   cgsnow*snowdd
          gamm(mg,k)= calgamm0
          hs1tm(mg,k+3)=gamm(mg,k)*tgsoil(mg,k)
          dtg=dt/gamm(mg,k)
          rhs(mg,k+3) = tgsoil(mg,k)
          c1(mg,k+3)= -dtg*coef(mg,k+3)
          c3(mg,k+3) = -dtg*coef(mg,k+4)
          c2(mg,k+3)= 1.-c1(mg,k+3)-c3(mg,k+3)
      endif 
  224 continue
       enddo

CSJP      do 225 mg=1,ln2
CSJP      if(tester2(mg))then
CSJP        rhs(mg,1) = rhs(mg,1)+tbhfl(mg)*dt/sgamm(mg,1)
CSJP      endif
CSJP  225 continue

      do 225 mg=1,ln2
      if(tester2(mg))then
c        rhs(mg,1) = rhs(mg,1)+tbhfl(mg)*dt/sgamm(mg,1) !old version
c                         modified by EAK to remove rare instability
        rhs(mg,1) = rhs(mg,1)+dt/sgamm(mg,1) 
     &              * (tbhfl(mg)-tggsn(mg,1,lg)*dgdtg(mg))
        c2(mg,1)=c2(mg,1)-dt/sgamm(mg,1)*dgdtg(mg)

      endif 
  225 continue

      npts = 0
      do mg=1,ln2
        if(tester2(mg))then
          npts = npts + 1
          INDEX1(npts) = mg
        endif
      enddo
c
      ms3 = ms+3
      do k=1,ms3
      do i=1,npts
        wc1(i,k)  = c1(INDEX1(i),k)
        wc2(i,k)  = c2(INDEX1(i),k)
        wc3(i,k)  = c3(INDEX1(i),k)
        wrhs(i,k)  = rhs(INDEX1(i),k)
      enddo
      enddo

      call trim3x(wc1,wc2,wc3,wrhs,ttgsoil,ms3,npts)

       do k=1,3
      do 226 i=1,npts
         tggsn(INDEX1(i),k,lg)=ttgsoil(i,k)
  226 continue
       enddo

       do k=1,ms
      do 227 i=1,npts
         tgsoil(INDEX1(i),k)=ttgsoil(i,k+3)
  227 continue
       enddo

      do 228 mg=1,ln2
      if(tester2(mg))then
        sgflux(mg,lg)=coef(mg,2)*(tggsn(mg,1,lg)
     &                      -tggsn(mg,2,lg))
        gflux(mg,lg)=coef(mg,4)*(tgsoil(mg,1)-tgsoil(mg,2))
      endif 
  228 continue
c
      endif ! if flag2

      return
      end
c
      function calgamm(mg,lg,k,wbliqp,wbicep,snowdd,zsmod)

      implicit none
C Global parameters
      include 'PARAMS.f'
      real cgsnow
      parameter (cgsnow=2105.)

C Argument list
      real calgamm
      integer mg,lg,k
      real wbliqp
      real wbicep
      real snowdd
      real zsmod

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'SURF1.f'

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
      real calgamm0

C Local data, functions etc
      real csice,cswat,rhowat
      data csice /2.100e3/, cswat /4.218e3/, rhowat /1000.0/

C Start code : ----------------------------------------------------------
   
      calgamm0 = max( (1.-pssat(mg,lg))*pcss(mg,lg)*
     & prhos(mg,lg) + wbliqp*cswat*rhowat + wbicep*csice*rhowat*0.9,
     & pcss(mg,lg)*prhos(mg,lg) ) * zsmod
      if(k.eq.1.and.isflag(mg,lg).eq.0) calgamm0 = calgamm0+
     &   cgsnow*snowdd
c    &   ssdnn(mg,lg)*cgsnow*snowdd/ssdnn(mg,lg)
      calgamm = calgamm0

      return
      end
c
      subroutine snowpr(lg,otggsn,dt,snowd,tester1)

      implicit none

!$OMP THREADPRIVATE ( /LOGIMSL/ )
!$OMP THREADPRIVATE ( /SOILPR/ )

C Global parameters
      include 'PARAMS.f'
      real snmin
      parameter (snmin=0.15)

C Argument list
      integer lg
      real otggsn(ln2)
      real dt
      real snowd(ln2)
      logical tester1(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'LOGIMSL.f'

      real tgsoil,ssatcur,gamm,coef,scond1s,cap1s,sgamm,sdepth
     & ,hs1tm,hs1tp
      common/soilpr/tgsoil(ln2,ms),ssatcur(ln2,ms+1),gamm(ln2,ms)
     &,coef(ln2,ms+4),scond1s(ln2,3),cap1s(ln2,3)
     &,sgamm(ln2,3),sdepth(ln2,3),hs1tm(ln2,9),hs1tp(ln2,9)

C Global data blocks
      include 'SURF1.f'

      real prhos,pcss,pcnsd,pbch,pssat,pswilt,phsbh,psfc,psucs
     & ,phyds,pcdr3,prlaim,prlais,pslveg,pscveg
      integer ip2bp3,ipbp2
      common/psoilvegpr/prhos(ln2,lat),pcss(ln2,lat),
     &pcnsd(ln2,lat),pbch(ln2,lat),pssat(ln2,lat),
     &pswilt(ln2,lat),phsbh(ln2,lat),ip2bp3(ln2,lat),
     &ipbp2(ln2,lat),psfc(ln2,lat),psucs(ln2,lat),
     &phyds(ln2,lat),pcdr3(ln2,lat),prlaim(ln2,lat),
     &prlais(ln2,lat),pslveg(ln2,lat),pscveg(ln2,lat)

C Local work arrays and variables
      real etac(ln2,3)

      integer k
      integer mg

      real ccoef
      real etat
      real excm
      real excd
      real pr
      real sdd
      real sd1
      real sm1
      real smass1n
      real smass2n
      real smass3n
      real smass1o
      real smass2o
      real smass3o
      real ssdn31n
      real ssdn32n
      real ssdn33n
      real ssdn31o
      real ssdn32o
      real ssdn33o
      real tggd
      real tr1
      real tr1d140
      real z1

C Local data, functions etc

C Start code : ----------------------------------------------------------

c
c        dt  - time step 
c        snowd - snow depth 
   
      do 3841 mg=1,ln2
      if(tester1(mg))then

      if( isflag(mg,lg).eq.0) then
          tggsn(mg,1,lg) = tggsl(mg,1,lg)
          tggsn(mg,2,lg) = tggsl(mg,1,lg)
          tggsn(mg,3,lg) = tggsl(mg,1,lg)
          ssdn3(mg,2,lg) = ssdn3(mg,1,lg)
          ssdn3(mg,3,lg) = ssdn3(mg,1,lg)
          sdepth(mg,1) = 0.07
          sdd=(snowd(mg)-0.07*ssdn3(mg,1,lg))/ssdn3(mg,1,lg)
          sdepth(mg,2)=max(0.02,0.45*sdd)
          sdepth(mg,3)=max(0.02,0.55*sdd)
          if(snowd(mg).gt.20.)sdepth(mg,2)=max(0.02,0.3*sdd)
          if(snowd(mg).gt.20.)sdepth(mg,3)=max(0.02,0.7*sdd) 
          smass(mg,1,lg)=        0.07*ssdn3(mg,1,lg)
          smass(mg,2,lg)=sdepth(mg,2)*ssdn3(mg,1,lg)
          smass(mg,3,lg)=sdepth(mg,3)*ssdn3(mg,1,lg)
      endif

      endif
 3841 continue

      do k=1,3
      do 3842 mg=1,ln2
      if(tester1(mg))then
        ccoef=0.0
        if(ssdn3(mg,k,lg) .ge.150.0) ccoef=4.6e-2
        tggd=min(273.1,tggsn(mg,k,lg))
        ssdn3(mg,k,lg) = ssdn3(mg,k,lg)+dt*ssdn3(mg,k,lg)
     &  *3.1e-6*exp(-3.e-2*(273.1-tggd)-
     &  ccoef*(ssdn3(mg,k,lg)-150.0))
c    &  *2.8e-6*exp(-3.e-2*(273.1-tggd)-
c       etat=3.7e7*exp(8.1e-2*(273.1-tggd))
        etat=3.0e7*exp(8.1e-2*(273.1-tggd))
        etac(mg,k)=etat*exp(0.021*ssdn3(mg,k,lg))
      endif
 3842 continue
      enddo

      do 3843 mg=1,ln2
      if(tester1(mg))then
        ssdn3(mg,1,lg)=ssdn3(mg,1,lg)+dt*9.81*ssdn3(mg,1,lg)
     &   /etac(mg,1)*0.5*ssdn3(mg,1,lg)*0.07
        ssdn3(mg,2,lg)=ssdn3(mg,2,lg)+dt*9.81*ssdn3(mg,2,lg)
     &   /etac(mg,2)*(ssdn3(mg,1,lg)*0.07+
     &    0.5*ssdn3(mg,2,lg)*smass(mg,2,lg)/ssdn3(mg,2,lg))
        ssdn3(mg,3,lg)=ssdn3(mg,3,lg)+dt*9.81*ssdn3(mg,3,lg)
     &   /etac(mg,3)*(ssdn3(mg,1,lg)*0.07+
     &    ssdn3(mg,2,lg)*smass(mg,2,lg)/ssdn3(mg,2,lg)+
     &    0.5*ssdn3(mg,3,lg)*smass(mg,3,lg)/ssdn3(mg,3,lg))
      endif
 3843 continue

      do 3844 mg=1,ln2
      if(tester1(mg))then
        tr1 =  snowd(mg)-osnowd(mg,lg)
        ssdn31o=ssdn3(mg,1,lg)
        ssdn32o=ssdn3(mg,2,lg)
        ssdn33o=ssdn3(mg,3,lg)
        smass1o=smass(mg,1,lg)
        smass2o=smass(mg,2,lg)
        smass3o=smass(mg,3,lg)

       if( tr1.ge.0.) then
        tr1d140 =  tr1/140.
        ssdn31n=max((smass1o+tr1)/(smass1o/ssdn31o+tr1d140),140.)
        smass1n=0.07*ssdn31n
        sdepth(mg,1)=0.07
        excm=smass1o+tr1-smass1n
        excd=excm/ssdn31n

        smass2n=max(0.01,smass2o+0.4*excm )
        ssdn32n = max(140.,min(500.,smass2n/
     &                    (smass2o/ssdn32o+0.4*excd)))
        sdepth(mg,2)=max(.02,smass2n/ssdn32n)

        smass3n=max(0.01,snowd(mg)-smass1n-smass2n)
        sdepth(mg,3) =max(0.02,smass3o/ssdn33o+0.6*excm/ssdn32n)
        ssdn33n =max(140.,min(500.,smass3n/sdepth(mg,3)))
        if(ssdn33n.lt.ssdn32n) then
          ssdn33n = ssdn32n
          sdepth(mg,3) =max(0.02,smass3n/ssdn33n)
        endif
        ssdn3(mg,1,lg)=ssdn31n
        ssdn3(mg,2,lg)=ssdn32n
        ssdn3(mg,3,lg)=ssdn33n
        smass(mg,1,lg)=smass1n
        smass(mg,2,lg)=smass2n
        smass(mg,3,lg)=smass3n
      else
c                                            snow melting
c       tr1 =  abs(tr1)/ssdn31o
        sdepth(mg,1)=0.07
c                                            current depth of 
c                                            the 1st layer
        sd1=max(0.005,smass1o/ssdn31o)
        sm1=max(0.01,smass1o) !current mass of the 1st layer
        excd=0.07-sd1 
        ssdn31n=max(140.,min(500.,(sd1*ssdn31o+
     &                        excd*ssdn32o)/0.07))
        smass1n=0.07*ssdn31n
        excm=smass1n-sm1
        excd=excm/ssdn31n

        pr=min(smass2o/(smass3o+smass2o),.9)
        smass2n=max(0.01,smass2o-pr*excm)
        sdepth(mg,2)=max(0.02,smass2o/ssdn32o-pr*excd)
        ssdn32n = max(140.,min(500.,smass2n/sdepth(mg,2)))
        if( ssdn32n .lt. ssdn32o ) then
            ssdn32n=ssdn32o
            smass2n=0.45*(snowd(mg)-smass1n)
            sdepth(mg,2)=max(0.02,smass2n/ssdn32n)
        endif

        smass3n=max(0.01,snowd(mg)-smass1n-smass2n)
        sdepth(mg,3)=max(0.02,smass3n/ssdn33o)
        ssdn3(mg,1,lg)=ssdn31n
        ssdn3(mg,2,lg)=ssdn32n

        smass(mg,1,lg)=smass1n
        smass(mg,2,lg)=smass2n
        smass(mg,3,lg)=smass3n

      endif

      endif
 3844 continue

      do 3845 mg=1,ln2
      if(tester1(mg))then
        isflag(mg,lg)=1
        z1=sdepth(mg,1)+sdepth(mg,2)+sdepth(mg,3)
        ssdnn(mg,lg)=(ssdn3(mg,1,lg)*sdepth(mg,1)+
     &                    ssdn3(mg,2,lg)*sdepth(mg,2)+
     &                    ssdn3(mg,3,lg)*sdepth(mg,3))/z1
        otggsn(mg)=tggsn(mg,1,lg)
      endif
 3845 continue

      do k=1,3
      do 3846 mg=1,ln2
      if(tester1(mg))then
        scond1s(mg,k) = max(0.2,min(2.576e-6*ssdn3(mg,k,lg)*
     &                  ssdn3(mg,k,lg)+0.074,1.0))
        cap1s(mg,k)=ssdn3(mg,k,lg)*2105.
        sgamm(mg,k)=cap1s(mg,k) * sdepth(mg,k)
      endif
 3846 continue
      enddo

      return
      end

      subroutine surfa_riv(ind,lg,zsmod)

      implicit none

!$OMP THREADPRIVATE ( /SURFBC/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer ind
      integer lg
      real zsmod(ln2,ms)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

      real egg,condxg,tsigmf,evapxf,ewwp,pmcmax
      common/surfbc/egg(ln2),condxg(ln2),tsigmf(ln2),evapxf(ln2)
     &             ,ewwp(ln2),pmcmax(ln2)

C Global data blocks
      include 'SURF1.f'

      logical lakesg,lakeind
      common/glakes/lakesg(ln2,lat),lakeind(lat)

      real resvr1,resrem,srunoff,sresvr1
      common/rivers/resvr1(0:lon+1,0:lat2),resrem(lon,lat2)
     & ,srunoff(lon,lat2),sresvr1(lon,lat2)

      real aftfh,aftfhg,ocondx,opreci,pmc
      integer isoilm
      common/surf/aftfh(ln2,lat),aftfhg(ln2,lat),ocondx(ln2,lat)
     &          ,isoilm(ln2,lat),opreci(ln2,lat)
     &          ,pmc(ln2,lat)

C Local work arrays and variables
      integer k
      integer ka
      integer k_unst
      integer k_unfr
      integer lgns
      integer ma
      integer mg
      integer ns

      real avgt
      real sumice
      real sumt
      real sumz

C Local data, functions etc

C Start code : ----------------------------------------------------------

      if(ind.eq.1)then
c-----------------------------------------------------------------------
c Set lakes properties. Fill up puddle water from reservoirs
c-----------------------------------------------------------------------
       do ns=1,2
         lgns=(ns-1)*lg+(lat2+1-lg)*(2-ns)
         do mg=1,lon
          ma=mg+(ns-1)*lon
c For T63 Great Lakes, large pmcmax.  Also keep pmc (puddle depth)
c   at pmcmax by topping up from river/reservoir flow.
c   Put water used to top up puddles (mms of water) into resrem (m of water)
          resrem(mg,lgns)=0.0
          if(lakesg(ma,lg))then
            tsigmf(ma)=0.0 ! No vegetation
            pmcmax(ma)=20.0
            if(pmc(ma,lg).lt.pmcmax(ma))then
              resrem(mg,lgns)=(pmcmax(ma)-pmc(ma,lg))*0.001 ! m water
              pmc(ma,lg)=pmcmax(ma)
            endif
c Set lake layer depths
            zsmod(ma,1)=3.0
            zsmod(ma,2)=5.0
            zsmod(ma,3)=7.0
            zsmod(ma,4)=8.0
            zsmod(ma,5)=8.0
            zsmod(ma,6)=8.0

          endif
         enddo
       enddo
     
       return

      endif ! (ind.eq.1)



      if(ind.eq.2)then
c-----------------------------------------------------------------------
c Make lake ice float at top. Then stabely stratify lake temps
c-----------------------------------------------------------------------

c Make sure that any ice in lakes floats at the surface
       do mg=1,ln2  
        if(lakesg(mg,lg))then

         sumice=0.0
         do k=1,ms
          sumice=sumice+wbice(mg,k,lg)*zsmod(mg,k)
         enddo

         if(sumice.gt.0.0)then
          do k=1,ms
           if(sumice.gt.zsmod(mg,k))then
            wbice(mg,k,lg)=1.0
            sumice=sumice-zsmod(mg,k)
           else
            wbice(mg,k,lg)=sumice/zsmod(mg,k)
            sumice=0.0
           endif
          enddo
         endif

c        write(6,7878)mg,lg,(tggsl(mg,k,lg)-273.0,k=1,ms)
c    &    ,(wbice(mg,k,lg),k=1,ms)
c7878 format(i3,1x,i2,6(1x,f5.1),6(1x,f6.2))

c Stratify unfrozen lake layers. First search for topmost unfrozen 
c  layer ( k = k_unfr ).
         k_unfr=ms
         do k=1,ms
           if(wbice(mg,k,lg).eq.0.0)then
             k_unfr=k
             go to 705
           endif
         enddo

c If k_unfr=ms, then all levels 1:ms-1 frozen. No stabilizing needed.
  705    if(k_unfr.eq.ms)go to 709

c 2 or more layers of unfrozen lake water found, starting at k=k_unfr

c---- search for level k=k_unst where temp(k) < temp(k+1) (unstable)
  706     continue

          do k=k_unfr,ms-1
           if(tggsl(mg,k,lg).lt.tggsl(mg,k+1,lg))then
             k_unst=k
             go to 707 ! Unstable layer found
           endif
          enddo

c All layers are stable if "go to 707" not invoked
          go to 709

c Mix layers downward until stable (i.e. average > temp below)
  707     ka=k_unst
          sumz=zsmod(mg,ka)
          sumt=tggsl(mg,ka,lg)*zsmod(mg,ka)
  708     ka=ka+1
          sumz=sumz+zsmod(mg,ka)
          sumt=sumt+tggsl(mg,ka,lg)*zsmod(mg,ka)
          avgt=sumt/sumz

c Mixing (averaging) has stabilized this part of the column when
c   avgt > temp(ka+1), or when mixing reached bottom level
          if( (avgt.gt.tggsl(mg,ka+1,lg)) .or. (ka.eq.ms) )then
           do k=k_unst,ka
            tggsl(mg,k,lg)=avgt
           enddo
           go to 706 ! Now re-check the column
          endif

          go to 708 ! Go back and mix in the next level down

  709     continue

        endif ! lakesg(mg,lg)
       enddo ! mg=1,ln2

       return

      endif ! (ind.eq.2)

      end
