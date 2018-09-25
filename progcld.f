c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c  $Log: progcld.f,v $
c  Revision 1.36  2001/11/07 06:56:28  rot032
c  Some further minor tuning from HBG
c
c  Revision 1.35  2001/10/12 02:28:09  rot032
c  Merge HBG and LDR changes.
c
c  Revision 1.34  2001/10/12 02:06:56  rot032
c  HBG changes, chiefly for new river-routing scheme and conservation
c
c  Revision 1.33.1.1  2001/10/12 02:13:44  rot032
c  LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c  Revision 1.33  2001/06/04 02:26:53  rot032
c  Changes from LDR to bring sulfur cycle to run P41
c
c  Revision 1.32  2001/03/07 04:28:57  rot032
c  Changes from LDR to bring sulfur cycle to run N63
c
c  Revision 1.31  2001/02/28 04:42:14  rot032
c  Merge LDR and HBG cnahges.
c
c  Revision 1.30  2001/02/28 04:36:36  rot032
c  Further tidy ups from HBG
c
c  Revision 1.29.1.1  2001/02/28 04:40:57  rot032
c  Bring sulfur cycle to run N58.
c
c  Revision 1.29  2001/02/22 06:46:47  rot032
c  Merge LDR and HBG changes.
c
c  Revision 1.28  2001/02/22 05:34:38  rot032
c  Changes from HBG to complete concatenation of NH/SH latitudes
c
c  Revision 1.27.1.1  2001/02/22 05:56:37  rot032
c  Changes from LDR to sulfur cycle, to bring it to run N52
c
c  Revision 1.27  2000/12/08 03:58:52  rot032
c  Changes to aerosol scheme (and small newrain fix) from LDR
c
c  Revision 1.26  2000/11/14 06:47:45  rot032
c  Changes from LDR for aerosol model
c
c  Revision 1.25  2000/06/20 02:27:09  rot032
c  Initial coupling of ECHAM chemical transport model (LDR)
c
c  Revision 1.24  1999/06/30 05:29:37  rot032
c  Pass in new args so that new autoconv treatment can be easily inserted.
c
c  Revision 1.23  1999/05/20 06:23:49  rot032
c  HBG changes to V5-2
c
c Revision 1.22  1998/12/10  00:55:31  ldr
c HBG changes to V5-1-21
c
c Revision 1.21  1997/11/27  05:35:00  ldr
c Flexible treatment of direct and indirect aerosol effects.
c
c Revision 1.20  1997/11/24  23:25:25  ldr
c Indirect sulfates changes to give version used for GRL paper
c
c Revision 1.19  1997/10/06  04:15:03  ldr
c Final corrections to V5-1.
c
c Revision 1.18  1997/10/03  05:45:46  ldr
c Changes for sulfates from LDR
c
c Revision 1.17  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.16  1997/06/11  02:21:29  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.15  1996/11/22  03:41:41  ldr
c Comments and tidy-ups from LDR; results slightly changed by extra
c calculation of qsg in progcld.
c
c******************************************************************************
c
c This is the interface to LDR's cloud scheme.
c
c INPUT/OUTPUT
c
c Input:
c
c parameters from include file PARAMS.f
c      lon - number of points around a latitude circle
c      ln2 - number of points for NH+SH latitude circles sequentially
c      nl - number of vertical levels
c
c from common/fewflags in FEWFLAGS.f
c      debug - namelist flag to control single column debugging
c      lgdebug - latitude index for single column debugging
c      insdebug - hemisphere index for single column debugging
c      mgdebug  - longitude index for single column debugging
c      ukconv - logical flag from UKMO convection scheme
c
c see also include files PHYSPARAMS.f (model physical constants)
c                        CPARAMS.f    (cloud scheme parameters)
c
c from arguments
c      lg - latitude index (ranges from 1 at poles to LAT at equator)
c      ns - hemisphere index
c      land - logical variable for surface type ( = T for land points)
c      tdt - leapfrog timestep (seconds)
c      fluxc - flux of convective rainfall in timestep (kg/m**2)
c      kbase - k index of convective cloud base 
c      ktop - k index of convective cloud top
c      rainx - convective rainfall rate at cloud base (mm/day)
c      CCW - UKMO CONVECTIVE CLOUD LIQUID WATER (KG/KG) ON MODEL LEVELS
c      prf - pressure at full levels (in hPa. NB: not SI units)
c      prh - pressure at half levels (in hPa. NB: not SI units)
c      dprf - layer pressure thickness (in hPa.   Ditto       )
c
c In/Out:
c
c from arguments
c      ttg - temperature (K)
c      qtg - water vapour mixing ratio (kg/kg)
c      qlg - cloud liquid water mixing ratio (kg/kg)
c      qfg - cloud ice mixing ratio (kg/kg)
c      precs - amount of stratiform precipitation in timestep (mm)
c      xtg - Tracer mixing ratio (kg/kg)
c
c Output:
c
c from arguments
c      cfrac - stratiform cloud fraction
c      ccov - stratiform cloud *cover* looking from above (currently = cfrac)
c      qlpath - cloud liquid water path in kg/m**2 (diagnostic only)
c      preci - amount of stratiform snowfall in timestep (mm)
c      cldcon - convective cloud cover in column
c      qlgsav - qlg saved for radiation (avg of values before/after precip.)
c      qfgsav - qfg saved for radiation (avg of values before/after precip.)
c      cfsav - ccov saved for radiation (avg of values before/after precip.)
c      rhg - relative humidity, approximated as q/qs (diagnostic only)
c      qccon - prescribed cloud water mixing ratio of convective clouds (kg/kg)
c      clcon - convective cloud fraction in individual layers
c      qevap - evaporation of rainfall in mixing ratio units (kg/kg)
c      qsubl - sublimation of snowfall in mixing ratio units (kg/kg)
c
c from common /radav in RADAV.f (diagnostics only)
c      qlz - zonal mean cloud liquid water mixing ratio (kg/kg)
c      qfz - zonal mean cloud ice mixing ratio          (kg/kg)
c      cfracz - zonal mean cloud fraction
c      qsublz - zonal mean sublimation of snowfall (kg/kg/s)
c      qevapz - zonal mean evaporation of rainfall (kg/kg/s)
c      qautoz - zonal mean autoconversion of cloud liquid water (kg/kg/s)
c      qcollz - zonal mean collection by rain of cloud liquid water (kg/kg/s)
c      qaccrz - zonal mean accretion by snow of cloud liquid water (kg/kg/s)
c 
c******************************************************************************

      subroutine progcld(lg,tdt,radstep,land,fluxc,kbase,ktop, !Inputs
     &                   rainx,CCW,sg,fscav,xtu,               !Inputs
     &                   ttg,qtg,qlg,qfg,precs,cfrac,ccov,xtg, !In and out
     &                   conwd,so2wd,so4wd,                    !In and out
     &                   so2oh,so2h2,so2o3,dmsoh,dmsn3,        !Outputs
     &                   qlpath,preci,cldcon,qlgsav,qfgsav,    !Outputs
     &                   cfsav,rhg,qccon,clcon,qevap,qsubl,cdso4,qfpath)

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'CPARAMS.f'
      include 'ECPARM.f'
      real TICEuk
      parameter (TICEuk = 258.15) !From UKPARAMS.f. Should be consistent with value there

C Argument list
      integer lg
      real tdt
      logical radstep
      logical land(ln2)
      real fluxc(ln2,nl)
      integer kbase(ln2)
      integer ktop(ln2)
      real rainx(ln2)
      real CCW(ln2,nl) ! UKMO CONVECTIVE CLOUD LIQUID WATER (KG/KG)
      real sg(ln2)
      real fscav(ln2,nl,ntrac)   !Convective tracer scavenging fraction (b/w 0 and 1)
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real precs(ln2)
      real cfrac(ln2,nl)
      real ccov(ln2,nl)
      real xtg(ln2,nl,ntrac)   !Tracer mixing ratio [kg/kg]
      real xtu(ln2,nl,ntrac)   !Tracer mixing ratio in convective updraught[kg/kg]; k=1 at TOA
      real conwd(ln2,ntrac)
      real so2wd(ln2)
      real so4wd(ln2)
      real so2oh(ln2)
      real so2h2(ln2)
      real so2o3(ln2)
      real dmsoh(ln2)
      real dmsn3(ln2)
      real qlpath(ln2)
      real preci(ln2)
      real cldcon(ln2)
      real qlgsav(ln2,nl)
      real qfgsav(ln2,nl)
      real cfsav(ln2,nl)
      real rhg(ln2,nl)
      real qccon(ln2,nl)
      real clcon(ln2,nl)
      real qevap(ln2,nl)
      real qsubl(ln2,nl)
      real cdso4(ln2,nl)
      real qfpath(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'

C Global data blocks
      include 'FEWFLAGS.f' !Input debug, mgdebug, lgdebug, insdebug.
      include 'FILES.f'      !Input so4_ind_rain_file, so4_ind_rad_file
      include 'RADAV.f'  !Output qsublz, qevapz etc. (diagnostics)
      include 'SULFATE1.f'   !Input so4rad,so4rain

      real surfalb,surfso4dir,surfso4rad,surfso4rain
      common /albedo/ surfalb, surfso4dir, surfso4rad, surfso4rain !Tidy this for SCM

C Local work arrays and variables
      real ccrain(ln2,nl)
      real cfa(ln2,nl)
      real dz(ln2,nl)
      real fluxi(ln2,nl)
      real fluxm(ln2,nl)
      real fluxr(ln2,nl)
      real pcfcover(ln2,nl)
      real pclcover(ln2,nl)
      real pcscav(ln2,nl)
      real pdpp1(ln2,nl)
      real pfevap(ln2,nl)
      real pfmelt(ln2,nl)
      real pfprec(ln2,nl)
      real pfconv(ln2,nl)
      real pfsnow(ln2,nl)
      real pfstay(ln2,nl)
      real pfsubl(ln2,nl)
      real pmaccr(ln2,nl)
      real pmiwc(ln2,nl)
      real pmlwc(ln2,nl)
      real pmratep(ln2,nl)
      real pqfsed(ln2,nl)
      real prhop1(ln2,nl)
      real ptp1(ln2,nl)
      real pxtm1(ln2,nl,ntrac)
      real qaccr(ln2,nl)
      real qauto(ln2,nl)
      real qca(ln2,nl)
      real qcl(ln2,nl)
      real qcoll(ln2,nl)
      real qenv(ln2,nl)
      real qsg(ln2,nl)
      real rhoa(ln2,nl)
      real so4m0(ln2)
      real tenv(ln2,nl)
      real wcon(ln2)
      real xte(ln2,nl,ntrac)
      real zk(ln2,0:nl)

      integer k
      integer ma
      integer mg
      integer ncl
      integer ns

      real acon
      real bcon
      real cstrat
      real cvcl
      real fl
      real hm
      real pk
      real so4mk
      real zcfrac
      real zmid
      real zqaccr
      real zqauto
      real zqcoll
      real zqevap
      real zqf
      real zql
      real zqsubl
      real slopes(ln2,nl)
      real plambs(ln2,nl)
      real prscav(ln2,nl)
      real pclcon(ln2,nl)
      real pccw(ln2,nl)
      real qint, rhodz
c Following code is for checking conservation
CD      real ttgold(ln2,nl),rain(ln2),snow(ln2),qlgo(ln2,nl),qfgo(ln2,nl),
CD     &     qtgo(ln2,nl)
CD      real sumt,sumq,conrev,snx,rnx
CD      real sumccwl,sumccwi,sumdq,sumdrs

C Local data, functions etc
      include 'ESTABL.f'

C Start code : ----------------------------------------------------------

c******************************************************************************

CD      ttgold(:,:)=ttg(:,:)
CD      qlgo(:,:)=qlg(:,:)
CD      qfgo(:,:)=qfg(:,:)
CD      qtgo(:,:)=qtg(:,:)

c Diagnostics for sulfur cycle

      if(coupled_aero.and.debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,3i3)')'IPASS=1, before progcld.'
          write(25,9)'DMS ',(xtg(mg,k,1),k=1,nl)
          write(25,9)'SO2 ',(xtg(mg,k,2),k=1,nl)
          write(25,9)'SO4 ',(xtg(mg,k,3),k=1,nl)
        endif
      endif

c Set up resolution dependent parameters

      acon=0.2    !Cloud fraction for non-precipitating convection
      if(nl.lt.24)then
        bcon=0.07  !Rate at which conv cloud frac increases with R.
      else
        bcon=0.05
      endif

c Set up convective cloud and anvil cirrus

      do mg=1,ln2
        if(ktop(mg).gt.0)then
          if(ukconv)then
            cvcl=rainx(mg) ! UKMO conv cloud amount = CCA() from convukmo.f
          else
            cvcl=min(acon + bcon*log(1.0+rainx(mg)),0.8) !NCAR
          endif
          cldcon(mg)=cvcl       !Use this to save total conv cloud
          wcon(mg)=wlc
        else
          cldcon(mg)=0.
          wcon(mg)=0.
        endif
      enddo

      do k=1,nl
        do mg=1,ln2
          rhoa(mg,k)=100.*prf(mg,k)/(rdry*ttg(mg,k))
          dz(mg,k)=100.*dprf(mg,k)/(rhoa(mg,k)*grav)
          pk=100.0*prf(mg,k)
          qsg(mg,k)=qsat(pk,ttg(mg,k))
          rhg(mg,k)=qtg(mg,k)/qsg(mg,k)

          if(k.le.ktop(mg).and.k.ge.kbase(mg))then
            ncl=ktop(mg)-kbase(mg)+1
            clcon(mg,k)=1.0-(1.0-cldcon(mg))**(1.0/ncl)
c            qcl(mg,k)=qsg(mg,k)+wcon(mg)/rhoa(mg,k) !dodgy
            qcl(mg,k)=qsg(mg,k)
            qenv(mg,k)=max(1.e-8,
     &                 qtg(mg,k)-clcon(mg,k)*qcl(mg,k))/(1-clcon(mg,k))
            qcl(mg,k)=(qtg(mg,k)-(1-clcon(mg,k))*qenv(mg,k))/clcon(mg,k)
          else
            clcon(mg,k)=0.
            qcl(mg,k)=0.
            qenv(mg,k)=qtg(mg,k)
          endif
          tenv(mg,k)=ttg(mg,k)
          ccrain(mg,k)=clcon(mg,k) !Assume raining conv cloud frac = clcon
c          pclcon(mg,nlp-k)=clcon(mg,k)
          if(clcon(mg,k).gt.1.e-20)then
            pclcon(mg,nlp-k)=0.1 !Prescribe this for chemistry and wet dep.
c            pclcon(mg,nlp-k)=0.05 !Prescribe this for chemistry and wet dep.
          else
            pclcon(mg,nlp-k)=0.
          endif
          if(ttg(mg,k).gt.TICEuk)then
            pccw(mg,nlp-k)=ccw(mg,k)
          else
            pccw(mg,nlp-k)=0.
          endif
        enddo
      enddo

      if(.not.ukconv)then
        do k=1,nl
          do mg=1,ln2
             CCW(mg,k)=wcon(mg)/rhoa(mg,k)
          enddo
        enddo
      endif
      
c Weight qlg and qfg so that they have the correct environmental values

      do k=1,nl
        do mg=1,ln2
          if(k.ge.kbase(mg).and.k.le.ktop(mg))then
            qlg(mg,k)=qlg(mg,k)/(1-clcon(mg,k))
            qfg(mg,k)=qfg(mg,k)/(1-clcon(mg,k))
          endif
        enddo
      enddo

c Calculate vertical distribution of sulfate mass and cloud droplet conc.
c following the Boucher and Lohmann parameterization. Apply a minimum droplet
c conc. of 20 per cm^3.

c Namelist options:
c naerosol_i=0 ignores sulfates and uses fixed droplet conc. from CPARAMS.f.
c naerosol_i=1 uses annual mean sulfate.
c naerosol_i=2 uses monthly mean sulfate.
c naerosol_i=3 uses interactive sulfate.
c Index in naerosol_i(:) refers to 1st or 2nd indirect effect
c Separate filenames so4_ind_rad_file and so4_ind_rain_file allow different 
c sulfate datasets to be used for radiation and rainfall respectively.
c 
      if(naerosol_i(2).gt.0.and.naerosol_i(2).lt.3)then
        do k=0,nl
          do mg=1,ln2
            zk(mg,k)=0.
          enddo
        enddo
        do k=1,2*nl/3
          do mg=1,ln2
            zk(mg,k)=zk(mg,k-1)+dz(mg,k) !At top of level
          enddo
        enddo

        Hm=1000.  !scale height

        if(SCM)then
          so4rain(1,lg)=surfso4rain
          print*,'surfso4rain= ',surfso4rain
        endif

c Exponentially decreasing with height...

        do mg=1,ln2
          so4m0(mg)=1.e6*so4rain(mg,lg)/Hm    !Convert from g/m**3 to ug/m**3
        enddo


c Or, distribute sulfate uniformly over about 10km in vertical...

C***        do mg=1,ln2
C***          if(land(mg))then
C***            so4m0(mg)=1.e6*so4rain(mg,lg)/Hm !Convert from g/m**3 to ug/m**3        
C***          else
C***            delz=zk(mg,2*nl/3)
C***            so4m0(mg)=1.e6*so4rain(mg,lg)/delz !Convert from g/m**3 to ug/m**3
C***          endif
C***        enddo

        do k=1,2*nl/3
          do mg=1,ln2
            zmid=0.5*(zk(mg,k)+zk(mg,k-1))
            so4mk=so4m0(mg)*exp(-zmid/Hm) !Exponential form
c            so4mk=so4m0(mg) !Uniform form
            if(land(mg))then
c              so4mk=so4m0(mg)*exp(-zmid/Hm) !Exponential form
              cdso4(mg,k)=max(20.e6, 173.8e6*so4mk**0.26)
            else
c              so4mk=so4m0(mg)   !Uniform form
              cdso4(mg,k)=max(20.e6, 114.8e6*so4mk**0.48)
            endif
c            cdso4(mg,k)=1.e6*(90.7*so4mk**0.45+23.) !Hegg 1994            
          enddo
        enddo
        do k=2*nl/3+1,nl
          do mg=1,ln2
            cdso4(mg,k)=20.e6
          enddo
        enddo

      elseif(naerosol_i(2).eq.3)then ! Use ECHAM SO4 to get cdso4.

        do k=1,nl
          do mg=1,ln2
            so4mk=max(1.e-5,3.e9*rhoa(mg,k)*xtg(mg,k,3)) ! x 3 to convert to ug/m3 SO4
            if(land(mg))then
              cdso4(mg,k)=max(20.e6, 173.8e6*so4mk**0.26)
            else
              cdso4(mg,k)=max(20.e6, 114.8e6*so4mk**0.48)
            endif
          enddo
        enddo


      endif

      call newcloud(tdt,lg,land,prf,kbase,ktop,rhoa,cdso4, !Inputs
     &     tenv,qenv,qlg,qfg,   !In and out
     &     cfrac,ccov,cfa,qca)   !Outputs

c Weight variables according to non-convective fraction of grid-box      
      
      do k=1,nl
        do mg=1,ln2
          ttg(mg,k)=clcon(mg,k)*ttg(mg,k)+(1-clcon(mg,k))*tenv(mg,k)
          qtg(mg,k)=clcon(mg,k)*qcl(mg,k)+(1-clcon(mg,k))*qenv(mg,k)
          if(k.ge.kbase(mg).and.k.le.ktop(mg))then
            cfrac(mg,k)=cfrac(mg,k)*(1-clcon(mg,k))
            ccov(mg,k)=ccov(mg,k)*(1-clcon(mg,k))              
            qlg(mg,k)=qlg(mg,k)*(1-clcon(mg,k))
            qfg(mg,k)=qfg(mg,k)*(1-clcon(mg,k))              
          endif
          qlgsav(mg,k)=qlg(mg,k)
          qfgsav(mg,k)=qfg(mg,k)
          cfsav(mg,k)=ccov(mg,k)
        enddo
      enddo

c Calculate precipitation and related processes
      
      call newrain(land,lg,tdt,fluxc,rhoa,dz,ccrain,prf,cdso4,  !Inputs
     &    cfa,qca,
     &    ttg,qlg,qfg,precs,qtg,cfrac,ccov,                     !In and Out
     &    preci,qevap,qsubl,qauto,qcoll,qaccr,fluxr,fluxi,fluxm,!Outputs
     &    pfstay,pqfsed,slopes,prscav)                          !Outputs
      
CD      print *,'Sum dT,dQ ',lg
CD      do 110 mg=1,ln2
CD        rain(mg)=fluxr(mg,1)/tdt
CD        snow(mg)=fluxi(mg,1)/tdt
CD        conrev=100.0/grav/tdt
CD        sumt=0.0
CD        sumdq=0.0
CD        sumccwi=0.0
CD        sumccwl=0.0
CD        do 100 k=1,nl
CD          sumt=sumt+(ttg(mg,k)-ttgold(mg,k))*dprf(mg,k)
CD          sumccwi=sumccwi+(qfg(mg,k)-qfgo(mg,k))*dprf(mg,k)
CD          sumccwl=sumccwl+(qlg(mg,k)-qlgo(mg,k))*dprf(mg,k)
CD          sumdq=sumdq+(qtgo(mg,k)-qtg(mg,k))*dprf(mg,k)
CD 100    continue
CD        sumdrs=(RAIN(mg)+SNOW(mg))/conrev+sumccwl+sumccwi-sumdq !Should = 0
CD        sumq=(hl*RAIN(mg)+(hl+hlf)*SNOW(mg))/conrev/CP
CD        sumccwl=sumccwl*hl/CP
CD        sumccwi=sumccwi*(hl+hlf)/CP
CD        sumq=sumq+sumccwl+sumccwi
CD        rnx=(0.5*tdt)*RAIN(mg)
CD        snx=(0.5*tdt)*SNOW(mg)
CD        if(sumq.gt.0.0)
CD     &       write (6,105) mg,lg,sumq,sumt,(sumq/sumt),rnx,snx
CD     &       ,sumdrs
CD 105    format (1x,2i3,' sumq,sumt ',2e13.6,2x,f6.3,2(1x,f6.2),
CD     &       4(1x,e13.6))
CD 110  continue
CD      if(lg.eq.lat)stop

      if(coupled_aero)then

        pfprec(:,1)=0. !At TOA
        pfmelt(:,1)=0. !At TOA
        pfsnow(:,1)=0. !At TOA
        pfconv(:,1)=0. !At TOA
        do k=1,nl-1
          pfprec(:,nlp-k)=(fluxr(:,k+1)+fluxm(:,k))/tdt !flux *entering* layer k
          pfmelt(:,nlp-k)=fluxm(:,k)/tdt !flux melting in layer k
          pfsnow(:,nlp-k)=(fluxi(:,k+1)-fluxm(:,k))/tdt !flux *entering* layer k
          pfconv(:,nlp-k)=fluxc(:,k)/tdt              !flux *leaving* layer k
        enddo

        do k=1,nl
          pdpp1(:,nlp-k)=100.*dprf(:,k)
          prhop1(:,nlp-k)=rhoa(:,k)
          ptp1(:,nlp-k)=ttg(:,k)
          pxtm1(:,nlp-k,:)=xtg(:,k,:)
          pfevap(:,nlp-k)=qevap(:,k)*rhoa(:,k)*dz(:,k)/tdt !flux evaporating in k
          pfsubl(:,nlp-k)=qsubl(:,k)*rhoa(:,k)*dz(:,k)/tdt !flux sublimating or staying in k
          plambs(:,nlp-k)=slopes(:,k)
c          pcscav(:,nlp-k)=fscav(:,k)
          pcscav(:,nlp-k)=0.
          do mg=1,ln2
            if(qlgsav(mg,k)+qfgsav(mg,k).gt.1.e-12)then
              cstrat=min(cfrac(mg,k), 1.-pclcon(mg,nlp-k))
              pclcover(mg,nlp-k)=cstrat*qlgsav(mg,k)
     &                           /(qlgsav(mg,k)+qfgsav(mg,k)) !Liquid-cloud fraction
              pcfcover(mg,nlp-k)=cstrat*qfgsav(mg,k)
     &                           /(qlgsav(mg,k)+qfgsav(mg,k)) !Ice-cloud fraction
              pmlwc(mg,nlp-k)=qlgsav(mg,k)
              pmiwc(mg,nlp-k)=qfgsav(mg,k)
              pmratep(mg,nlp-k)=(qauto(mg,k)+qcoll(mg,k))/tdt
              pmaccr(mg,nlp-k)=qaccr(mg,k)/tdt
            else
              pclcover(mg,nlp-k)=0.
              pcfcover(mg,nlp-k)=0.
              pmlwc(mg,nlp-k)=0.
              pmiwc(mg,nlp-k)=0.
              pmratep(mg,nlp-k)=0.
              pmaccr(mg,nlp-k)=0.
            endif
          enddo
        enddo

        call xtchemie (1, lg, tdt, pdpp1, pmratep, pfprec,      !Inputs
     &        pclcover, pmlwc, prhop1, ptp1, sg, pxtm1, pfevap, !Inputs
     &        pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,pfstay,pqfsed,  !Inputs
     &        plambs,prscav,pclcon,pccw,pfconv,pcscav,xtu,   !Input
     &        conwd,so2wd, so4wd,                  !In and Out
     &        xte, so2oh, so2h2, so2o3, dmsoh, dmsn3)  !Output

        do k=1,nl
          xtg(:,k,:)=xtg(:,k,:)+xte(:,nlp-k,:)*tdt
        enddo

C***      do k=1,nl
C***        do mg=1,ln2
C***          if(xtg(mg,k,3).gt.1.e-3)then
C***            print*,'so4,lg,mg,k ',xtg(mg,k,3),lg,mg,k
C***          endif
C***        enddo
C***      enddo

      endif

c Recalculate vertical distribution of sulfate mass and cloud droplet conc.
c ** This pass is to give values to go into radiation scheme only. **

      if(radstep.and.naerosol_i(1).gt.0.and.naerosol_i(1).lt.3)then

        if(SCM)then
          so4rad(1,lg)=surfso4rad
          print*,'surfso4rad= ',surfso4rad
        endif

c Repeat this calculation here...

        do k=0,nl
          do mg=1,ln2
            zk(mg,k)=0.
          enddo
        enddo
        do k=1,2*nl/3
          do mg=1,ln2
            zk(mg,k)=zk(mg,k-1)+dz(mg,k) !At top of level
          enddo
        enddo

c Exponentially decreasing with height...

        do mg=1,ln2
          so4m0(mg)=1.e6*so4rad(mg,lg)/Hm   !Convert from g/m**3 to ug/m**3
        enddo

c Or, distribute sulfate uniformly over about 10km in vertical...

C***        do mg=1,ln2
C***          delz=zk(mg,2*nl/3)
C***          so4m0(mg)=1.e6*so4rad(mg,lg)/delz    !Convert from g/m**3 to ug/m**3
C***        enddo

        do k=1,2*nl/3
          do mg=1,ln2
            zmid=0.5*(zk(mg,k)+zk(mg,k-1))
            so4mk=so4m0(mg)*exp(-zmid/Hm) !Exponential form
c            so4mk=so4m0(mg) !Uniform form
            if(land(mg))then
              cdso4(mg,k)=max(20.e6, 173.8e6*so4mk**0.26)
            else
              cdso4(mg,k)=max(20.e6, 114.8e6*so4mk**0.48)
            endif
c            cdso4(mg,k)=1.e6*(90.7*so4mk**0.45+23.) !Hegg 1994            
          enddo
        enddo

      elseif(naerosol_i(1).eq.3)then ! Use ECHAM SO4 to get cdso4.

        do k=1,nl
          do mg=1,ln2
            so4mk=max(1.e-5,3.e9*rhoa(mg,k)*xtg(mg,k,3)) ! x 3 to convert to ug/m3 SO4
            if(land(mg))then
              cdso4(mg,k)=max(20.e6, 173.8e6*so4mk**0.26)
            else
              cdso4(mg,k)=max(20.e6, 114.8e6*so4mk**0.48)
            endif
          enddo
        enddo

      endif

c If using interactive sulfate for direct effect (naerosol_d.eq.3), save it here.
c Factor 1.e3 to convert to g/m2, x 100 for P, x 3 to get sulfate from sulfur

      if(radstep.and.naerosol_d.eq.3)then
        so4dir(:,lg)=0.
        do k=1,nl
          do mg=1,ln2
            so4dir(mg,lg)=so4dir(mg,lg)+3.e5*xtg(mg,k,3)*dprf(mg,k)/grav
          enddo
        enddo
      endif

c Set up cloud variables to be passed into radiation scheme
      
      do k=1,nl
        do mg=1,ln2
          qlgsav(mg,k)=0.5*(qlg(mg,k)+qlgsav(mg,k))
          qfgsav(mg,k)=0.5*(qfg(mg,k)+qfgsav(mg,k))
          cfsav(mg,k)=0.5*(cfsav(mg,k)+ccov(mg,k))
          rhg(mg,k)=qtg(mg,k)/qsg(mg,k) !Not strictly RH!!
        enddo
      enddo

c The convective cloud water amounts are reduced, argument being that values supplied
c by convukmo are for updraft, not cloud as a whole...

c     if(lw.eq.22)then
        qccon(:,:)=clcon(:,:)*CCW(:,:)*0.25
c     else
c       qccon(:,:)=clcon(:,:)*CCW(:,:)*0.5
c     endif


c Some purely diagnostic stuff
c NB: ql, qf and cf exclude convective part.

      if(.not.SCM)then
        do ns=1,2
          ma=(ns-1)*lon

          do k=1,nl
            zql=0.0
            zqf=0.0
            zcfrac=0.0

            zqsubl=0.0
            zqevap=0.0
            zqauto=0.0
            zqcoll=0.0
            zqaccr=0.0
            do mg=1+ma,lon+ma

              zql=zql+qlgsav(mg,k) !Excluding convective part
              zqf=zqf+qfgsav(mg,k)
              zcfrac=zcfrac+cfsav(mg,k)

              zqsubl=zqsubl+qsubl(mg,k)
              zqevap=zqevap+qevap(mg,k)
              zqauto=zqauto+qauto(mg,k)
              zqcoll=zqcoll+qcoll(mg,k)
              zqaccr=zqaccr+qaccr(mg,k)
            enddo

            qlz(lg,ns,k)=qlz(lg,ns,k)+zql/lon ! cloud liquid water
            qfz(lg,ns,k)=qfz(lg,ns,k)+zqf/lon ! cloud ice
            cfracz(lg,ns,k)=cfracz(lg,ns,k)+zcfrac/lon ! cloud fraction

            qsublz(lg,ns,k)=qsublz(lg,ns,k)+zqsubl/tdt/lon !Now in kg/kg/s
            qevapz(lg,ns,k)=qevapz(lg,ns,k)+zqevap/tdt/lon
            qautoz(lg,ns,k)=qautoz(lg,ns,k)+zqauto/tdt/lon
            qcollz(lg,ns,k)=qcollz(lg,ns,k)+zqcoll/tdt/lon
            qaccrz(lg,ns,k)=qaccrz(lg,ns,k)+zqaccr/tdt/lon
          enddo

        enddo                   ! ns=1,2
      endif    !.not.SCM
      
c Add convective cloud water into fields for radiation

      do k=1,nl-1
        do mg=1,ln2
          cfsav(mg,k)=min(1.,cfsav(mg,k)+clcon(mg,k))
          fl=max(0.0,min(1.0,(ttg(mg,k)-ticon)/(273.15-ticon)))
          qlgsav(mg,k)=qlgsav(mg,k)+fl*qccon(mg,k)
          qfgsav(mg,k)=qfgsav(mg,k)+(1.-fl)*qccon(mg,k)
        enddo
      enddo


c Add up liquid water paths 

      do mg=1,ln2
        qlpath(mg)=0.
        qfpath(mg)=0.
      enddo

      do k=1,2*nl/3+1
        do mg=1,ln2
          qlpath(mg)=qlpath(mg)+qlgsav(mg,k)*rhoa(mg,k)*dz(mg,k)
        enddo
      enddo

      do k=1,nl
        do mg=1,ln2
          qfpath(mg)=qfpath(mg)+qfgsav(mg,k)*rhoa(mg,k)*dz(mg,k)
        enddo
      enddo

      if(debug)then

        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,3i3)')'IPASS=1, after progcld.'
          write(25,9)'DMS ',(xtg(mg,k,1),k=1,nl)
          write(25,9)'SO2 ',(xtg(mg,k,2),k=1,nl)
          write(25,9)'SO4 ',(xtg(mg,k,3),k=1,nl)
          write(25,91)'clcon ',(clcon(mg,k),k=1,nl)
          write(25,9)'qccon ',(qccon(mg,k),k=1,nl)
          write(25,91)'rhoa ',(rhoa(mg,k),k=1,nl)
          write(25,91)'dz ',(dz(mg,k),k=1,nl)
          write(25,9)'cfsav ',(cfsav(mg,k),k=1,nl)
          write(25,9)'cfrac ',(cfrac(mg,k),k=1,nl)
          write(25,9)'pclcov',(pclcover(mg,k),k=nl,1,-1)
          write(25,9)'qlgsav ',(qlgsav(mg,k),k=1,nl)
          write(25,9)'qfgsav ',(qfgsav(mg,k),k=1,nl)
          write(25,1)' rainx ',rainx(mg),' cldcon ',cldcon(mg),
     &         'qlpath ',qlpath(mg)
          write(25,*)

          qint=0.
          do k=1,nl
            rhodz=100.*dprf(mg,k)/9.80616
            qint=qint+(qtg(mg,k)+qlg(mg,k)+qfg(mg,k))*rhodz
          enddo
          write(25,1)'integrated q (mm) = ',qint

        endif
      endif
 1    format(3(a,g12.6))
 91   format(a,30f10.3)
 9    format(a,30g10.3)

      
      return
      end
