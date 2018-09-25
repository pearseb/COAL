c (1) Fixing syntax errors identified by the g95 Fortran compiler.
c (2) Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Calculation of fcol removed - this value is not used, and the calculation
c was generating FPEs.
c SJP 2003/04/17
c
c $Log: newrain.f,v $
c Revision 1.52  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.51  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.50.1.1  2001/10/12 02:13:44  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.50  2001/02/22 06:46:47  rot032
c Merge LDR and HBG changes.
c
c Revision 1.49  2001/02/22 05:34:36  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.48.1.1  2001/02/22 05:56:36  rot032
c Changes from LDR to sulfur cycle, to bring it to run N52
c
c Revision 1.48  2000/12/08 03:58:51  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.47  2000/11/14 06:47:45  rot032
c Changes from LDR for aerosol model
c
c Revision 1.46  2000/06/20 02:27:08  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.45  1999/06/30 06:40:57  rot032
c Remove assignment of Ecol
c
c Revision 1.44  1999/06/30 05:29:37  rot032
c Pass in new args so that new autoconv treatment can be easily inserted.
c
c Revision 1.43  1999/05/20 06:23:48  rot032
c HBG changes to V5-2
c
c Revision 1.42  1998/12/10  00:55:28  ldr
c HBG changes to V5-1-21
c
c Revision 1.41  1998/05/26  06:14:02  ldr
c New mixed-phase cloud scheme from LDR (use plates option).
c
c Revision 1.40  1997/12/23  04:10:02  ldr
c Remove tau1, tau2 scaling for evaporation of rain.
c
c Revision 1.39  1997/11/24  23:25:21  ldr
c Indirect sulfates changes to give version used for GRL paper
c
c Revision 1.38  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.37  1997/08/21  02:52:50  ldr
c Reduce cloudiness after rainfall.
c
c Revision 1.36  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.35  1997/06/11  02:21:27  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.34  1996/11/22  03:41:41  ldr
c Comments and tidy-ups from LDR; results slightly changed by extra
c calculation of qsg in progcld.
c
***************************************************************************
c
c This routine is part of the prognostic cloud scheme. It calculates rainfall
c and the evaporation of rain, and also calls icefall, which does the frozen 
c precipitation. It is called by progcld.
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
c
c see also include files PHYSPARAMS.f (model physical constants)
c                        CPARAMS.f    (cloud scheme parameters)
c
c from arguments
c      land - logical variable for surface type ( = T for land points)
c      tdt - leapfrog timestep (seconds)
c      lg - latitude index (ranges from 1 at poles to LAT at equator)
c      fluxc - flux of convective rain in timestep (kg/m**2)
c      rhoa - air density (kg/m**3)
c      dz - layer thicknes (m)
c      ccrain - raining convective cloud fraction at level k
c      prf - pressure at full levels (in hPa. NB: not SI units)
c
c In/Out:
c
c from arguments
c      ttg - temperature (K)
c      qlg - cloud liquid water mixing ratio (kg/kg)
c      qfg - cloud ice mixing ratio (kg/kg)
c      precs - amount of stratiform precipitation in timestep (mm)
c      qtg - water vapour mixing ratio (kg/kg)
c      cfrac - stratiform cloud fraction
c      ccov - stratiform cloud *cover* looking from above (currently = cfrac)
c
c Output:
c
c from arguments
c      preci - amount of stratiform snowfall in timestep (mm)
c      qevap - evaporation of rainfall (kg/kg)
c      qsubl - sublimation of snowfall (kg/kg)
c      qauto - autoconversion of cloud liquid water (kg/kg)
c      qcoll - collection by rain of cloud liquid water (kg/kg)
c      qaccr - accretion by snow of cloud liquid water (kg/kg)
c
***************************************************************************

      subroutine newrain(land,lg,tdt,fluxc,rhoa,dz,ccrain,prf,cdso4,!Inputs
     &                  cfa,qca,
     &                  ttg,qlg,qfg,precs,qtg,cfrac,ccov,   !In and Out
     &                  preci,qevap,qsubl,qauto,qcoll,qaccr,fluxr,fluxi,
     &                  fluxm,pfstay,pqfsed,slopes,prscav)     !Outputs

      implicit none
C Global parameters
      include 'PARAMS.f'     !Input model grid dimensions
      include 'PHYSPARAMS.f' !Input physical constants
      include 'CPARAMS.f'    !Input cloud scheme parameters

C Argument list
      logical land(ln2)
      integer lg
      real tdt
      real fluxc(ln2,nl)
      real rhoa(ln2,nl)
      real dz(ln2,nl)
      real ccrain(ln2,nl)
      real prf(ln2,nl)
      real cdso4(ln2,nl)
      real ttg(ln2,nl)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real precs(ln2)
      real qtg(ln2,nl)
      real cfrac(ln2,nl)
      real ccov(ln2,nl)
      real preci(ln2)
      real qevap(ln2,nl)
      real qsubl(ln2,nl)
      real qauto(ln2,nl)
      real qcoll(ln2,nl)
      real qaccr(ln2,nl)
      real cfa(ln2,nl)
      real qca(ln2,nl)
      real pqfsed(ln2,nl)
      real pfstay(ln2,nl)
      real slopes(ln2,nl)
      real prscav(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f'   !Input debug, lgdebug etc.
      include 'PRINTT.f'
      include 'TIMEX.f'

      character*1 cnmlp(ln2,lat),rainp(ln2,lat)
     & ,snmap(ln2,lat),snwice(ln2,lat),snicel(ln2,lat)
      common/cloud1/cnmlp,rainp,snmap,snwice,snicel

      integer jclb,nlow,nmid,ktied,kbgel,klcmc,klowrn,kmidrn
      real aftsea
      common/levdata/jclb,nlow,nmid,aftsea,ktied,kbgel
     &,klcmc,klowrn,kmidrn

C Local work arrays and variables
      real clfr(ln2,nl)
      real qsg(ln2,nl),fluxr(ln2,nl)
      real cifr(ln2,nl)
      real fluxm(ln2,nl)
      real frclr(ln2,nl)
      real fluxi(ln2,nl),qsl(ln2,nl)
      real fluxa(ln2,nl)
      real clfra(ln2)
      real ccra(ln2)
      real cfrain(ln2,nl)
      real cfmelt(ln2,nl)
      real fluxrain(ln2)
      real Cdrop(ln2,nl)
      real fracr(ln2,nl)
      integer irain(ln2,nl)
c     logical qlgtest,qfgtest
c     double precision Frpr, Frprb, Frprc

      integer irn
      integer jh
      integer k
      integer mb
      integer mg
      integer njumps
      integer ns
      integer nt

      real apr
      real bpr
      real bl
      real cdt
      real cev
      real cftemp
      real clrevap
      real coll
      real crate
      real delt
      real dql
      real dqsdt
      real es
      real evap
CSJP      real fcol
      real fr
      real frb
      real frc
      real pk
      real qcic
      real qcrit
      real ql
      real ql1
      real ql2
      real qpf
      real rhodz
      real satevap
      real selfcoll
      real tk

C Local data, functions etc
      real esdiff(-40:0)  !Table of es(liq) - es(ice) (MKS), -40 <= t <= 0 C
      data esdiff / 
     & 6.22, 6.76, 7.32, 7.92, 8.56, 9.23, 9.94,10.68,11.46,12.27,
     & 13.11,13.99,14.89,15.82,16.76,17.73,18.70,19.68,20.65,21.61,
     & 22.55,23.45,24.30,25.08,25.78,26.38,26.86,27.18,27.33,27.27,
     & 26.96,26.38,25.47,24.20,22.51,20.34,17.64,14.34,10.37, 5.65,
     & 0.08 /
      character*1 cont(8)
      data cont/' ','L','M','4','H','5','6','7'/
      include 'ESTABL.f'
      real pow75,x
      pow75(x)=sqrt(x*sqrt(x))

C Start code : ----------------------------------------------------------

      delt=tdt
      njumps=nint(tdt/delt)
      delt=tdt/njumps !To make sure tdt it a multiple of delt

      do k=1,nl
        do mg=1,ln2
          fracr(mg,k)=0.
          fluxr(mg,k)=0.
          frclr(mg,k)=0.
          fluxm(mg,k)=0.
          fluxi(mg,k)=0.
          fluxa(mg,k)=0.
          qevap(mg,k)=0.
          qauto(mg,k)=0.
          qcoll(mg,k)=0.
          cfrain(mg,k)=0.
          pk=100.0*prf(mg,k)
          qsg(mg,k)=qsati(pk,ttg(mg,k))
          if(qfg(mg,k).gt.0.)then
            cifr(mg,k)=cfrac(mg,k)*qfg(mg,k)/(qlg(mg,k)+qfg(mg,k))
          else
            cifr(mg,k)=0.
          endif
c          clfr(mg,k)=max(cfrac(mg,k)-cifr(mg,k),0.) 
c Previous line not good for roundoff; next few lines are better
          if(qlg(mg,k).gt.0.)then
            clfr(mg,k)=cfrac(mg,k)*qlg(mg,k)/(qlg(mg,k)+qfg(mg,k))
          else
            clfr(mg,k)=0.
          endif
        enddo
      enddo

c Define Cdrop

      if(naerosol_i(2).gt.0)then
        do k=1,nl-1
        do mg=1,ln2
          Cdrop(mg,k)=cdso4(mg,k)
        enddo
        enddo
      else
        do mg=1,ln2
          if(land(mg))then
            Cdrop(mg,1)=Cdropl
          else
            Cdrop(mg,1)=Cdrops
          endif
        enddo
        do k=2,nl-1
        do mg=1,ln2
          Cdrop(mg,k)=Cdrop(mg,1)
        enddo
        enddo
      endif


***************** Cut here to insert new auto scheme ********************            

      do k=nl-1,1,-1
        do mg=1,ln2
            cfrain(mg,k)=0.0
            rhodz=rhoa(mg,k)*dz(mg,k)

            if(clfr(mg,k).gt.0.)then

c Using new autoconv scheme...            

C***              ql=qlg(mg,k)
C***              cfla=0.
C***              dqla=0.
C***              if(cfa(mg,k).gt.0.)then
C***                qcrit=(4*pi/3)*rhow*Rcm**3*Cdrop(mg,k)/rhoa(mg,k)
C***                cfla=cfa(mg,k)*clfr(mg,k)/(clfr(mg,k)+cifr(mg,k))
C***                qla=qca(mg,k)/cfa(mg,k)
C***                if(qla.gt.qcrit)then
C***                  Crate=Aurate*
C***     &                rhoa(mg,k)*(rhoa(mg,k)/(Cdrop(mg,k)*rhow))**(1./3)
C***c                  Crate=0.1818*Crate !TEMPORARY
C***                  ql1=1/pow75(qla**(-4./3)+(4./3)*Crate*delt)
C***c                  ql1=qla-crate*delt !Simple explicit
C***                  ql1=max(ql1, qcrit) !Intermediate qlg after auto
C***                  Frb=dz(mg,k)*rhoa(mg,k)*(qla-ql1)/delt
C***                  cdt=delt*0.5*Ecol*0.24*pow75(Frb)
C***                  selfcoll=min(ql1,ql1*cdt)
C***c                  selfcoll=min(ql1,ql1*cdt/(1+0.5*cdt))
C***                  ql2=ql1-selfcoll
C***                  cfrain(mg,k)=cfla
C***                else
C***                  ql2=qla
C***                endif
C***                dqla=cfla*(qla-ql2)
C***                ql=qlg(mg,k)-dqla
C***                dqla=cfla*(qla-ql1)
C***              endif
C***              dql=qlg(mg,k)-ql

c Or, using old autoconv scheme...

              qcrit=(4*pi/3)*rhow*Rcm**3*Cdrop(mg,k)/rhoa(mg,k)
              qcic=qlg(mg,k)/clfr(mg,k) !In cloud value

              if(qcic.lt.qcrit)then
                ql=qlg(mg,k)
              else
                Crate=Aurate*
     &                rhoa(mg,k)*(rhoa(mg,k)/(Cdrop(mg,k)*rhow))**(1./3)
                ql1=1./pow75(qcic**(-4./3)+(4./3)*Crate*delt)
                ql1=max(ql1, qcrit) !Intermediate qlg after auto
                Frb=dz(mg,k)*rhoa(mg,k)*(qcic-ql1)/delt
                cdt=delt*0.5*Ecol*0.24*pow75(Frb) !old
                  selfcoll=min(ql1,ql1*cdt)
c                selfcoll=min(ql1,ql1*cdt/(1+0.5*cdt))
                ql2=ql1-selfcoll
                ql=clfr(mg,k)*ql2
                cfrain(mg,k)=clfr(mg,k)
              endif
              dql=qlg(mg,k)-ql

              qauto(mg,k)=qauto(mg,k)+dql
              qlg(mg,k)=qlg(mg,k)-dql
              fluxa(mg,k)=dql*rhodz
            endif
          enddo
        enddo


c Call frozen precipitation routine

      call icefall (lg,tdt,rhoa,dz,prf, !Inputs
     &              ttg,qsg,qlg,qfg,qtg,cfrac,cfmelt,            !In and Out
     &              fluxi,fluxm,clfr,cifr,qsubl,qaccr,pfstay,pqfsed,
     &              slopes)            !Outputs




      do nt=1,njumps
        
        do mg=1,ln2
          clfra(mg)=0.
          ccra(mg)=0.
          fluxrain(mg)=0.
          prscav(mg,1)=0.
        enddo
        
c Now work down through the levels...
        
        do k=nl-1,1,-1
          do mg=1,ln2
            rhodz=rhoa(mg,k)*dz(mg,k)
            evap=0.

c Add flux of melted snow and flux of rain due to autoconversion to fluxrain

            fluxrain(mg)=fluxrain(mg)+(fluxm(mg,k)+fluxa(mg,k+1))/njumps
            
c Evaporation of rain
            
            if(fluxrain(mg).gt.0.and.cfrac(mg,k).lt.1.0)then
              pk=100.0*prf(mg,k)
              qsg(mg,k)=qsati(pk,ttg(mg,k))
              if(ttg(mg,k).lt.273.15.and.ttg(mg,k).ge.tice)then
                qsl(mg,k)=qsg(mg,k)+epsil
     &             *esdiff(nint(ttg(mg,k)-273.15))/(100.0*prf(mg,k))
              else
                qsl(mg,k)=qsg(mg,k)
              endif             !qsl is qs value over liquid surface
              qpf=fluxrain(mg)/rhodz !Mix ratio of rain which falls into layer
              Tk=ttg(mg,k)
              es=qsl(mg,k)*pk/epsil 
              Apr=(hl/(rKa*Tk))*(hl/(rvap*Tk)-1)
              Bpr=rvap*Tk/((Dva/pk)*es)
              Fr=fluxrain(mg)/delt/clfra(mg)
              Cev=clfra(mg)
     &       *3.8e2*sqrt(Fr/rhoa(mg,k))/(qsl(mg,k)*(Apr+Bpr))
              dqsdt=hl*qsl(mg,k)/(rvap*ttg(mg,k)**2)
              bl=1+0.5*Cev*delt*(1+hlcp*dqsdt)
              evap=delt*(Cev/bl)*(qsl(mg,k)-qtg(mg,k))
              satevap=(qsl(mg,k)-qtg(mg,k))/(1+hlcp*dqsdt) !Evap to saturate
c              Vr=11.3*Fr**(1./9.)/sqrt(rhoa(mg,k)) !Actual fall speed
c              Vr=5./sqrt(rhoa(mg,k))                !Nominal fall speed

              clrevap=(1-clfr(mg,k))*qpf
              evap=max(0., min(evap,satevap,qpf,clrevap))
              frclr(mg,k)=rhodz*(clrevap-evap)          !over delt
              qevap(mg,k)=qevap(mg,k)+evap
              qtg(mg,k)=qtg(mg,k)+evap
              ttg(mg,k)=ttg(mg,k) - hlcp*evap
            endif
            

c Now do the collection term.

c            if(qlg(mg,k).gt.1.e-10)then
              if(fluxrain(mg).gt.0.)then
                Fr=fluxrain(mg)/clfra(mg)/delt
                cfrain(mg,k)= min(1.,cfrain(mg,k)+clfr(mg,k)*clfra(mg))
c                cfrain(mg,k)= min(clfra(mg),clfr(mg,k)) !max overlap
              else
                Fr=0.
              endif


              if(fluxc(mg,k+1).gt.0.)then
                Frc=max(0.,fluxc(mg,k+1)/max(ccra(mg),0.01)/tdt) ! over tdt
                cfrain(mg,k)= min(1.,cfrain(mg,k)+clfr(mg,k)*ccra(mg))
              else
                Frc=0.
              endif

c The collection term comprises collection by stratiform rain falling from
c above (Fr), stratiform rain released in this grid box (Frb), and
c convective rain (Frc).
c Frb term now done above.

CSJP              fcol=min(1.,clfra(mg)/clfr(mg,k))
c              cdt=delt*Ecol*0.24*(fcol*pow75(Fr)       !Max o'lap
              cdt=delt*Ecol*0.24*(clfra(mg)*pow75(Fr) !Random o'lap
     &                           +ccra(mg)*pow75(Frc))
c              prscav(mg,nlp-k)=cdt/Ecol !Inc conv part
              prscav(mg,nlp-k)=delt*0.24*clfra(mg)*pow75(Fr) !Strat only

              coll=min(qlg(mg,k),qlg(mg,k)*cdt/(1+0.5*cdt))
              qcoll(mg,k)=qcoll(mg,k)+coll
              qlg(mg,k)=qlg(mg,k)-coll
              fluxrain(mg)=fluxrain(mg)+coll*rhodz
c            endif

c subtract evaporated rain

            fluxrain(mg)=fluxrain(mg)-rhodz*evap
            fluxrain(mg)=max(fluxrain(mg),0.) !To avoid roundoff -ve's

c Calculate the raining cloud cover down to this level, for stratiform (clfra)
c and convective (ccra).

            cfrain(mg,k)=min(1.,
     &             cfrain(mg,k)+cfmelt(mg,k)-cfrain(mg,k)*cfmelt(mg,k))
            if(frclr(mg,k).lt.1.e-15)then
              clfra(mg)=max(1.e-15,cfrain(mg,k))
            else
              clfra(mg)=max(1.e-15,
     &                  clfra(mg)+cfrain(mg,k)-clfra(mg)*cfrain(mg,k))
c              clfra(mg)=max(clfra(mg),cfrain(mg,k)) !max overlap
            endif
            ccra(mg)=ccra(mg)+ccrain(mg,k)-ccra(mg)*ccrain(mg,k)
            fracr(mg,k)=clfra(mg)
            fluxr(mg,k)=fluxa(mg,k)+fluxrain(mg) !Flux leaving layer k
          enddo
        enddo
        
      enddo

c Factor 0.5 here accounts for leapfrog scheme

      do mg=1,ln2
        precs(mg)=precs(mg)+0.5*(fluxr(mg,1)+fluxi(mg,1))
c preci() may already have contribution from ukconv
        preci(mg)=preci(mg)+0.5*fluxi(mg,1)
      enddo

      if ((mod(mins+int(mstep, 8),1440_8).eq.0_8).and.cvrnm) then
c**** update the rainfall map at end of day only
c**** Low rainfall (up to sig=0.9), Middle Rain (up to sig=0.45), High above
c**** (see initax.f for klowrn,kmidrn)
        do k=1,nl
          do mg=1,ln2
            if((fluxr(mg,k)+fluxi(mg,k)).gt.0.0)then
              irain(mg,k)=1
            else
              irain(mg,k)=0
            endif
          enddo
        enddo
        do mg=1,ln2
          irn=1
          jh=0
          do k=kmidrn+1,nl
            jh=jh+irain(mg,k)
          enddo
          if (jh.ge.1) irn=5
          jh=0
          do k=klowrn+1,kmidrn
            jh=jh+irain(mg,k)
          enddo
          if (jh.ge.1) irn=irn+2
          jh=0
          do k=1,klowrn
            jh=jh+irain(mg,k)
          enddo
          if (jh.ge.1) irn=irn+1
          rainp(mg,lg)=cont(irn)
        enddo
      end if

c Remove small amounts of cloud

      do k=1,nl
        do mg=1,ln2
          if(qlg(mg,k).lt.1e-10.or.clfr(mg,k).lt.1e-5)then
            qtg(mg,k)=qtg(mg,k)+qlg(mg,k)
            ttg(mg,k)=ttg(mg,k)-hlcp*qlg(mg,k)
            qlg(mg,k)=0.
            clfr(mg,k)=0.
          endif
          if(qfg(mg,k).lt.1e-10.or.cifr(mg,k).lt.1e-5)then
            qtg(mg,k)=qtg(mg,k)+qfg(mg,k)
            ttg(mg,k)=ttg(mg,k)-hlscp*qfg(mg,k)
            qfg(mg,k)=0.
            cifr(mg,k)=0.
          endif
        enddo
      enddo

c Adjust cloud fraction (and cloud cover) after precipitation

      do k=1,nl-1
        do mg=1,ln2
c Next 3 lines commented for m35-m38 runs.
C***          if(qlgsav(mg,k).gt.0.)then
C***            clfr(mg,k)=clfr(mg,k)*qlg(mg,k)/qlgsav(mg,k)
C***          endif
          cftemp=min(cifr(mg,k)+clfr(mg,k), 1.)
          if(cfrac(mg,k).gt.0.)then
            ccov(mg,k)=min(1., ccov(mg,k)*cftemp/cfrac(mg,k))
          else
            ccov(mg,k)=cftemp
          endif
          cfrac(mg,k)=cftemp
        enddo
      enddo

      if(debug)then
        do k=1,nl
          do mg=1,ln2
            if(cfrac(mg,k).gt.0)then
              if(qlg(mg,k).le.0.and.qfg(mg,k).le.0)then
                print*,'end newrain cloud with zero/neg water: cfrac=',
     &               cfrac(mg,k),' ccov= ',ccov(mg,k)
                print*,'cifr ,clfr ',cifr(mg,k),clfr(mg,k)
                print*,'qfg, qlg=', qfg(mg,k),qlg(mg,k)
                ns=(mg-1+lon)/lon
                mb=mg-(ns-1)*lon
                print*,'lg,ns,mg,k ',lg,ns,mb,k
                print*
              endif
            else
              if(qlg(mg,k).gt.0.or.qfg(mg,k).gt.0)then
               print*,'end newrain cloud water with no cloud: qfg,qlg=',
     &               qfg(mg,k),qlg(mg,k)
               print*,'cifr ,clfr ',cifr(mg,k),clfr(mg,k)
               ns=(mg-1+lon)/lon
               mb=mg-(ns-1)*lon
               print*,'ns,lg,k,mg ',ns,lg,k,mb
               print*
              endif
            endif
          enddo
        enddo

c Diagnostics for debugging

        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,3i3)')'IPASS=1, after newrain.'
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,9)'qtg ',(qtg(mg,k),k=1,nl)
          write(25,9)'qsg ',(qsg(mg,k),k=1,nl)
c          write(25,9)'qsl ',(qsl(mg,k),k=1,nl)
          write(25,1)'precs ',precs(mg)
          write(25,9)'cfrac ',(cfrac(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,9)'cdso4 ',(cdso4(mg,k),k=1,nl)
          write(25,9)'cfa ',(cfa(mg,k),k=1,nl)
          write(25,9)'qca ',(qca(mg,k),k=1,nl)
          write(25,9)'fluxr ',(fluxr(mg,k),k=1,nl)
          write(25,9)'frclr ',(frclr(mg,k),k=1,nl)
          write(25,9)'fluxi ',(fluxi(mg,k),k=1,nl)
          write(25,9)'fluxc ',(fluxc(mg,k),k=1,nl)
          write(25,9)'qcoll ',(qcoll(mg,k),k=1,nl)
          write(25,9)'qauto ',(qauto(mg,k),k=1,nl)
          write(25,9)'qaccr ',(qaccr(mg,k),k=1,nl)
          write(25,9)'qevap ',(qevap(mg,k),k=1,nl)
          write(25,9)'ccrain ',(ccrain(mg,k),k=1,nl-1)
          write(25,9)'cfrain ',(cfrain(mg,k),k=1,nl-1)
          write(25,9)'fracr ',(fracr(mg,k),k=1,nl-1)
          write(25,*)
        endif
      endif
 1    format(3(a,f10.5))
 91   format(a,30f10.3)
 9    format(a,30g10.3)


      return

      end
