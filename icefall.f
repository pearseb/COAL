c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c $Log: icefall.f,v $
c Revision 1.35  2001/02/22 05:56:37  rot032
c Changes from LDR to sulfur cycle, to bring it to run N52
c
c Revision 1.34  2000/12/08 03:58:51  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.33  2000/11/14 06:47:45  rot032
c Changes from LDR for aerosol model
c
c Revision 1.32  2000/06/20 02:54:07  rot032
c Merge of HBG and LDR changes
c
c Revision 1.31  2000/06/20 02:08:31  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.30.1.1  2000/06/20 02:27:09  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.30  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.29  1998/12/10  00:55:29  ldr
c HBG changes to V5-1-21
c
c Revision 1.28  1998/05/26  05:42:26  ldr
c Qcloud changes from LDR: new icefall Kessler option, Platt's (1996)
c optical properties and emissivity "fix" for low clouds.
c
c Revision 1.27  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.26  1997/08/21  02:52:50  ldr
c Reduce cloudiness after rainfall.
c
c Revision 1.25  1997/07/24  05:42:48  ldr
c Tuning of ukconv version and T63 fixes from HBG
c
c Revision 1.24  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.23  1997/05/30  07:25:32  ldr
c Reinstate the RCS header!
c
c******************************************************************************
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
c
c      lg - latitude index (ranges from 1 at poles to LAT at equator)
c      tdt - leapfrog timestep (seconds)
c      rhoa - air density (kg/m**3)
c      dz - layer thicknes (m)
c      prf - pressure at full levels (in hPa. NB: not SI units)
c
c In/Out:
c
c from arguments
c
c      ttg - temperature (K)
c      qsg - saturation mixing ratio (kg/kg)
c      qlg - cloud liquid water mixing ratio (kg/kg)
c      qfg - cloud ice mixing ratio (kg/kg)
c      qtg - water vapour mixing ratio (kg/kg)
c      cfrac - stratiform cloud fraction
c      cfmelt - fraction of grid box occupied by falling ice that melts
c
c Output:
c
c from arguments
c
c      fluxi - flux of falling ice in timestep (kg/m**2)
c      fluxm - flux of falling ice that melts in timestep (kg/m**2)
c      clfr - liquid cloud fraction
c      qsubl - sublimation of snowfall (kg/kg)
c      qaccr - accretion by snow of cloud liquid water (kg/kg)
c      pqfsed - (dqf/qf) due to ice falling out of layer (fraction)
c      pfstay - incoming flux of ice staying in layer  (kg/m**2/s)
c
c******************************************************************************

      subroutine icefall (lg,tdt,rhoa,dz,prf,               !Inputs
     &                   ttg,qsg,qlg,qfg,qtg,cfrac,cfmelt,  !In and Out
     &                   fluxi,fluxm,clfr,cifr,qsubl,qaccr, !Outputs
     &                   pfstay,pqfsed,slopes)              !Outputs

C Global parameters
      include 'PARAMS.f'     !Input model grid dimensions
      include 'PHYSPARAMS.f' !Input physical constants
      include 'CPARAMS.f'    !Input cloud scheme parameters

C Argument list
      integer lg
      real tdt
      real rhoa(ln2,nl)
      real dz(ln2,nl)
      real prf(ln2,nl)
      real ttg(ln2,nl)
      real qsg(ln2,nl)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real qtg(ln2,nl)
      real cfrac(ln2,nl)
      real cfmelt(ln2,nl)
      real fluxi(ln2,nl)
      real fluxm(ln2,nl)
      real clfr(ln2,nl)
      real cifr(ln2,nl)
      real qsubl(ln2,nl)
      real qaccr(ln2,nl)
      real pqfsed(ln2,nl)
      real pfstay(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f'   !Input debug, lgdebug etc.

C Local work arrays and variables
      real fluxice(ln2)
      real rhoi(ln2,nl)
      real cifra(ln2)
      real rica(ln2)
      real vi2(ln2,nl),Csbsav(ln2,nl-1),gam(ln2,nl-1),slopes(ln2,nl)
c      real qfdiv(ln2,nl)
c      real qaggi(ln2,nl)
c      real qacci(ln2,nl)
      real fout(ln2,nl-1)
      real fthru(ln2,nl-1)
c      real alpha1(ln2,nl)

C Local data, functions etc
      real esdiff(-40:1)  !Table of es(liq) - es(ice) (MKS), -40 <= t <= 1 C
      data esdiff / 
     & 6.22, 6.76, 7.32, 7.92, 8.56, 9.23, 9.94,10.68,11.46,12.27,
     & 13.11,13.99,14.89,15.82,16.76,17.73,18.70,19.68,20.65,21.61,
     & 22.55,23.45,24.30,25.08,25.78,26.38,26.86,27.18,27.33,27.27,
     & 26.96,26.38,25.47,24.20,22.51,20.34,17.64,14.34,10.37, 5.65,
     & 0.08, 0. /

C Start code : ----------------------------------------------------------

c Set up timestep for ice processes

      delt=tdt
c      delt=360.
      njumps=nint(tdt/delt)
      delt=tdt/njumps !To make sure tdt it a multiple of delt


c Convert from mixing ratio to density of ice, and work out ice cloud fraction

      do k=1,nl
        do mg=1,ln2
          fluxi(mg,k)=0.
          rhoi(mg,k)=rhoa(mg,k)*qfg(mg,k)
          if(qfg(mg,k).gt.0.)then
            cifr(mg,k)=cfrac(mg,k)*qfg(mg,k)/(qlg(mg,k)+qfg(mg,k))
          else
            cifr(mg,k)=0.
          endif
          clfr(mg,k)=max(cfrac(mg,k)-cifr(mg,k),0.)
          qsubl(mg,k)=0.
          qaccr(mg,k)=0.
c          qaggi(mg,k)=0.
c          qacci(mg,k)=0.
c          qfdiv(mg,k)=0.
          cfmelt(mg,k)=0.
        enddo
      enddo

      if(debug)then
        if(lg.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, before icefall.'
          write(25,91)'cfrac ',(cfrac(mg,k),k=1,nl)
          write(25,91)'cifr ',(cifr(mg,k),k=1,nl)
          write(25,91)'clfr ',(clfr(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,*)
        endif

        do k=1,nl
          do mg=1,ln2
            if(cfrac(mg,k).gt.0)then
              if(qlg(mg,k).eq.0.and.qfg(mg,k).eq.0)then
                print*,
     &         'start icefall cloud with zero water: cfrac=',cfrac(mg,k)
              endif
            else
              if(qlg(mg,k).gt.0.or.qfg(mg,k).gt.0)then
             print*,'start icefall cloud water with no cloud: qfg,qlg=',
     &               qfg(mg,k),qlg(mg,k)
              endif
            endif
          enddo
        enddo
      endif

c Set up ice fall speed field and other arrays

      do mg=1,ln2
        vi2(mg,nl)=0.1
        slopes(mg,nl)=0.
      enddo

      do k=nl-1,1,-1
        do mg=1,ln2
          if(cifr(mg,k).gt.0.)then
c           vi2(mg,k)=3.23*(rhoi(mg,k)/cifr(mg,k))**0.17
c           vi2(mg,k)=0.9*3.23*(rhoi(mg,k)/cifr(mg,k))**0.17
            vi2(mg,k)=max(0.1,2.05+0.35*log10(qfg(mg,k)/cifr(mg,k)))
          else
            vi2(mg,k)=vi2(mg,k+1)
          endif
          tc=ttg(mg,k)-273.15
          slopes(mg,k)=1.6e3*10**(-0.023*tc)
          alphaf=hls*qsg(mg,k)/(rvap*ttg(mg,k)**2)
          gam(mg,k)=hlscp*alphaf !(L/cp)*dqsdt (HBG notation)

c Set up the Rate constant for snow sublimation
          
          Tk=ttg(mg,k)
          pk=100*prf(mg,k)
          es=qsg(mg,k)*pk/epsil
          Aprpr=(hls/(rKa*Tk))*(hls/(rvap*Tk)-1)
          Bprpr=rvap*Tk/((Dva/pk)*es)
         curly=0.65*slopes(mg,k)**2+0.493*slopes(mg,k)!Factor in curly brackets
     &         *sqrt(slopes(mg,k)*vi2(mg,k+1)*rhoa(mg,k)/um)

c Define the rate constant for sublimation of snow, omitting factor rhoi

          Csbsav(mg,k)=4*curly/
     &          (rhoa(mg,k)*qsg(mg,k)*(Aprpr+Bprpr)*pi*vi2(mg,k+1)*rhos)

c Set up the parameters for the flux-divergence calculation

          alph=delt*vi2(mg,k)/dz(mg,k)
          fout(mg,k)=1-exp(-alph) !analytical
          fthru(mg,k)=1-fout(mg,k)/alph !analytical
c          alpha1(mg,k)=1.e-3*exp(0.025*(ttg(mg,k)-273.15))
        enddo
      enddo

c Save sedimentation rate for aerosol scheme

      do k=1,nl-1
        do mg=1,ln2
c          pqfsed(mg,nlp-k)=fout(mg,k)*qfg(mg,k)/tdt
          pqfsed(mg,nlp-k)=fout(mg,k)
        enddo
      enddo
      do mg=1,ln2
        pqfsed(mg,1)=0.
      enddo


      do nt=1,njumps !Need njumps=1 for new code

c Assume no cloud at top level

        do mg=1,ln2
          fluxice(mg)=0.
          cifra(mg)=0.
          rica(mg)=0.
          pfstay(mg,1)=0.
        enddo

c Now work down through the levels...

        do k=nl-1,1,-1
          do mg=1,ln2

            sublflux=0.
            fsclr=0.
            caccr=0.
            dqf=0.

c Melt falling ice if > 0 deg C

            if(ttg(mg,k).gt.273.15.and.fluxice(mg).gt.0.)then
              qif=fluxice(mg)/(dz(mg,k)*rhoa(mg,k))      !Mixing ratio of ice
              dttg=-hlfcp*qif
              ttg(mg,k)=ttg(mg,k)+dttg
              qsg(mg,k)=qsg(mg,k) + gam(mg,k)*dttg/hlscp
              fluxm(mg,k)=fluxm(mg,k)+fluxice(mg)
              cfmelt(mg,k)=cifra(mg)
              fluxice(mg)=0.
              cifra(mg)=0.
            endif


c Compute the sublimation of ice falling from level k+1 into level k

            if(fluxice(mg).gt.0.and.qtg(mg,k).lt.qsg(mg,k))then ! sublime snow

              Csb=Csbsav(mg,k)*fluxice(mg)/delt !LDR
              bf=1+0.5*Csb*delt*(1+gam(mg,k))
              dqs=max(0.,delt*(Csb/bf)*(qsg(mg,k)-qtg(mg,k)))
              dqs=min(dqs,(qsg(mg,k)-qtg(mg,k))/(1+gam(mg,k))) !Don't supersat.

              fsclr=(1-cifr(mg,k)-clfr(mg,k))*fluxice(mg)
              sublflux=min(dqs*rhoa(mg,k)*dz(mg,k),fsclr)
              fluxice(mg)=fluxice(mg)-sublflux
              fsclr=fsclr-sublflux
              dqs=sublflux/(rhoa(mg,k)*dz(mg,k))
              qsubl(mg,k)=qsubl(mg,k)+dqs
              qtg(mg,k)=qtg(mg,k)+dqs
              dttg=-hlscp*dqs
              ttg(mg,k)=ttg(mg,k)+dttg
              qsg(mg,k)=qsg(mg,k) + gam(mg,k)*dttg/hlscp
            endif

c Save this for the wet deposition scheme.

            pfstay(mg,nlp-k)=fluxice(mg)*(1-fthru(mg,k))/tdt !Flux staying in layer k

c Accretion of cloud water by falling ice
c This calculation uses the incoming fluxice without subtracting sublimation
c (since subl occurs only outside cloud), so add sublflux back to fluxice.
            
            if(fluxice(mg)+sublflux.gt.0.and.qlg(mg,k).gt.0.)then
              ql=qlg(mg,k)
              cdt=Eac*slopes(mg,k)*(fluxice(mg)+sublflux)/(2*rhos)
              dqf=min(ql,cifra(mg)*ql,ql*cdt/(1+0.5*cdt))

              clfr(mg,k)=clfr(mg,k)*(1.-dqf/qlg(mg,k))
              caccr=clfr(mg,k)*dqf/qlg(mg,k)
              qlg(mg,k)=qlg(mg,k)-dqf
              qaccr(mg,k)=qaccr(mg,k)+dqf
              fluxice(mg)=fluxice(mg)+rhoa(mg,k)*dz(mg,k)*dqf
              dttg=hlfcp*dqf
              ttg(mg,k)=ttg(mg,k)+dttg
              qsg(mg,k)=qsg(mg,k) + gam(mg,k)*dttg/hlscp

            endif


            if(fsclr.eq.0.)then
              cifra(mg)=max(0.01,cifr(mg,k)+caccr)
            else
              ci=cifr(mg,k)+caccr
              cifra(mg)=max(0.01,ci+cifra(mg)-ci*cifra(mg))
            endif


c Compute fluxes into the box
            
            if(fluxice(mg).gt.0.)then
              rhoiin=fluxice(mg)/dz(mg,k)
              cffluxin=min(1.,rhoiin/rica(mg))
            else
              rhoiin=0.
              cffluxin=0.
            endif


c Compute the fluxes of ice and cloud amount leaving the box
            
            if(cifr(mg,k).gt.0.)then

c Use the flux-divergent form as in Rotstayn (QJRMS, 1997)
              rhoiout=rhoi(mg,k)*fout(mg,k)
              cffluxout=cifr(mg,k)*fout(mg,k)

c Or, use a Kessler type formulation with a fixed qcrit option and 
c     accretion of ice
C***              riic=rhoi(mg,k)/cifr(mg,k)
C***c              rcrit=1.e-6*rhoa(mg,k) !i.e. qcrit=const
C***c or use qcrit = f(T) following Sinha & Shine, J. Climate (1994)
C***              rcrit=0.0007e-3*exp(0.0987*(ttg(mg,k)-213.15))
C***              cdt=delt*alpha1(mg,k)
C***c              dra=min(rhoi(mg,k), cifr(mg,k)*dim(riic,rcrit)*cdt) !Explicit
C***              dra=cifr(mg,k)*dim(riic,rcrit)*min(1., cdt/(1+0.5*cdt)) !T-centred
C***c              dra=cifr(mg,k)*dim(riic,rcrit)*(1.-exp(-cdt))      !Analytical
C***c              qaggi(mg,k)=qaggi(mg,k)+dra   !diagnostic only
C***              rhoi2=rhoi(mg,k)-dra !Intermediate rhoi after aggregation
C***
C***c Accretion of ice by falling snow
C***c Exclude accreted LW (graupel) from flxi, but include 50% of accreted ice.
C***
C***              Esi=1.e3*alpha1(mg,k) !Efficiency for accretion of ice by snow
C***              flxi=fluxice(mg)+sublflux+rhoa(mg,k)*dz(mg,k)*(.5*dra-dqf)
C***              cdt=Esi*slopes(mg,k)*flxi/(2*rhos)
C***c              drs=min(rhoi2,rhoi2*cdt) !Explicit
C***              drs=rhoi2*cdt/(1+0.5*cdt) !Time-centered
C***c              drs=rhoi2*(1.-exp(-cdt))   !analytical
C***c              drs=0. !!!!!
C***c              qacci(mg,k)=qacci(mg,k)+drs !diagnostic only
C***
C***              rhoiout=min(rhoi(mg,k), dra+drs)
C***              cffluxout=cifr(mg,k)*max(0.,rhoiout/rhoi(mg,k))

c End of Kessler-type formulation 

              rica(mg)=rhoi(mg,k)/cifr(mg,k) !in cloud rhoi above
            else !Keep value of rica from above
              rhoiout=0.
              cffluxout=0.
            endif
            
c Update the rhoi and cifr fields

            cifr(mg,k)=min(1.-clfr(mg,k),
     &                (cifr(mg,k)-cffluxout) + cffluxin*(1-fthru(mg,k)))
            rhoi(mg,k)=(rhoi(mg,k)-rhoiout) + rhoiin*(1-fthru(mg,k))
            fluxice(mg)=rhoiout*dz(mg,k)+fluxice(mg)*fthru(mg,k) 
c Now fluxice is flux leaving layer k
            fluxi(mg,k)=fluxi(mg,k)+fluxice(mg)


          enddo
        enddo

      enddo

c End of small timestep loop

c Re-create qfg field

      do k=1,nl
        do mg=1,ln2
          qfg(mg,k)=rhoi(mg,k)/rhoa(mg,k)
        enddo
      enddo

c Diagnostics for debugging

      if(debug)then
        if(lg.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, after icefall.'
          write(25,91)'cfrac ',(cfrac(mg,k),k=1,nl)
          write(25,91)'cifr ',(cifr(mg,k),k=1,nl)
          write(25,91)'clfr ',(clfr(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,9)'qsubl ',(qsubl(mg,k),k=1,nl)
          write(25,9)'vi2 ',(vi2(mg,k),k=1,nl)
          write(25,9)'dz ',(dz(mg,k),k=1,nl)
          write(25,9)'rhoa ',(rhoa(mg,k),k=1,nl)
          write(25,9)'fluxi ',(fluxi(mg,k),k=1,nl)
c          write(25,9)'qfdiv ',(qfdiv(mg,k),k=1,nl)
          write(25,*)

C***          ma=(ins-1)*lon
C***          do k=1,nl
C***            do mg=1+ma,lon+ma
C***              mb=mg-ma
C***              if(qfg(mg,k).gt.0.1)then
C***                print*,'end icefall qfg = ',qfg(mg,k)
C***                print*,'lg,ins,mg,k ',lg,ins,mb,k
C***              endif
C***              
C***              if(cfrac(mg,k).gt.0)then
C***                if(qlg(mg,k).le.0.and.qfg(mg,k).le.0)then
C***                  print*,'end icefall cloud with zero/neg water:cfrac=',
C***     &                 cfrac(mg,k)
C***                  print*,'cifr ,clfr ',cifr(mg,k),clfr(mg,k)
C***                  print*,'qfg, qlg=', qfg(mg,k),qlg(mg,k)
C***                  print*,'lg,ins,mg,k ',lg,ins,mb,k
C***                  print*
C***                endif
C***              else
C***                if(qlg(mg,k).gt.0.or.qfg(mg,k).gt.0)then
C***                  print*,'end icefall, cloud water/no cloud: qfg,qlg=',
C***     &                 qfg(mg,k),qlg(mg,k)
C***                  print*,'cifr ,clfr ',cifr(mg,k),clfr(mg,k)
C***                  print*,'ins,lg,k,mg ',ins,lg,k,mb
C***                  print*
C***                endif
C***              endif
C***            enddo
C***          enddo
        endif
      endif

 91   format(a,30f10.3)
 9    format(a,30g10.3)

      return

      end
