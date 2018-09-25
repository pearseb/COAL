c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Moved RCRITL and RCRITS to /CLOUDPAR/, so that they can be read from the
c namelist input. The default values are set in READNML1 to 0.75 and 0.85
c respectively.
c SJP 2004/09/24
c
c Value of rcrits increased from 0.85 to 0.9 in the case of 18 vertical levels,
c as this value has been diagnosed, in conjunction with modifications to the
c values of refac1 and refac2 in cloud2.f, to give surface energy balance.
c SJP 2003/05/24
c
c $Log: newcloud.f,v $
c Revision 1.52  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.51  2001/10/12 02:06:57  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.50.1.1  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.50  2001/02/22 05:56:37  rot032
c Changes from LDR to sulfur cycle, to bring it to run N52
c
c Revision 1.49  1999/06/30 05:29:37  rot032
c Framework for new autoconv treatment so it can be applied easily.
c
c Revision 1.48  1999/05/20 06:23:54  rot032
c HBG changes to V5-2
c
c Revision 1.47  1998/12/10  00:55:46  ldr
c HBG changes to V5-1-21
c
c Revision 1.46  1998/05/27  02:01:02  ldr
c Fix up previous version.
c
c Revision 1.45  1998/05/26  06:08:06  ldr
c New mixed-phase cloud scheme from LDR (use plates option).
c
c Revision 1.44  1997/10/03  05:51:08  ldr
c Merge of LDR sulfate stuff with HBG/MRD NEC stuff.
c
c Revision 1.43  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.42.1.1  1997/10/03  05:45:48  ldr
c Changes for sulfates from LDR
c
c Revision 1.42  1997/08/21  02:52:50  ldr
c Vertically sub-grid cloud for 18L version too.
c
c Revision 1.41  1997/07/24  06:02:14  ldr
c Merge HBG and LDR changes.
c
c Revision 1.40  1997/07/24  05:42:51  ldr
c Tuning of ukconv version and T63 fixes from HBG
c
c Revision 1.39.1.1  1997/07/24  05:59:17  ldr
c Mods to re-tune qcloud for 9L model from LDR.
c
c Revision 1.39  1997/06/19  01:56:13  ldr
c Changes from HBG to ukconv version to improve SWCF, LWCF in tropics.
c
c Revision 1.37  1997/02/21  00:26:10  ldr
c Go back to rcrits/l=0.85/0.75 for 18L version, tidy cloud2.f and make
c clddia.f, cloud.f general for NL.
c
c Revision 1.36  1996/11/22  03:41:41  ldr
c Comments and tidy-ups from LDR; results slightly changed by extra
c calculation of qsg in progcld.
c
c*****************************************************************************
c
c This routine is part of the prognostic cloud scheme. It calculates the
c formation and dissipation of cloud, and the liquid fraction in mixed-phase
c clouds. It is called by progcld.
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
c      tdt - leapfrog timestep (seconds)
c      lg - latitude index (ranges from 1 at poles to LAT at equator)
c      
c      land - logical variable for surface type ( = T for land points)
c      prf - pressure at full levels (in hPa. NB: not SI units)
c      kbase - k index of convective cloud base 
c      ktop - k index of convective cloud top
c
c In/Out:
c
c from arguments
c      ttg - temperature (K)
c      qtg - water vapour mixing ratio (kg/kg)
c      qlg - cloud liquid water mixing ratio (kg/kg)
c      qfg - cloud ice mixing ratio (kg/kg)
c
c Output:
c
c from arguments
c      cfrac - cloudy fraction of grid box
c      ccov - cloud cover looking from above (currently = cloud fraction)
c 
c******************************************************************************

      subroutine newcloud(tdt,lg,land,prf,kbase,ktop,rhoa,cdso4,  !Inputs
     &                    ttg,qtg,qlg,qfg,    !In and out
     &                    cfrac,ccov,cfa,qca)         !Outputs

c This routine is part of the prognostic cloud water scheme

C Global parameters
      include 'PARAMS.f'     !Input model grid dimensions
      include 'PHYSPARAMS.f' !Input physical constants
      include 'CPARAMS.f'    !Input cloud scheme parameters

C Argument list
      real tdt
      integer lg
      logical land(ln2)
      real prf(ln2,nl)
      integer kbase(ln2)
      integer ktop(ln2)
      real rhoa(ln2,nl)
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real cfrac(ln2,nl)
      real ccov(ln2,nl)
      real cdso4(ln2,nl)
      real cfa(ln2,nl)
      real qca(ln2,nl)


C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CLOUDPAR.f'
      include 'CNSTA.f'
      include 'FEWFLAGS.f'   !Input debug, lgdebug etc.
c     include 'TIMEX.f'

C Local work arrays and variables
      real qcg(ln2,nl)
      real qtot(ln2,nl),tliq(ln2,nl),qsg(ln2,nl)
      real fice(ln2,nl)
      real rcritb(ln2)
      real rcrcvb(ln2)
      real Cdrop(ln2,nl)

C Local data, functions etc
      real esdiff(-40:0)  !Table of es(liq) - es(ice) (MKS), -40 <= t <= 0 C
      data esdiff / 
     & 6.22, 6.76, 7.32, 7.92, 8.56, 9.23, 9.94,10.68,11.46,12.27,
     & 13.11,13.99,14.89,15.82,16.76,17.73,18.70,19.68,20.65,21.61,
     & 22.55,23.45,24.30,25.08,25.78,26.38,26.86,27.18,27.33,27.27,
     & 26.96,26.38,25.47,24.20,22.51,20.34,17.64,14.34,10.37, 5.65,
     & 0.08 /
      include 'ESTABL.f'  !Contains arithmetic statement function qsat(p,T)

C Start code : ----------------------------------------------------------

c Set up resolution dependent parameters

CSJP      if(nl.lt.24)then
CSJP        rcritl=0.75
CSJP        rcrits=0.9
CSJP        rcrits=0.85
CSJP      else
CSJP        rcritl=0.75
CSJP        rcrits=0.9
CSJP      endif
        
      if(debug)then
        if(lg.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, before newcloud'
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,9)'qtg ',(qtg(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,*)
        endif
      endif

c Define Cdrop

      if(naerosol_i(2).gt.0)then
        do k=1,nl
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
        do k=2,nl
          do mg=1,ln2
            Cdrop(mg,k)=Cdrop(mg,1)
          enddo
        enddo
      endif

c First melt cloud ice or freeze cloud water to give correct ice fraction fice.
c Then calculate the cloud conserved variables qtot and tliq.
c Note that qcg is the total cloud water (liquid+frozen)



      do k=1,nl
        do mg=1,ln2
          if(ttg(mg,k).ge.273.15)then
            fice(mg,k)=0.
          elseif(ttg(mg,k).ge.tice)then
C***            tc=ttg(mg,k)-273.15
C***            fice(mg,k) = 1. - (max (0., 0.0664*tc + 0.621)) !Bower et al, QJ, 1996
            if(qfg(mg,k).gt.1.0e-12)then
              fice(mg,k)=min(qfg(mg,k)/(qfg(mg,k)+qlg(mg,k)), 1.)
            else
              fice(mg,k)=0.
            endif
          else
            fice(mg,k)=1.
          endif
          qcg(mg,k)=qlg(mg,k)+qfg(mg,k)
          qfnew=fice(mg,k)*qcg(mg,k)
          ttg(mg,k)=ttg(mg,k)+hlfcp*(qfnew-qfg(mg,k)) !Release L.H. of fusion
          qfg(mg,k)=fice(mg,k)*qcg(mg,k)
          qlg(mg,k)=max(0.,qcg(mg,k)-qfg(mg,k))
        enddo
      enddo

      if(ukconv)then
        cvsrcr=0.2
        cvlrcr=0.1
      else
        cvsrcr=0.0
        cvlrcr=0.0
      endif

c---- Set up latitude row of rcrit values :
      do mg=1,ln2
          if(land(mg))then
            rcritb(mg)=rcritl
            rcrcvb(mg)=cvlrcr
          else
            rcritb(mg)=rcrits
            rcrcvb(mg)=cvsrcr
          endif
      enddo

      do k=1,nl
        do mg=1,ln2
          fi=0.
          if(ttg(mg,k).le.Tice)fi=1.
c          fi=fice(mg,k)
          hlrvap=(hl+fi*hlf)/rvap
          qtot(mg,k)=qtg(mg,k)+qcg(mg,k)
          tliq(mg,k)=ttg(mg,k)-hlcp*qcg(mg,k)-hlfcp*qfg(mg,k)

c Calculate qs and gam=(L/cp)*dqsdt,  at temperature tliq
          pk=100.0*prf(mg,k)
          qsg(mg,k)=qsati(pk,tliq(mg,k)) !Ice value
          deles=esdiff(min(max(-40,(nint(tliq(mg,k)-273.15))),0))
          qsl=qsg(mg,k)+epsil*deles/pk !qs over liquid
          qs=qsg(mg,k)
          if(ttg(mg,k).lt.273.15.and.ttg(mg,k).gt.Tice)qs=qsl
          dqsdt=qs*hlrvap/tliq(mg,k)**2

c Calculate cloudy fraction of grid box (cfrac) and gridbox-mean cloud water
c using the triangular PDF of Smith (1990)

          rcrit=max(rcritb(mg),min(0.99,sig(k)))
          if((k.ge.kbase(mg)).and.(k.le.ktop(mg)))then
            rcrit=rcrit-rcrcvb(mg)
          endif

          al=1/(1+(hlcp+fi*hlfcp)*dqsdt)    !Smith's notation
          qc=qtot(mg,k)-qs

          delq=(1-rcrit)*qs      !UKMO style (equivalent to above)
          cfrac(mg,k)=1.
          qcg(mg,k)=al*qc
          if(qc.lt.delq)then
            cfrac(mg,k)=1-0.5*((qc-delq)/delq)**2
            qcg(mg,k)=al*(qc-(qc-delq)**3/(6*delq**2))
          endif
          if(qc.le.0.)then
            cfrac(mg,k)=0.5*((qc+delq)/delq)**2
            qcg(mg,k)=al*(qc+delq)**3/(6*delq**2)
          endif
          if(qc.le.-delq)then
            cfrac(mg,k)=0.
            qcg(mg,k)=0.
          endif

c Calculate the cloud fraction (cfa) in which ql exceeds qcrit, and
c the corresponding gridbox-mean cloud water mixing ratio qca. 
c This (qca) is the cloud-water mixing ratio inside cfa divided by cfa.
c The new variable qc2 is like qc above, but is used for integration limits
c only, not the integrand

C***          qcrit=(4*pi/3)*rhow*Rcm**3*Cdrop(mg,k)/rhoa(mg,k)
C***          qc2=qtot(mg,k)-qs-qcrit/al 
C***          cfa(mg,k)=1.
C***          qca(mg,k)=al*qc
C***          if(qc2.lt.delq)then
C***            cfa(mg,k)=1-0.5*((qc2-delq)/delq)**2
C***            qto=(qtot(mg,k)-delq+2*(qs+qcrit/al))/3.
C***            qca(mg,k)=al*(qtot(mg,k) - qto + cfa(mg,k)*(qto-qs))
C***          endif
C***          if(qc2.le.0.)then
C***            cfa(mg,k)=0.5*((qc2+delq)/delq)**2
C***            qca(mg,k)=cfa(mg,k)*(al/3.)*(2*qcrit/al + qc+delq)
C***          endif
C***          if(qc2.le.-delq)then
C***            cfa(mg,k)=0.
C***            qca(mg,k)=0.
C***          endif

        enddo
      enddo

c Use prescribed ice fraction
c Bower et al...
        
C***        do k=1,nl
C***          do mg=1,ln2
C***            qfg(mg,k)=qcg(mg,k)*fice(mg,k)
C***            qlg(mg,k)=qcg(mg,k)-qfg(mg,k)
C***          enddo
C***        enddo

c OR...
c Condense or evaporate ql first.

      do k=1,nl
        do mg=1,ln2
          if(ttg(mg,k).gt.Tice)then
            qlg(mg,k) = max (0., qcg(mg,k)-qfg(mg,k))
            qfg(mg,k) = qcg(mg,k) - qlg(mg,k)
          else
            qlg(mg,k) = 0.
            qfg(mg,k) = qcg(mg,k)
          endif
        enddo
      enddo


c Do the vapour deposition calculation in the liquid part of mixed-phase clouds

c Calculate deposition on cloud ice

      do k=1,nl-1
        do mg=1,ln2
          if(cfrac(mg,k).gt.0.)then
            Tk=tliq(mg,k)+hlcp*(qlg(mg,k)+qfg(mg,k))/cfrac(mg,k) !T in liq cloud
c            fl=qlg(mg,k)/(qfg(mg,k)+qlg(mg,k))
            if(Tk.lt.273.15.and.Tk.ge.Tice.and.qlg(mg,k).gt.0.)then
              pk=100*prf(mg,k)
              es=100*exp(23.33086-6111.72784/Tk+0.15215*alog(Tk)) !ice value
              qs=0.622*es/pk !ice value
c              qs=qsati(pk,Tk)
c              es=qs*pk/0.622 !ice value
              Aprpr=hl/(rKa*Tk)*(hls/(rvap*Tk)-1)
              Bprpr=rvap*Tk/((Dva/pk)*es)
c              Cice=1.0e-2*exp(0.6*(273.15-ttg(mg,k))) !Fletcher 1962
c              deles=esdiff (nint(Tk-273.15))
              esl=100*exp(53.67957-6743.769/Tk-4.8451*alog(Tk))
              deles=esl-es
              Cice=1.e3*exp(12.96*deles/es - 0.639) !Meyers et al 1992
c              Cice=Cice/100.  !Sensitivity experiment
c              write(30,'(f12.3,g12.3)')Tk-273.15,cice

              cm0=1.e-12 !Initial crystal mass
              qi0=cm0*Cice/rhoa(mg,k) !Initial ice mixing ratio

c Next 2 lines are for assumption of fully mixed ql and qf (also a line further down).

              qi0=max(qi0, qfg(mg,k)/cfrac(mg,k)) !Assume all qf and ql are mixed
              fd=1.   !Fraction of cloud in which deposition occurs

c              fd=fl   !Or, use option of adjacent ql,qf


C***              alf=1./2
C***              Crate=65.2*sqrt(Cice/rhoa(mg,k))*deles/((Aprpr+Bprpr)*es) !plates
C***              qfdep=fd*cfrac(mg,k)*0.25*(Crate*tdt+2*sqrt(qi0))**2  !Plates

C***              alf=0.339
C***              Crate=0.839*(Cice/rhoa(mg,k))**0.661
C***     &             *deles/((Aprpr+Bprpr)*es)                 !Columns

C***              alf=1./3
              rhoic=700.
              Crate=7.8*((Cice/rhoa(mg,k))**2/rhoic)**(1./3) !Spheres
     &              *deles/((Aprpr+Bprpr)*es)
              qfdep=fd*cfrac(mg,k)*sqrt                         !Spheres
     &             (((2./3)*Crate*tdt+qi0**(2./3))**3)
C***c              qrate=qfdep/fd/cfrac(mg,k)/tdt

C***              qfdep=fd*cfrac(mg,k)*                             !General
C***     &              ( (1-alf)*Crate*tdt + qi0**(1-alf) )**(1./(1-alf))

c Also need this line for fully-mixed option...
              qfdep = qfdep - qfg(mg,k)

              qfdep=min(qfdep,qlg(mg,k))
              qlg(mg,k)=qlg(mg,k)-qfdep
              qfg(mg,k)=qfg(mg,k)+qfdep
            endif

            fice(mg,k)=qfg(mg,k)/(qfg(mg,k)+qlg(mg,k))

          endif
        enddo
      enddo  

c Recalculate qcg and cfrac for mixed-phase clouds that are fully glaciated
c Use qs=qsg.

      do k=1,nl
        do mg=1,ln2
          if(ttg(mg,k).lt.273.15.and.ttg(mg,k).ge.Tice.and.
     &         fice(mg,k).gt.0.999999)then
            qs=qsg(mg,k)
            fi=1.
            hlrvap=(hl+fi*hlf)/rvap

            dqsdt=qs*hlrvap/tliq(mg,k)**2

c Calculate cloudy fraction of grid box (cfrac) and gridbox-mean cloud water
c using the triangular PDF of Smith (1990)

            rcrit=max(rcritb(mg),min(0.99,sig(k)))
            if((k.ge.kbase(mg)).and.(k.le.ktop(mg)))then
              rcrit=rcrit-rcrcvb(mg)
            endif

            al=1/(1+(hlcp+fi*hlfcp)*dqsdt) !Smith's notation
            qc=qtot(mg,k)-qs

            delq=(1-rcrit)*qs   !UKMO style (equivalent to above)
            cfrac(mg,k)=1.
            qcg(mg,k)=al*qc
            if(qc.lt.delq)then
              cfrac(mg,k)=1-0.5*((qc-delq)/delq)**2
              qcg(mg,k)=al*(qc-(qc-delq)**3/(6*delq**2))
            endif
            if(qc.le.0.)then
              cfrac(mg,k)=0.5*((qc+delq)/delq)**2
              qcg(mg,k)=al*(qc+delq)**3/(6*delq**2)
            endif
            if(qc.le.-delq)then
              cfrac(mg,k)=0.
              qcg(mg,k)=0.
            endif

            qfg(mg,k)=qcg(mg,k)
            qlg(mg,k)=0.
          endif

        enddo
      enddo


c Calculate new values of vapour mixing ratio and temperature

      do k=1,nl
        do mg=1,ln2
          qtg(mg,k)=qtot(mg,k)-qcg(mg,k)
          ttg(mg,k)=tliq(mg,k)+hlcp*qcg(mg,k)+hlfcp*qfg(mg,k)
          pk=100.0*prf(mg,k)
          qsg(mg,k)=qsati(pk,ttg(mg,k)) !Ice value
          ccov(mg,k)=cfrac(mg,k)      !Do this for now
        enddo
      enddo

c Vertically sub-grid cloud

      if(nl.lt.18)then
        do k=1,nl-1
          do mg=1,ln2
            if(cfrac(mg,k).gt.1.0e-2)then
              ccov(mg,k)=cfrac(mg,k)**(2./3)
            endif
          enddo
        enddo
      else
        do k=2,nl-1
          do mg=1,ln2
            if(cfrac(mg,k).gt.1.0e-2
     &           .and.cfrac(mg,k+1).eq.0..and.cfrac(mg,k-1).eq.0.)then
c            ccov(mg,k)=cfrac(mg,k)**(2./3)
              ccov(mg,k)=sqrt(cfrac(mg,k))
            endif
          enddo
        enddo
      endif


      
      if(debug)then
        if(lg.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, after newcloud'
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,9)'qtg ',(qtg(mg,k),k=1,nl)
          write(25,9)'qsg ',(qsg(mg,k),k=1,nl)
          write(25,9)'cfrac ',(cfrac(mg,k),k=1,nl)
          write(25,9)'cdso4 ',(cdso4(mg,k),k=1,nl)
          write(25,9)'cfa ',(cfa(mg,k),k=1,nl)
          write(25,9)'qca ',(qca(mg,k),k=1,nl)
          write(25,9)'qtot ',(qtot(mg,k),k=1,nl)
          write(25,9)'qcg ',(qcg(mg,k),k=1,nl)
          write(25,9)'qlg ',(qlg(mg,k),k=1,nl)
          write(25,9)'qfg ',(qfg(mg,k),k=1,nl)
          write(25,*)
        endif
      endif
 91   format(a,30f10.3)
 9    format(a,30g10.3)

      return
      end
