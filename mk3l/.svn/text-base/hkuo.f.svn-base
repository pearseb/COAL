c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c Revision 1.22  1999/05/20 06:23:58  rot032
c HBG changes to V5-2
c
c Revision 1.21  1998/12/10  00:56:02  ldr
c HBG changes to V5-1-21
c
c Revision 1.20  1997/12/17  23:23:07  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.19  1996/10/24  01:02:51  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.18  1996/03/21  03:18:46  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.17  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.16  1995/06/30  02:44:41  ldr
c Changes to qcloud (run F88) and put qsat formula into ESTABL.f
c
c Revision 1.15  1995/02/27  01:40:34  ldr
c Bugfixes received from HBG in response to problems on Fujitsu.
c
c Revision 1.14  1994/08/08  17:21:26  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.13  94/03/22  15:30:54  ldr
c Double precision needed for Seca.
c 
c Revision 1.12  93/12/17  15:32:43  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.11  93/09/06  12:28:13  ldr
c Use correct formula for qs and dqsdt (i.e. don't subtract es from p).
c 
c Revision 1.10  93/08/19  15:08:13  ldr
c Minor cosmetic changes.
c 
c Revision 1.9  93/07/06  16:09:01  ldr
c      New map for ice depth averaged over leads area. This requires extra array
c      in common block cloudm1. The mapping in surfset.f reqires pl(mg).
c 
c Revision 1.8  92/12/19  19:17:26  ldr
c Changes to hkuo to enable new diagnostic cloud scheme as in conv (hbg).
c 
c Revision 1.7  92/12/09  14:43:30  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.6  92/06/16  12:02:40  ldr
c Rewritten to pass physical variables as arguments rather than in common,
c in order to implement sea-ice model.
c
c     INPUT/OUTPUT
c     Input:   from common/fewflags in FEWFLAGS.f
c                  hybrid - if true, use new hybrid vertical coordinate
c
c              from common/hybrpr in HYBRPR.f
c                  dprf - pressure thickness 
c                  prf  - pressure at full levels
c
c              from common/timex in TIMEX.f
c                  mstep - time step in minutes
c
c              from arguments
c                  dqgdt - mixing ratio at current time - dt
c                  ipass - model pass counter   tdt - 2 X timestep (dt)
c
c     Output:  from common/cloud1 in this subroutine
c                  cnmlp - convection map   rainp - rain map 
c
c              from arguments
c                  rhg - relative humidity
c
c     In/Out:  from common/kuocom in this subroutine
c                  kbconv - level from which moisture convergence is
c                           calculated. Set to min( kbconv, kb )
c
c              from common/newcl in this subroutine
c                  avcvrn - average convective rain, ice present
c                  icln   - convective cloud base/top indicator, ice present
c                  ipcln  - convective cloud base/top indicator, no ice
c                  pavcvrn- average convective rain, no ice
c
c              from common/timex in TIMEX.f
c                  mins - current model time in mins
c
c              from arguments
c                  precc - convective precipitation
c                  qtg   - mixing ratio      ttg - physical temperature
c 
c 
      subroutine hkuo(ipass,lg,tdt,pg,dqgdt,
     &                ttg,qtg,precs,precc,rhg,
     &            fmp,kta,kba,fluxc,rainx,qevc)

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer ipass,lg
      real tdt
      real pg(ln2)
      real dqgdt(ln2,nl)
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real precs(ln2)
      real precc(ln2)
      real rhg(ln2,nl)
      real fmp(ln2,nl)
      integer kta(ln2)
      integer kba(ln2)
      real fluxc(ln2,nl)
      real rainx(ln2)
      real qevc(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'

C Global data blocks
      include 'FEWFLAGS.f'
      include 'PRINTT.f'
      include 'TIMEX.f'

      character*1 cnmlp(ln2,lat),rainp(ln2,lat)
     & ,snmap(ln2,lat),snwice(ln2,lat),snicel(ln2,lat)
      common/cloud1/cnmlp,rainp,snmap,snwice,snicel

      integer kuocb,kbconv
      real rhcut,rhkuo,rhsat
      common/kuocom/kuocb,kbconv,rhcut,rhkuo,rhsat

      integer jclb,nlow,nmid,ktied,kbgel,klcmc,klowrn,kmidrn
      real aftsea
      common/levdata/jclb,nlow,nmid,aftsea,ktied,kbgel
     &,klcmc,klowrn,kmidrn

      integer icln,ipcln
      real avcvrn,pavcvrn
      common/newcl/icln(ln2,2,lat),avcvrn(ln2,lat)
     &           ,ipcln(ln2,2,lat),pavcvrn(ln2,lat)

C Local work arrays and variables
      real h(nl),p(nl),phi(nl),q(nl),qc(nl),qold(nl),qs(nl)
      real sdq(nl)
      real bet(ln2,nl)
      real dsig(ln2,nl) ! d(p)/P* = d(sigma) if not hybrid (i.e. dsk)

      double precision tt(nl),tc(nl),tint

      integer irain(ln2,nl)
c     qls, tls, qcn, tcn are only used for optional diagnostic printouts
c     real qls(nl),tls(nl),qcn(nl),tcn(nl)
      character*1 chcon

      integer irn
      integer iter
      integer jh
      integer k
      integer kb
      integer kbcon
      integer kt
      integer mg

      real b
      real c
      real convi
      real delq
      real delt
      real dqeff
      real dqsdt
      real dspart
      real hlrvap
      real qconst
      real qint
      real qnew
      real rainc
      real rains
      real rhmean
      real rhpart
      real sumdq
      real tconst
      real xrn

C Local data, functions etc
      character*1 cont(8)
      data cont/' ','L','M','4','H','5','6','7'/
      include 'ESTABL.f'  !Contains statement function form of establ

C Start code : ----------------------------------------------------------

c     Now in terms of moist static energy 30/8/90 from nested kuohc
C     and Q weighting functions based on QC not QS
c     the variable kbcon is the level from which moisture convergence
c     calculated. Set to min( kbconv, kb )
C     A VERSION OF KUO CONVECTION set up for GCM model 4/9 JMcG  27/3/90
C     It should be called as the final physics routine affecting QTG.
c     This version has set KUOCB=.8, and uses Anthes-type B

c
c     TDT is 2*timestep size (secs)
c     TTG and QTG are current T and mixing ratio (kg/kg).
c     PG is surface pressure (in mbs)
c     dQGdt is mixing ratio tendency
c     PRECS() is array for accumulated large scale precip
c     PRECC() is array for accumulated convective precip
c     SIG() contains sigma values
c     DSK() contains sigma thicknesses

c
c     set various parameters used in kuo scheme; also physical constants
c
      hlrvap=hl/rvap
      c=grav/6.5e-3
      do mg=1,ln2
        kba(mg)=0
        kta(mg)=0
        bet(mg,1)=c *((prf(mg,1)/pg(mg))**(-rdry/c)-1)
      enddo
      do 12 k=1,nl-1
      do 12 mg=1,ln2
   12   bet(mg,k+1)=rdry*log(prf(mg,k)/prf(mg,k+1))*.5

      do k=1,nl
        do mg=1,ln2
          fluxc(mg,k)=0.
          qevc(mg,k)=0.
          fmp(mg,k)=0.0
          dsig(mg,k)=dprf(mg,k)/pg(mg)
          irain(mg,k)=0
        enddo
      enddo


c  reset convective cloud stats each day for pure restarts

      if((mod(mins,1440_8).eq.0_8).and.(ipass.eq.1))then
        do 22 mg=1,ln2
        icln(mg,1,lg)=0
        icln(mg,2,lg)=0
   22   avcvrn(mg,lg)=0.0
        do 24 mg=1,ln2
        ipcln(mg,1,lg)=0
        ipcln(mg,2,lg)=0
   24   pavcvrn(mg,lg)=0.0
      end if
C
C     LOOP OVER ALL POINTS
C
      do 88 mg=1,ln2
c
c     set some input arrays for the kuo scheme; sigma level pressures
      do 14 k=1,nl
      p(k)=100.0*prf(mg,k)
c     qls(k)=0.
c     tls(k)=0.
c     qcn(k)=0.
c     tcn(k)=0.
      qold(k)=qtg(mg,k)-tdt*dqgdt(mg,k)
      tc(k)=ttg(mg,k)
      qc(k)=qtg(mg,k)
      tt(k)=ttg(mg,k)
   14 q(k)=qtg(mg,k)
c
      rainc=0.
c
C     CALCULATE  HEIGHT AND SAVE IN phi
      phi(1)=bet(mg,1)*tt(1)
      do 15 k=2,nl
   15 phi(k)=phi(k-1)+bet(mg,k)*(tt(k)+tt(k-1))
c
c     calculate saturation values of specific humidity
      do 17 k=1,nl
      qs(k)=qsat(p(k),tt(k))
      qs(k)=max(qs(k),0.)
c
c     calculate the moist static energy
   17 h(k)=cp*tt(k)+phi(k)+hl*qs(k)
      kb=kuocb
      kt=0
      iter=0
      go to 21
c
    2 kb=kb+1
c
c     skip adjustment if cloud base too high
c
   21 if(kb.gt.nl-2)go to 7
      iter=iter+1
      do 3 k=kb,nl-2
c
c     check base relative humidity
c
      if(q(k).gt.rhkuo*qs(k))go to 35
    3 continue
c
c     if the cloud base sp. humidity does not exceed a certain
c     value, also skip the kuo adjustment
c
      go to 7
c----------------------------------------------------------------------
c
c     kuo adjustment procedure
c
c----------------------------------------------------------------------
c
c     set kb to level where sp. humidity first exceeds kuo criterion
c
   35 kb=k
c
c     find cloud top level by checking that all layers are unstable
c     and calculate integrated moisture conv. etc. (dsig is +ive)
c
      kbcon=min(kbconv,kb)
      convi=0.
      do 36 k=kbcon,kb
   36 convi=convi+dsig(mg,k)*(q(k)-qold(k))
      do 37 k=kb+1,nl-1
      if(h(kb).lt.h(k))go to 4
c
c     may later add test for int(q.dsig)/int(qs.dsig)
c
   37 convi=convi+dsig(mg,k)*(q(k)-qold(k))
c
c     if moisture is divergent start process again
c
    4 if(convi.lt.0.)go to 2
c
c     successful cloud base has been found; cloud top
c     set to nl-1 or at level where stable layer forms
c
      kt=k-1
c
c     ensure that convective clouds occupy at least 2 layers;
c     otherwise start again
c
      if(kt.eq.kb)go to 2
c
c     set values of humidity and temp.
c
      qc(kb)=qs(kb)
      tc(kb)=tt(kb)
cZZZZ trial generation of convective cloud
c          icln(mg,1,lg)=kb+1
c          icln(mg,2,lg)=kt
c.... save the max kt and min kb over the radiation step (2 Hrs)
c.... for the Slingo 87 cloud formulation.
c.... icln will be reset to zero following the call to CLOUD
      if(ipass.eq.1)then
         if(icln(mg,1,lg).eq.0)icln(mg,1,lg)=kb
         icln(mg,1,lg)=min(icln(mg,1,lg),kb)
         icln(mg,2,lg)=max(icln(mg,2,lg),kt)
       elseif(ipass.eq.2)then
         if(ipcln(mg,1,lg).eq.0)ipcln(mg,1,lg)=kb
         ipcln(mg,1,lg)=min(ipcln(mg,1,lg),kb)
         ipcln(mg,2,lg)=max(ipcln(mg,2,lg),kt)
       endif

c
c     iterate over all cloud levels
c
      do 45 k=kb+1,kt
c
c     calculate tc and qc together
c
      dqsdt=qs(k)*hlrvap/tt(k)**2
      tc(k)=(h(kb)-phi(k)-hl*(qs(k)-dqsdt*tt(k)))/(cp+hl*dqsdt)
   45 qc(k)=(h(kb)-phi(k)-cp*tc(k))/hl
c
         kba(mg)=kb
         kta(mg)=kt
c
c     calculation of b parameter using anthes scheme-------
c
      dspart=0.
      rhpart=0.
      do 50 k=kb,kt
      dspart=dspart+dsig(mg,k)
   50 rhpart=rhpart+dsig(mg,k)*(q(k)/qs(k))
c     currently use b such that b=1. for rh<rhcut (0.5,say); b=0. for rh
c     this assumes rhcut <= 1.; to make b=0., can set rhcut=1.E20 say

      rhmean=min(rhpart/dspart,1.)
      b=min((1.-rhmean)**2/(1.-rhcut)**2,1.)
c
c     calculate intermediate temp. and hum. values for precip.
c     calculation
c
      tint=0.
      qint=0.
      do 56 k=kbcon,kt
      qint=qint+dsig(mg,k)*max(qc(k)-q(k),1.e-7)
   56 tint=tint+dsig(mg,k)*(tc(k)-tt(k))
      tconst=(1.-b)*convi/tint
      qconst=b*convi/qint
c
c     calculate precip. and adjusted humidities
      do 6 k=kbcon,kt
      delq=qconst*max(qc(k)-q(k),1.e-7)
c     qcn(k)=qcn(k)+delq
    6 q(k)=qold(k)+delq
c
c     calculate adjusted temps.
      do 61 k=kb,kt
      dqeff=tconst*(tt(k)-tc(k))
c
c     prec will be in mm in mks model because density=1000.
      rainc=rainc-100.0*dprf(mg,k)*dqeff/grav
      delt=-hlcp*dqeff
c     tcn(k)=tcn(k)+delt
   61 tt(k)=tt(k)+delt
c
c     now recalculate qs for cloud levels
      do 65 k=kb,kt
   65 qs(k)=qsat(p(k),tt(k))
c
c     set new cloud base at next level above previous cloud top
c     and repeat the adjustment process
C***      kb=kt+1
C***      go to 21
c

    7 if(.not.qcloud)then

c      do a supersaturation check; recalculate outputs if necessary
      rains=0.0
      do 75 k=1,nl
      if(q(k).gt.rhsat*qs(k))then
        irain(mg,k)=1
        dqsdt=qs(k)*hlrvap/tt(k)**2
        qnew=rhsat*(qs(k)+hlcp*q(k)*dqsdt)/(1.+hlcp*rhsat*dqsdt)
        delq=qnew-q(k)
        tt(k)=tt(k)-hlcp*delq
c       qls(k)=qls(k)+delq
c       tls(k)=tls(k)-hlcp*delq
        rains=rains-100.0*dprf(mg,k)*delq/grav
        q(k)=qnew
      endif
   75 continue
      precs(mg)=precs(mg)+0.5*rains

      endif ! .not.qcloud
C****
      do 76 k=1,nl
      qs(k)=qsat(p(k),tt(k))
      rhg(mg,k)=q(k)/qs(k)
   76 continue
C
C ACCUMULATE PRECIP. - multiply by .5 because of leapfrog
C  time differencing.
C
      xrn=0.5*rainc
      precc(mg)=precc(mg)+xrn
c Stuff for qcloud (approx!)
      rainx(mg)=xrn*(1440.0/mstep)
      if(kt.gt.kb)then
        sumdq=0.0
        do k=kt,kb,-1
          sumdq=sumdq+max(0.0,qtg(mg,k)-q(k))*prf(mg,k)
          sdq(k)=sumdq
        enddo
        do k=kb,kt
          fluxc(mg,k)=2.0*xrn*sdq(k)/sumdq
        enddo
      endif
C
C TTG and QTG are updated in split fashion
C TO BE PASSED BACK TO THE MAIN ROUTINE
C
      do 81 k=1,nl
          ttg(mg,k)=tt(k)
   81     qtg(mg,k)=q(k)

c....
c.... add convective rainfall rate at cloud base for conv cloud
c.... determination (Slingo,87) (rate each step in units mms/day)
c.... => divide by number of steps when using
      if(ipass.eq.1)then
        avcvrn(mg,lg)=avcvrn(mg,lg) + xrn*(1440.0/mstep)
      elseif(ipass.eq.2)then
        pavcvrn(mg,lg)=pavcvrn(mg,lg) + xrn*(1440.0/mstep)
      endif
      
   88 continue 

      if ((mod(mins+int(mstep, 8),1440_8).eq.0).and.cvrnm) then
c**** update convection map at end of day only
c**** (see initax.f for klcmc)
      do 370 mg=1,ln2
         chcon=' '
         if (kta(mg).eq.0) go to 370
         chcon='L'
         if (kta(mg).le.klcmc) go to 370
         chcon='M'
         if (kba(mg).ge.klcmc) go to 370
         chcon='P'
 370     cnmlp(mg,lg)=chcon

        if (.not.qcloud) then
c**** update the rainfall map at end of day only
c**** Low rainfall (up to sig=0.9), Middle Rain (up to sig=0.45), High above
c**** (see initax.f for klowrn,kmidrn)
         do 350 mg=1,ln2
            irn=1
            jh=0
            do 352 k=kmidrn+1,nl
  352       jh=jh+irain(mg,k)
            if (jh.ge.1) irn=5
            jh=0
            do 354 k=klowrn+1,kmidrn
  354       jh=jh+irain(mg,k)
            if (jh.ge.1) irn=irn+2
            jh=0
            do 356 k=1,klowrn
  356       jh=jh+irain(mg,k)
            if (jh.ge.1) irn=irn+1
  350       rainp(mg,lg)=cont(irn)
        endif
      end if

      return
      end
