c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/10
c
c Changes to add machine type 'ALPH'.
c SJP 2001/11/22
c
c $Log: conv.f,v $
c Revision 1.59  2001/02/22 05:34:40  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.58  2000/11/14 03:11:36  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.57  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.56  1998/12/10  00:55:38  ldr
c HBG changes to V5-1-21
c
c Revision 1.55  1997/12/23  04:10:02  ldr
c Remove tau1, tau2 scaling for evaporation of rain.
c
c Revision 1.54  1997/12/17  23:22:48  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.53  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.52  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.51  1996/10/24  01:27:40  ldr
c Merge of TIE changes with LDR/MRD changes since V5-0-12.
c
c Revision 1.50  1996/08/08  02:40:39  ldr
c Update qcloud to 24P: No anvils, usual vadv and tightened-up conservation.
c
c Revision 1.49.1.1  1996/10/24  01:02:34  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.49  1996/06/25  04:56:58  ldr
c Only add condensation to qfg if qcloud=T, and sneak this into V5-0-12.
c
c Revision 1.48  1996/06/13  02:23:00  ldr
c Merge of TIE and LDR changes.
c
c Revision 1.46.1.1  1996/06/13  02:05:55  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.47  1996/06/13  01:50:52  ldr
c Update qcloud to run 24E (works for 24L and 18L)
c
c Revision 1.46  1996/03/25  04:31:29  ldr
c Replace clhx by hlcp in conv, rainda and radin (needed for V5-0-7).
c
c Revision 1.45  1996/03/21  03:18:32  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.44  1996/02/19  04:09:42  ldr
c Generalize for 24 levels.
c
c Revision 1.43  1996/01/09  06:18:16  ldr
c This is the version used for the long run H06 and write-up.
c
c Revision 1.42  1995/11/30  04:40:16  ldr
c Final updates to qcloud scheme to bring it to long run h63 used for
c BMRC workshop and write-up.
c
c Revision 1.41  1995/11/23  06:03:27  ldr
c Update qcloud to run h32 (used for BMRC & long run)
c
c Revision 1.40  1995/08/31  04:30:40  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.39  1995/08/29  02:29:42  ldr
c Merge of HBG's corrections to V4-7-13h with LDR's changes for run g61.
c
c Revision 1.38  1995/08/29  01:48:50  ldr
c MErge of HBG's corrections to V4-7-13h with LDR changes.
c
c Revision 1.36.1.1  1995/08/29  01:30:09  ldr
c HBG's corrections to his V4-7-13h, to make results with hybrid=F agree
c with previous version.
c
c Revision 1.35.1.2  1995/08/29  01:50:06  ldr
c Update qcloud to run g61.
c
c Revision 1.36  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.35  1995/06/30  02:44:39  ldr
c Changes to qcloud (run F88) and put qsat formula into ESTABL.f
c
c Revision 1.34  1995/05/02  06:10:36  ldr
c Changes to V4-6-16l from LDR - mainly affecting qcloud scheme.
c This version used for run F35.
c
c Revision 1.33  1994/12/15  06:29:45  ldr
c Further mods to LDR's cloud scheme - this version is as used in E58 run
c and as described in seminar of 27/10/94.
c
c Revision 1.32  94/09/12  12:50:04  ldr
c Changes to qcloud scheme of V4-6 to create run E32, except that E32 has
c enhanced collection of cloud water by falling ice.
c 
c Revision 1.31  94/08/09  12:01:51  ldr
c Set RHMOIS back to 0 for qcloud version.
c 
c Revision 1.30  94/08/08  17:20:56  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.29  94/08/08  13:16:11  ldr
c Push the debug prints down into individual routines.
c 
c Revision 1.28  94/08/04  16:54:38  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.27  94/04/29  15:18:10  ldr
c Minor changes for 18L: start deep and shallow conv from k=3, increase Cd.
c 
c Revision 1.26  94/03/30  12:58:28  ldr
c Merge of HBG and LDR changes.
c 
c Revision 1.25  94/03/30  10:15:08  ldr
c Merge of 1.23.1.1 and 1.24 changes.
c 
c Revision 1.24  94/03/22  15:30:32  ldr
c Double precision needed for Seca (18 levels).
c 
c Revision 1.23.1.1  94/03/30  10:14:12  ldr
c Go back to V4-3 heating profile.
c 
c Revision 1.23.2.1  94/03/30  12:34:08  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.23  93/12/17  15:31:57  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.22  93/10/15  14:16:49  ldr
c Changes to create special 'FUJI' version of model with fast FFTs and
c Legendre transform.
c 
c Revision 1.21  93/10/07  12:09:34  ldr
c Move machine parameter statement into new include file MACHINE.f.
c 
c Revision 1.20  93/10/05  13:05:36  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.19  93/09/07  10:26:21  ldr
c Move dry adjustment from rainda to conv.
c 
c Revision 1.18  93/08/17  16:26:23  ldr
c Reduced RHMOIS from 60% to 40% (HBG).
c 
c Revision 1.17  93/07/06  16:08:58  ldr
c      New map for ice depth averaged over leads area. This requires extra array
c      in common block cloudm1. The mapping in surfset.f reqires pl(mg).
c 
c Revision 1.16  93/06/16  14:07:36  ldr
c   In conv.f, if needed following 2nd iteration , heating reset to 0 for
c    cloud base, and heating made >= 0 for all levels above cloud base (HBG).
c 
c Revision 1.15  93/03/12  12:10:32  ldr
c Removed /johnb and put its contents into /fewflags.
c 
c Revision 1.14  93/01/28  13:01:32  ldr
c Modified by HBG to implement ECMWF evap of rain below cloud base but use old
c moistening scheme inside towers.
c 
c Revision 1.13  93/01/26  16:59:26  ldr
c Generalized slightly by Hal to allow return (almost) to old (V4-0) 
c conv scheme by commenting out a couple of lines.
c
c     INPUT/OUTPUT:
c     Input:   from common/cnsta in CNSTA.f
c                  dsk - sigma thicknesses (1 to nl)
c
c              from common/const in CONST.f
c                  cappa - specific gas const dry air/spec heat dry air at 
c                          const pressure
c
c              from common/fewflags in FEWFLAGS.f
c                  debug    - flag to control column debugging
c                  insdebug - flag to control hemisphere debugging
c                  lgdebug  - flag to control latitude debugging
c                  mgdebug  - flag to control longitude debugging
c                  qcloud   - true if using prognostic cloud scheme
c
c
c              from common/hybrpr in HYBRPR.f
c                  dprf - pressure thickness at each sigma level
c                  muf  - pressure derivative with respect to coordinates, at
c                         full levels
c                  prf  - pressure at full levels
c
c              from common/levdata in this subroutine
c                  jclb  - convection cloud base
c                  klcmc - convective mapping indicator
c
c              from common/printt in PRINTT.f
c                  cvrnm - convection/rain mapping indicator
c
c              from common/timex in TIMEX.f
c                  mstep - timestep in minutes
c
c              from arguments
c                  gam - (L/cp)dqsdt
c                  lg  - latitude index        tdt - 2 X timestep (dt)
c                  pg  - surface pressure(mbs) qsg - ice value
c                  ipass - model pass counter
c
c     In/Out:  from common/newcl in this subroutine
c                  avcvrn - average convective rain, ice present
c                  icln   - convective cloud base/top indicator, ice present
c                  ipcln  - convective cloud base/top indicator, no ice
c                  pavcvrn- average convective rain, no ice
c
c              from common/timex in TIMEX.f
c                  mins -  current model time in mins
c
c              from arguments
c                  kba - level of cloud base  kta - level of cloud top
c                  qtg - mixing ratio         rhg - relative humidity
c                  precc - convective precipitation
c                  rainx - rainfall rate at cloud base
c                  ttg - current temperature
c
c     Output:  from common/cloud1 in this subroutine
c                  cnmlp - convection map
c
c              from arguments
c                  fluxc - flux of convective rain 
c                  fmp   - convective mass flux
c                  qevc  - evaporation of convective rainfall(kg/kg) 
c
c
 
      subroutine conv(ipass,lg,tdt,pg,qsg,gam,
     &                ttg,qtg,rhg,precc,qfg,
     &                fmp,kta,kba,fluxc,rainx,qevc)

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      real Dva,rKa
      parameter (Dva=2.21) !Diffusivity of qv in air (0 deg. and 1 Pa)
      parameter (rKa=2.4e-2) !Thermal conductivity of air (0 deg)

C Argument list
      integer ipass,lg
      real tdt
      real pg(ln2)
      real qsg(ln2,nl)
      real gam(ln2,nl)
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real rhg(ln2,nl)
      real precc(ln2)
      real qfg(ln2,nl)
      real fmp(ln2,nl)
      integer kta(ln2)
      integer kba(ln2)
      real fluxc(ln2,nl)
      real rainx(ln2)
      real qevc(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'

C Global data blocks
      include 'CNSTA.f'
      include 'FEWFLAGS.f'
      include 'PRINTT.f'
      include 'TIMEX.f'

      character*1 cnmlp(ln2,lat),rainp(ln2,lat)
     & ,snmap(ln2,lat),snwice(ln2,lat),snicel(ln2,lat)
      common/cloud1/cnmlp,rainp,snmap,snwice,snicel

      integer jclb,nlow,nmid,ktied,kbgel,klcmc,klowrn,kmidrn
      real aftsea
      common/levdata/jclb,nlow,nmid,aftsea,ktied,kbgel
     &,klcmc,klowrn,kmidrn

      integer icln,ipcln
      real avcvrn,pavcvrn
      common/newcl/icln(ln2,2,lat),avcvrn(ln2,lat)
     &           ,ipcln(ln2,2,lat),pavcvrn(ln2,lat)

C Local work arrays and variables
      real dskm(ln2,nl),deltm(ln2,nl)
      real dq(nl),ss(ln2,nl),ug(ln2,nl),revq(nl)
      real dryttg(ln2,nl)         ! Amount of dry adjustment
      integer icvb(ln2)

      double precision u,cam,ds
      dimension u(nl),cam(nl-1,nl),ds(nl)
      double precision heatpc, htpold, beta, x, y, sumd, alxx,
     &     smin, htpolx, temp, convp, xb, sum, bcnv

      character*1 chcon
      logical test

      integer i
      integer iter
      integer j
      integer jclt
      integer jj
      integer k
      integer kb
      integer kt
      integer m
      integer mg
      integer ns

      real alpha
      real apr
      real btcnv
      real bpr
      real bl
      real cev
      real clfra
      real conrev
      real dqsdt
      real dta
      real ecmwa1
      real ecmwa2
      real ecmwcc
      real ecmwev
      real es
      real evap
      real fr
      real pk
      real qclt
      real ratx
      real rhmois
      real rhmx
      real rhoa
      real rnrt
      real ssrat
      real sumqs
      real sumqt
      real sumss
      real temt
      real temx
      real temxk
      real temxkp

C Local data, functions etc
      real rhcv
      data rhcv/0.75/
      include 'ESTABL.f'

C Start code : ----------------------------------------------------------

c   convective cloud

c Fluxc will be the flux of conv rain leaving layer k in kg/m**2,

      do k=1,nl
        do mg=1,ln2
          fluxc(mg,k)=0.
          qevc(mg,k)=0.
        enddo
      enddo
      do k=2,nl
        do mg=1,ln2
          deltm(mg,k)=-0.5*cappa*log(prf(mg,k)/prf(mg,k-1))
        enddo
      enddo

      if(qcloud)then
        do k=1,nl
          do mg=1,ln2
            pk=100.0*prf(mg,k)
            qsg(mg,k)=qsat(pk,ttg(mg,k))
            rhg(mg,k)=qtg(mg,k)/qsg(mg,k)
            dqsdt=qsg(mg,k)*hl/(rvap*ttg(mg,k)**2)
            gam(mg,k)=(hl/cp)*dqsdt
          enddo
        enddo
      endif

c---- Keep the initial temperatures so that the amount of dry adjustment
c---- (if any) can be computed. Dry instabilities are removed before
c---- any moist convection is computed to avoid numerical instablities.
c---- The original dry instabilities can then be replaced following moist
c---- convection adjustments if so desired.
      do k=1,nl
      do mg=1,ln2
        dryttg(mg,k)=ttg(mg,k)
      enddo
      enddo

c**** check for dry instability. (adjust with overkill).

CSJP  Former machine dependence at this point
         k=0
 6       k=k+1
         if (k.eq.nl) goto 9
 7       test=.true.
         do 8 mg=1,ln2
c**   dry adjustment of temps if (s(k)-s(k+1))/cp>0 .
c**   if adjustment occurs, then skip back down to
c**   level underneath to check if adjustment has
c**   created instability at that level.
            temxk=1.0-deltm(mg,k+1)
            temxkp=1.0+deltm(mg,k+1)
            temx=temxk*ttg(mg,k)-temxkp*ttg(mg,k+1)
            if (temx.le.0.0) go to 8
            test=.false.
            dta=0.1+temx*dprf(mg,k+1)/
     &      (dprf(mg,k+1)*temxk+dprf(mg,k)*temxkp)
            ttg(mg,k)=ttg(mg,k)-dta
            ttg(mg,k+1)=ttg(mg,k+1)+(dprf(mg,k)/dprf(mg,k+1))*dta
 8        continue
c--   check if any points adjusted.
c--   if so, and k.ne.1 , skip back 1 level.
         if (k.eq.1) go to 6
         if (test) go to 6
c--   adjustment has occured : skip back 1 level
         k=k-1
         go to 7

 9    continue 

      do k=1,nl
      do mg=1,ln2
c---- compute the temp adjustments due to dry instability removal
        dryttg(mg,k)=dryttg(mg,k)-ttg(mg,k)
      enddo
      enddo

c****
c**** convection
c****

ccccc New scheme : ECMWF evaporation of rain below cloud base only.
ccccc Using old method with moistening parameter inside towers

      rhmois=0.4
      ECMWcc=sqrt(0.05)
c         ECMWcc=0.

c---- convection cloud base set at jclb (see initax.f)
c----  select maximum level for cloud top (kuo uses nl-2)
c---- jclt=nl-2
      jclt=nl 
         btcnv=1.0
         bcnv=tdt/(btcnv*3600.0)
c**** preset conv mass flux history array to 0.0
      do 10 k=1,nl
         do 10 mg=1,ln2
 10      fmp(mg,k)=0.0
      do 20 mg=1,ln2
         icvb(mg)=nl
         kba(mg)=0
 20      kta(mg)=0

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

c      PRECOMPUTE INSTABILITY PARAMETERS FOR ENTIRE ROW 
      do 41 k=1,jclt
         do 41 mg=1,ln2
can use  dskm(mg,k)=dprf(mg,k)/pg(mg)
   41    dskm(mg,k)=dsk(k)*muf(mg,k)/pg(mg)
      do 40 k=jclb+1,jclt
         do 40 mg=1,ln2
c--   dry static energy ss()
         ss(mg,k)=ttg(mg,k)*(1.0+deltm(mg,k))
     &           -ttg(mg,k-1)*(1.0-deltm(mg,k))
c--   moist instability parameters
         ug(mg,k)=max(hlcp*(qtg(mg,k-1)-qsg(mg,k))-ss(mg,k),0.0)
c--   if the rh value is insufficient, reset ug() value to 0.0
         if (rhg(mg,k-1).lt.rhcv) ug(mg,k)=0.0
c--   create an array indicating cloud base (if any)
         if (ug(mg,k).gt.0.0) icvb(mg)=min(icvb(mg),k)
 40   continue

c**** check for conv at each point from level=jclb up
      do 360 mg=1,ln2
c**** test for convection (u(k+1)=(h(k)-h*(k+1))/cp>0)
c**** from precomputed instability parameters ug(mg,k)=u(k)
c**** set min cloud base for conv to be from level jclb
c**** convection control by rh>rhcv
c----    stability parameter : u(j+1)=(h(j)-h*(j+1))/cp
c     u(j+1)=max(hlcp*(qtg(mg,j)-qsg(mg,j+1))-ss(mg,j+1),0.0)
c---- instability shown by icvb()<nl
         if (icvb(mg).eq.nl) go to 360
         kb=icvb(mg)-1
         u(kb+1)=ug(mg,kb+1)
         sumss=ss(mg,kb+1)
         kt=kb+1
         if (kt.eq.jclt) go to 60
c---- run up the levels stopping when no instability
         do 50 k=kb+2,jclt
c--   u(k)=(h(kb)-h*(k))/cp
            sumss=sumss+ss(mg,k)
            u(k)=hlcp*(qtg(mg,kb)-qsg(mg,k))-sumss
c---- if unstable go to next level up
            kt=jclt
            if (u(k).gt.0.0) go to 50
c---- u(k)<=0.0 : stable level reached : jump out of loop
            kt=k-1
            sumss=sumss-ss(mg,k)
            go to 60
 50      continue
c---- if end of (do 50) reached then conv has risen to
c----  maximum allowable level = jclt
c---- m=number of levels involved in convection
 60      m=kt-kb+1


c---- set up moisture changes
         temt=u(kt)/(1.0+gam(mg,kt))
         qclt=qsg(mg,kt)+(gam(mg,kt)*temt/hlcp)
         sumqs=0.0
         sumqt=0.0
c---- this next do loop is much faster if not vectorised.
c---- (it is incrementing at steps of 64 => memory conflicts)
CDIR$ NEXTSCALAR
         do 70 j=kb,kt
 70         sumqt=sumqt+qtg(mg,j)*dskm(mg,j)

CDIR$ NEXTSCALAR
         do 71 j=kb+1,kt
 71         sumqs=sumqs+qsg(mg,j)*dskm(mg,j)
	 rhmx=min(rhmois,0.90*sumqt/sumqs)
         alpha=(qtg(mg,kb)-qclt)/(sumqt-rhmx*sumqs)
         dq(1)=alpha*qtg(mg,kb)
         do 80 j=kb+1,kt
            i=j-kb+1
 80         dq(i)=alpha*(qtg(mg,j)-rhmx*qsg(mg,j))

c---- set up total heating
         heatpc=temt+sumss

cZZZZ trial generation of convective cloud
c        icln(mg,1,lg)=kb+1
c        icln(mg,2,lg)=kt
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


c---- To determine heating profile : 2 iterations
c---- first preset heating at cloud base (ds(1)) to be zero.
c---- This determines the approx shape of the profile from
c---- kb+1 to kt. If ds(2) is found to be less than zero, then
c---- recompute with new ds(1) which is less than the ds(2)
c---- value just computed.
c---- (if 2nd iteration needed, current  method uses ds1=3*ds2)

         if (m.gt.2) go to 90
c---- special case for m=2 (no iterations)
c---- m=2 : convection from 1 level to next only - special coding
c        ds(1)=0.0
c        ds(2)=heatpc/dskm(mg,kt)
c.... use ratio of lowest 2 levels of heating
         ratx=ss(mg,kb+1)/(ss(mg,kb+2)+ss(mg,kb+1))
         ds(2)=(heatpc/dskm(mg,kt))/(1.0+ratx)
         ds(1)=ratx*ds(2)*dskm(mg,kt)/dskm(mg,kb)
         go to 260
c**** m>=3 : convection for 3 or more levels
 90      ds(2)=0.0
         htpold=heatpc
         iter=0
 100     iter=iter+1
         ds(1)=3.0*ds(2)
         beta=ds(1)*dskm(mg,kb)/htpold
         heatpc=(1.0-beta)*htpold
         do 110 j=1,m-1
            do 110 i=1,m-1
 110        cam(i,j)=0.0
c---- set up matrix cam(,) for solution by back substitution
c--
         xb=(1.0-deltm(mg,kb+1))*ds(1)-hlcp*dq(1)
         do 120 k=kb+1,kt-1
            i=k-kb
 120        cam(i,m)=xb*(u(k+1)-u(kb+1))
c--   'tr' elimination terms
         do 140 k=kb+1,kt-1
            i=k-kb
            cam(i,1)=u(k+1)*(1.0+gam(mg,kb+1)+deltm(mg,kb+1))-u(kb+1)*
     &       (deltm(mg,kb+1)+deltm(mg,kb+2))
            cam(i,i+1)=-u(kb+1)*(1.0+gam(mg,k+1)+deltm(mg,k+1))
            if (m.eq.3) go to 140
            jj=kb+1
            do 130 j=2,i
               jj=jj+1
 130           cam(i,j)=-u(kb+1)*(deltm(mg,jj)+deltm(mg,jj+1))
 140     continue
c---- energy conservation terms
         i=i+1
         j=0
         do 150 k=kb+1,kt
            j=j+1
 150        cam(i,j)=dskm(mg,k)
         cam(i,m)=heatpc
c---- solve this matrix set of equations
         do 180 i=m-1,2,-1
            x=cam(i,i)
            y=cam(i-1,i)
            do 160 j=1,m
               cam(i,j)=cam(i,j)*y
 160           cam(i-1,j)=cam(i-1,j)*x
            do 170 j=1,m
 170           cam(i-1,j)=cam(i-1,j)-cam(i,j)
 180     continue
         ds(2)=cam(1,m)/cam(1,1)
         do 200 i=2,m-1
            sum=cam(i,m)
            do 190 j=1,i-1
 190           sum=sum-cam(i,j)*ds(j+1)
 200        ds(i+1)=sum/cam(i,i)
c--   check for appropriate lower gradient after 1 iteration
         if (ds(2).ge.0.0) go to 210
         if (iter.eq.1) go to 100

c---- reform some profiles for >= 0 heating rates.
c---- (following second iteration etc)
 210     continue
c..... next bit if >=0.0 heating at lowest level required
         if(ds(1).lt.0.0)then
         sumd=0.0
         do 220 k=kb,kt
 220        sumd=sumd+dskm(mg,k)
            alxx=htpold/(htpold-ds(1)*sumd)
            do 225 j=2,m
 225           ds(j)=alxx*(ds(j)-ds(1))
            ds(1)=0.0
         end if
c.... next bit ensures >=0.0 heating at all levels above cloud base
         smin=0.0
         do 230 j=2,m
 230        smin=min(smin,ds(j))
         if (smin.lt.0.0) then
         sumd=0.0
         do 235 k=kb+1,kt
 235        sumd=sumd+dskm(mg,k)
            htpolx=htpold-dskm(mg,kb)*ds(1)
            alxx=htpolx/(htpolx-smin*sumd)
            do 240 j=2,m
 240           ds(j)=alxx*(ds(j)-smin)
         end if
c.... adjust ratio of lowest 2 levels of heating
         ssrat=(ss(mg,kb+1)/(ss(mg,kb+2)+ss(mg,kb+1)))
         ratx=ssrat*dskm(mg,kb+1)/dskm(mg,kb)
         ds(2)=ds(2)/(1.0+ssrat)
         ds(1)=ratx*ds(2)

c---- determine the mass flux : use the u(kb+1) exp decay equation
 260   temp=u(kb+1)/(hlcp*dq(1)+(1.0+gam(mg,kb+1)+deltm(mg,kb+1))*ds(2)
     &   -(1.0-deltm(mg,kb+1))*ds(1))
         convp=bcnv*temp
         do 270 j=1,m
            ds(j)=convp*ds(j)
 270        dq(j)=convp*dq(j)

c---- add the convection changes to the temp and moisture
c---- fields. store the conv mass flux - for mmtm mixing.
c---- update the rh values. calc total precip.

         Rnrt=0.0
ch       conrev=1000.0*(100.0*pg(mg)/grav/tdt)
         conrev=1000.0*(100.0/grav/tdt)
c        sumt=0.0
c---- Determine additive rainfall from all layers
CDIR$ NEXTSCALAR
         do 340 k=kt,kb,-1
         i=k-kb+1
c        Rnrt=SUM[dq(i)*dsk(k)]*100.0*pg(mg)/grav/tdt  !kgm/m**2/sec
c        Rnrt=1000.0*Rnrt                     !gm/m**2/sec
cq         if(ttg(mg,k).gt.253.15.or..not.qcloud.or.ipass.eq.2)then
           Rnrt=Rnrt+dq(i)*dprf(mg,k)
cq         else
cq           qfg(mg,k)=qfg(mg,k)+dq(i)
cq         endif
c           sumt=sumt+ds(i)*dskm(mg,k)
        fluxc(mg,k)=Rnrt*conrev*tdt/1000 !Flux of conv rain leaving k in kg/m^2
  340 continue
         Rnrt=Rnrt*conrev
         do 3401 k=kb,kt
         i=k-kb+1
            qtg(mg,k)=qtg(mg,k)-dq(i)
            ttg(mg,k)=ttg(mg,k)+ds(i)
            fmp(mg,k)=convp/tdt
            rhg(mg,k)=qtg(mg,k)/(qsg(mg,k)+gam(mg,k)*ds(i)/hlcp)
 3401 continue
c....
c.... add convective rainfall rate at cloud base for conv cloud
c.... determination (Slingo,87) (rate each step in units mms/day)
c.... => divide by number of steps when using
        rainx(mg)=0.5*(Rnrt/1000.0)*tdt*(1440.0/mstep)
        avcvrn(mg,lg)=avcvrn(mg,lg) + rainx(mg)*(2.0-ipass)
        pavcvrn(mg,lg)=pavcvrn(mg,lg) + rainx(mg)*(ipass-1.0)

c Reevaporation of convective rainfall

       if(qcloud)then !Use LDR style evap of rainfall
         do k=kb-1,1,-1
           es=qsg(mg,k)*100*prf(mg,k)/epsil
           Apr=(hl/(rKa*ttg(mg,k)))*(hl/(rvap*ttg(mg,k))-1)
           Bpr=rvap*ttg(mg,k)/((Dva/(100*prf(mg,k)))*es)
           clfra=0.1 !Rainy fraction of box
           Fr=Rnrt/1000.
           rhoa=100.*prf(mg,k)/(rdry*ttg(mg,k))
ck           rhorav=Rnrt/5. !g/m^3 (Kessler)
ck           Cev=5.44e-4*rhorav**0.65 !Kessler
           Cev=3.8e2*sqrt(clfra*Fr/rhoa)/(qsg(mg,k)*(Apr+Bpr))
           bl=1+0.5*Cev*tdt*(1+gam(mg,k))          !Time centred scheme
           evap=tdt*(Cev/bl)*(qsg(mg,k)-qtg(mg,k)) !In mixing ratio units
           evap=min(evap, (qsg(mg,k)-qtg(mg,k))/(1+gam(mg,k))) !Don't supersat
           if(evap.gt.0.)then
             revq(k)=min(evap,Rnrt/(conrev*dprf(mg,k)))
           else
             revq(k)=0.
           endif
          Rnrt=max(0.0,Rnrt-revq(k)*conrev*dprf(mg,k))
c          sumt=sumt-(revq(k)*hlcp)*dskm(mg,k)
          fluxc(mg,k)=Rnrt*tdt/1000 !Flux of conv rain leaving k in kg/m^2
        enddo
      else   !Use ECMWF style evap
c---- ECMWF based evaporation of convective rainfall
c---- passing through lower layers. Convective rain has an assumed
c---- grid coverage of Cc = 5%. The evap rate is controlled by
c---- E = Cc.a1.(Qsat-Q) ( sqrt(sigma) * Rnrt / a2 / Cc)**a3
c---- where a1=5.44E-04, a2=5.09E-03, and a3 = 0.5777
c---- Here Rnrt is the rainfall rate in gm/m**2/sec.
c---- See "Research manual 3" from ECMWF.
         ECMWa1=5.44e-04
c????    ECMWa2=5.09e-03         ! This is ECMWF value : ??????
         ECMWa2=5.09             ! This is the Kessler value.
c        ECMWa3=0.5777           ! This ia approximated by sqrt
c---- Calculate evap from top down : The evap may become
c---- equal to the total rain (integrated down to that level) before
c---- reaching the next level (or the surface).
         do 341 k=kb-1,1,-1
c        Vo=(ECMWa2*Rnrt**(1.0/8.0)/sqrt(sig(k)))**(8.0/9.0)
c     ECMWev is in sec**(-1)
            ECMWev=ECMWcc*ECMWa1*max((qsg(mg,k)-qtg(mg,k)),0.0)*
     &       sqrt(sqrt(prf(mg,k)/pg(mg))*Rnrt/ECMWa2)
c     revq is the change in q due to evap over a timestep of 2dt
c     Total evap is not allowed to exceed 100% of available rainfall
            revq(k)=min(tdt*ECMWev,Rnrt/(conrev*dprf(mg,k)))
            Rnrt=max(0.0,Rnrt-revq(k)*conrev*dprf(mg,k))
            fluxc(mg,k)=Rnrt*tdt/1000 !Flux of conv rain leaving k in kg/m^2
c           sumt=sumt-(revq(k)*hlcp)*dskm(mg,k)
  341 continue
      endif
         do 3411 k=kb-1,1,-1
            qtg(mg,k)=qtg(mg,k)+revq(k)
            qevc(mg,k)=revq(k)
            ttg(mg,k)=ttg(mg,k)-revq(k)*hlcp
            rhg(mg,k)=qtg(mg,k)/(qsg(mg,k)-gam(mg,k)*revq(k))
 3411 continue

         kba(mg)=kb
         kta(mg)=kt

ch       sumq=Rnrt/conrev
c        sumq=Rnrt/conrev/pg(mg)
c        ns=(mg-1+lon)/lon
c        ma=mg-(ns-1)*lon
c        write (6,350) ma,lg,ns,sumq,sumt/hlcp
c350     format (1x,3i3,' sumq,sumt ',2e15.8)

         precc(mg)=precc(mg)+0.5*(Rnrt/1000.0)*tdt
 360  continue

      do k=1,nl
      do mg=1,ln2
c---- Code to allow the vertical temp profile to retain the original 
c---- well mixed dry instabilities (i.e. replace dryttg)
        ttg(mg,k)=ttg(mg,k)  +  dryttg(mg,k)
      enddo
      enddo

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,i1)')'After conv. IPASS = ',ipass
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,99)'qtg ',(qtg(mg,k),k=1,nl)
          write(25,91)'rhg ',(rhg(mg,k),k=1,nl)
          write(25,1)'precc ',precc(mg)
          write(25,*)
        endif
      endif
 1    format(3(a,g10.3))
 91   format(a,30f10.3)
 99   format(a,30g10.3)

      if ((mod(mins+int(mstep, 8),1440_8).eq.0_8).and.cvrnm) then
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
      endif

      return
      end
