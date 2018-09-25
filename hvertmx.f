c (1) Fixing syntax errors identified by the g95 Fortran compiler.
c (2) Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c Removed unnecessary "include 'MACHINE.f'"
c SJP 2001/11/22
c
c $Log: hvertmx.f,v $
c Revision 1.66  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.65  2001/10/12 02:06:58  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.64.1.1  2001/10/12 02:13:46  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.64  2001/06/04 02:27:03  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.63  2001/02/28 04:36:40  rot032
c Further tidy ups from HBG
c
c Revision 1.62  2001/02/22 06:46:47  rot032
c Merge LDR and HBG changes.
c
c Revision 1.61  2001/02/22 05:34:44  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.60.1.1  2001/02/22 05:56:38  rot032
c Changes from LDR to sulfur cycle, to bring it to run N52
c
c Revision 1.60  2000/06/20 02:54:07  rot032
c Merge of HBG and LDR changes
c
c Revision 1.59  2000/06/20 02:08:34  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.58.1.1  2000/06/20 02:27:11  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.58  1999/05/20 06:23:57  rot032
c HBG changes to V5-2
c
c Revision 1.57  1998/12/10  00:55:55  ldr
c HBG changes to V5-1-21
c
c Revision 1.56  1997/12/17  23:22:58  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.55  1997/10/06  07:57:30  ldr
c Final corrections to V5-1.
c
c Revision 1.54  1997/07/24  05:59:17  ldr
c Mods to re-tune qcloud for 9L model from LDR.
c
c Revision 1.53  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.52  1996/10/24  01:02:52  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.51  1996/06/13  02:06:45  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.50  1996/03/21  03:18:48  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.49  1996/02/19  04:09:50  ldr
c Generalize for 24 levels.
c
c Revision 1.48  1996/01/09  06:18:16  ldr
c This is the version used for the long run H06 and write-up.
c
c Revision 1.47  1995/11/30  04:40:16  ldr
c Final updates to qcloud scheme to bring it to long run h63 used for
c BMRC workshop and write-up.
c
c Revision 1.46  1995/11/23  06:03:29  ldr
c Update qcloud to run h32 (used for BMRC & long run)
c
c Revision 1.45  1995/08/31  05:15:06  ldr
c Update qcloud to run g69, conbined with HBG tidy-ups to V4-7-16l.
c
c Revision 1.44  1995/08/31  04:30:43  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.43  1995/08/29  02:29:42  ldr
c Merge of HBG's corrections to V4-7-13h with LDR's changes for run g61.
c
c Revision 1.42  1995/08/29  01:48:50  ldr
c MErge of HBG's corrections to V4-7-13h with LDR changes.
c
c Revision 1.40.1.1  1995/08/29  01:30:09  ldr
c HBG's corrections to his V4-7-13h, to make results with hybrid=F agree
c with previous version.
c
c Revision 1.39.1.2  1995/08/29  01:50:06  ldr
c Update qcloud to run g61.
c
c Revision 1.39.1.1  1995/08/18  06:14:46  ldr
c LDR's changes to V4-7-12l to bring cloud scheme to run g57.
c
c Revision 1.39  1995/08/08  01:58:09  ldr
c Corrections to V4-5-30mic to make it work on cherax, merged into V4-7-2l.
c
c Revision 1.30.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.38  1995/06/30  02:44:42  ldr
c Changes to qcloud (run F88) and put qsat formula into ESTABL.f
c
c Revision 1.37  1995/05/02  06:10:36  ldr
c Changes to V4-6-16l from LDR - mainly affecting qcloud scheme.
c This version used for run F35.
c
c Revision 1.36  1994/12/23  02:56:31  ldr
c Little fix to avoid divzero on Cray.
c
c Revision 1.35  1994/12/22  23:51:06  ldr
c Mods to cloud scheme of V4-6-13l to produce faster version for HBG and
c to allow vertically subgrid cloud below inversions in BL. Set slowice=T
c for accurate (but slow) ice fall calculation. Also, Ri now includes hlfusionn
c
c Revision 1.34  94/09/12  12:50:11  ldr
c Changes to qcloud scheme of V4-6 to create run E32, except that E32 has
c enhanced collection of cloud water by falling ice.
c 
c Revision 1.33  94/08/08  17:21:30  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.32  94/08/08  13:13:55  ldr
c Fix Tiedtke SC scheme to stop epart going -ve if q is -ve at k=1.
c 
c Revision 1.31  94/08/04  11:09:46  ldr
c Use equal weights for Ri calculation in qcloud and modify Tiedtke SC scheme
c to use LCL instead of fixed cloud base.
c 
c Revision 1.30  94/06/17  13:44:02  ldr
c Pass in T, q rather than Tliq, qtot for qcloud scheme.
c 
c Revision 1.29  94/05/13  14:55:50  ldr
c Changes to introduce separate prognostic variable for cloud ice.
c 
c Revision 1.28  94/04/29  15:18:07  ldr
c Minor changes for 18L: start deep and shallow conv from k=3, increase Cd.
c 
c Revision 1.27  94/03/30  12:58:39  ldr
c Merge of HBG and LDR changes.
c 
c Revision 1.26  94/03/30  10:24:00  ldr
c Do Geleyn SC from k=2 only.
c 
c Revision 1.25.1.1  94/03/30  12:34:29  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c
c     INPUT/OUTPUT
c     Input    from common/fewflags in FEWFLAGS.f
c                  debug    - flag to control column debugging
c                  mgdebug  - flag to control longitude debugging
c                  lgdebug  - flag to control latitude debugging
c                  insdebug - flag to control hemisphere debugging
c                  hybrid   - if true, use new hybrid vertical coordinate
c                  qcloud   - T if using prognostic cloud scheme
c                  sltrace  - if T compute tracer transport by semi-
c                             Langrangian advection and by mixing and
c                             convection in the vertical
c
c              from common/hybrpr in HYBRPR.f
c                  dprf  - pressure thickness at each sigma level
c                  pdpsk - (pressure/surface pressure)**cappa 
c                  prf   - pressure at full levels
c                  prh   - pressure at half levels
c
c              from common/levdata in this subroutine
c                  kbgel - level at which shallow convection begins
c                  ktied - highest shallow convection cloud base
c
c              from common/logimsl in LOGIMSL.f
c                  land - logical variable for surface type
c
c              from common/traceblk in TRACEBLK.f
c                  con - concentration of atmospheric tracers at t
c
c              from arguments
c                  fg    - sensible heat flux   eg - latent heat flux
c                  gam - (L/cp)dqsdt            ns - hemisphere index
c                  ipass - model pass counter   lg - latitude index
c                  pg - surface pressure(mbs)   tg - surface temperature
c                  tdt - 2 X timestep (dt)
c                  cduv - momentum drag coefficient
c                  cfrac - cloudy fraction of grid box 
c
c     Output:  from arguments
c                  cten - rate of change of tracer concentration due to 
c                         convection
c
c     In/Out:  from common/timex in TIMEX.f
c                  mins - current model time in mins
c
c              from arguments
c                  qtg - mixing ratio    tscrn - temperature at screen level
c                  u   - zonal wind      v     - meridional wind
c                  dqgdt - mixing ratio at current time - dt
c                  ktop  - index of (data level) pressure of cloud top,
c                          used in the longwave program
c                  ttg   - physical temperature
c                  qfg   - frozen H2O cloud/ice mixing ratio (kg/kg)
c                  qlg   - liquid H2O cloud/ice mixing ratio (kg/kg)
c                  uten, vten  - effective time derivative of u, v
c                  dkevm - change in KE  
c                  kbase - bottom of shallow convective cloud
c                  rich  - Richardson number
c 
      subroutine hvertmx(ipass,lg,tdt,cduv,fg,eg,pg,tg,
     &                   qlg,qfg,cfrac,gam,tscrn,ustar,xtem, !inputs
     &                   qtg,dqgdt,u,v,ttg,uten,vten,xtg,     !input+ output
     &                   dkevm,cten,rich,kbase,ktop)      !output

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /LOGIMSL/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'CPARAMS.f'
      include 'ECPARM.f'

C Argument list
      integer ipass,lg
      real tdt
      real cduv(ln2)
      real fg(ln2)
      real eg(ln2)
      real pg(ln2)
      real tg(ln2)
      real qlg(ln2,nl)
      real qfg(ln2,nl)
      real cfrac(ln2,nl)
      real gam(ln2,nl)
      real tscrn(ln2)
      real ustar(ln2)
      real xtg(ln2,nl,ntrac)
      real xtem(ln2,ntrac)
      real qtg(ln2,nl)
      real dqgdt(ln2,nl)
      real u(ln2,nl)
      real v(ln2,nl)
      real ttg(ln2,nl)
      real uten(ln2,nl)
      real vten(ln2,nl)
      real dkevm(ln2,nl)
c     real cten(ln2,nl,ntrace) ! see below (after TRACEBLK.f)
      real rich(ln2,nl)
      integer kbase(ln2)
      integer ktop(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'
      include 'LOGIMSL.f'

C Global data blocks
      include 'FEWFLAGS.f'
      include 'TIMEX.f'
      include 'TRACEBLK.f'
      real cten(ln2,nl,ntrace) ! Must be after TRACEBLK.f

      real roncp,rong,rlogs1,rlogs2,rlogh1,rlog12
      common/hverdat/roncp,rong,rlogs1,rlogs2,rlogh1,rlog12

      integer jclb,nlow,nmid,ktied,kbgel,klcmc,klowrn,kmidrn
      real aftsea
      common/levdata/jclb,nlow,nmid,aftsea,ktied,kbgel
     &,klcmc,klowrn,kmidrn

      real ttmsav,pttmsav
      common/ttmarr/ttmsav(ln2,2,lat),pttmsav(ln2,2,lat)

C Local work arrays and variables
      real sigsp(ln2),Tsp(ln2)
      real cons(ln2,nl)
c     real dw(ln2)
      real betatt(ln2,nl),betaqt(ln2,nl)
      real thetam(ln2,nl)
      real au(ln2,nl),at(ln2,nl),cu(ln2,nl),ct(ln2,nl)
      real rhs(ln2,nl),dz(ln2),delthet(ln2),qs(ln2,nl)
      real sqmxl(ln2)
      real dskx(ln2,nl)
      real fm(ln2),fh(ln2),dzr(ln2),delsig(ln2)
      real theta(ln2,nl),x(ln2)
      real zh(ln2,0:nl),tmnht(ln2,nl),gt(ln2,nl),gu(ln2,nl)
      real rkm(ln2,nl),rkh(ln2,nl)
      real wke(ln2,nl),wkt(ln2,nl),wkb(ln2,nl)
      real xt(ln2,nl)
c  Optional ncarpbl stuff
        real dzk(ln2,nl)

      integer k
      integer kb
      integer kt
      integer lgn
      integer ma
      integer mid
      integer mg
      integer ns
      integer nt

      real adkedt
      real al
      real arg
      real betac
      real betaq
      real betat
      real condrag
      real conflux
      real dqtot
      real dvmod
      real epart
      real fice
      real hlcpfac
      real pk
      real psp
      real qc
      real qsatb
      real qsk
      real q1
      real ri
      real sfrich
      real sighkx
      real sqmxlav
      real the
      real the1
      real vkona
      real w1
      real w2

C Local data, functions etc
      real delta
      parameter (delta=1/epsil-1 ) !Used by qcloud scheme

      real amxlsq
c     data amxlsq/100./
      data amxlsq/30./

      real bprm,cm,ch,vkar,rdiv
c     can set coefficients for Louis scheme based on Kansas data:
cc    data bprm/4.7/,cm/7.4/,ch/5.3/,vkar/.35/,rdiv/.74/
c     or use the following improved coefficients with vkar=.4:
      data bprm/5./,cm/5./,ch/2.6/,vkar/.4/,rdiv/1./
      include 'ESTABL.f'

C Start code : ----------------------------------------------------------

c     Stability dependent vertical mixing routine
c     n.b. cduv, cdtq are drag coeffs mult by vmod
c
c     ttg, u, v, qtg are current values
c     pg is surface pressure(mbs), tg is surface temperature
c     dw is soil wetness   -  not used
c     sig gives sigma values
c     uten, vten are returned as effective time derivatives of u, v
c     dt is time step size

c---- Optional NCAR non-local boundary layer formulation.
c---- On if ncarpbl=.true.  Calls NCAR subroutine pbldif().
c---- (Original coding provided by Kim Nguyen, 20/07/98,
c----  and reworked for CSIRO spectral AGCM by hbg.)

      if((.not.SCM.and.(mod(mins,1440_8).eq.0_8)).or.
     &                   (SCM.and.(nsteps.eq.0)))then
        if(ipass.eq.1)then
        do 12 k=1,2
        do 12 mg=1,ln2
        pttmsav(mg,k,lg)=ttg(mg,k)
   12   ttmsav(mg,k,lg)=ttg(mg,k)
        endif
      endif

c       Asymptotic mixing length set above to 100.
      vkona=vkar/amxlsq
c Following values now preset in initax.f
c     roncp=rdry/cp
c     rong=rdry/grav
c     rlogs1=log(sig(1))
c     rlogs2=log(sig(2))
c     rlogh1=log(sigh(2))
c     rlog12=1./(rlogs1-rlogs2)
c     if(ukconv)then
c       ksc=0 ! No shallow conv if UKMO convection
c     else
c       if(qcloud)then  !Use Tiedtke shallow convection scheme with ktop=ksc
c         if(nl.gt.9)then
c           ksc=nl/3
c         else
c           ksc=-1
c         endif
c       else            !Use Geleyn shallow convection scheme
c         ksc=-1  
c       endif
c     endif

      do 148 mg=1,ln2
      tmnht(mg,1)=(ttg(mg,2)*rlogs1-ttg(mg,1)*rlogs2+
     &               (ttg(mg,1)-ttg(mg,2))*rlogh1)*rlog12
c     N.B. an approximate zh is quite adequate for this routine
      zh(mg,0)=0.
  148 zh(mg,1)=ttg(mg,1)*rong*dprf(mg,1)/prf(mg,1)
      do 15 k=2,nl
      do 15 mg=1,ln2
      zh(mg,k)=zh(mg,k-1)+
     &  ttg(mg,k)*rong*dprf(mg,k)/prf(mg,k)
 15   continue
      do 16 k=2,nl-1
      do 16 mg=1,ln2
   16 tmnht(mg,k) =(ttg(mg,k)+ttg(mg,k+1))*.5
      do 17 k=1,nl
      do 17 mg=1,ln2
      dskx(mg,k)=dprf(mg,k)/pg(mg)
      rkh(mg,k)=0.
      rkm(mg,k)=0.
      gt(mg,k)=0.
   17 gu(mg,k)=0.
      if(ipass.eq.1)then
      do 171 k=1,2
      do 171 mg=1,ln2
      thetam(mg,k)=.5*(ttmsav(mg,k,lg)+ttg(mg,k))/pdpsk(mg,k)
  171 theta(mg,k)=ttg(mg,k)/pdpsk(mg,k)
      else
      do 172 k=1,2
      do 172 mg=1,ln2
      thetam(mg,k)=.5*(pttmsav(mg,k,lg)+ttg(mg,k))/pdpsk(mg,k)
  172 theta(mg,k)=ttg(mg,k)/pdpsk(mg,k)
      endif
      do 173 k=3,nl
      do 173 mg=1,ln2
      theta(mg,k)=ttg(mg,k)/pdpsk(mg,k)
  173 thetam(mg,k)=theta(mg,k)

c     ************ section for Tiedtke shallow convection ***********
c     Effect of Tiedtke differs from bmrc's in mixing momentum

      do mg=1,ln2
        kbase(mg)=nl
        ktop(mg)=0
      enddo

      if(ksc.gt.1)then

c Calculate LCL for near surface air
c Assume qstar=qtg(mg,1) but allow some sub-grid variability
c The formula for Tsp is eqn (21) of Bolton (1980), MWR 108, 1046-1053.

        do mg=1,ln2
          q1=max(qtg(mg,1),1.0e-8)
          if(land(mg))then !Assume qstar=qtg(mg,1)
            epart=1.05*q1*pg(mg)/epsil !in hPa
          else !Could use tsurf from hsflux for sea and MLO points
            epart=1.01*q1*pg(mg)/epsil !in hPa
          endif
          Tsp(mg)=2840/(3.5*log(tscrn(mg))-log(epart)-4.805) + 55 
          sigsp(mg)=(Tsp(mg)/tscrn(mg))**(1/roncp)
c.... sigsp in hybrid coordinates is Psp/P*
        enddo

c Look for the lowest layer s.t. the top of the layer is above the LCL,
c and call that layer cloud base. ktied set up in initax.f

        do k=ktied,1,-1
          do mg=1,ln2
            if(sigsp(mg).gt.(prh(mg,k+1)/pg(mg)))kbase(mg)=k
          enddo
        enddo

c Increase Kh, Km up to cloud top, found by using thetaes profile
c Formula for thetaes is eqn (3a) of Betts (1982) JAS 39, 1484-1505.

        do 23 mg=1,ln2
        if(kbase(mg).lt.nl)then
          psp=100.0*pg(mg)*sigsp(mg)
          qsatb=qsat(psp,Tsp(mg))
          if(tscrn(mg).gt.273.15)then
            hlcpfac=2670. !Not exactly hl/cp; see Betts (1982)
          else
            hlcpfac=3000. !Approx value including hlfusion
          endif
c          the1=thetam(mg,1)*exp(hlcpfac*qsatb/Tsp(mg)) !thetaes at LCL
          arg=hlcpfac*qsatb/Tsp(mg)
          the1=thetam(mg,1)*(1+arg+0.5*arg**2) !approx thetaes at LCL
          do 22 k=kbase(mg)+1,ksc
          pk=100.0*prf(mg,k)
          qsk=qsat(pk,ttg(mg,k))
c          the=thetam(mg,k)*exp(hlcpfac*qsk/ttg(mg,k))
          arg=hlcpfac*qsk/ttg(mg,k)
          the=thetam(mg,k)*(1+arg+0.5*arg**2)
          if(the.lt.the1)then
            rkm(mg,k-1)=6.
            rkh(mg,k-1)=6.
            rkm(mg,k)=2.
            rkh(mg,k)=2.
            ktop(mg)=k
          else
            go to 23
          endif
   22     continue
        endif
   23   continue
      endif
c     *********** end of Tiedtke shallow convection section ****

c---- If Geleyn shallow convection precompute relevant qs(mg,k) values
      if(ksc.lt.0)then
c       S/C acting between levels kbgel and (nl+1)/2
c        (see initax.f for kbgel)
        kb=kbgel
        kt=(nl+1)/2+1
        do 2 k=kb,kt
        do 2 mg=1,ln2
        pk=100.0*prf(mg,k)
 2      qs(mg,k)=qsat(pk,ttg(mg,k))
      end if

c Pre-calculate the buoyancy parameters if using qcloud scheme.
c Follow Smith's (1990) notation; gam() is HBG's notation for (L/cp)dqsdt.
c The factor of (1/sigkap)=T/theta in betatt differs from Smith's formulation
c because we mix theta rather than T.

      if ( qcloud ) then
        do k=1,nl
          do mg=1,ln2
            thetam(mg,k)=thetam(mg,k) !Convert to thetal - used only to calc Ri
     &                  -(hlcp*qlg(mg,k)+hlscp*qfg(mg,k))/pdpsk(mg,k)
            betat=1/ttg(mg,k)
            qc=qlg(mg,k)+qfg(mg,k)
            if(qc.gt.1.0e-12)then
              fice=qfg(mg,k)/qc
            else
              fice=0.
            endif
            betaq=delta/(1+delta*qtg(mg,k)-qc)
            al=1/(1+gam(mg,k))
            betac=cfrac(mg,k)
     &           *al * ((hlcp+fice*hlfcp)*betat - betaq/(1-epsil) )
            betatt(mg,k)=(betat-(gam(mg,k)/hlcp)*betac)
     &                                               *pdpsk(mg,k) !Beta_t_tilde
            betaqt(mg,k)=betaq+betac                              !Beta_q_tilde
          enddo
        enddo
      endif
c----

      do 4 k=1,nl-1
      do 3 mg=1,ln2
      dz(mg) =tmnht(mg,k)*rong*(prf(mg,k)-prf(mg,k+1))/prh(mg,k+1)
      dzk(mg,k) =dz(mg)
    3 delthet(mg)=thetam(mg,k+1)-thetam(mg,k)

c     ************ section for Geleyn shallow convection *******
      if(ksc.lt.0)then
        if(k.ge.kb.and.k.le.kt-1)then
        do 31 mg=1,ln2
   31   delthet(mg)=delthet(mg)-hlcp*
     &  max(0. , qs(mg,k+1)-qtg(mg,k+1)-qs(mg,k)+qtg(mg,k)   )
        endif
      endif
c     *********** end of Geleyn shallow convection section ******

      if(hybrid)then
        do mg=1,ln2
          sqmxl(mg)=(vkar*zh(mg,k)/(1.+vkona*zh(mg,k)))**2
c         alc=min(30.0,vkar*zh(mg,k)/(1.+vkona*zh(mg,k)))
c         sqmxl(mg)=alc*alc
        enddo
      else
c     An average csq and sqmxl for the latitude circle will be fine
c     Could similarly do: dz, zh, and maybe tmnht
        do ns=1,2
        ma=(ns-1)*lon
        mid=ma+(lon+1)/2 !Works OK for GCM, and also SCM (when lon=1)
        sqmxlav=(vkar*zh(mid,k)/(1.+vkona*zh(mid,k)))**2
c       alc=min(30.0,vkar*zh(mid,k)/(1.+vkona*zh(mid,k)))
c       sqmxlav=alc*alc
        do mg=1+ma,lon+ma
          sqmxl(mg)=sqmxlav
        enddo
        enddo
      endif

c     x is bulk ri *(dvmod **2)

      if ( qcloud ) then !follow Smith (1990) scheme and notation
        do mg=1,ln2
c          dz1=zh(mg,k)-zh(mg,k-1) !thickness of layer k
c          dz2=zh(mg,k+1)-zh(mg,k) !thickness of layer k+1
c          w1=dz1/(dz1+dz2)        !weight for lower level
          w1=0.5
          w2=1.0-w1                 !weight for upper level
          dqtot=qtg(mg,k+1)+qlg(mg,k+1)+qfg(mg,k+1)
     &         -  (qtg(mg,k)+qlg(mg,k)+qfg(mg,k))
          x(mg)=grav*dz(mg)*(
     &          (w1*betatt(mg,k)+w2*betatt(mg,k+1))*delthet(mg) +
     &          (w1*betaqt(mg,k)+w2*betaqt(mg,k+1))*dqtot )
        enddo
      else
        do mg=1,ln2
          sighkx=(prh(mg,k+1)/pg(mg))**roncp
          x(mg)=grav*dz(mg)*(delthet(mg)/
     &    (tmnht(mg,k)/sighkx)+.61*(qtg(mg,k+1)-qtg(mg,k)))
        enddo
      endif

c      Smith (1990) Fh,Fm

        do mg=1,ln2
        dvmod=max(sqrt((u(mg,k+1)-u(mg,k))**2  +
     &                 (v(mg,k+1)-v(mg,k))**2) , 1.)
        ri=x(mg)*dvmod**(-2)
        rich(mg,k)=ri
        if(ri.lt. 0.)then
c       Unstable case
c       first do momentum
          fm(mg)=dvmod*(1.0-10.0*ri/(1.0+2.5*sqrt(abs(ri))))
c       Now heat
          fh(mg)=dvmod*(1.0-10.0*ri/(1.0+0.4*sqrt(abs(ri))))
        else
c       Stable case
          if(ri.lt.0.05)then
            fm(mg)=dvmod*(1.0-5.0*ri)**2
          else
            fm(mg)=dvmod*1.6875/(1.0+40.0*ri)
          endif
          fh(mg)=fm(mg)
        endif
        enddo
c
c     Calculate k's, gu and gt
c
c----
      do 393 mg=1,ln2
      dzr(mg)=1./dz(mg)
  393 delsig(mg)=(prf(mg,k+1)-prf(mg,k))/pg(mg)
      do 395 mg=1,ln2
          rkm(mg,k)= rkm(mg,k)+fm(mg)*sqmxl(mg)*dzr(mg)
          rkh(mg,k)= rkh(mg,k)+fh(mg)*sqmxl(mg)*dzr(mg)
  395 continue
    4 continue

C     NCAR
      if(ncarpbl) call pbldif(lg,tdt,u,v,ttg,tg,pg,fg,eg,ustar,xtem, !Inputs
     &                  rkm,rkh,theta,qtg,xtg)                    !In & Out
C
C     Computing guv, gt term
C
      do k=1,nl-1
         do mg=1,ln2
            dzr(mg)=1./dzk(mg,k)
            delsig(mg)=(prf(mg,k+1)-prf(mg,k))/pg(mg)
            gu(mg,k)=rkm(mg,k)*tdt*delsig(mg)*dzr(mg)**2
            gt(mg,k)=rkh(mg,k)*tdt*delsig(mg)*dzr(mg)**2
         enddo
      enddo
c
      condrag=-grav*tdt/rdry
      conflux=-grav*tdt/100

      do 45 mg=1,ln2
      at(mg,1)=0.0
   45 au(mg,1) =cduv(mg)*(condrag/dskx(mg,1))/tg(mg)

      do 55 k=2,nl
      do 55 mg=1,ln2
      at(mg,k) =gt(mg,k-1)/dskx(mg,k)
   55 au(mg,k) =gu(mg,k-1)/dskx(mg,k)
c----
      do 56 k=1,nl
      do 56 mg=1,ln2
      cu(mg,k) =gu(mg,k)/dskx(mg,k)
   56 ct(mg,k) =gt(mg,k)/dskx(mg,k)
c
c     First do u ; pass back time tendency only, through uten array
      do 58 k=1,nl
      do 58 mg=1,ln2
   58 rhs(mg,k)=u(mg,k)
      call trim(rhs,au,cu,rhs,0,wke,wkt,wkb)
      do 59 k=1,nl
      do 59 mg=1,ln2
      uten(mg,k)=(rhs(mg,k)-u(mg,k))/tdt
   59 u(mg,k)=rhs(mg,k)

c     Now do v ; pass back time tendency only, through vten array
      do 61 k=1,nl
      do 61 mg=1,ln2
   61 rhs(mg,k)=v(mg,k)
      call trim(rhs,au,cu,rhs,1,wke,wkt,wkb)
      do 62 k=1,nl
      do 62 mg=1,ln2
      vten(mg,k)=(rhs(mg,k)-v(mg,k))/tdt
   62 v(mg,k)=rhs(mg,k)

c---- Compute the KE change per pairs of levels for calculating heating
c---- caused by vertical mixing
      do 621 k=1,nl-1
      do 621 mg=1,ln2
      adkedt=(dprf(mg,k)*cu(mg,k)/tdt)*
     & ((u(mg,k)-u(mg,k+1))**2+(v(mg,k)-v(mg,k+1))**2)
     &  /(dprf(mg,k)+dprf(mg,k+1))
      dkevm(mg,k  )=dkevm(mg,k  )+adkedt
  621 dkevm(mg,k+1)=dkevm(mg,k+1)+adkedt
c---- Compute the loss of KE at level 1 caused by surface drag.
      do 622 mg=1,ln2
c     if(land(mg))then
       sfrich=(au(mg,1)/tdt)*(u(mg,1)**2+v(mg,1)**2) ! Loss of KE
       dkevm(mg,1)=dkevm(mg,1)+sfrich
c     endif
  622 continue

c     Now do theta (then convert back to ttg)
      do 63 k=1,nl
      do 63 mg=1,ln2
      rhs(mg,k)=theta(mg,k)
   63 continue
      do 64 mg=1,ln2
      rhs(mg,1)=rhs(mg,1)-(conflux/dprf(mg,1)/cp)*fg(mg)
   64 continue
      call trim(rhs,at,ct,theta,0,wke,wkt,wkb)
      do 66 k=1,nl
      do 66 mg=1,ln2
   66 ttg(mg,k)=theta(mg,k)*pdpsk(mg,k)
c
c     Now do moisture
c     Put input value of qtg into theta array temporarily
      do 67 k=1,nl
      do 67 mg=1,ln2
      theta(mg,k)=qtg(mg,k)
      rhs(mg,k)=qtg(mg,k)
   67 continue
      do 68 mg=1,ln2
   68 rhs(mg,1)=rhs(mg,1)-(conflux/dprf(mg,1)/hl)*eg(mg)
      call trim(rhs,at,ct,qtg,1,wke,wkt,wkb)

c Do cloud liquid water and cloud ice

      if(qcloud)then
        do k=1,nl
          do mg=1,ln2
            rhs(mg,k)=qlg(mg,k)
          enddo
        enddo
        call trim(rhs,at,ct,qlg,1,wke,wkt,wkb)
        do k=1,nl
          do mg=1,ln2
            rhs(mg,k)=qfg(mg,k)
          enddo
        enddo
        call trim(rhs,at,ct,qfg,1,wke,wkt,wkb)
      endif

c     Provide a correction to dqgdt for use by hkuo later this timestep
c     Also update current ttmsav for use next time around
      do 69 k=1,nl
      do 69 mg=1,ln2
   69 dqgdt(mg,k)=dqgdt(mg,k) + (qtg(mg,k)-theta(mg,k))/tdt
      if(ipass.eq.1)then
      do 70 k=1,2
      do 70 mg=1,ln2
   70 ttmsav(mg,k,lg)=ttg(mg,k)
      else
      do 71 k=1,2
      do 71 mg=1,ln2
   71 pttmsav(mg,k,lg)=ttg(mg,k)
      endif

c Vertical mixing of tracers.
c Pass back tendency cten: add to con later in radin.
      
      if(sltrace)then
        lgn=lat2p-lg
        do 80 nt=1,ntrace

        do 73 k=1,nl
        do 73 mg=1,lon
        cons(mg,k)=con(mg,lgn,k,nt)
        rhs(mg,k)=con(mg,lgn,k,nt)
        cons(mg+lon,k)=con(mg,lg,k,nt)
        rhs(mg+lon,k)=con(mg,lg,k,nt)
   73   continue
        call trim(rhs,at,ct,cons,1,wke,wkt,wkb)
        do 75 k=1,nl
        do 75 mg=1,lon
        cten(mg,k,nt)=(cons(mg,k)-con(mg,lgn,k,nt))/tdt
        cten(mg+lon,k,nt)=(cons(mg+lon,k)-con(mg,lg,k,nt))/tdt
   75   continue

 80     continue
      endif

c Pass back the updated tracer mixing ratios xtg.

      if(coupled_aero)then
        do nt=1,ntrac
          xt(:,:)=xtg(:,:,nt)
          rhs(:,:)=xtg(:,:,nt)
          rhs(:,1)=rhs(:,1)-(conflux/dprf(:,1))*xtem(:,nt)
          call trim(rhs,at,ct,xt,1,wke,wkt,wkb)              
          xtg(:,:,nt)=xt(:,:)
        enddo
      endif          

      if(debug)then
        if(lg.eq.lgdebug)then
          ns=insdebug
          mg=mgdebug+(ns-1)*lon
          write(25,'(a,i1)')'After hvertmx. IPASS = ',ipass
          write(25,91)'ttg ',(ttg(mg,k),k=1,nl)
          write(25,9)'qtg ',(qtg(mg,k),k=1,nl)
          write(25,9)'rich ',(rich(mg,k),k=1,nl)
          write(25,*)' kbase ',kbase(mg),' ktop ',ktop(mg)
          write(25,1)'Tsp ',Tsp(mg),' sigsp ',sigsp(mg)
          write(25,9)'uten ',(uten(mg,k),k=1,nl)
          write(25,9)'vten ',(vten(mg,k),k=1,nl)
          write(25,9)'vten ',(vten(mg,k),k=1,nl)
          if(coupled_aero)then
            write(25,9)'DMS ',(xtg(mg,k,1),k=1,nl)
            write(25,9)'SO2 ',(xtg(mg,k,2),k=1,nl)
            write(25,9)'SO4 ',(xtg(mg,k,3),k=1,nl)
          endif
        endif
      endif
 1    format(3(a,f10.3))
 9    format(a,30g10.3)
 91   format(a,30f10.3)

      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
CCCCCC NCAR computing rkh & rkm, theta and qg due to counter gradient 
      subroutine pbldif(lg,tdt,u,v,ttg,tg,pg,fg,eg,ustar,xtem, !Inputs
     &                  rkm,rkh,theta,qtg,xtg)              !In & Out
C------------------------------------------------------------------------
C 
C Atmospheric boundary layer computation.
C
C Nonlocal scheme that determines eddy diffusivities based on a
C diagnosed boundary layer height and a turbulent velocity scale;
C also, countergradient effects for heat and moisture, and constituents
C are included, along with temperature and humidity perturbations which 
C measure the strength of convective thermals in the lower part of the 
C atmospheric boundary layer.
C
C For more information, see Holtslag, A.A.M., and B.A. Boville, 1993:
C Local versus Nonlocal Boundary-Layer Diffusion in a Global Climate
C Model. J. Clim., vol. 6., p. 1825--1842.
c
c Updated by Holtslag and Hack to exclude the surface layer from the
c definition of the boundary layer Richardson number. Ri is now defined
c across the outer layer of the pbl (between the top of the surface
c layer and the pbl top) instead of the full pbl (between the surface and
c the pbl top). For simiplicity, the surface layer is assumed to be the
c region below the first model level (otherwise the boundary layer depth 
c determination would require iteration).
C
C------------------------------Code history--------------------------------
C
C Original version:  B. Boville
C Standardized:      J. Rosinski, June 1992
C Reviewed:          B. Boville, P. Rasch, August 1992
C Reviewed:          B. Boville, P. Rasch, April 1996
C
C Modified for boundary layer height diagnosis: Bert Holtslag, june 1994
C >>>>>>>>>  (Use ricr = 0.3 in this formulation)
C
C-----------------------------------------------------------------------
c
C
C

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'FEWFLAGS.f'
      include 'ECPARM.f'        ! Input NTRAC (no. of tracers) and NE (.eq.1 for tracers)
      real tiny                 ! lower bound for wind magnitude
      parameter (tiny=1.e-36)

C Argument list
C
C Input arguments
C
      integer lg
      real tdt
      real u(ln2,nl)            ! zonal velocity        [m/sec]
      real v(ln2,nl)            ! meridional velocity   [m/sec]
      real ttg(ln2,nl)          ! temperature           [K]
      real tg(ln2)              ! surface temperature   [K]
      real pg(ln2)              ! surface pressure      [mbs]
      real fg(ln2)              ! surface heat flux     [W/m2 = J/m2/s]
      real eg(ln2)              ! surface evaporation   [W/m2 = J/m2/s]
                                ! note that in NCAR (cflx = eg/hl) already
                                ! done in srfoce.F before pass to pbldif.F
      real ustar(ln2)           ! friction velocity     [m/sec]
      real xtem(ln2,ntrac)      ! surface tracer flux   [kg/m2/s]

C     Input and Output arrays
      real theta(ln2,nl)        ! potential temperature [K]
      real qtg(ln2,nl)          ! water vapour          [Kgm/Kgm]
      real rkm(ln2,nl)          ! Km from Louis
      real rkh(ln2,nl)          ! Kh from Louis
      real xtg(ln2,nl,ntrac)    ! tracer m.r.           [Kgm/Kgm]
C
C     working arrays
C
      real cgh(ln2,nl)              ! counter-gradient term for heat [K/m]
      real cgq(ln2,nl)              ! counter-gradient term for moisture 
      real cgx(ln2,nl,ntrac)        ! counter-gradient term for tracers
      real cgs(ln2,nl)              ! counter-gradient star (cg/flux) [s/m2]
      real pblh(ln2)                ! boundary-layer height [m]

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'

C Local work arrays and variables
      real heatv(ln2)            ! surface virtual heat flux
      real thvref(ln2)           ! reference level virtual temperature
      real z(ln2,nl)             ! height above surface  [m]
      real phiminv(ln2)          ! inverse phi function for momentum
      real phihinv(ln2)          ! inverse phi function for heat 
      real wm(ln2)               ! turbulent velocity scale for momentum
      real zm(ln2)               ! current level height
      real zp(ln2)               ! current level height + one level up
      real khfs(ln2)             ! surface kinematic heat flux [mK/s]
      real kqfs(ln2)             ! sfc kinematic constituent flux [m/s]
      real kxfs(ln2,ntrac)       ! sfc kinematic tracer flux [m/s]
      real rino(ln2,nl)          ! bulk Richardson no. from level to ref lev
      real tlv(ln2)              ! ref. level pot tmp + tmp excess
      real fak1(ln2)             ! k*ustar*pblh
      real fak2(ln2)             ! k*wm*pblh
      real fak3(ln2)             ! fakn*wstr/wm 
      real pblk(ln2)             ! level eddy diffusivity for momentum
      real pr(ln2)               ! Prandtl number for eddy diffusivities
      real zl(ln2)               ! zmzp / Obukhov length
      real zh(ln2)               ! zmzp / pblh      at half levels
      real zzh(ln2)              ! (1-(zmzp/pblh))**2      at half levels
      real wstr(ln2)             ! w*, convective velocity scale
      real obklen(ln2)           ! Obukhov length

      logical unstbl(ln2)        ! pts w/unstbl pbl (positive virtual ht flx)
      logical stblev(ln2)        ! stable pbl with levels within pbl
      logical unslev(ln2)        ! unstbl pbl with levels within pbl
      logical unssrf(ln2)        ! unstb pbl w/lvls within srf pbl lyr
      logical unsout(ln2)        ! unstb pbl w/lvls in outer pbl lyr
      logical check(ln2)         ! True=>chk if Richardson no.>critcal

C Local variables
      integer mg                ! longitude index
      integer k                 ! level index
      integer nt                ! tracer index

C extra variables
C
      real tkv                  ! model level potential temperature
      real therm                ! thermal virtual temperature excess
      real pmid                 ! midpoint pressures
      real vvk                  ! velocity magnitude squared
      real zmzp                 ! level height halfway between zm and zp
      real rrho                 ! 1./bottom level density (temporary)
      real term                 ! intermediate calculation
      real fac                  ! interpolation factor
C      real pblmin               ! min pbl height due to mechanical mixing

      real betam   ! Constant in wind gradient expression
      real betas   ! Constant in surface layer gradient expression
      real betah   ! Constant in temperature gradient expression
      real fak     ! Constant b (eq 4.d.19) in surface temperature excess
      real onet    ! 1/3 power in wind gradient expression
      real fakn    ! Constant in turbulent prandtl number
      real ricr    ! Critical richardson number
      real sffrac  ! Surface layer fraction of boundary layer
      real vk      ! Von Karman's constant
      real ccon    ! fak * sffrac * vk
      real binm    ! betam * sffrac
      real binh    ! betah * sffrac
      real zkmin   ! Minimum kneutral*f(ri)
      real bet1    ! 0.5*rdry/grav
C
      real qgx(ln2,nl),xtx(ln2,nl,ntrac)
      real tmp1,khfv,sigmh,sigmhk,sigotbk,sigotbkm1
      real ztodtgor,delsig
      integer kk


      data betam, betas, betah, sffrac / 15.0, 5.0, 15.0, 0.1 /
      data fac,fak,fakn,ricr,vk,zkmin / 100., 8.5, 7.2, .3, 0.4,0.01/
C
      onet   = 1./3.
      ccon   = fak*sffrac*vk
      binm   = betam*sffrac
      binh   = betah*sffrac
ccccccc
C
C     Save theta,qtg and initialize background Kh and Km using local scheme
C
      do k=1,nl
         do mg=1,ln2
            qgx(mg,k)=qtg(mg,k)
            cgh(mg,k)=0.
            cgq(mg,k)=0.
            cgx(mg,k,:)=0.
            cgs(mg,k)=0.
         enddo
       enddo
       if(ne.eq.1)xtx(:,:,:)=xtg(:,:,:)

C Heights (FULL level) of model levels
       bet1=0.5*rdry/grav
       do mg=1,ln2
         z(mg,1)= bet1*(ttg(mg,1)+tg(mg))*alog(pg(mg)/prf(mg,1))
         do k=2,nl
            z(mg,k)= z(mg,k-1) +
     &      bet1*(ttg(mg,k)+ttg(mg,k-1))*alog(prf(mg,k-1)/prf(mg,k))
         enddo
       end do       !mg=1,nl

C
C      Compute kinematic surface fluxes
       do mg=1,ln2
           pmid=prh(mg,1)*100.0
           rrho = rdry*ttg(mg,1)/pmid
           ustar(mg) = max(ustar(mg),0.01)
           khfs(mg) = fg(mg)*rrho/cp    !khfs=w'theta', mK/s
           kqfs(mg) = eg(mg)*rrho/hl    !(m/s) moisture is dimensionless (kg/kg)
           if(ne.eq.1)then
             kxfs(mg,:) = xtem(mg,:)*rrho !(m/s) tracer m.r. is     "
           endif
C
C          Compute various arrays for use later:
C
C thvref(mg) is the quantity calculated from the temp and moisture of the
C first model level and of the surface (???) by applying Geleyn (1988)
C what does this mean ??
           thvref(mg) = theta(mg,1)*(1.0 + 0.61*qtg(mg,1))

C heatv here is the virtual heat flux at the surface as said in
C Ncar tech report page 88.
C should I take the temperature at the surface instead 
C of theta(mg,1) for heatv ??
           heatv(mg)  = khfs(mg) + 0.61*theta(mg,1)*kqfs(mg)
           wm(mg)     = 0.
           therm      = 0.
           fak3(mg)   = 0.  
           zh(mg)     = 0.  
c          obklen at t point
           obklen(mg) = -thvref(mg)*ustar(mg)**3/
     &             (grav*vk*(heatv(mg) + sign(1.e-10,heatv(mg))))
C
C          Calculate virtual potential temperature first level
C          and initialize pbl height to z1 i.e  1st full level
C
           pblh(mg) = z(mg,1)     
           check(mg) = .true.
           rino(mg,1) = 0.0
       enddo !mg=1,ln2
C
C
C PBL height calculation:
C Search for level of pbl. Scan upward until the Richardson number between
C the first level and the current level exceeds the "critical" value.
C Richardson no. is computed using eq. (4.d.18) NCAR technical report, CCM3)
C
       do mg=1,ln2
           do k=2,nl
             if (check(mg)) then
              vvk = (u(mg,k) - u(mg,1))**2 
     &           + (v(mg,k) - v(mg,1))**2
     &           + fac*ustar(mg)**2
               vvk = max(vvk,tiny)
               tkv = theta(mg,k)*(1. + 0.61*qtg(mg,k))
               rino(mg,k) = grav*(tkv - thvref(mg))*(z(mg,k)-z(mg,1))
     &                 /(thvref(mg)*vvk)
               if(rino(mg,k).ge.ricr) then
                 pblh(mg) = z(mg,k-1) + (ricr - rino(mg,k-1))/
     &                 (rino(mg,k) - rino(mg,k-1))*(z(mg,k) - z(mg,k-1))
                 check(mg) = .false.
                 goto 100
               end if
             end if
           end do
100        continue
       enddo        !mg=1,ln2

C
C Set pbl height to maximum value where computation exceeds number of
C layers allowed
C
       do mg=1,ln2
           if (check(mg)) pblh(mg) = z(mg,nl)
       end do
C

C Improve estimate of pbl height for the unstable points.
C Find unstable points (virtual heat flux is positive):
C
       do mg=1,ln2
          if (heatv(mg) .gt. 0.) then
            unstbl(mg) = .true.
            check(mg) = .true.
          else
            unstbl(mg) = .false.
            check(mg) = .false.
          end if   
       end do
C
C For the unstable case, compute velocity scale and the
C convective temperature excess:
C
       do mg=1,ln2
          if (check(mg)) then
            phiminv(mg) = (1. - binm*pblh(mg)/obklen(mg))**onet
            wm(mg)= ustar(mg)*phiminv(mg)   !eqn (4.d.27)
C           therm: 2nd term in eq. (4.d.19): 
C           temperature excess due to convective thermal
            therm = heatv(mg)*fak/wm(mg)      
            rino(mg,1) = 0.0
C           eq. (4.d.19) : tlv then used in eq. (4.d.18) to improve
C           pblh
            tlv(mg) = thvref(mg) + therm
          end if
       end do
C
C Improve pblh estimate for unstable conditions using the
C convective temperature excess:
C
       do mg=1,ln2
          do k=2,nl
            if (check(mg)) then
              vvk = (u(mg,k) - u(mg,1))**2 
     &           + (v(mg,k) - v(mg,1))**2
     &           + fac*ustar(mg)**2
              vvk = max(vvk,tiny)
              tkv = theta(mg,k)*(1. + 0.61*qtg(mg,k))
              rino(mg,k) = grav*(tkv - tlv(mg))*(z(mg,k)-z(mg,1))
     &                 /(thvref(mg)*vvk)     ! mistake ???? (see (4.d.18)
              if(rino(mg,k).ge.ricr) then
                pblh(mg) = z(mg,k-1) + (ricr - rino(mg,k-1))/
     &                 (rino(mg,k) - rino(mg,k-1))*(z(mg,k) - z(mg,k-1))
                check(mg) = .false.
              end if
            end if
          end do
       end do

C
C Points for which pblh exceeds number of pbl layers allowed;
C set to maximum
C
       do mg=1,ln2
          if (check(mg)) pblh(mg) = z(mg,nl)
       end do
C
C PBL height must be greater than some minimum mechanical mixing depth
C Several investigators have proposed minimum mechanical mixing depth
C relationships as a function of the local friction velocity, u*.  We 
C make use of a linear relationship of the form h = c u* where c=700.
C The scaling arguments that give rise to this relationship most often 
C represent the coefficient c as some constant over the local coriolis
C parameter.  Here we make use of the experimental results of Koracin 
C and Berkowicz (1988) [BLM, Vol 43] for wich they recommend 0.07/f
C where f was evaluated at 39.5 N and 52 N.  Thus we use a typical mid
C latitude value for f so that c = 0.07/f = 700.
C
C      do mg=1,ln2
C         pblmin  = 700.0*ustar(mg)
C         pblh(mg) = max(pblh(mg),pblmin)
C      end do
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C pblh is now available; do preparation for diffusivity calculation:
C
       do mg=1,ln2
          fak1(mg) = ustar(mg)*pblh(mg)*vk
C
C Do additional preparation for unstable cases only, set temperature
C and moisture perturbations depending on stability.
C
          if (unstbl(mg)) then
            phiminv(mg) = (1. - binm*pblh(mg)/obklen(mg))**onet
            phihinv(mg) = sqrt(1. - binh*pblh(mg)/obklen(mg))
            wm(mg)      = ustar(mg)*phiminv(mg)
            wstr(mg)    = (heatv(mg)*grav*pblh(mg)/thvref(mg))**onet 
            fak3(mg)    = fakn*wstr(mg)/wm(mg)
            fak2(mg)    = wm(mg)*pblh(mg)*vk
          end if
       end do
C
C Main level loop to compute the diffusivities and 
C counter-gradient terms:
C
       do 1000 k=1,nl-1
C
C Find levels within boundary layer:
C
           do mg=1,ln2
              unslev(mg) = .false.
              stblev(mg) = .false.
              zm(mg) = z(mg,k)
              zp(mg) = z(mg,k+1)
              if ( zp(mg).gt.pblh(mg))zp(mg) = pblh(mg)
              if (zm(mg) .lt. pblh(mg)) then
                zmzp = 0.5*(zm(mg) + zp(mg))
                zh(mg) = zmzp/pblh(mg)
                zl(mg) = zmzp/obklen(mg)
                zzh(mg) = 0.
                if (zh(mg).le.1.0) zzh(mg) = (1. - zh(mg))**2
C
C stblev for points zm < plbh and stable and neutral
C unslev for points zm < plbh and unstable
C
                if (unstbl(mg)) then
                  unslev(mg) = .true.
                else
                  stblev(mg) = .true.
                end if
              end if
           end do
C
C Stable and neutral points; set diffusivities; counter-gradient
C terms zero for stable case:
C term: pblk is Kc in eq. (4.d.16)
C
           do mg=1,ln2
             if (stblev(mg)) then
               if (zl(mg).le.1.) then   ! 0 < z/L < 1.
                 pblk(mg) = fak1(mg)*zh(mg)*zzh(mg)/(1. + betas*zl(mg))
               else                    !  z/L > 1.
                 pblk(mg) = fak1(mg)*zh(mg)*zzh(mg)/(betas + zl(mg))
               end if
               rkm(mg,k) = max(pblk(mg),rkm(mg,k))
               rkh(mg,k) = rkm(mg,k)
             end if
           end do
C
C unssrf, unstable within surface layer of pbl
C unsout, unstable within outer   layer of pbl
C
           do mg=1,ln2
             unssrf(mg) = .false.
             unsout(mg) = .false.
             if (unslev(mg)) then
               if (zh(mg).lt.sffrac) then
                 unssrf(mg) = .true.
               else
                 unsout(mg) = .true.
               end if
             end if
           end do
C
C Unstable for surface layer; counter-gradient terms zero
C

           do mg=1,ln2
            if (unssrf(mg)) then
              term = (1. - betam*zl(mg))**onet
              pblk(mg) = fak1(mg)*zh(mg)*zzh(mg)*term
              pr(mg) = term/sqrt(1. - betah*zl(mg))
            end if
           end do
C
C Unstable for outer layer; counter-gradient terms non-zero:
C
           do mg=1,ln2
             if (unsout(mg)) then
               pblk(mg) = fak2(mg)*zh(mg)*zzh(mg)
               cgs(mg,k) = fak3(mg)/(pblh(mg)*wm(mg))
               khfv=khfs(mg) + 0.61*theta(mg,1)*kqfs(mg)
               cgh(mg,k) = khfv*cgs(mg,k)               !eq. (4.d.17)
               cgq(mg,k) = kqfs(mg)*cgs(mg,k)           !eq. (4.d.17)
               pr(mg) = phiminv(mg)/phihinv(mg) 
     &                + ccon*fak3(mg)/fak
             end if
           end do
           if(ne.eq.1)then
             do mg=1,ln2
               if (unsout(mg)) then
                 cgx(mg,k,:)=kxfs(mg,:)*cgs(mg,k)
               endif
             enddo
           endif
C
C For all unstable layers, set diffusivities
C
           do mg=1,ln2
             if (unslev(mg)) then
               rkm(mg,k) = max(pblk(mg),rkm(mg,k))
               rkh(mg,k) = max(pblk(mg)/pr(mg),rkh(mg,k))
             end if
           end do
1000   continue      !end of k loop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     update theta and moisture (qtg) due to counter gradient 
c     NCAR techique note p92. Recommended that should not
c     include counter gradient term in implicit solution

c      tdt = 2 x dt
c
       ztodtgor = tdt*grav/rdry

C      update theta and qtg due to counter gradient
       do mg=1,ln2
          do k=2,nl-1
            delsig=(prh(mg,k+1)-prh(mg,k))/pg(mg)
            tmp1 = ztodtgor/delsig

            sigmh=prh(mg,k+1)/pg(mg)
            sigmhk=prh(mg,k)/pg(mg)

            sigotbk=sigmh/(0.5*(ttg(mg,k+1) + ttg(mg,k)))
            sigotbkm1=sigmhk/(0.5*(ttg(mg,k-1) + ttg(mg,k)))

            theta(mg,k)=theta(mg,k) + tmp1*
     $           (sigotbk*rkh(mg,k)*cgh(mg,k) -
     $            sigotbkm1*rkh(mg,k-1)*cgh(mg,k-1))
            qtg(mg,k) = qgx(mg,k) + tmp1 *
     $           (sigotbk*rkh(mg,k)*cgq(mg,k) -
     $            sigotbkm1*rkh(mg,k-1)*cgq(mg,k-1))

            if(ne.eq.1)then
              xtg(mg,k,:) = xtx(mg,k,:) + tmp1 *
     $           (sigotbk*rkh(mg,k)*cgx(mg,k,:) -
     $            sigotbkm1*rkh(mg,k-1)*cgx(mg,k-1,:))
            endif              
          end do
          k=1
          delsig = (prh(mg,k+1)-prh(mg,k))/pg(mg)
          tmp1 = ztodtgor/delsig
          sigotbk=(prh(mg,k+1)/pg(mg))/(0.5*(ttg(mg,k+1) + ttg(mg,k)))
          theta(mg,k)=theta(mg,k) + tmp1*
     $           (sigotbk*rkh(mg,k)*cgh(mg,k))
          qtg(mg,k) = qgx(mg,k) + tmp1 *
     $           (sigotbk*rkh(mg,k)*cgq(mg,k))
          if(ne.eq.1)then
              xtg(mg,k,:) = xtx(mg,k,:) + tmp1 *
     $           (sigotbk*rkh(mg,k)*cgx(mg,k,:))
          endif
       end do

C
C      Check for neg qtg's and put the original vertical
C      profile back if a neg value is found. A neg value implies that the
C      quasi-equilibrium conditions assumed for the countergradient term are
C      strongly violated.
C      Original code rewritten by Rosinski 7/8/91 to vectorize in longitude.

       do mg=1,ln2
         do k=1,nl
            if (qtg(mg,k) .lt. 0) then
             do kk=1,nl
               qtg(mg,k)=qgx(mg,k)
             enddo
             goto 300
            endif
300         continue
         enddo
       enddo

c Do the same for tracer mixing ratios.

       if(ne.eq.1)then
         do nt=1,ntrac
           do mg=1,ln2
             do k=1,nl
               if (xtg(mg,k,nt) .lt. 0)then
                 do kk=1,nl
                   xtg(mg,k,nt)=xtx(mg,k,nt)
                 enddo
                 goto 301
               endif
301            continue
             enddo
           enddo
         enddo
       endif

      if(debug)then
        if(lg.eq.lgdebug)then
          mg=mgdebug+(insdebug-1)*lon
          write(25,'(a,i1)')'After ncarpbl'
          write(25,1)'kxfs(2)',kxfs(mg,2),' xtem(2) ',xtem(mg,2)
          write(25,9)'rkh ',(rkh(mg,k),k=1,nl)
          write(25,9)'cgx ',(cgx(mg,k,2),k=1,nl)
          write(25,9)'cgs ',(cgs(mg,k),k=1,nl)
        endif
      endif
 1    format(3(a,g10.3))
 9    format(a,30g10.3)

       return
       end
