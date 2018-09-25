c Extracted the COMMON block /giant4/ to a header file.
c SJP 2009/03/12
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Removed DEFAULT(NONE) clauses from the OpenMP directives, in order to fix
c problems on both the SC and blizzard.
c SJP 2004/03/16
c
c Modified to make use of DYNAMIC scheduling, which leads to a considerable
c improvement in performance
c SJP 2003/03/29
c
c Loops re-rolled so that code passes array-bounds checking.
c SJP 2003/03/08
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c Changes to add machine type 'ALPH', which makes use of the super-fast FFTW
c FFT library.
c SJP 2001/11/22
c
c $Log: dynm.f,v $
c Revision 1.34  2001/02/28 04:36:40  rot032
c Further tidy ups from HBG
c
c Revision 1.33  2001/02/22 05:34:44  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.32  2000/08/16 02:59:26  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.31  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.30  1999/05/20 06:23:56  rot032
c HBG changes to V5-2
c
c Revision 1.29  1997/12/23  00:23:38  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.28  1997/12/17  23:22:57  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.27.1.1  1997/12/19  02:03:18  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.27  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.26  1997/03/06  04:35:28  mrd
c Force appropriate R or T truncation in call to ftospec, no matter how
c LLMAX is defined.
c
c Revision 1.25  1996/11/19  06:05:53  ldr
c Parallel Legendre transforms for kaos.
c
c Revision 1.24  1996/10/24  01:27:40  ldr
c Merge of TIE changes with LDR/MRD changes since V5-0-12.
c
c Revision 1.23  1996/08/12  01:51:22  mrd
c Generalise for triangular truncations other than T63
c
c Revision 1.22.1.1  1996/10/24  01:02:38  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.22  1996/06/13  02:06:17  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.21  1996/03/21  03:18:36  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.20  1995/09/07  06:38:04  ldr
c Corrected padding on /legnd/.
c
c Revision 1.19  1995/07/05  07:47:32  ldr
c Removed spurious SGI mp_ stuff, to avoid unresolved externals on Cray.
c
c Revision 1.18  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c     INPUT/OUTPUT:
c     Input:   from common/fewflags in FEWFLAGS.f
c                  ifwds - default F, only T when doing forward step from
c                          new restart file
c                  impvor- if true, use implicit treatment of vorticity
c                          equation, default T
c
c              from common/gausl in GAUSL.f
c                  w    - normal guass weights
c                  wocs - gauss weights/cos(lat)**2
c
c              from common/glegnd in this subroutine
c                  plmg  - global array of Legendre polynomials
c                  cplmg - global array of gradients of above
c
c              from common/glegndr in this subroutine
c                  cplmr) rotated versions of plmg and cplmg, used for
c                  plmr ) vectorisation on the Cray
c
c              from common/masiv1 in this subroutine
c                  ronmx - global array of measure of east-west surface stress
c                  sonmx - global array of measure of north-sth surface stress
c                  iphys ) - current phys time step and previous phys time
c                  iphysm)   step, these indices swap around
c
c              from common/masiv4 in MASIV4.f
c                  savegrid - global surface variables
c
c              from common/ramp in this subroutine
c                  rampm - ramp in m for all k levels
c
c              from common/rmgrid in RMGRID.f
c                  rmg - pressure weighted moistures at current timestep
c
c              from common/timex in TIMEX.f
c                  mstep - time step in minutes
c
c              from common/uvpgd in UVPGD.f
c                  pgd - global array of pressure
c                  ugd - global array of zonal velocity
c                  vgd - global array of meridional velocity
c
c              from common/worknsd in this subroutine
c                  real and imaginary Fourier coefficients, n=northern
c                  hemisphere, s=southern hemisphere
c
c                  afni, afnr, - vorticity and divergence
c                  bfni, bfnr, - vorticity and divergence
c                  efni, efnr, - energy
c                  ffni, ffnr, - horizontal advection of temp
c                  gfni, gfnr,
c                  hfni, hfnr, - vertical motion/heating term
c                                           for temperature
c
c              from common/worknsvo in this subroutine
c                  real and imaginary Fourier coefficients for implicit
c                  vorticity version, n=northern hemisphere, s=southern hem.
c
c                  apni, apnr, bpni, bpnr
c
c     Output:  from common/bgrnd in BGRND.f
c                  kdynm - counts no of 1/4 days of dynamical statistical
c                          data that have been saved
c
c              from common/legnd in this subroutine
c                  cplm - gradients of plm, where plm holds the Legendre
c                         polys for (lw,mw)
c
c              from common/giant4 in this subroutine
c                  z4 - grid elevations [non-spectral]
c                  ron - measure of east-west stress at surface
c                  son - measure of north-south stress at surface
c
c              from common/work1x in this subroutine
c                  pn - surface pressure
c                  un - zonal wind velocity m/s
c                  vn - meridional wind velocity m/s
c
c              from common/glmean in GLEAN.f
c                  pslbar - global mean sea level pressure
c
c     In/Out:  from common/giant1 in this subroutine
c                  aon, bon, eon, fon, gon, hon - grid values of A, B, E,
c                  F, G, H  before Fast Fourier Transforms
c
c              from common/legnd in this subroutine
c                  plm - Legendre polynomials for (lw,mw)
c
c              from common/timex in TIMEX.f
c                  mins - current model time in mins
c                  minw - minutes between dynamical field prints in DYNFP=T
c
c              from common/work1x in this subroutine
c                  rmn - pressure weighted grid point moisture
c
c              from common/workrid in WORKRID.f
c                  deli, delr, dipi, dipr, diti, ditr, dixi, dixr - rotated
c                  versions of eli, elr, ipi, ipr, iti, itr, ixi, ixr, used
c                  for efficient vectorisation on the Cray
c
c              from arguments
c                  tdt - 2 X timestep (dt)
c

*VOCL TOTAL,REPEAT(999999)

      subroutine dynm(tdt)

      implicit none

!$OMP THREADPRIVATE ( /GIANT1/ )
!$OMP THREADPRIVATE ( /GIANT4/ )
!$OMP THREADPRIVATE ( /LEGND/ )
!$OMP THREADPRIVATE ( /WORK1X/ )
!$OMP THREADPRIVATE ( /WORKNSD/ )

C Global parameters
      include 'PARAMS.f'

C Argument list
      real tdt

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      real aon,bon,eon,fon,gon,hon
      common/giant1/aon(ln2,nl)
     & ,bon(ln2,nl),eon(ln2,nl),fon(ln2,nl),gon(ln2,nl),hon(ln2,nl)

      include 'GIANT4.f'

      real plmx,plm,cplm,pad
      common/legnd/plmx(lw1,mw),plm(lw,mw),cplm(lw,mw),pad(2,mw)

      real un,vn,pn,dvn,von,ten,cpn,pln,rmn
      common/work1x/un(ln2,nl),vn(ln2,nl),pn(ln2),dvn(ln2,nl)
     & ,von(ln2,nl),ten(ln2,nl),cpn(ln2),pln(ln2),rmn(ln2,nl)

      real afnr,bfnr,efnr,ffnr,gfnr,hfnr,afni,bfni,efni,ffni
     & ,gfni,hfni,apnr,bpnr,apni,bpni
      common/worknsd/
     &  afnr(mw,nl,2),bfnr(mw,nl,2),efnr(mw,nl,2),ffnr(mw,nl,2)
     & ,gfnr(mw,nl,2),hfnr(mw,nl,2)
     & ,afni(mw,nl,2),bfni(mw,nl,2),efni(mw,nl,2),ffni(mw,nl,2)
     & ,gfni(mw,nl,2),hfni(mw,nl,2)
     & ,apnr(mw,nl,2),bpnr(mw,nl,2)
     & ,apni(mw,nl,2),bpni(mw,nl,2)

C Global data blocks
      include 'BGRND.f'
      include 'FEWFLAGS.f'
      include 'GAUSL.f'
      include 'GLMEAN.f'
      include 'MASIV4.f'
      include 'RMGRID.f'
      include 'TIMEX.f'
      include 'UVPGD.f'
      include 'WORKRID.f'
      include 'ZMEAN.f'

      integer lwh
      parameter (lwh=(lw+1)/2)
      real plmg,cplmg,plmgo,plmge,cplmgo,cplmge,rampm
      common/glegnd/plmg(lw1,mw,lat),cplmg(lw,mw,lat)
     &           ,plmgo(lat,lwh,mw),plmge(lat,lwh,mw)
     &           ,cplmgo(lat,lwh,mw),cplmge(lat,lwh,mw)
     &           ,rampm(mw,nl)

      real ronmx,sonmx
      integer iphysm,iphys
      common/masiv1/ronmx(ln2,nl,lat,2),sonmx(ln2,nl,lat,2)
     &             ,iphysm,iphys

C Local work arrays and variables
c-- Dimension following at mw+1 to prevent Bank Conflicts at T63
      real efor(mw+1,nl,lat),efoi(mw+1,nl,lat),efer(mw+1,nl,lat),
     &     efei(mw+1,nl,lat),gfor(mw+1,nl,lat),gfoi(mw+1,nl,lat),
     &     gfer(mw+1,nl,lat),gfei(mw+1,nl,lat),ffor(mw+1,nl,lat),
     &     ffer(mw+1,nl,lat),hfoi(mw+1,nl,lat),hfei(mw+1,nl,lat),
     &     pbor(mw+1,nl,lat),pboi(mw+1,nl,lat),pber(mw+1,nl,lat),
     &     pbei(mw+1,nl,lat),xaor(mw+1,nl,lat),xaoi(mw+1,nl,lat),
     &     xaer(mw+1,nl,lat),xaei(mw+1,nl,lat),paor(mw+1,nl,lat),
     &     paoi(mw+1,nl,lat),paer(mw+1,nl,lat),paei(mw+1,nl,lat),
     &     xbor(mw+1,nl,lat),xboi(mw+1,nl,lat),xber(mw+1,nl,lat),
     &     xbei(mw+1,nl,lat) 
      real zonpsl(2)
      logical lstat

      integer iphyst
      integer istday
      integer istd4m1
      integer k
      integer lg
      integer lgn
      integer ll
      integer ll_lim
      integer mg
c     integer mm ! defined in LLMAX.f
      integer nex
      integer i, j

      real x2
      real x2mm1
      real x3

C Local data, functions etc
      include 'LLMAX.f'

C Start code : ----------------------------------------------------------

c Multiprocessor version of DYNM

*PDIR SERIAL
c       lstat=mod(mins+mstep,360).eq.0 ! stats each 1/4 day
c Following gets stats more often than 1/4 day (avoids 1/4 sampling bias)
      istday=1440/mstep
      istd4m1=int(istday/4.0-1.0)
        lstat=mod(nsteps+1,istd4m1).eq.0
      if(lstat)kdynm=kdynm+1
      iphyst=iphysm
      if(ifwds)iphyst=iphys
*PDIR ENDSERIAL

!$OMP  PARALLEL
!$OMP& PRIVATE (i, j, k, lg, lgn, ll, ll_lim, mg, mm, nex, x2, x2mm1,
!$OMP&          x3, zonpsl)
!$OMP& SHARED  (cplmg, cplmge, cplmgo, efei, efer, efoi, efor, eli,
!$OMP&          elr, ffer, ffor, gfei, gfer, gfoi, gfor, hfei, hfoi,
!$OMP&          impvor, iphyst, ipi, ipr, iti, itr, ixi, ixr, lstat,
!$OMP&          paei, paer, paoi, paor, pbei, pber, pboi,
!$OMP&          pbor, pgd, plmg, plmge, plmgo, rampm, rmg, ronmx,
!$OMP&          savegrid, sonmx, tdt, ugd, vgd, w, wocs, xaei, xaer,
!$OMP&          xaoi, xaor, xbei, xber, xboi, xbor, zpslk)

C**** 1111 : GAUSS INTEGRAL LOOP FOR EACH G-LAT CIRCLE
C     NORTH G-LAT FIRST,CORRESPONDING SOUTH G-LAT SECOND

!$OMP DO SCHEDULE(DYNAMIC)

      do 1111 lg=1,lat
C***  THE NORTH & SOUTH GRID POINT VALUES ARE COMPUTED TOGETHER
c----
c     extract the appropriate legndre polys and store in plm()
c     for (lw,mw).
      do 33 mm=1,mw 
      do 33 ll=1,llmax(mm)
   33 plm(ll,mm) = plmg(ll,mm,lg)
      if(trunc.eq.'T')then
         do mm=2,mw
            plm(mm,lw+2-mm)=0.0
         end do
      end if
c     extract cos(lat)d(plm)/d(lat) and store in cplm()
      do 34 mm=1,mw
      do 34 ll=1,llmax(mm)
   34 cplm(ll,mm)=cplmg(ll,mm,lg)

C**** READ T-1 FRICTION & Z* FROM STORAGE ARRAYS
      do 37 k=1,nl
      do 37 mg=1,ln2
      ron(mg,k)=ronmx(mg,k,lg,iphyst)
   37 son(mg,k)=sonmx(mg,k,lg,iphyst)
      do 38 mg=1,ln2
   38 z4(mg)=savegrid(mg,1,lg)

C**** TRANSFORM SPECTRAL QUANTITIES (& SOME GRADIENTS) INTO FOURIER
C**   SPACE THEN ONTO LATITUDE GRID BY FFT,DONE IN SUB DTOG
C     CALL TIMER('DTOG    ',3)
      if(mw.eq.64)then
        call dtogcr63    !T63 version (optimized for Cray)
      else
        call dtogcray    !R21/R42 version (optimized for Cray)
      endif
C     CALL TIMER('DTOG    ',4)

C**** GET UN,VN,PN
      lgn=lat2p-lg
      do 40 k=1,nl
      do 40 mg=1,lon
      un(mg,k)=ugd(mg,lgn,k)
      vn(mg,k)=vgd(mg,lgn,k)
      un(mg+lon,k)=ugd(mg,lg,k)
   40 vn(mg+lon,k)=vgd(mg,lg,k)
      do 42 mg=1,lon
      pn(mg)=pgd(mg,lgn)
   42 pn(mg+lon)=pgd(mg,lg)

      do 44 k=1,nl
      do 44 mg=1,ln2
   44 rmn(mg,k)=rmg(mg,k,lg)

C**   FOR EACH LAT GRID POINT COMPUTE THE NON-LINEAR TERMS,STORE IN
C     AON,BON, ETC. TRANSFORM INTO FOURIER SPACE,THEN PERFORM THE
C     SEQUENTIAL ADDITIONS TO THE GAUSS INTEGRALS IP,IX,IT ETC.
C**
C**** FORM THE NON-LINEAR PRODUCTS AT THE GRID POINTS.
C     CALL TIMER('DYNMNL  ',3)
      call dynmnl(lg,lstat,zonpsl,lgn)
C     CALL TIMER('DYNMNL  ',4)
      if(lstat)call dynmst(lg)
      do k=1,2
        zpslk(lg,k)=zonpsl(k)
      enddo

C*    TRANSFORM THE GRID POINT NON-LINEAR TERMS INTO FOURIER SPACE
C     CALL TIMER('DY MFFTM',3)
c   Dynamics FFT transforming 6 variables (6*nl) and for NH,SH at
c    same time : (6*nl)*2
c   Next argument (2) => Dynamics variables
c
c   The input data is stacked (aon,bon etc) in a common block (giant1)
c    to process as many "levels" of data at one time.
c   The return data (afnr,afni etc) is also stacked in a common block
c    (worknsd). However, to facilitate vectorizing below, a change of
c    position of the last two arguments in (afnr,afni) relative to (aon)
c    etc is required. The change from [aon(ln2,nl) to afnr(mw,nl,2) etc]
c    is achieved through precomputed pointer "igtos_d" : 
c    see gauleg.f and mfftm.f

      nex=(6*nl)*2

CSJP  Former machine dependence at this point
      call mfftma(aon, afnr, afni, nex, 2)

C     CALL TIMER('DY MFFTM',4)

      if(impvor)call dynmvo(lg,tdt)

C*    PERFORM THE SEQUENTIAL ADDITIONS TO THE GAUSS INTEGRAL ARRAYS
C*    (SPLIT INTO REAL & IMAGINARY PARTS).
C**   MULTIPLY BY THE CORRESPONDING GAUSS INTEGRAL WEIGHTS
C*    W() ARE NORMAL GAUSS WEIGHTS,WOCS() ARE GAUSS WEIGHTS/COS(LAT)**2
C     CALL TIMER('DY SPADD',3)
      x2=wocs(lg)
      x3=w(lg)
      do 68 k=1,nl
      do 68 mm=1,mw
        x2mm1=x2*rampm(mm,k)
          efor(mm,k,lg)=(efnr(mm,k,1)+efnr(mm,k,2))*x3 
          efoi(mm,k,lg)=(efni(mm,k,1)+efni(mm,k,2))*x3 
          efer(mm,k,lg)=(efnr(mm,k,1)-efnr(mm,k,2))*x3 
          efei(mm,k,lg)=(efni(mm,k,1)-efni(mm,k,2))*x3 
          gfor(mm,k,lg)=(gfnr(mm,k,1)-gfnr(mm,k,2))*x2 
          gfoi(mm,k,lg)=(gfni(mm,k,1)-gfni(mm,k,2))*x2 
          gfer(mm,k,lg)=(gfnr(mm,k,1)+gfnr(mm,k,2))*x2 
          gfei(mm,k,lg)=(gfni(mm,k,1)+gfni(mm,k,2))*x2 
          ffor(mm,k,lg)=(ffni(mm,k,1)+ffni(mm,k,2))*x2mm1
     &                 +(hfnr(mm,k,1)+hfnr(mm,k,2))*x3 
          ffer(mm,k,lg)=(ffni(mm,k,1)-ffni(mm,k,2))*x2mm1
     &                 +(hfnr(mm,k,1)-hfnr(mm,k,2))*x3 
          hfoi(mm,k,lg)=(hfni(mm,k,1)+hfni(mm,k,2))*x3 
     &                 -(ffnr(mm,k,1)+ffnr(mm,k,2))*x2mm1
          hfei(mm,k,lg)=(hfni(mm,k,1)-hfni(mm,k,2))*x3 
     &                 -(ffnr(mm,k,1)-ffnr(mm,k,2))*x2mm1
          xaor(mm,k,lg)=(afnr(mm,k,1)-afnr(mm,k,2))*x2 
          xaoi(mm,k,lg)=(afni(mm,k,1)-afni(mm,k,2))*x2 
          xaer(mm,k,lg)=(afnr(mm,k,1)+afnr(mm,k,2))*x2 
          xaei(mm,k,lg)=(afni(mm,k,1)+afni(mm,k,2))*x2 
          xbor(mm,k,lg)=(bfni(mm,k,1)+bfni(mm,k,2))*x2mm1
          xboi(mm,k,lg)=(bfnr(mm,k,1)+bfnr(mm,k,2))*x2mm1
          xber(mm,k,lg)=(bfni(mm,k,1)-bfni(mm,k,2))*x2mm1
   68     xbei(mm,k,lg)=(bfnr(mm,k,1)-bfnr(mm,k,2))*x2mm1

      If(impvor)Then

        do 682 k=1,nl
        do 682 mm=1,mw
          x2mm1=x2*rampm(mm,k)
          pbor(mm,k,lg)=(bpnr(mm,k,1)-bpnr(mm,k,2))*x2 
          pboi(mm,k,lg)=(bpni(mm,k,1)-bpni(mm,k,2))*x2 
          pber(mm,k,lg)=(bpnr(mm,k,1)+bpnr(mm,k,2))*x2 
          pbei(mm,k,lg)=(bpni(mm,k,1)+bpni(mm,k,2))*x2 
          paor(mm,k,lg)=(apni(mm,k,1)+apni(mm,k,2))*x2mm1
          paoi(mm,k,lg)=(apnr(mm,k,1)+apnr(mm,k,2))*x2mm1
          paer(mm,k,lg)=(apni(mm,k,1)-apni(mm,k,2))*x2mm1
  682     paei(mm,k,lg)=(apnr(mm,k,1)-apnr(mm,k,2))*x2mm1

      Else

        do 684 k=1,nl
        do 684 mm=1,mw
          x2mm1=x2*rampm(mm,k)
          pbor(mm,k,lg)=(bfnr(mm,k,1)-bfnr(mm,k,2))*x2 
          pboi(mm,k,lg)=(bfni(mm,k,1)-bfni(mm,k,2))*x2 
          pber(mm,k,lg)=(bfnr(mm,k,1)+bfnr(mm,k,2))*x2 
          pbei(mm,k,lg)=(bfni(mm,k,1)+bfni(mm,k,2))*x2 
          paor(mm,k,lg)=(afni(mm,k,1)+afni(mm,k,2))*x2mm1
          paoi(mm,k,lg)=(afnr(mm,k,1)+afnr(mm,k,2))*x2mm1
          paer(mm,k,lg)=(afni(mm,k,1)-afni(mm,k,2))*x2mm1
  684     paei(mm,k,lg)=(afnr(mm,k,1)-afnr(mm,k,2))*x2mm1
 
      Endif
 1111 continue

!$OMP END DO

C**** END OF DYNAMICS LOOP.

c Zero elr,eli,itr,iti,ipr,ipi,ixr,ixi arrays

!$OMP DO SCHEDULE(DYNAMIC)

      do 82 k=1,nl
CSJP        do ll=1,lw*mw
CSJP          elr(ll,1,k)=0.
CSJP          eli(ll,1,k)=0.
CSJP          itr(ll,1,k)=0.
CSJP          iti(ll,1,k)=0.
CSJP          ipr(ll,1,k)=0.
CSJP          ipi(ll,1,k)=0.
CSJP          ixr(ll,1,k)=0.
CSJP          ixi(ll,1,k)=0.
CSJP        enddo
        do 82 j = 1, mw
          do 82 i = 1, lw
            elr(i, j, k) = 0.0
            eli(i, j, k) = 0.0
            itr(i, j, k) = 0.0
            iti(i, j, k) = 0.0
            ipr(i, j, k) = 0.0
            ipi(i, j, k) = 0.0
            ixr(i, j, k) = 0.0
            ixi(i, j, k) = 0.0
 82   continue

!$OMP END DO

c Do the Legendre transform, parallelizing over mm
c
c In the following subroutine calls, the last two arguments are the output, 
c being the returned spectral variables, real and imaginary.  All the other 
c arguments are input.

c In a triangular model the truncation is enforced here, no matter how
c LLMAX is defined. This is necessary to maintain the truncation because
c plmg has an extra meridional wave

!$OMP DO SCHEDULE(DYNAMIC)

      do mm=1,mw

         ll_lim = lw1-mm
         if ( trunc .eq. 'R' ) ll_lim = lw

         call ftospec(ll_lim,mm,
     &            plmgo(1,1,mm),plmge(1,1,mm),
     &            efor,efoi,efer,efei,elr,eli)
         call ftospec3(ll_lim,mm,
     &            plmgo(1,1,mm),plmge(1,1,mm),
     &            cplmgo(1,1,mm),cplmge(1,1,mm),
     &            ffor,gfor,hfoi,gfoi,ffer,gfer,hfei,gfei,itr,iti,
     &            paor,pbor,paoi,pboi,paer,pber,paei,pbei,ipr,ipi,
     &            xbor,xaor,xboi,xaoi,xber,xaer,xbei,xaei,ixr,ixi)
      enddo

!$OMP END DO

!$OMP END PARALLEL

*PDIR SERIAL
      pslbar=0.0
      do lg = 1, lat
        x3=0.5*w(lg)/lon
	pslbar = pslbar+(zpslk(lg,1)+zpslk(lg,2))*x3
      enddo
      iphyst=iphysm
      iphysm=iphys
      iphys=iphyst
*PDIR ENDSERIAL

      return
      end
