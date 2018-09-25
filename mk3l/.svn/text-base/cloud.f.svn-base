c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/10
c
c $Log: cloud.f,v $
c Revision 1.34  2001/02/28 04:36:36  rot032
c Further tidy ups from HBG
c
c Revision 1.33  2001/02/22 05:34:37  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.32  2001/02/12 05:39:44  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.31  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.30  1998/12/10  00:55:30  ldr
c HBG changes to V5-1-21
c
c Revision 1.29  1997/12/17  23:22:44  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.28  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.27  1997/05/19  07:35:21  ldr
c Implement zmean cfrac diagnostic for dcloud scheme.
c
c Revision 1.26  1997/02/21  00:26:10  ldr
c Go back to rcrits/l=0.85/0.75 for 18L version, tidy cloud2.f and make
c clddia.f, cloud.f general for NL.
c
c Revision 1.25  1996/10/24  01:02:31  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.24  1996/06/13  02:05:47  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.23  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.22  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.19.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.21  1994/08/08  17:20:53  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.20  94/08/04  16:54:22  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.19  94/04/29  15:10:28  ldr
c Changes from HBG to reduce high cirrus, increase convective cirrus to
c previous value (2*cvx), to improve jet again. Also change low cloud and
c Cd over sea to achieve balance.
c 
c Revision 1.18  93/12/17  15:31:52  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.17  93/10/11  14:22:19  ldr
c Pass clcon (conv cloud) back up to radin to add to diagnosed cloud for
c qcloud scheme.
c 
c Revision 1.16  93/10/05  13:05:28  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.15  93/08/17  16:27:03  ldr
c New albedos etc. to balance R21 version (HBG).
c 
c Revision 1.14  93/08/03  11:25:05  ldr
c Replace intrinsic functions with generic forms where possible.
c 
c Revision 1.13  93/07/06  16:42:54  ldr
c New cloud albedos to go woth changes in clddia (HBG).
c 
c Revision 1.12  93/06/16  14:02:38  ldr
c HBG: Cloud tuning (changed IR absorptivity).
c
c     INPUT/OUTPUT
c     Input:   from common/radnml in this subroutine
c                  dcflag - if T use diagnostic clouds, if F use
c                           fixed clouds, default F
c                  nrad - no timesteps between calls to the radiation scheme,
c                         set to zero for no radiation, default 4
c
c              from common/uvpgd in UVPGD.f
c                  sdot - d(sigma)/d(t)
c
c              from arguments
c                  cldoff - no cloud for clear sky radiation calculations
c                  ipass  - model pass counter
c                  rhg    - relative humidity
c
c     In/Out:  from arguments
c                  alatn/s  - latitude (from -pi/2 to pi/2)
c                  bvnf  - Brunt-Vaisala frequency, full levels
c                  mins  - current model time in mins
c                  clat  - cloud amount diagnostic
c                  clcon - convective cloud amount
c                  clh   - high level cloud    cll - low level cloud
c                  clm   - mid level cloud     ins - hemisphere index
c                  lg    - latitude index  
c
c 
      subroutine cloud(ipass,lg,alatn,alats
     & ,cll,clm,clh,cldoff,rhg,clcon,bvnf,clat)

      implicit none

!$OMP THREADPRIVATE ( /RADISW/ )

C Global parameters
      include 'HCON.f'
      include 'PARAMS.f'
      include 'RDPARM.f'

C Argument list
      integer ipass
      integer lg
      real alatn
      real alats
      real cll(imax)
      real clm(imax)
      real clh(imax)
      logical cldoff
      real rhg(imax,l)
      real clcon(imax,l)
      real bvnf(imax,nl)
      real clat(imax,l)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'RADISW.f'

C Global data blocks
      include 'FEWFLAGS.f'
      include 'TIMEX.f'
      include 'UVPGD.f'

      integer icln,ipcln
      real avcvrn,pavcvrn
      common/newcl/icln(ln2,2,lat),avcvrn(ln2,lat)
     &           ,ipcln(ln2,2,lat),pavcvrn(ln2,lat)

C Local work arrays and variables
      real sigdt(imax,l)
      real rhg1(imax,l)

      integer ich(imax), icm(imax), ict(imax), icb(imax)

      integer i
      integer ip
      integer is
      integer k
      integer kb
      integer kt
      integer lgn
      integer nc
      integer ncl
      integer nclds0
      integer nradcv

      real cvrn
      real pgx
      real refac

C Local data, functions etc
      real emch,emcm,emcl
      save emch,emcm,emcl
      data emch,emcm,emcl/1.,1.,1./
CXXXX NOTE : coca,cwca,cwcb defined for number of distinct layers
CXXXX  of cloud allowed. Currently three layers of cloud :
CXXXX    1 unused layer + 3 cloud layers (high,mid,low) + 1 extra = 5)
      real coca(5),cwca(5),cwcb(5)
      save coca,cwca,cwcb
c     data coca/0.,.210,.540,.660,.100/ original
c     data cwca/0.,.190,.460,.500,0./ original
c     data coca/0.,.200,.510,.615,.100/ to Dec 30,1990
c     data cwca/0.,.180,.430,.460,0./ to Dec 30 1990
c     data coca/0.,.205,.525,.640,.100/ to Feb 1 1991 
c     data cwca/0.,.185,.445,.480,0./ to Feb1 1991  
c     data coca/0.,.207,.534,.653,.100/ as han run, to 2907.hbd
c     data coca/0.,.202,.529,.647,.100/                        
c     data cwca/0.,.183,.451,.490,0./          
c     data coca/0.,.2048,.5265,.6435,.100/                        
c     data cwca/0.,.1853,.4485,.4875,0./          
CH    data coca/0.,.2030,.5230,.6380,.100/                        
CH    data cwca/0.,.1840,.4450,.4840,0./          
CX    data coca/0.,.2030,.5000,.6000,.100/                        
CX    data cwca/0.,.1840,.4200,.5000,0./          
CX    data coca/0.,.2030,.5500,.6600,.100/                        
CX    data cwca/0.,.1840,.4700,.5200,0./          
c     data coca/0.,.2030,.5500,.6400,.100/                        
C     data coca/0.,.2030,.5500,.6000,.100/                        
C     data cwca/0.,.1840,.4700,.5000,0./          
cx    data coca/0.,.2030,.5100,.5700,.100/                        
cx    data cwca/0.,.1840,.4300,.4700,0./          
      data coca/0.,.2030,.5500,.6400,.100/                        
      data cwca/0.,.1840,.4700,.5000,0./          
      data cwcb/0.,.040,.150,.300,0./

C Start code : ----------------------------------------------------------

      If (cldoff) Then

         do 10 ip=1,imax
            cll(ip)=0.
            clm(ip)=0.
            clh(ip)=0.
            nclds(ip)=0
 10      continue

      Else

c----
c---- Either compute cloud (dcflag=.true.) or use obs
c----
            do 20 k=1,l
               do 20 ip=1,imax
               clat(ip,k)=0.0
 20            clcon(ip,k)=0.0
         IF (dcflag) THEN

c---- create convection cloud vertical extent array from
c---- kb,kt parameters stored in icln() array (set in conv)
c Modified to handle leads in ice (using ipcln and pavcvrn arrays - LDR)
c.... amount of conv cloud from average precip over radiation
c.... step. From Slingo 1987. (set minimum of 4% conv cloud)
c.... step. From Slingo 1987. (set minimum of 10% conv cloud)
      if(mod(mins,1440_8).eq.0_8)then
        nradcv=1
      else
        nradcv=nrad
      endif
            if (ipass.eq.1) then
               do 40 ip=1,imax
                  kb=icln(ip,1,lg)
                  if (kb.gt.0) then
                     kt=icln(ip,2,lg)
                     cvrn=avcvrn(ip,lg)/nradcv
                     IF(ukconv)THEN
                       cvrn=min(0.8,max(0.10,cvrn))
                     ELSE
c                      cvrn=max(0.19966564,cvrn)
                       cvrn=max(0.32144763,cvrn)
                       cvrn=min(0.243+0.126*log(cvrn),0.8)
                     ENDIF
                     clcon(ip,1)=cvrn
c.... Adjust the convective cloud amount cvrn such that when
c.... randomly overlapped it will reproduce the cloud amount cvrn.
c.... (See CCM2 model : Cl=1 - (1 - Cvrn)**(1/n) )
                     ncl=kt-kb+1
                     cvrn=1.0-(1.0-cvrn)**(1.0/ncl)
                     do 30 k=kb,kt
 30                     clcon(ip,l+1-k)=cvrn
                  end if
 40            continue
            else if (ipass.eq.2) then
               do 60 ip=1,imax
                  kb=ipcln(ip,1,lg)
                  if (kb.gt.0) then
                     kt=ipcln(ip,2,lg)
                     cvrn=pavcvrn(ip,lg)/nradcv
                     IF(ukconv)THEN
                       cvrn=min(0.8,max(0.10,cvrn))
                     ELSE
c                      cvrn=max(0.19966564,cvrn)
                       cvrn=max(0.32144763,cvrn)
                       cvrn=min(0.243+0.126*log(cvrn),0.8)
                     ENDIF
                     clcon(ip,1)=cvrn
                     ncl=kt-kb+1
                     cvrn=1.0-(1.0-cvrn)**(1.0/ncl)
                     do 50 k=kb,kt
 50                     clcon(ip,l+1-k)=cvrn
                  end if
 60            continue
            end if

c---- set up diagnostic cloud
            lgn=lat2+1-lg
            do 78 k=1,l
               do 78 ip=1,lon
c.... pgx is surface pressure in mbs
c.... sdot*pgx gives mbs/sec : sigdt needed in Pa/sec
c.... 1mb = 100 Pa
               pgx=press(ip,lp1)/1000.0
               sigdt(ip,k)=sdot(ip,l+1-k,lgn)*(pgx*100.0)
               sigdt(ip+lon,k)=sdot(ip,l+1-k,lg)*(pgx*100.0)
 78         continue
            do 80 k=1,l
               do 80 ip=1,imax
               rhg1(ip,k)=min(100.*rhg(ip,l+1-k),100.)
 80         continue
               do 81 ip=1,imax
c              rhg1(ip,l-1)=max(rhg1(ip,l-1),rhg1(ip,l))
               if(rhg1(ip,l).gt.rhg1(ip,l-1))
     &          rhg1(ip,l-1)=0.5*(rhg1(ip,l)+rhg1(ip,l-1))
 81         continue
c.... clddia now vectorized for a latitude row
            call clddia (rhg1,press,temp,clh,clm,cll,ich,icm,ict,
     &       icb,clcon,bvnf,clat,sigdt,lg)

         ELSE

            do 90 ip=1,imax/2
               call cldset (alatn,mins,
     &         clh(ip),clm(ip),cll(ip),ich(ip),icm(ip),ict(ip),icb(ip))
               is=ip+imax/2
               call cldset (alats,mins,
     &         clh(is),clm(is),cll(is),ich(is),icm(is),ict(is),icb(is))
 90         continue

         END IF

         do 100 ip=1,imax
 100        nclds(ip)=0

         if (ipass.eq.1) then
            do 110 i=1,imax
               avcvrn(i,lg)=0.0
               icln(i,1,lg)=0
 110           icln(i,2,lg)=0
         else if (ipass.eq.2) then
            do 120 i=1,imax
               pavcvrn(i,lg)=0.0
               ipcln(i,1,lg)=0
 120           ipcln(i,2,lg)=0
         end if

      End If
c***INITIALIZE THE CLOUD AND CLOUD INDEX FIELDS
c     Except for the ground layer (nc=1) the assumption is that
c     no cloud exists. also, for this purpose, the cloud index is set
c     at one (p=0)
c     Don't set cirrf, cuvrf, cirab at the surface because these are set
c     the albedo by radfs.
      do 130 i=1,imax
         camt(i,1)=zero
         emcld(i,1)=one
         ktop(i,1)=1
         kbtm(i,1)=1
         ktopsw(i,1)=1
         kbtmsw(i,1)=1
 130  continue
      do 140 k=2,lp1
         do 140 i=1,imax
         camt(i,k)=zero
         emcld(i,k)=one
         cirrf(i,k)=0.
         cuvrf(i,k)=0.
         cirab(i,k)=0.
         ktop(i,k)=1
         kbtm(i,k)=1
         ktopsw(i,k)=1
         kbtmsw(i,k)=1
 140  continue
c***NOW SET CLOUD AND CLOUD INDEX FIELDS DEPENDING ON THE NO. OF CLOUDS
      nc=1
      do 150 ip=1,imax
c---FIRST, THE ground layer (nc=1)
         emcld(ip,nc)=one
         camt(ip,nc)=one
         ktop(ip,nc)=lp1
         kbtm(ip,nc)=lp1
         ktopsw(ip,nc)=lp1
         kbtmsw(ip,nc)=lp1
 150  continue

      do k=1,l
        do ip=1,imax
          clat(ip,k)=0. !Overwrite clat diagnostic so that it contains the 
        enddo           !cloud seen by the radiation.
      enddo

      if(nl.ge.18)then
        refac=0.9
        emch=0.7  !Overwrite default value of 1.
      else
        refac=1.0
      endif

c---THEN,VALUES FOR THE ACTUAL CLOUDS; IN THIS CASE, UP TO 3 CLOUDS
      do 160 ip=1,imax
c     Low cloud
         nc=2
         if (cll(ip).gt.0) then
            nclds(ip)=nclds(ip)+1
            camt(ip,nc)=cll(ip)
            clat(ip,ict(ip))=cll(ip)
            ktop(ip,nc)=ict(ip)
            kbtm(ip,nc)=icb(ip)
            emcld(ip,nc)=emcl
            ktopsw(ip,nc)=ict(ip)
            kbtmsw(ip,nc)=icb(ip)+1
            cuvrf(ip,nc)=refac*coca(4)
            cirrf(ip,nc)=refac*cwca(4)
            cirab(ip,nc)=refac*cwcb(4)
            nc=nc+1
         end if
c     Medium cloud
         if (clm(ip).gt.0) then
            nclds(ip)=nclds(ip)+1
            camt(ip,nc)=clm(ip)
            clat(ip,icm(ip))=clm(ip)
            ktop(ip,nc)=icm(ip)
            kbtm(ip,nc)=icm(ip)
            emcld(ip,nc)=emcm
            ktopsw(ip,nc)=icm(ip)
            kbtmsw(ip,nc)=icm(ip)+1
            cuvrf(ip,nc)=refac*coca(3)
            cirrf(ip,nc)=refac*cwca(3)
            cirab(ip,nc)=refac*cwcb(3)
            nc=nc+1
         end if
c     High cloud
         if (clh(ip).gt.0) then
            nclds(ip)=nclds(ip)+1
            camt(ip,nc)=clh(ip)
            clat(ip,ich(ip))=clh(ip)
            ktop(ip,nc)=ich(ip)
            kbtm(ip,nc)=ich(ip)
            emcld(ip,nc)=emch
            ktopsw(ip,nc)=ich(ip)
            kbtmsw(ip,nc)=ich(ip)+1
            cuvrf(ip,nc)=refac*coca(2)
            cirrf(ip,nc)=refac*cwca(2)
            cirab(ip,nc)=refac*cwcb(2)
         end if
 160  continue
c***THE FOLLOWING HANDLES OVERLAPPING SHORTWAVE INDICES.NOTE THAT
c  THE BOTTOM SW INDEX IS 1 GREATER THAN THE CORRESPONDING LW INDEX.
c  THEN FURTHER ADJUSTMENT IS MADE,BELOW,TO AVOID OVERLAPPING CLOUDS.
      do 180 ip=1,imax
         if (nclds(ip).le.1) go to 180
         nclds0=nclds(ip)
         do 170 i=2,nclds0
            if (ktopsw(ip,i).lt.kbtmsw(ip,i+1)) then
c    CASE 1- IF LOWER CLOUD IS THICK,ADD 1 TO TOP INDEX OF LOWER CLOUD
               if (ktopsw(ip,i).lt.kbtmsw(ip,i)) then
                  ktopsw(ip,i)=ktopsw(ip,i)+1
               else
c    CASE 2- THICK UPPER CLOUD.(IT IS NOT POSSIBLE TO HAVE THIN UPPER
c    AND LOWER CLOUDS AT ONCE AT THE SAME INDEX LEVEL.) SUBTRACT
c    1 FROM BOTTOM INDEX OF UPPER CLOUD.
                  kbtmsw(ip,i+1)=kbtmsw(ip,i+1)-1
               end if
            end if
 170     continue
 180  continue
      return
      end
