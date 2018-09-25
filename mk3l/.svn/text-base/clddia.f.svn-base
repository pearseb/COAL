c Added IMPLICIT NONE statement to CLDDBLK.
c SJP 2007/05/28
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/10
c
c $Log: clddia.f,v $
c Revision 1.34  2001/02/22 05:34:46  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.33  1999/05/20 06:23:58  rot032
c HBG changes to V5-2
c
c Revision 1.32  1998/12/10  00:56:03  ldr
c HBG changes to V5-1-21
c
c Revision 1.31  1997/12/17  23:23:10  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.30  1997/02/21  00:26:10  ldr
c Go back to rcrits/l=0.85/0.75 for 18L version, tidy cloud2.f and make
c clddia.f, cloud.f general for NL.
c
c Revision 1.29  1996/10/24  01:02:30  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.28  1996/06/13  02:05:37  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.27  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.26  1995/07/18  04:10:12  ldr
c Corrected lg.eq.22 to read lw.eq.22. Caution: 1.23 to 1.25 were wrong!
c
c Revision 1.25  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.24  1995/06/30  03:56:17  ldr
c Corrected lq.eq.22 to read lg.eq.22.
c
c Revision 1.23  1995/05/29  04:59:00  ldr
c Tuning for T63 version for INS/BGH run.
c
c Revision 1.22  94/08/08  17:20:50  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.21  94/08/04  16:54:07  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.20  94/04/29  15:08:46  ldr
c Changes from HBG to reduce high cirrus, increase convective cirrus to
c previous value (2*cvx), to improve jet again. Also change low cloud and
c Cd over sea to achieve balance.
c 
c Revision 1.19  94/03/31  10:22:53  ldr
c Slightly increase high cloud and decrease low cloud to balance.
c 
c Revision 1.18  94/03/30  12:58:08  ldr
c Merge of HBG and LDR changes.
c 
c Revision 1.17  94/03/30  10:24:22  ldr
c Decrease low cloud thresholds to V43 levels, to go with shallow conv
c starting at k=2.
c 
c Revision 1.16.1.1  94/03/30  12:33:55  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.16  93/12/17  15:31:48  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.15  93/10/06  11:52:42  ldr
c Final tuning to high cloud amounts for R42: increase clhmax to .825.
c 
c Revision 1.14  93/10/05  17:52:22  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.13  93/09/24  10:38:47  ldr
c Final tunings to modify V4-4: Modify cloud amounts.
c 
c Revision 1.12  93/07/23  16:02:12  ldr
c Put data statements into block data (for VP).
c 
c Revision 1.11  93/07/06  16:26:49  ldr
c      Using correct check for thick low cloud (25% not 0.25% difference in RH)
c      Modified code for > 2 levels available for thick low cloud
c      New cloud albedoes needed in cloud.f
c 
c Revision 1.10  93/06/16  15:56:26  ldr
c HBG changes to cloud scheme:
c    LS low cloud for RH>70%, but if land RH>60%.   LS high cloud *70%
c    Modified inversion cloud: based on max base RH in inversion levels,
c     and is graded between 50%-70% (was 60%-80%)
c
c     INPUT/OUTPUT
c     Input:   from common/clddiab in this subroutine
c                  icld
c
c              from common/hybrpr in HYBRPR.f
c                  pdpsk - (pressure/surface pressure)**cappa 
c
c              from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c              from common/physical in PHYSICAL.f
c                  he - variance of sub-grid scale topography
c
c              from arguments
c                  bvnf - Brunt-Vaisala frequency, full levels
c                  clcon - convective cloud amount
c                  lg - latitude index
c                  press - pressure (CGS units) at data levels of model
c                  temp - temperature (K) at data levels of model 
c                  sigdt - vertical velocity in Pa/s
c
c     Output:  from arguments
c                  icb - bottom of low cloud    ich - level of high cloud  
c                  icm - level of medium cloud  ict - top of low cloud
c
c     In/Out:  from arguments
c                  ch, cl - total cloud  clat - cloud amount diagnostic
c                  cm, rhum - relative humidity


      block data clddblk

      implicit none

      integer icld
      real cldlev
      common/clddiab/ icld(2,3), cldlev(2,3)  !Communicates with initfs
C
C DATA BLOCKS TO DEFINE THE RELATIVE HUMIDITY LIMITS
C          AND OTHER FACTORS
C
      data cldlev /0.15,0.43,0.43,0.8,0.8,1.0/
      data icld /6*0/
      end

c******************************************************************************

      subroutine clddia(rhum,press,temp,ch,cm,cl,ich,icm,ict,icb
     & ,clcon,bvnf,clat,sigdt,lg)

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /PHYSICAL/ )

C Global parameters
      include 'PARAMS.f'
      include 'RDPARM.f'

C Argument list
      real rhum(imax,l)
      real press(imax,lp1)
      real temp(imax,lp1)
      real ch(imax)
      real cm(imax)
      real cl(imax)
      integer ich(imax)
      integer icm(imax)
      integer ict(imax)
      integer icb(imax)
      real clcon(imax,l)
      real bvnf(imax,l)
      real clat(imax,l)
      real sigdt(imax,l)
      integer lg

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'
      include 'PHYSICAL.f'

C Global data blocks
      include 'LSMI.f'

      integer icld
      real cldlev
      common/clddiab/ icld(2,3), cldlev(2,3)  !Communicates with initfs

C Local work arrays and variables
c     real cd(imax,l)
c     integer kbh(imax,l),kth(imax,l)
      real cinv(imax)
      real rx(imax)

      integer irx(imax)

      integer i
      integer ic1
      integer ic2
      integer ip
      integer j

      real bvfact
      real clhmax
      real clh2
      real cll2
      real cll3
      real clm2
      real clx
      real gradl
      real pb1
      real pm1
      real rhbase
      real rhcr
      real rhl
      real sqn2
      real tb1
      real tmm1
      real z4std

C Local data, functions etc

C Start code : ----------------------------------------------------------

C
C SUBROUTINE TO CALCULATE CLOUD AMOUNT AND LEVEL FROM
C THE RELATIVE HUMIDITY BASED STRATIFORM CLOUD WITH WALKER(1978)
C MODIFICATION FOR LOW CLOUD.
C---------------------------------------------------------------------*
C   VERSION 1.8                                                       *
C---------------------------------------------------------------------*
C     COMMON /TRACKR/ VERNO
C     DATA VERNO /1.8 /
C
C INPUTS INCLUDE TEMP, PRESS  AND CAN
C        BE FOUND IN THE COMMON BLOCKS BELOW.
C RELATIVE HUMIDITY IS CALCULATED OUTSIDE AND INPUT AS RHUM
C
C OUTPUT INCLUDES CLOUDY (CLOUD AMOUNT),
C                 KTH2 - CLOUD TOP LEVEL
C                 KBH2 - CLOUD BOTTOM LEVEL,
C AND THE OTHER REPRESENTATIVE VALUES.  ALL ARE STORED IN THE
C COMMON BLOCKS BELOW.
C
C
C---------------------------------------------------------------------*
C   CHANGES FROM V0.0 ARE:                                            *
C   V0.1 - RHCRIT CALCULATED ONLY ONCE.                               *
C        - RHUM CALCULATED EXTERNALLY.                                *
C   V1.0 - UPPER LIMIT FOR HIGH CLOUD MOVED TO 0.15                   *
C        - LOW CLOUD LOWER LIMIT SET TO 1.0 BUT BOTTOM LEVEL DROPPED. *
C        - BOUNDARY BETWEEN HIGH AND LOW MOVED TO 0.43.               *
C   V1.1 - HIGH CLOUD TOP LEVEL WITH LEVEL ABOVE IF CLOUD THERE.      *
C        - DITTO FOR MIDDLE CLOUD (HAVE TO JUSTIFY THIS).             *
C        - STABILITY CLOUD RESTRICTED TO 2ND BOTTOM LAYER AND MAXIMUM *
C          OVERLAP ASSUMED IF REL HUM CLOUD ALSO PRESENT.             *
C          IF THICK RH CLOUD DIAGNOSED THEN CHECK THAT IT EXISTS IN   *
C          SC LAYER. IF SO MODIFY CD(4) TO MAX OVERLAP. IF NOT IGNORE *
C          SC CLOUD.                                                  *
C   V1.2 - REMOVED AVERAGING WITH LAYER ABOVE FOR BOTH HIGH AND MIDDLE*
C        - TEMPERATURE CONSTRAINT FOR STATIC STABILITY CLOUD (273).   *
C        - CRITICAL RH RESET FOR ALL LEVELS.                          *
C   V1.7 -       -RHC : HIGH 60 MIDDLE 50  LOW 60                     *
C   V1.8 -TEMPERATURE CONSTRAINT ON STATIC STABILITY CLOUD REMOVED    *
C---------------------------------------------------------------------*
C
C---------------------------------------------------------------------*

C     PARAMETER SETTINGS FOR THE LONGWAVE AND SHORTWAVE RADIATION CODE:
C          IMAX   =  NO. POINTS ALONG THE LAT. CIRCLE USED IN CALCS.
C          L      =  NO. VERTICAL LEVELS (ALSO LAYERS) IN MODEL
C***NOTE: THE USER NORMALLY WILL MODIFY ONLY THE IMAX AND L PARAMETERS
C          NBLW   =  NO. FREQ. BANDS FOR APPROX COMPUTATIONS. SEE
C                      BANDTA FOR DEFINITION
C          NBLX   =  NO. FREQ BANDS FOR APPROX CTS COMPUTATIONS
C          NBLY   =  NO. FREQ. BANDS FOR EXACT CTS COMPUTATIONS. SEE
C                      BDCOMB FOR DEFINITION
C          INLTE  =  NO. LEVELS USED FOR NLTE CALCS.
C          NNLTE  =  INDEX NO. OF FREQ. BAND IN NLTE CALCS.
C          NB,KO2 ARE SHORTWAVE PARAMETERS; OTHER QUANTITIES ARE DERIVED
C                    FROM THE ABOVE PARAMETERS.

c Initialization at start of run now done in initfs (LDR 12/92)
c 
c------ Begin cloud calculation
c 
c---------------------------------------------------------------------*
c                                                                     *
c HIGH CLOUD PARAMETER CALCULATION                                    *
c                                                                     *
c---------------------------------------------------------------------*
      if(lw.eq.22)then
        clhmax=0.90
        rhcr=90.0
        bvfact=2.5e-04
      else
        clhmax=1.0
        rhcr=70.0
        bvfact=7.5e-04
      endif
      j=1
      ic1=icld(1,j)
      ic2=icld(2,j)
      do 10 ip=1,imax
         ch(ip)=0.0
         rx(ip)=-4.
 10      ich(ip)=ic2
      do 20 i=ic1,ic2
         do 20 ip=1,imax
         if (rhum(ip,i).gt.rx(ip)) then
           rx(ip)=rhum(ip,i)
           ich(ip)=i
         end if
         sqn2=min(bvnf(ip,lp1-i)**2/bvfact,1.0)
         rhl=99.0-(100.0-rhcr)*(1.0-sqn2)
         rhum(ip,i)=(rhum(ip,i)-100.0*clcon(ip,i))/(1.0-clcon(ip,i))
         clh2=clhmax*(max((rhum(ip,i)-rhl)/(100.0-rhl),0.0))**2
c combine LS cloud and Convective cloud per level : clx
         clx=clcon(ip,i)*(1.0-clh2)+clh2
         clat(ip,i)=clx
c  CALCULATE CONVECTIVE AND RELATIVE HUMIDITY DEPENDENT CLOUD.
c  random overlap down through levels
         ch(ip)=ch(ip)+clx-clx*ch(ip)
 20   continue
c     do 30 ip=1,imax
c      cd(ip,j+1)=ch(ip)
c     kth(ip,j+1) = ich(ip)
c30   kbh(ip,j+1) = ich(ip)
c---------------------------------------------------------------------*
c                                                                     *
c MIDDLE CLOUD PARAMETER CALCULATION                                  *
c                                                                     *
c---------------------------------------------------------------------*
      j=2
      ic1=icld(1,j)
      ic2=icld(2,j)
      if(lw.eq.22)then
        if(nl.ge.18)then
          rhcr=85.0
        else
          rhcr=80.0
        endif
      else
        rhcr=70.0
      endif
      do 40 ip=1,imax
         cm(ip)=0.0
         rx(ip)=-4.
 40      icm(ip)=ic2
      do 50 i=ic1,ic2
         do 50 ip=1,imax
         if (rhum(ip,i).gt.rx(ip)) then
           rx(ip)=rhum(ip,i)
           icm(ip)=i
         end if
         sqn2=min(bvnf(ip,lp1-i)**2/2.5e-04,1.0)
         rhl=99.0-(100.0-rhcr)*(1.0-sqn2)
         rhum(ip,i)=(rhum(ip,i)-100.0*clcon(ip,i))/(1.0-clcon(ip,i))
         clm2=(max((rhum(ip,i)-rhl)/(100.0-rhl),0.0))**2
c combine LS cloud and Convective cloud per level : clx
         clx=clcon(ip,i)*(1.0-clm2)+clm2
         clat(ip,i)=clx
c  CALCULATE CONVECTIVE AND RELATIVE HUMIDITY DEPENDENT CLOUD.
c  random overlap down through levels
         cm(ip)=cm(ip)+clx-clx*cm(ip)
 50   continue
c     do 60 ip=1,imax
c      cd(ip,j+1)=cm(ip)
c     kth(ip,j+1) = icm(ip)
c60   kbh(ip,j+1) = icm(ip)
c---------------------------------------------------------------------*
c                                                                     *
c LOW CLOUD PARAMETER CALCULATION                                     *
c                                                                     *
c---------------------------------------------------------------------*
      ic1=icld(1,3)
      ic2=icld(2,3)

c---------------------------------------------------------------------*
c  CALCULATE THE THERMAL GRADIENT FOR DETECTION OF STABLE AIRMASS     *
c---------------------------------------------------------------------*
      do 140 ip=1,imax
         cinv(ip)=0.0
 140     irx(ip)=0
c  FIND THE MOST STABLE LAYER
      do 150 i=ic1,ic2
         do 150 ip=1,imax
         pb1=press(ip,i+1)*0.001
         tb1=temp(ip,i+1)/pdpsk(ip,lp1-i-1)
         pm1=press(ip,i)*0.001
         tmm1=temp(ip,i)/pdpsk(ip,lp1-i)
         gradl=-10.0*(tmm1-tb1)/(pm1-pb1)
c  IF A LEVEL, WHICH IS STABLE ENOUGH TO PRODUCE CLOUD, HAS BEEN FOUND
c  THEN THE RH AT BASE OF INVERSION MUST BE > 50% FOR CLOUD TO
c  FORM (ECMWF RESEARCH MANUAL 3)
c        rhbase=rhum(ip,i+1)
         rhbase=max(rhum(ip,i),rhum(ip,i+1))
         if((gradl.gt.0.5).and.(rhbase.gt.50.0)) then
            cll3=min(gradl-0.500,1.0)
            if (rhbase.lt.70.0) cll3=cll3*(1.0-(70.0-rhbase)/20.0)
            if(cll3.gt.cinv(ip))then
              cinv(ip)=cll3
              irx(ip)=min(l-1,i+1)
            endif
         endif
 150  continue

c---------------------------------------------------------------------*
c  CALCULATE LARGE SCALE LOW LEVEL CLOUD                              *
c---------------------------------------------------------------------*

      do 70 ip=1,imax
 70      cl(ip)=0.0
      do 80 i=ic1,ic2
         do 80 ip=1,imax
         rhum(ip,i)=(rhum(ip,i)-100.0*clcon(ip,i))/(1.0-clcon(ip,i))
         if(nl.ge.18)then
           rhl=85.0
         else
           rhl=75.0
         endif
         if(imsl(ip,lg).eq.4)then
c.... if land, then allow for more cloud over mountainous regions
c.... (subgridscale topographic cloud) by reducing the threshold relative
c.... humidity. This reduction depends upon the stability (via N**2) and
c.... the roughness of the surface via the std deviation.
c.... The basal threshold is also lower than over sea to allow for extra CCN
           sqn2=min(bvnf(ip,lp1-i)**2/2.5e-04,1.0)
           z4std=he(ip)*0.5
           if(nl.ge.18)then
             rhl=85.0-20.0*(1.0-sqn2)*min(1.0,z4std/1000.0)
           else
             rhl=70.0-20.0*(1.0-sqn2)*min(1.0,z4std/1000.0)
           endif
         end if
         cll2=(max((rhum(ip,i)-rhl)/(100.0-rhl),0.0))**2
c VERTICAL VELOCITY MAY BE USED TO CONSTRAIN LOW CLOUD
c        if((sigdt(ip,i).gt.0.0).and.(rhum(ip,i).lt.99.0)) cll2=0.0
         if(sigdt(ip,i).gt.0.0) cll2=0.0
c combine LS cloud and inversion cloud (if any)
         if(irx(ip).eq.i)cll2=max(cll2,cinv(ip))
c combine LS cloud and Convective cloud per level : clx
         clx=clcon(ip,i)*(1.0-cll2)+cll2
         clat(ip,i)=clx
c  CALCULATE CONVECTIVE AND RELATIVE HUMIDITY DEPENDENT CLOUD.
c  random overlap down through levels
         cl(ip)=cl(ip)+clx-clx*cl(ip)
 80   continue

c set low cloud position (No thick cloud in this version)
c (Needs new coding for setting low cloud when >2 levels 
c available for low cloud)
c              ----------------
c              ......CCCCCCCCCC      = level ic1
c              ----------------
c              ----------------
c              CCCCCCCCC.......      = level ic2
c              ----------------
c----
      do 110 ip=1,imax
         ict(ip)=ic2
 110     icb(ip)=ic2
c     do 115 ip=1,imax
c        if(clat(ip,ic1).gt.2.0*clat(ip,ic2))then
c          icb(ip)=ic1
c          ict(ip)=ic1
c        endif
c        if(clat(ip,ic2).gt.2.0*clat(ip,ic1))then
c          icb(ip)=ic2
c          ict(ip)=ic2
c        endif
c115  continue
ch set low cloud, focussing on lowest level
      do 115 i=ic1,ic2
      do 115 ip=1,imax
         if(clat(ip,i).gt.0.0)then
           ict(ip)=i
           icb(ip)=i
         endif
 115  continue

      do 170 ip=1,imax
c     cd(ip,4)=cl(ip)
c     cd(ip,5)=cinv(ip)
c 
c CALCULATE THE RELEVANT SHORT-WAVE CLOUD PARAMETERS
c 
c     kth(ip,4) = ict(ip)
c     kbh(ip,4) = icb(ip)
         if (cl(ip).lt.0.02) then
            ict(ip)=l
            icb(ip)=l
         end if
 170  continue
     
      return
      end
