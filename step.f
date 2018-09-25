c Added linear interpolation in time between monthly values of heat and
c  freshwater fluxes used for bulk forcing of OGCM
c Pearse J Buchanan 16/07/2018

c 2015 Dec 9    Pearse J Buchanan
c     - general clean up of biogeochemical loops

c Further modifications to the ocean model output routines to improve memory
c usage. The features that were causing stack overflows have now been resolved.
c SJP 2009/05/11
c
c Modifying the ocean model output routines for improved memory usage.
c SJP 2009/05/06
c
c Adding density as an ocean model statistic.
c SJP 2009/04/21
c
c Enabling user control over the ocean model statistics that are saved to file.
c SJP 2009/04/19
c
c Implementing the equation of state of McDougall et al (2003).
c SJP 2009/04/08
c
c (1) Implementing a Robert time filter.
c (2) Adding code to check for conservation of tracers within the ocean.
c SJP 2008/12/17
c
c (1) Modified so that flux adjustments are only applied when FLUXADJ=T.
c (2) Removed one array bounds violation.
c SJP 2008/03/08
c
c Modified for coupling of the new, high-resolution version of the ocean model
c to the atmosphere.
c SJP 2007/12/20
c
c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Rationalise the diagnostic information written to standard output.
c SJP 2007/11/21
c
c Major tidy-up of the ocean model source code.
c SJP 2007/06/17
c
c Removed the redundant calls to STEP6, TRACRI and DENCAL, plus the redundant
c references to the arrays HMIX and HSUM, and the header files ICEOC.f and
c RHOJMIX.f.
c SJP 2007/06/02
c
c Added IMPLICIT NONE statement, plus variable declarations as necessary.
c SJP 2007/05/30
c
c Transferred COMMON block declarations to separate header files.
c SJP 2007/05/29
c
c Modified to enable relaxation of the coupled model SSTs and SSSs towards
c prescribed values.
c SJP 2006/01/05
c
c Call to STEP1 removed, which is redundant.
c SJP 2005/03/19
c
c (1) Calls to STEP2, STEP3 and STEP5 removed, as these are redundant.
c (2) Remove unnecessary variable declarations.
c SJP 2004/01/06
c
c Subroutines OCDATRA, OCDATRO, R2FCFLD and STEP0 split off into separate
c files.
c SJP 2003/12/30
c
c Adjacent loops with identical bounds fused, where possible.
c SJP 2003/12/20
c
c (1) Modified so that the Gent-McWilliams eddy parameterisation scheme is
c always used, enabling some optimisation. The parameter IGM is thus rendered
c irrelevant.
c (2) Some tidying up.
c SJP 2003/12/18
c
c (1) Removed the adjustments made by the stand-alone OGCM to the
c     climatological SSTs and SSSs at two locations: the tips of South America
c     and the Antarctic Peninsula. This is unnecessary if the climatological
c     data has been interpolated to provide values at these points.
c (2) Renamed the files read in by the stand-alone OGCM as follows:
c       hr4seas.dat     -> stress.dat (surface stresses)
c       sstd            -> sst.dat    (SSTs)
c       stsall_fix6.dat -> sss.dat    (SSSs)
c (3) Writes to Fortran unit 81 commented out.
c SJP 2003/09/07
c
c (1) Modified for changes to /orestart/.
c (2) Removed arrays ATMK and BTMK from STEP0, as these are filled, but never
c     used.
c (3) Commented out the reading of the wind stress reduction factors by the
c     stand-alone OGCM, as I intend to force the model with climatological
c     momentum fluxes calculated by the AGCM.
c (4) Commented out the reading of the wind stress adjustments by the coupled
c     model, as I do not intend to use these.
c SJP 2003/09/04
c
c Replaced calls to OGET/OPUT with data access via COMMON block /orestart/.
c SJP 2003/09/02
c
c Some loops restructured to avoid array bounds violations.
c SJP 2003/05/02
c
c Fix up calls to OGET and OPUT.
c SJP 2003/05/01
c
c Calls to OFIND removed, as they did nothing. Parameter definitions moved
c to include file OPARAMS.f.
c SJP 2003/04/29
c
c Writes to fort.76 in STEP6 commented out, as this data is not required.
c SJP 2002/06/27
c
c Commented out read from cwice.dat12 in OCDATRA, as this data is only
c used when the ocean model is run in stand-alone mode.
c SJP 2002/02/15
c
C $Log: step.f,v $
C Revision 1.23  2001/10/12 02:13:44  rot032
C LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
C
C Revision 1.22  2000/08/16 02:59:24  rot032
C NCEP initialization option and CGCM changes from HBG
C
C Revision 1.21  1999/05/20 06:23:47  rot032
C HBG changes to V5-2
C
c Revision 1.20  1998/05/27  03:01:38  ldr
c Tidy up last version.
c
c Revision 1.19  1998/05/27  02:10:17  ldr
c Merge TIE and ACH changes.
c
c Revision 1.18  1998/05/26  05:10:49  ldr
c Final 10k run changes (mainly ACH salinity stuff) merged into V5-1-9.
c
c Revision 1.17.1.1  1998/05/27  02:07:33  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.17  1997/12/23  00:23:32  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.16.1.1  1998/05/26  04:48:54  ldr
c Final changes (inc. ACH salinity stuff) for 10k run.
c
c Revision 1.16  1997/12/19  01:25:30  ldr
c Changes from ACH for 21 OGCM levels and optinal GM mixing scheme...
c
c Include array for storing depth of convective mixed layer
c Include code to save monthly averages of uedd, vedd and wedd
c Also include code to print certain diffusivity profiles every six months
c Change to print out slabs of T and S near antarctica
c at half yearly intervals
c Change to 21 levels in ocean model, and insertion of
c eddy-induced transport (major changes delineated)
c
c Revision 1.15.1.1  1997/12/19  02:03:08  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.15  1997/03/06  23:33:59  ldr
c Corrections from HBG to the recent tidy-ups to ocean routines.
c
c Revision 1.14  1996/03/21  03:19:05  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.13  1994/08/08  17:22:35  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.12  94/03/30  12:35:06  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.11  93/12/17  15:33:49  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.10  93/12/17  11:51:45  ldr
c Changes to V4-4-45l from HBG for coupled model
c 
c Revision 1.9  93/11/29  11:38:39  ldr
c Changes to V4-4-32l from HBG for coupled model
c 
c Revision 1.8  93/11/03  11:44:31  ldr
c Changes to V4-4-24l from HBG for coupled model
c 
c Revision 1.7  93/10/05  13:07:33  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.6  93/08/19  15:10:49  ldr
c Minor cosmetic changes.
c 
c Revision 1.5  93/08/13  14:43:43  ldr
c Minor changes to get coupled model to run on SGI.
c 
c Revision 1.4  93/08/10  15:28:01  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c 
      SUBROUTINE STEP
C
C=======================================================================
C                                                                    ===
C  STEP IS CALLED ONCE PER TIMESTEP.  IT INITIALIZES VARIOUS         ===
C       QUANTITIES, BOOTSTRAPS THE BASIC ROW BY ROW COMPUTATION      ===
C       OF PROGNOSTIC VARIABLES, MANAGES THE I/O FOR THE LATTER,     ===
C       AND PERFORMS VARIOUS ANALYSIS PROCEDURES ON THE PROGRESSING  ===
C       SOLUTION.                                                    ===
C                                                                    ===
C=======================================================================
C                                                                    ===
C      Note : step.f calls                                           ===
C                  ocdatro (data if uncoupled) or                    ===
C                      ocdatra (data if coupled)                     ===
C                  STINIT/M2003                                      ===
C                  STATE/M2003                                       ===
C --- DO 380 J=2,JMTM1 --------------------------------------------  ===
C                  CLINIC                                            ===
C                  TRACER                                            ===
C 380 CONTINUE ----------------------------------------------------  ===
C                  RELAX                                             ===
C                                                                    ===
C=======================================================================
C
      implicit none
C
C---------------------------------------------------------------------
C  DEFINE GLOBAL DATA
C---------------------------------------------------------------------
C
      include 'OPARAMS.f'
      include 'OCEAN_NML.f'
      include 'ORESTART.f'
      include 'FULLWD.f'
      include 'SCALAR.f'
      include 'ONEDIM.f'
      include 'WORKSP.f'
      include 'TIME1.f'
      include 'ACHEXT.f'
      include 'SFC4.f'
      include 'TTFC.f'  !agcm
      include 'ETRANS1.f'
      include 'ETRANS2.f'
      include 'KTSIN.f'
      include 'CHANGEM.f'
      include 'CONVECTD.f'
      include 'MDAY.f'
      include 'LEVD.f'
      include 'SSTSAL.f'
      include 'ATM2OC.f'        !agcm
      include 'FLXBLOK.f'
      include 'A2O.f'
      include 'O2A.f'
      include 'FEWFLAGS.f'
      include 'CRELAX.f'
      include 'SFC1.f'
      include 'COEFFS.f'
      include 'OHIST.f'
      include 'bulkf_variables.f' !pjb
C
C---------------------------------------------------------------------
C  DEMENSION AND EQUIVALENCE LOCAL DATA
C---------------------------------------------------------------------
C
      integer ii, jj, kk, ll, i, j, k, m, nk, nj, ni, iprt, istrt,
     &        istop, kdep, ndiskc, ndiskx, l
      real vbr, tbrn, tbrs, ttn, tmt, rdate, fx, tglobe, sglobe,
     &     vglobe, fsglobe, fhglobe, diag1, diag2, scl, totdx,
     &     totdz, vbrz, tbrz, daysyr, ttyear, ttday, darea, dsf,
     &     dtf, dvol, hint, sbar1, sbar2, sint, tbar1, tbar2, tint,
     &     vint
      DIMENSION VBR(KM),TBRN(KM,NT),TBRS(KM,NT),TTN(8,JMT,NTMIN2),
     &          TMT(JMT,KM)
      save ttn
      real tlevel(km,4),alevel(km),tmin(km)
      integer itmin(km),jtmin(km)
      real vbredd(km),tmtedd(jmt,km)
      real sden(imt), tden(imt), rden(imt)

c rjm....
      include "bio.h"
      include "extra.h"
      integer ncount,ncday
      real xcount,xicount
      save nsum1,itlast,ncount
      data nsum1 /0/, itlast/0/, ncount/0/
      integer lmon1_t, lmon2_t
      save lmon1_t, lmon2_t
c rjm
C
C---------------------------------------------------------------------
C  BEGIN EXECUTABLE CODE
C---------------------------------------------------------------------
C
      entry ocstep
C
C=======================================================================
C  BEGIN SECTION FOR THE INITIALIZATION OF  ============================
C  VARIOUS QUANTITIES ON EACH TIMESTEP      ============================
C=======================================================================
C
C---------------------------------------------------------------------
C  UPDATE TIMESTEP COUNTER AND TOTAL ELAPSED TIME
C---------------------------------------------------------------------
C
      ITT=ITT+1
      ittx=ITT
      TTSEC=TTSEC+DTTSF

c..........Calculate weights for Levitus Climatology........
c
      if(itt.eq.itfs+1)then

        do m = 1, ntmin2
          do j = 1, jmt
            do ll = 1, 8
              TTN(LL,J,M)=0.0
            end do
          end do
        end do

        lmon1=nmth-1
        if(lmon1.eq.0)lmon1=12
        lmon1n=lmon1
        lmon2=nmth
        date=0.5*float(mdays(lmon1)*knitd)
          rdate=date/float(knitd)
          print *,'starting values: lmon1=',lmon1,' lmon2=',
     &     lmon2,'  date=', date,' rdate=',rdate
        wtint=0.5*float((mdays(lmon1)+mdays(lmon2))*knitd)
      end if
      date=date+1.0
      if(date.gt.wtint)then
        lmon1=nmth
        lmon1n=lmon1
        lmon2=nmth+1
        if(lmon2.eq.13)lmon2=1
        date=date-wtint
        wtint=0.5*float((mdays(lmon1)+mdays(lmon2))*knitd)
      end if

c---- having set up month indicators, check if data read required
      if(lmon1n.ne.lmon1o)then
c-- data read is required upon starting model (lmon1o=-1,lmon1n=0) or
c-- when there is a changeover at the middle of a month (lmon1n.ne.lmon1o).
c-- Upon entry to ocdatra/ocdatro, lmon1 and lmon2 have values 1 to 12
         If(lcouple)Then
           if (fluxadj) call ocdatra
           if (crelax_flag) call ocdatro
         Else
           call ocdatro
         End If
         lmon1o=lmon1n
c-- reset lmon1,lmon2 to 1,2 for data storage arrays holding only 2 months
c-- (not 12 months).

c rjm ...
c What is the reason to copy everything to an array with 2 times when the memory
c required is trival!  I will store the index for my calculations
         lmon1_t=lmon1
         lmon2_t=lmon2
c rjm 
         lmon1=1
         lmon2=2
      end if
c----
        w2=date/wtint
        w1=1.-w2
        wt2=w2
        wt1=w1
      If (.not. lcouple .or. crelax_flag) then
c........compute new levitus sst and salinity values for this timestep..

         do 151 j=2,jmtm1
         do 151 i=1,imt

c rjm...
         if (nt.gt.2) then
* ice correction include when reading file
* 1 - windspeed^2 *(1-ice), 2- incident solar rad, 3- sea ice, 4- atmos Fe deposition
* 5 - reactive NO3 deposition *(1-ice), 6 - atmospheric pressure (atm)           
            sbcbio(i,j,1) = obc(i,j,lmon1_t,1)*w1+obc(i,j,lmon2_t,1)*w2
            sbcbio(i,j,2) = obc(i,j,lmon1_t,2)*w1+obc(i,j,lmon2_t,2)*w2
            sbcbio(i,j,3) = obc(i,j,lmon1_t,3)*w1+obc(i,j,lmon2_t,3)*w2
            sbcbio(i,j,4) = obc(i,j,lmon1_t,4)*w1+obc(i,j,lmon2_t,4)*w2
            sbcbio(i,j,5) = obc(i,j,lmon1_t,5)*w1+obc(i,j,lmon2_t,5)*w2
            sbcbio(i,j,6) = obc(i,j,lmon1_t,6)*w1+obc(i,j,lmon2_t,6)*w2
         endif
c rjm

         if (bulk_force) then
            spechum(i,j) = w1*spechumm(i,j,lmon1)+w2*spechumm(i,j,lmon2)
            airtem(i,j) = w1*airtemm(i,j,lmon1)+w2*airtemm(i,j,lmon2)
            slp(i,j) = w1*slpm(i,j,lmon1)+w2*slpm(i,j,lmon2)
            swdown(i,j) = w1*swdownm(i,j,lmon1)+w2*swdownm(i,j,lmon2)
            lwdown(i,j) = w1*lwdownm(i,j,lmon1)+w2*lwdownm(i,j,lmon2)
            fclt(i,j) = w1*fcltm(i,j,lmon1)+w2*fcltm(i,j,lmon2)
            rain(i,j) = w1*rainm(i,j,lmon1)+w2*rainm(i,j,lmon2)
            snow(i,j) = w1*snowm(i,j,lmon1)+w2*snowm(i,j,lmon2)
            roff(i,j) = w1*roffm(i,j,lmon1)+w2*roffm(i,j,lmon2)
            fice(i,j) = w1*ficem(i,j,lmon1)+w2*ficem(i,j,lmon2)
            uwnd(i,j) = w1*uwndm(i,j,lmon1)+w2*uwndm(i,j,lmon2)
            vwnd(i,j) = w1*vwndm(i,j,lmon1)+w2*vwndm(i,j,lmon2)
         endif
         if (nt.gt.2) then
            dust(i,j) = w1*dustm(i,j,lmon1)+w2*dustm(i,j,lmon2)
            nox(i,j) = w1*noxm(i,j,lmon1)+w2*noxm(i,j,lmon2)
            if (.not.bulk_force) then ! get fields for BGC when restoring to SST and SSS
              slp(i,j) = w1*slpm(i,j,lmon1)+w2*slpm(i,j,lmon2)
              swdown(i,j) = w1*swdownm(i,j,lmon1)+w2*swdownm(i,j,lmon2)
              roff(i,j) = w1*roffm(i,j,lmon1)+w2*roffm(i,j,lmon2)
              fice(i,j) = w1*ficem(i,j,lmon1)+w2*ficem(i,j,lmon2)
              uwnd(i,j) = w1*uwndm(i,j,lmon1)+w2*uwndm(i,j,lmon2)
              vwnd(i,j) = w1*vwndm(i,j,lmon1)+w2*vwndm(i,j,lmon2)
            endif
         endif

         sst(i,j)=w1*sstm(i,j,lmon1)+w2*sstm(i,j,lmon2)
 151     sal(i,j)=w1*salm(i,j,lmon1)+w2*salm(i,j,lmon2)

      End If
! for dust and NOx deposition use existing climatology for couple run
      do j=2,jmtm1
         do i=1,imt
            dust(i,j) = dustm(i,j,lmon1_t)*w1+dustm(i,j,lmon2_t)*w2
            nox(i,j) = noxm(i,j,lmon1_t)*w1+noxm(i,j,lmon2_t)*w2
         enddo
      enddo
C
C---------------------------------------------------------------------
C  UPDATE PERMUTING DISC I/O UNITS
C---------------------------------------------------------------------
C
      NDISKB=MOD(ITT+1,2)+1
      NDISK =MOD(ITT  ,2)+1
      NDISKA=NDISKB
C
C---------------------------------------------------------------------
C  ADJUST VARIOUS QUANTITIES FOR MIXING TIMESTEP
C---------------------------------------------------------------------
C
      MIX=0
      MXP=0
      C2DTTS=2.0*DTTSF
      C2DTUV=2.0*DTUVF
      C2DTSF=2.0*DTSFF
      IF ((.not. robert_time_filter) .and. (MOD(ITT,NMIX).EQ.1)) THEN
        MIX=1
        C2DTTS=DTTSF
        C2DTUV=DTUVF
        C2DTSF=DTSFF

        DO 170 J=1,JMT
        DO 170 I=1,IMT
          PB(I,J)=P(I,J)
 170    CONTINUE

      ENDIF
 182  CONTINUE

      if (lcouple) then

C******COPY OTX AND OTY AGCM WINDS STRESSES INTO STRX AND STRY

      do j = 2, jmtm1
        do i = 2, imtm1
          strx(i, j, 1) = otx(i-1, j-1)
          stry(i, j, 1) = oty(i-1, j-1)
        end do
      end do
c rjm  copy atmospheric fields to the obgc field
      if (nt.ge. 3) then

         ! fix for getting an approximate daily average
         ncount=ncount +1  
         ncday = 86400/DTTSF !number of timesteps per day
         xcount=min(ncday,ncount)*1.
         xicount=1./xcount

         do j = 2, jmtm1
            do i = 2, imtm1
               sbcbio(i,j,1) = oswnd(i-1,j-1)**2 *(1.-osice(i-1,j-1))  

               ! average over a day for consistency with ocean only formulation of EP
               if (ncount.eq.1) then 
                  sbcbio(i,j,2) = max(0.0, osrad(i-1,j-1)) * 
     +                            (1.-osice(i-1,j-1))
               endif
               sbcbio(i,j,2) = ( sbcbio(i,j,2) * (xcount-1) + 
     +                           max( 0.0,osrad(i-1,j-1) )
     +                           * (1.-osice(i-1,j-1)) ) * xicount

               sbcbio(i,j,3) = osice(i-1,j-1)  
               !sbcbio(i,j,6) = pslgrid(i-1,j-1)*(1./1013.2501) ! convert hPa --> atm
                  !	  sbcbio(i,j,4) = 0  ! need dust from atmospheric model
               ! use the climatological value
            end do
         end do

!      print*,'rjm ',xcount,sbcbio(40,40,2),osrad(39,39),osice(39,39)
      endif
c rjm

      DO 6001 J=2,JMTM1
      STRX(1,J,1)=STRX(IMT-1,J,1)
      STRY(1,J,1)=STRY(IMT-1,J,1)
      STRX(IMT,J,1)=STRX(2,J,1)
 6001 STRY(IMT,J,1)=STRY(2,J,1)
      DO 6002 I=1,IMT
      STRX(I,1,1)=STRX(I,2,1)
      STRY(I,1,1)=STRY(I,2,1)
      STRX(I,JMT,1)=STRX(I,JMTM1,1)
 6002 STRY(I,JMT,1)=STRY(I,JMTM1,1)

      end if
C
C---------------------------------------------------------------------
C  ESTABLISH OVER DIMENSIONED ARRAYS FOR VECTORIZATION
C---------------------------------------------------------------------
C

      DO 184 K=1,KM
      DO 184 I=1,IMT
        DXTQ  (I,K)=DXT  (I)
        DXT4RQ(I,K)=DXT4R(I)
        DXUQ  (I,K)=DXU  (I)
        DXU2RQ(I,K)=DXU2R(I)
        DZZQ  (I,K)=DZZ  (K)
        DZ2RQ (I,K)=DZ2R (K)
        DZZ2RQ(I,K)=DZZ2R(K)
        C2DZQ (I,K)=C2DZ (K)
        AHIQ  (I,K)=AHI  (K)
        AHHQ  (I,K)=AHH  (K)
        aheq  (i,k)=ahe  (k)
        EEMQ  (I,K)=EEM  (K)
        FFMQ  (I,K)=FFM  (K)
        DTXQ   (I,K)=DTXF(K)
 184  CONTINUE

C
C---------------------------------------------------------------------
C  LOAD COEFFICIENT ARRAYS FOR SUBSEQUENT CALLS TO "STATE" AND "STATEC"
C---------------------------------------------------------------------
C
      if (.not. m2003_eos) CALL STINIT
C
C---------------------------------------------------------------------
C  RESET KM+1 BOUNDARY VELOCITY TO ZERO
C---------------------------------------------------------------------
C

      DO 122 I=1,IMT
        UUNDER(I)=0.
        VUNDER(I)=0.
 122  CONTINUE

C
C---------------------------------------------------------------------
C  INITIALIZE VARIOUS QUANTITIES USED FOR ANALYSIS OF THE SOLUTION
C---------------------------------------------------------------------
C
      EKTOT=0.0

      DO 130 M=1,NT
        DTABS(M)=0.0
        TVAR(M)=0.0
 130  CONTINUE

      NERGY=0
      IF(MOD(ITT,NNERGY).EQ.0) NERGY=1
      IF(NERGY.EQ.1 .AND. MXP.EQ.0) THEN
        BUOY=0.0

        DO 190 LL=1,8
          ENGINT(LL)=0.0
          ENGEXT(LL)=0.0
        DO 190 I=1,IMT
          ZUSENG(I,LL)=0.0
          ZVSENG(I,LL)=0.0
 190    CONTINUE

        DO 192 M=1,NT
        DO 192 LL=1,6
          TTDTOT(LL,M)=0.0
 192    CONTINUE

        DO 194 J=1,JMT
          DO 193 M=1,NTMIN2
          DO 193 LL=1,8
            TTN(LL,J,M)=0.0
 193      CONTINUE
        DO 194 K=1,KM
          TMT(J,K)=0.0
          tmtedd(j,k)=0.0
 194    CONTINUE

      ENDIF

C                 initialize deep ocean temperature sums

      tglobe = 0.
      sglobe = 0.
      vglobe = 0.
      fsglobe = 0.
      fhglobe = 0.

      do m = 1,4
      do k = 1,km
        tlevel(k,m) = 0.  
        alevel(k)  = 0.  
        tmin(k)  = 30.  
        itmin(k)  = 0. 
        jtmin(k)  = 0.
      enddo
      enddo

C
C=======================================================================
C  END OF SECTION FOR INITIALIZATION  ==================================
C=======================================================================
C

c...  If checking conservation of tracers, calculate mean temperature and
c...  salinity of ocean at t-1 time level
      if (check_conservation) then

        ndiskc = ndiskb
        if (mix .eq. 1) ndiskc = ndisk

        tint = 0.0
        sint = 0.0
        vint = 0.0
        do j = 2, jmtm1
          do i = 2, imtm1
            do k = 1, kmt(i, j)
              dvol = cst(j) * dxt(i) * dyt(j) * dz(k)
              tint = tint + dvol * odam_t(i, j, k, 1, ndiskc)
              sint = sint + dvol * odam_t(i, j, k, 2, ndiskc)
              vint = vint + dvol
            end do
          end do
        end do
        tbar1 = tint / vint
        sbar1 = 35.0 + 1000.0 * (sint / vint)

      end if

C=======================================================================
C  BEGIN A BOOTSTRAP PROCEDURE TO PREPARE FOR THE  =====================
C  ROW-BY-ROW COMPUTATION OF PROGNOSTIC VARIABLES  =====================
C=======================================================================
C
C---------------------------------------------------------------------
C  FETCH DATA FOR ROW 2 FROM THE DISC
C---------------------------------------------------------------------
C

      do kk = 1, nt
        do jj = 1, km
          do ii = 1, imt
            tbp(ii, jj, kk) = odam_t(ii, 2, jj, kk, ndiskb)
            tp(ii, jj, kk) = odam_t(ii, 2, jj, kk, ndisk)
          end do
        end do
      end do

      do jj = 1, km
        do ii = 1, imt
          ubp(ii, jj) = odam_u(ii, 2, jj, ndiskb)
          vbp(ii, jj) = odam_v(ii, 2, jj, ndiskb)
          up(ii, jj) = odam_u(ii, 2, jj, ndisk)
          vp(ii, jj) = odam_v(ii, 2, jj, ndisk)
        end do
      end do

C
C---------------------------------------------------------------------
C  MOVE TAU-1 DATA TO TAU LEVEL ON A MIXING TIMESTEP
C---------------------------------------------------------------------
C
      IF(MIX.EQ.1) THEN

        DO 224 M=1,NT
        DO 224 K=1,KM
        DO 224 I=1,IMT
          TBP(I,K,M)=TP(I,K,M)
 224    CONTINUE

        DO 226 K=1,KM
        DO 226 I=1,IMT
          UBP(I,K)=UP(I,K)
          VBP(I,K)=VP(I,K)
 226    CONTINUE

      ENDIF
C
C---------------------------------------------------------------------
C  INITIALIZE ARRAYS FOR FIRST CALLS TO CLINIC AND TRACER
C---------------------------------------------------------------------
C
      FX=0.0

      do k = 1, nt
        do j = 1, km
          do i = 1, imt
            tb(i, j, k) = fx
            t(i, j, k) = fx
          end do
        end do
      end do
      do j = 1, km
        do i = 1, imt
          ub(i, j) = fx
          vb(i, j) = fx
          u(i, j) = fx
          v(i, j) = fx
        end do
      end do

      DO 250 K=1,KM
      DO 250 I=1,IMT
        FVST(I,K)=FX
        RHOS(I,K)=FX
        FMM (I,K)=FX
        FM  (I,K)=FX
C
C---------------------------------------------------------------------
C  CONSTRUCT MASK ARRAY FOR ROW 2 TRACERS
C---------------------------------------------------------------------
C
        IF(KMT(I,2).GE.KAR(K)) THEN
          FMP(I,K)=1.0
        ELSE
          FMP(I,K)=0.0
        ENDIF
 250  CONTINUE

C
C---------------------------------------------------------------------
C  SET VORTICITY COMPUTATION ARRAYS AT SOUTHERN WALL
C---------------------------------------------------------------------
C
      FX=0.0

      DO 258 I=1,IMT
        ZUS(I)=FX
        ZVS(I)=FX
 258  CONTINUE

C
C---------------------------------------------------------------------
C  SAVE INTERNAL MODE VELOCITIES FOR ROW 2
C  AND COMPUTE ADVECTIVE COEFFICIENT FOR SOUTH FACE OF ROW 2 U,V BOXES
C---------------------------------------------------------------------
C
      FX=DYU2R(2)*CSR(2)*CST(2)*0.5

      DO 260 K=1,KM
      DO 260 I=1,IMT
        UCLIN(I,K)=UP(I,K)
        VCLIN(I,K)=VP(I,K)
        FVSU(I,K)=(VP(I,K)+V(I,K))*FX
 260  CONTINUE

C
C---------------------------------------------------------------------
C  COMPUTE EXTERNAL MODE VELOCITIES FOR ROW 2
C---------------------------------------------------------------------
C
C  1ST, COMPUTE FOR TAU-1 TIME LEVEL
C
C  2ND, COMPUTE FOR TAU TIME LEVEL
C
      J=1

      DO 270 I=1,IMTM1
        DIAG1=PB(I+1,J+2)-PB(I  ,J+1)
        DIAG2=PB(I  ,J+2)-PB(I+1,J+1)
        SFUB(I)=-(DIAG1+DIAG2)*DYU2R(J+1)*HR(I,J+1)
        SFVB(I)= (DIAG1-DIAG2)*DXU2R(I  )*HR(I,J+1)*CSR(J+1)
        DIAG1=P (I+1,J+2)-P (I  ,J+1)
        DIAG2=P (I  ,J+2)-P (I+1,J+1)
        SFU (I)=-(DIAG1+DIAG2)*DYU2R(J+1)*HR(I,J+1)
        SFV (I)= (DIAG1-DIAG2)*DXU2R(I  )*HR(I,J+1)*CSR(J+1)
 270  CONTINUE

C
C  3RD, SET CYCLIC BOUNDARY CONDITIONS
C
      SFUB(IMT)=SFUB(2)
      SFVB(IMT)=SFVB(2)
      SFU (IMT)=SFU (2)
      SFV (IMT)=SFV (2)

c...  If using Robert time filter, save external mode
      if (robert_time_filter) then
        do i = 1, imt
          ssfubp(i) = sfub(i)
          ssfvbp(i) = sfvb(i)
        end do
      end if

C
C---------------------------------------------------------------------
C  ADD EXTERNAL MODE TO INTERNAL MODE FOR ROW 2  (OCEAN PTS. ONLY)
C---------------------------------------------------------------------
C

      DO 300 K=1,KM
      DO 300 I=1,IMU
        IF(KMU(I,2).GE.KAR(K)) THEN
          UBP(I,K)=UBP(I,K)+SFUB(I)
          VBP(I,K)=VBP(I,K)+SFVB(I)
          UP (I,K)=UP (I,K)+SFU (I)
          VP (I,K)=VP (I,K)+SFV (I)
        ENDIF
 300  CONTINUE

C
C---------------------------------------------------------------------
C  ACCUMULATE KINETIC ENERGY FROM ROW 2 EVERY NTSI TIMESTEPS
C---------------------------------------------------------------------
C
      IF(MOD(ITT,NTSI).EQ.0) THEN

        DO 305 K=1,KM
          FX=0.5*CS(J+1)*DYU(J+1)*DZ(K)
        DO 305 I=2,IMUM1
          EKTOT=EKTOT+(UP(I,K)*UP(I,K)+VP(I,K)*VP(I,K))*FX*DXU(I)
 305    CONTINUE

      ENDIF
C
C---------------------------------------------------------------------
C  COMPUTE DENSITY OF ROW 2
C---------------------------------------------------------------------
C
      if (m2003_eos) then

        do k = 1, km
          if (omask(2, k)) then
            do i = 1, imt
              sden(i) = 1000.0*tp(i, k, 2)+35.0
              tden(i) = tp(i, k, 1)
            end do
            call m2003(sden, tden, zdb(k), rden)
            do i = 1, imt
              rhos(i, k) = rden(i)
            end do
          else
            do i = 1, imt
              rhos(i, k) = 1.0
            end do
          end if
        end do

c.....  Save density if required
        if (save_rho .and. (mxp .eq. 0)) then
          do k = 1, km
            do i = 2, imtm1
              ohist_rho(i, 2, k) = ohist_rho(i, 2, k) + rhos(i, k)
            end do
          end do
        end if

      else

        CALL STATE(TP(1,1,1),TP(1,1,2),RHOS,TDIF(1,1,1),TDIF(1,1,2))

      end if
C
C  SET CYCLIC BOUNDARY CONDITIONS
C

      DO 310 K=1,KM
        RHOS(IMT,K)=RHOS(2,K)
 310  CONTINUE

C
C=======================================================================
C  END OF BOOTSTRAP PROCEDURE  =========================================
C=======================================================================
C
C=======================================================================
C  BEGIN ROW-BY-ROW COMPUTATION OF PROGNOSTIC VARIABLES  ===============
C=======================================================================
C

C-----------------------------------------------------------------------
C  Prepare some biogeochemical diagnostics for the next timestep

      if (vary_stoich) then
         ! Calculate average C:P ratio for next timestep
         if (count(carb2P.gt.0.0).gt.0.0) then
            ave_CtoP = sum(carb2P)/count(carb2P.gt.0.0)
         else
            ave_CtoP = 106.0
         endif
         EPmax = maxval(EP)
         print*, "C:P = ",ave_CtoP
      endif
      
      if (ReminTemp_Marsay .or. ReminTemp_Q10 .or. ReminPico) then
         ! Calculate Export Production (EP) and max EP
         EPmax = maxval(EP)
         print*, "Max C flux = ", EPmax
      endif
C-----------------------------------------------------------------------

      DO 380 J=2,JMTM1
C

        DO 774 K = 1,KM
        DO 774 I = 1,IMT
          AHIQ  (I,K)=AHI  (K)  *  AHIFAC(J,K)
          AHHQ  (I,K)=AHH  (K)  *  AHHFAC(J,K)
          aheq  (i,k)=ahe  (k)  *  ahefac(j,k)
  774  CONTINUE

C---------------------------------------------------------------------
C  MOVE ALL SLAB DATA DOWN ONE ROW
C---------------------------------------------------------------------
C

      do nk = 1, nt
        do nj = 1, km
          do ni = 1, imt
            tbm(ni, nj, nk) = tb(ni, nj, nk)
            tm(ni, nj, nk) = t(ni, nj, nk)
            tb(ni, nj, nk) = tbp(ni, nj, nk)
            t(ni, nj, nk) = tp(ni, nj, nk)
          end do
        end do
      end do

      do nj = 1, km
        do ni = 1, imt
          ubm(ni, nj) = ub(ni, nj)
          um(ni, nj) = u(ni, nj)
          ub(ni, nj) = ubp(ni, nj)
          u(ni, nj) = up(ni, nj)
          vbm(ni, nj) = vb(ni, nj)
          vm(ni, nj) = v(ni, nj)
          vb(ni, nj) = vbp(ni, nj)
          v(ni, nj) = vp(ni, nj)
        end do
      end do

      if (robert_time_filter) then
        do i = 1, imt
          ssfub(i) = ssfubp(i)
          ssfvb(i) = ssfvbp(i)
        end do
      end if

C---------------------------------------------------------------------
C  COMPLETE READIN OF J+1 SLAB (EXCEPT LAST ROW)
C---------------------------------------------------------------------
C
      IF(J.NE.JMTM1) THEN

        do kk = 1, nt
          do jj = 1, km
            do ii = 1, imt
              tbp(ii, jj, kk) = odam_t(ii, j+1, jj, kk, ndiskb)
              tp(ii, jj, kk) = odam_t(ii, j+1, jj, kk, ndisk)
            end do
          end do
        end do

        do jj = 1, km
          do ii = 1, imt
            ubp(ii, jj) = odam_u(ii, j+1, jj, ndiskb)
            vbp(ii, jj) = odam_v(ii, j+1, jj, ndiskb)
            up(ii, jj) = odam_u(ii, j+1, jj, ndisk)
            vp(ii, jj) = odam_v(ii, j+1, jj, ndisk)
          end do
        end do

      ENDIF
C
C---------------------------------------------------------------------
C  INITIATE WRITEOUT OF NEWLY COMPUTED DATA FROM PREVIOUS ROW
C---------------------------------------------------------------------
C
      IF(J.GT.2) then

        do kk = 1, nt
          do jj = 1, km
            do ii = 1, imt
              odam_t(ii, j-1, jj, kk, ndiska) = ta(ii, jj, kk)
            end do
          end do
        end do

        do jj = 1, km
          do ii = 1, imt
            odam_u(ii, j-1, jj, ndiska) = ua(ii, jj)
            odam_v(ii, j-1, jj, ndiska) = va(ii, jj)
          end do
        end do

      end if
C
C---------------------------------------------------------------------
C  MOVE TAU-1 DATA TO TAU LEVEL ON A MIXING TIMESTEP
C---------------------------------------------------------------------
C
      IF(MIX.EQ.1) THEN

        DO 337 M=1,NT
        DO 337 K=1,KM
        DO 337 I=1,IMT
          TBP(I,K,M)=TP(I,K,M)
 337    CONTINUE

        DO 338 K=1,KM
        DO 338 I=1,IMT
          UBP(I,K)=UP(I,K)
          VBP(I,K)=VP(I,K)
 338    CONTINUE

      ENDIF
C
C---------------------------------------------------------------------
C  SHIFT MASKS DOWN ONE ROW AND COMPUTE NEW MASKS
C---------------------------------------------------------------------
C

      DO 345 K=1,KM
      DO 345 I=1,IMT
        FMM(I,K)=FM (I,K)
        FM (I,K)=FMP(I,K)
        IF(KMT(I,J+1).GE.KAR(K)) THEN
           FMP(I,K)=1.0
        ELSE
           FMP(I,K)=0.0
        ENDIF
        IF(KMU(I,J).GE.KAR(K)) THEN
           GM(I,K)=1.0
        ELSE
           GM(I,K)=0.0
        ENDIF
 345  CONTINUE

C
C---------------------------------------------------------------------
C  CALL THE MAIN COMPUTATIONAL ROUTINES TO UPDATE THE ROW
C---------------------------------------------------------------------
C
      IF(J.NE.JMTM1) CALL CLINIC(J)
C
         CALL TRACER(J, EPmax)

c...  If using Robert time filter, filter the tracers and baroclinic velocities
c...  to produce new values at time level tau. The unfiltered time level tau
c...  values are then overwritten.

      if (robert_time_filter) then
        
        do m = 1, nt
          do k = 1, km
            do i = 1, imt
               tf(i, k, m) = pnu2m * t(i, k, m) + 
     &                       pnu * (ta(i, k, m) + tb(i, k, m))

               if (N15diagnostics) then
               if (m.eq.n_n15) then
               if (j.eq.latj .and. i.eq.loni .and. k.le.kmt(i,j)) then
                  print*, " "
                  print*, " Robert Time Filtering "
                  print*,"depth=",k,"lat=",j*1.59,"lon=",i*2.8125
                  print*, "pnu2m = ",pnu2m, "pnu = ", pnu
                  print*, "NO3 (tf) = ", tf(i,k,n_no3)
                  print*, "15N (tf) = ", tf(i,k,n_n15)
                  print*, "delta15N = ",(((tf(i,k,n_n15)/(tf(i,k,n_no3)
     +                                   -tf(i,k,n_n15)))/1.0)-1.)*1000
                  print*, " "
               endif
               endif
               endif

               
            end do
          end do
        end do

        if (j .ne. jmtm1) then

          do k = 1, km
            do i = 1, imt
              uf(i, k) = gm(i, k) * ssfub(i)
              uf(i, k) = pnu * (ub(i, k) - uf(i, k) + ua(i, k))
              uf(i, k) = pnu2m * usav(i, k) + uf(i, k)
              vf(i, k) = gm(i, k) * ssfvb(i)
              vf(i, k) = pnu * (vb(i, k) - vf(i, k) + va(i, k))
              vf(i, k) = pnu2m * vsav(i, k) + vf(i, k)
            end do
          end do

        else

          do k = 1, km
            do i = 1, imt
              uf(i, k) = 0.0
              vf(i, k) = 0.0
            end do
          end do

        end if

        do kk = 1, nt
          do jj = 1, km
            do ii = 1, imt
              odam_t(ii, j, jj, kk, ndisk) = tf(ii, jj, kk)
            end do
          end do
        end do

        do jj = 1, km
          do ii = 1, imt
            odam_u(ii, j, jj, ndisk) = uf(ii, jj)
            odam_v(ii, j, jj, ndisk) = vf(ii, jj)
          end do
        end do

      end if

C
C---------------------------------------------------------------------
C  PRINT THE PROGRESSING SOLUTION AT SPECIFIED ROWS ON ENERGY TSTEP
C---------------------------------------------------------------------
C

c...  Compute the following statistics every energy timestep:
c...
c...    1. The mean temperature of the ocean
c...    2. The mean salinity of the ocean
c...    3. The global-mean surface heat flux
c...    4. The global-mean surface salinity tendency
c...    5. The mean temperature of each model level
c...    6. The mean salinity of each model level
c...    7. The RMS horizontal velocity of each model level
c...    8. The RMS vertical velocity of each model level
c...    9. The minimum temperature of each model level

       if (nergy .eq. 1) then

         do k = 1,km
         do i = 2,imtm1
           if(t(i,k,1).ne.0.)then
             tglobe = tglobe + cst(j)*dz(k)*t(i,k,1)
             if(t(i,k,1).lt.tmin(k))then
               tmin(k) = t(i,k,1)
               itmin(k) = i
               jtmin(k) = j
             endif
             sglobe = sglobe + cst(j)*dz(k)*t(i,k,2)
             vglobe = vglobe + cst(j)*dz(k)
             if(k.eq.1)fhglobe = fhglobe + cst(j)*flux(i,j,1)    
             if(k.eq.1)fsglobe = fsglobe + cst(j)*flux(i,j,2)
             tlevel(k,1) = tlevel(k,1) + cst(j)*t(i,k,1)
             tlevel(k,2) = tlevel(k,2) + cst(j)*t(i,k,2)
             tlevel(k,3) = tlevel(k,3) + cst(j)*(u(i,k)**2+v(i,k)**2)
             tlevel(k,4) = tlevel(k,4) + cst(j)*(w(i,k)**2)
             alevel(k) = alevel(k) + cst(j)
           endif
         enddo
         enddo

       end if

C
C  *********    PRINT SOME VALUES AT FREQUENT INTERVALS      ***********

C  *********   SET PARAMETERS FOR PRINTOUT OF FULL SOLUTION   **********
      IF(NERGY.EQ.0.OR.MXP.EQ.1) GO TO 339
      if (j .ne. jplot) go to 8090
      IPRT = IMT
C
C  DETERMINE INDEX OF FIRST T OCEAN POINT
C

      DO 430 I=1,IMT
        ISTRT=I
        IF(KMT(I,J).NE.0) GO TO 431
 430  CONTINUE

 431  CONTINUE
      ISTOP=ISTRT+IPRT-1
      IF(ISTOP.GT.IMT) ISTOP=IMT

      DO 8015 M = 1,2
        IF(M.EQ.1) PRINT 8001,J,ITT
        IF(M.EQ.2) PRINT 8002,J,ITT
 8001   FORMAT(20H TEMPERATURE FOR J =,I4,12H AT TIMESTEP,I7)
 8002   FORMAT(20H SALINITY    FOR J =,I4,12H AT TIMESTEP,I7)
      SCL = 1.0
        IF(M.EQ.2) SCL=1.E-3
        CALL MATRIX(T(1,1,M),IMT,ISTRT,ISTOP,0,KM,SCL)
 8015 CONTINUE

      PRINT 8011,J,ITT
 8011 FORMAT(20H W VELOCITY  FOR J =,I4,12H AT TIMESTEP,I7)
C
C  SET CYCLIC BOUNDARY CONDITION ON W BEFORE PRINTING
C

      DO 433 K=1,KMP1
        W(1  ,K)=W(IMTM1,K)
        W(IMT,K)=W(2    ,K)
 433  CONTINUE

      SCL = 1.E-3
      CALL MATRIX(W,IMT,ISTRT,ISTOP,0,KMP1,SCL)
C
C  DETERMINE INDEX OF FIRST U,V OCEAN POINT
C

      DO 440 I=1,IMTM1
        ISTRT=I
        IF(KMU(I+1,J).NE.0) GO TO 441
 440  CONTINUE

 441  CONTINUE
      ISTOP=ISTRT+IPRT-1
      IF(ISTOP.GT.IMT) ISTOP=IMT
      PRINT 8021, J,ITT
 8021 FORMAT(20H U VELOCITY  FOR J =,I4,12H AT TIMESTEP,I7)
      SCL = 1.0
      CALL MATRIX(U,IMT,ISTRT,ISTOP,0,KM,SCL)
      PRINT 8022, J,ITT
 8022 FORMAT(20H V VELOCITY  FOR J =,I4,12H AT TIMESTEP,I7)
      CALL MATRIX(V,IMT,ISTRT,ISTOP,0,KM,SCL)
C
C---------------------------------------------------------------------
C  COMPUTE THE NORTHWARD TRANSPORT OF EACH TRACER QUANTITY
C  AS WELL AS THE ZONALLY INTEGRATED MERIDIONAL MASS TRANSPORT
C---------------------------------------------------------------------
C
 8090 IF(J.EQ.JMTM1) GO TO 8190

      DO 8092 K=1,KM
        VBR(K)=0.0
        vbredd(k) = 0.
      DO 8092 M=1,NT
        TBRS(K,M)=TBRN(K,M)
        TBRN(K,M)=0.0
 8092 CONTINUE

      IF(J.GT.2) GO TO 8110

      DO 8094 M=1,NT
      DO 8094 K=1,KM
        TBRS(K,M)=0.0
 8094 CONTINUE

      DO 8102 K=1,KM
        TOTDX=0.0
        DO 8100 I=2,IMTM1
          TOTDX=TOTDX+DXT(I)*(FM(I,K))
        DO 8100 M=1,NT
          TBRS(K,M)=TBRS(K,M)+T(I,K,M)*FM(I,K)*DXT(I)
 8100   CONTINUE
        IF(TOTDX.NE.0.0) THEN
          DO 8101 M=1,NT
            TBRS(K,M)=TBRS(K,M)/TOTDX
 8101     CONTINUE
        ENDIF
 8102 CONTINUE

 8110 CONTINUE

      DO 8130 K=1,KM
        TOTDX=0.0
        DO 8120 I=2,IMTM1
          TOTDX=TOTDX+DXT(I)*(FMP(I,K))
          VBR(K)=VBR(K)+V(I,K)*DXU(I)*CS(J)
          vbredd(k) = vbredd(k) + vedd(i,k)*dxu(i)*cs(j)
        DO 8120 M=1,NT
          TBRN(K,M)=TBRN(K,M)+TP(I,K,M)*FMP(I,K)*DXT(I)
 8120   CONTINUE
        IF(TOTDX.NE.0.0) THEN
          DO 8122 M=1,NT
            TBRN(K,M)=TBRN(K,M)/TOTDX
 8122     CONTINUE
        ENDIF
        IF(K.EQ.1) TMT(J,1)=VBR(1)*DZ(1)
        IF(K.GT.1) TMT(J,K)=TMT(J,K-1)+VBR(K)*DZ(K)
      DO 8130 M=1,NT
        TTN(1,J,M)=TTN(1,J,M)+VBR(K)*(TBRN(K,M)+TBRS(K,M))*0.5*DZ(K)
      DO 8130 I=2,IMTM1
        TTN(6,J,M)=TTN(6,J,M)+(V(I,K)*DXU(I)+V(I-1,K)*DXU(I-1))*
     &             (T(I,K,M)+TP(I,K,M))*CS(J)*0.25*DZ(K)
        TTN(7,J,M)=TTN(7,J,M)-ESAV(I,K,M)*FM(I,K)*DXT(I)*CS(J)*DZ(K)
 8130 CONTINUE

c                   sum eddy transport vertically from bottom to top

      tmtedd(j,km) = -vbredd(km)*dz(km)

      do k =kmm1,1,-1
        tmtedd(j,k) = tmtedd(j,k+1)-vbredd(k)*dz(k)
      enddo

      do k = 2,km
        if(tmtedd(j,k).eq.0.)tmtedd(j,k) = -999.9e12
      enddo

      DO 8140 M=1,NT
      DO 8140 I=2,IMTM1
        TOTDZ=0.0
        VBRZ=0.0
        TBRZ=0.0
        DO 8136 K=1,KM
          IF(KMT(I,J).GE.K.AND.KMT(I,J+1).GE.K) THEN
            VBRZ=VBRZ+(V(I,K)*DXU(I)+V(I-1,K)*DXU(I-1))*DZ(K)
            TBRZ=TBRZ+(T(I,K,M)+TP(I,K,M))*DZ(K)
            TOTDZ=TOTDZ+DZ(K)
          ENDIF
 8136   CONTINUE
        IF(TOTDZ.EQ.0.) GO TO 8140
        TBRZ=TBRZ/TOTDZ
        TTN(3,J,M)=TTN(3,J,M)+VBRZ*TBRZ*CS(J)*0.25
        TTN(5,J,M)=TTN(5,J,M)-(WSX(I,J)*DXU(I)+WSX(I-1,J)*DXU(I-1))*
     &             (T(I,1,M)+TP(I,1,M)-TBRZ)*CS(J)/(8.0*OMEGA*SINE(J))
 8140 CONTINUE

      DO 8150 M=1,NT
        TTN(2,J,M)=TTN(6,J,M)-TTN(1,J,M)
        TTN(4,J,M)=TTN(6,J,M)-TTN(3,J,M)-TTN(5,J,M)
        TTN(8,J,M)=TTN(6,J,M)+TTN(7,J,M)
 8150 CONTINUE

 8190 CONTINUE
 339  CONTINUE

c...  Add the model statistics for this timestep to the output arrays
      if (mxp .ne. 0) goto 7312

      if (save_temp) then
        do k = 1, km
          do i = 2, imtm1 
            ohist_temp(i, j, k) = ohist_temp(i, j, k) + t(i, k, 1)
          end do
        end do
       end if

      if (save_sal) then
        do k = 1, km
          do i = 2, imtm1
            ohist_sal(i, j, k) = ohist_sal(i, j, k) + t(i, k, 2)
          end do
        end do
      end if

      if (save_u) then
        do k = 1, km
          do i = 2, imtm1
            ohist_u(i, j, k) = ohist_u(i, j, k) + u(i, k)
          end do
        end do
      end if

      if (save_v .or. save_over) then
        do k = 1, km
          do i = 2, imtm1
            ohist_v(i, j, k) = ohist_v(i, j, k) + v(i, k)
          end do
        end do
      end if

      if (save_w) then
        do k = 1, km
          do i = 2, imtm1
            ohist_w(i, j, k) = ohist_w(i, j, k) + w(i, k)
          end do
        end do
      end if

      if (save_uedd) then
        do k = 1, km
          do i = 2, imtm1
            ohist_uedd(i, j, k) = ohist_uedd(i, j, k) + uedd(i, k)
          end do
        end do
      end if

      if (save_vedd .or. save_over) then
        do k = 1, km
          do i = 2, imtm1
            ohist_vedd(i, j, k) = ohist_vedd(i, j, k) + vedd(i, k)
          end do
        end do
      end if

      if (save_wedd) then
        do k = 1, km
          do i = 2, imtm1
            ohist_wedd(i, j, k) = ohist_wedd(i, j, k) + wedd(i, k)
          end do
        end do
      end if


* rjm - calculate biogeochemical fields each timestep

      if (nt.gt. 2 ) then   ! Do only if biogeochemistry is active
         
         ! Calculate ave_tr (3D tracer fields)
         do m = 1,nt 
            do k = 1, km
               do i = 2, imtm1
                  
                  ave_tr(i,j,k,m) = ave_tr(i,j,k,m) + fm(i,k) *T(I,K,m)

               enddo
            enddo
         enddo

! may want to change what is saved
! 4 extra fields set by nbio2d where first nt are reserved obgc tracer
         trname(nt+1) = "pco2 "
         trname(nt+2) = "POC Export "
         trname(nt+3) = "POP Export "
         trname(nt+4) = "PIC Export "

         trname(nt+5) = "wind squared "
         trname(nt+6) = "solar"
         trname(nt+7) = "sea ice con"
         trname(nt+8) = "Fe dust"
         trname(nt+9) = "Nr depos"
         trname(nt+10) = "Shortwave H"
         trname(nt+11) = "Longwave H"
         trname(nt+12) = "Longwave Up H"
         trname(nt+13) = "Sensible H"
         trname(nt+14) = "Latent H"

         do i = 2,imtm1    ! lat coords 2 - 113
            
            do l = 1,nt
               if (l.le.2) then     ! For T and S fluxes
                  ave_flux(i,j,l) = ave_flux(i,j,l) 
     +                              + fm(i,1)*flux(i,j,l)
               else                 ! For other obgc tracers
                  ave_flux(i,j,l) = ave_flux(i,j,l) 
     +                              + fm(i,1)*fluxgas(i,j,l-2)
               endif
            enddo

            ! Extra 2d fields
            ave_flux(i,j,nt+1) = ave_flux(i,j,nt+1) +
     +                           fm(i,1)*pco2o(i,j,n_dic-2)
            ave_flux(i,j,nt+2) = ave_flux(i,j,nt+2) 
     +                           + fm(i,1) * poc(i,1,j)
            ave_flux(i,j,nt+3) = ave_flux(i,j,nt+3) 
     +                           + fm(i,1) * pop(i,1,j)
            ave_flux(i,j,nt+4) = ave_flux(i,j,nt+4) 
     +                           + fm(i,1) * pic(i,1,j)

*            if (lcouple) then
            ! checks but may be useful for coupled run
               ave_flux(i,j,nt+5) = ave_flux(i,j,nt+5) +
     +                              sqrt(uwnd(i,j)**2+vwnd(i,j)**2)
               ave_flux(i,j,nt+6) = ave_flux(i,j,nt+6) + swdown(i,j)
               ave_flux(i,j,nt+7) = ave_flux(i,j,nt+7) + fice(i,j)
               ave_flux(i,j,nt+8) = ave_flux(i,j,nt+8) + sbcbio(i,j,4)
               ave_flux(i,j,nt+9) = ave_flux(i,j,nt+9) + sbcbio(i,j,5)
               ave_flux(i,j,nt+10) = ave_flux(i,j,nt+10) + swflx(i,j)
               ave_flux(i,j,nt+11) = ave_flux(i,j,nt+11) + lwflx(i,j)
               ave_flux(i,j,nt+12) = ave_flux(i,j,nt+12) + lwup(i,j)
               ave_flux(i,j,nt+13) = ave_flux(i,j,nt+13) + senflx(i,j)
               ave_flux(i,j,nt+14) = ave_flux(i,j,nt+14) + latflx(i,j)
*            endif
            !            
            if (vary_stoich) then
               ave_NtoP(i,j,1) = ave_NtoP(i,j,1) + NtoP(i,j)
               ave_carb2P(i,j,1) = ave_carb2P(i,j,1)+carb2P(i,j)
               ave_o2rem(i,j,1) = ave_o2rem(i,j,1) + o2_rem(i,j)
               ave_no3rem(i,j,1) = ave_no3rem(i,j,1)+no3_rem(i,j)
            endif

            if (vary_stoich .or. ReminPico .or. ReminTemp_Marsay .or.
     +          ReminTemp_Q10) then
               EP(i,j) = poc(i,1,j) * (1e-2*12.0*(c2dtts/2.0))  ! gives mg C m-2 time-1
     +                   / ((c2dtts/2.0)/3600.0)                ! gives mg C m-2 hr-1
            endif

            if (ReminPico .or. ReminTemp_Marsay .or. ReminTemp_Q10) then
               ave_MCb(i,j,1) = ave_MCb(i,j,1) + Mcurve_b(i,j)*fm(i,1)
            endif
     
            do k = 1,km
            
               ave_Pden(i,j,k,1) = ave_Pden(i,j,k,1) + 
     +                             denP(i,k,j)*fm(i,k)
               ave_fix(i,j,k,1) = ave_fix(i,j,k,1)+Nfix(i,k,j)*fm(i,k)
               ave_exp(i,j,k,1) = ave_exp(i,j,k,1)+Nexp(i,k,j)*fm(i,k)
               ave_oxyneg(i,j,k,1) = ave_oxyneg(i,j,k,1) - 
     +                               oxyneg(i,k,j) * fm(i,k)

               if (diff_out) then
                  ave_Udif(i,j,k,1) = ave_Udif(i,j,k,1) + 
     +                                Xdif(i,j,k)*fm(i,k)
                  ave_Vdif(i,j,k,1) = ave_Vdif(i,j,k,1) + 
     +                                Ydif(i,j,k)*fm(i,k)
                  ave_Wdif(i,j,k,1) = ave_Wdif(i,j,k,1) + 
     +                                Wdif(i,j,k)*fm(i,k)
                  ave_Gar84(i,j,k,1) = ave_Gar84(i,j,k,1) + 
     +                                 Gar84(i,j,k)*fm(i,k)
               endif
            
               if (sedfluxes) then
                  ave_opsed(i,j,k,1) = ave_opsed(i,j,k,1) - 
     +                                 op_sed(i,k,j)*fm(i,k)*dizt(k)
                  ave_opsed_f(i,j,k,1) = ave_opsed_f(i,j,k,1) -
     +                                   op_sed_f(i,k,j)*fm(i,k)*dizt(k)
                  ave_icsed(i,j,k,1) = ave_icsed(i,j,k,1) - 
     +                                 ic_sed(i,k,j)*fm(i,k)*dizt(k)
                  ave_sised(i,j,k,1) = ave_sised(i,j,k,1) - 
     +                                 si_sed(i,k,j)*fm(i,k)*dizt(k)
                  ave_Sden(i,j,k,1) = ave_Sden(i,j,k,1) + 
     +                                denS(i,k,j)*fm(i,k)
                  ave_Ssul(i,j,k,1) = ave_Ssul(i,j,k,1) + 
     +                                sulS(i,k,j)*fm(i,k)
               endif

               if (ReminDiagnostics) then
                  ! Organic matter (units are mmol P m-3 timestep-1)
                  ave_pop(i,j,k,1) = ave_pop(i,j,k,1) + 
     +                            (pop(i,k,j)*dizt(k)*(c2dtts/2.0)*1e-2)
                  ave_optot(i,j,k,1) = ave_optot(i,j,k,1) + 
     +                         (op_tot(i,k,j)*dizt(k)*(c2dtts/2.0)*1e-2)
                  ave_oprem(i,j,k,1) = ave_oprem(i,j,k,1) + 
     +                         (op_rem(i,k,j)*dizt(k)*(c2dtts/2.0)*1e-2)
                  ave_opden(i,j,k,1) = ave_opden(i,j,k,1) + 
     +                         (op_den(i,k,j)*dizt(k)*(c2dtts/2.0)*1e-2)
               endif
            
               if (n_n15.gt.0) then
                  ave_d15N(i,j,k,1) = ave_d15N(i,j,k,1) + 
     +                                 d15Norg(i,k,j)*fm(i,1)
               endif
            
            enddo ! k loop

         enddo ! i loop
      
      endif
* rjm - end biogeochemical field computations

 7312 continue

c...  Calculate the depth of convection for this timestep and update the output
c...  array. The depth of convection is indicated by an isopycnal slope in
c...  excess of the maximum value allowed for mixing surfaces. Note that
c...  SLOPE(I, K) is defined at the top of each TS gridbox.
      if (save_cdepthm) then
        do i = 2, imtm1
          kdep = 1
          do k = 2, km
            if (slope(i, k) .le. 1.0/slmxrf) then
              kdep = k - 1
              exit
            end if
          end do
          if (ohist_cdepthm(i, j) .lt. zdz(kdep))
     &      ohist_cdepthm(i, j) = zdz(kdep)
        end do
      end if

C collect uocn and vocn arrays for use by ice dynamics

      do 7440 i=2,imtm1
      uco(i,j)=U(i,1)*gm(i,1)
 7440 vco(i,j)=V(i,1)*gm(i,1)

      uco(1,j)=uco(imtm1,j)
      vco(1,j)=vco(imtm1,j)
      uco(imt,j)=uco(2,j)
      vco(imt,j)=vco(2,j)

 380  CONTINUE
     
      if (NO3diagnostics) then
      do j = 2,jmtm1
         do i = 1,imtm1
            do k = 1,km

               vol = (dxt(i)*dyt(j)*dzt(k)/1e6)
               
               o2loss(i,k,j) = -oxyneg(i,k,j)*vol
               Nl_Pden(i,k,j) = denP(i,k,j)*vol
               Nl_Sden(i,k,j) = denS(i,k,j)*vol
               no3fix(i,k,j) = Nfix(i,k,j)*vol
               no3dif(i,k,j) = no3change(i,k,j)*vol
               budget(i,k,j) = ( (no3fix(i,k,j)) -
     +                           Nl_Pden(i,k,j) - Nl_Sden(i,k,j) )

            enddo
         enddo
      enddo

      write(6,100) 'O2 Loss = ', sum(o2loss)
      write(6,100) 'NO3 Pden = ', sum(Nl_Pden)
      write(6,100) 'NO3 Sden = ', sum(Nl_Sden)
      write(6,100) 'NO3 NFix = ', sum(no3fix)
      write(6,100) 'NO3 change = ', sum(no3dif)
      write(6,100) 'Budget = ', sum(budget)
  100 format(A,E14.7)
      open(unit=2, file='NO3budget.dat')
      ttyear = (itt*dttsf)/(86400.*365.)
      write(2,101) ttyear, sum(budget) 
  101 format(2(E14.7,1x))
      endif ! NO3 diagnostics

      do j = 2,jmtm1
         do i = 1,imtm1
            OM(i,j) = 0.0
            do k = 1,km
               OM(i,j) = OM(i,j) + op_tot(i,k,j)+op_sed(i,k,j)
            enddo
            if (ABS(OM(i,j)).gt.epsilon(1e-10)) then ! epsilon value = 2.220446049250313E-016
               print*, ' '
               print*, 'Total Organic Matter not conserved!'
               print*, ' '
               print*, 'Sum of exchanges through water = ',OM(i,j)
               print*, 'op_tot = ',op_tot(i,k,j),
     +                 'op_sed = ',op_sed(i,k,j)
               print*, 'i = ',i,' j = ',j
               print*, 'lon = ',i*2.8125,' lat = ',j*1.59-90
               print*, ' '
            endif
            
            OM(i,j) = 0.0
            do k = 1,km
               OM(i,j) = OM(i,j) + op_tot_f(i,k,j)+op_sed_f(i,k,j)
            enddo
            if (ABS(OM(i,j)).gt.epsilon(1e-10)) then ! epsilon value = 2.220446049250313E-016
               print*, ' '
               print*, 'Total Organic Matter not conserved!'
               print*, ' '
               print*, 'Sum of exchanges through water = ',OM(i,j)
               print*, 'op_tot = ',op_tot_f(i,k,j),
     +                 'op_sed = ',op_sed_f(i,k,j)
               print*, 'i = ',i,' j = ',j
               print*, 'lon = ',i*2.8125,' lat = ',j*1.59-90
               print*, ' '
            endif
         enddo
      enddo

C
C=======================================================================
C  END ROW-BY-ROW COMPUTATION  =========================================
C=======================================================================
C
C---------------------------------------------------------------------
C  PRINT ONE LINE OF TIMESTEP INFORMATION ON SPECIFIED TIMESTEPS
C---------------------------------------------------------------------
C
C  *******    WRITE OUT INFO AT FREQUENT INTERVALS     *************
      IF(EB.AND.MIX.EQ.1) GO TO 390
      IF(MOD(ITT,NTSI).EQ.0) THEN
        EKTOT=EKTOT/VOLUME

        DO 381 M=1,NT
          DTABS(M)=DTABS(M)/VOLUME
 381    CONTINUE

        DAYSYR=365.00
        TTYEAR=TTSEC/(3600.*24.*DAYSYR)
        TTDAY=TTSEC/(3600.*24.)
        TTDAY=MOD(TTDAY,DAYSYR)
        print 912, itt, ttyear, ttday, ektot, dtabs(1), dtabs(2), mscan
  912   format (" TIMESTEP=", i10, " YEAR=", f10.3, " DAY=", f8.3,
     &          " ENERGY=", 1pe13.6, " DT=", 1pe13.6, " DS=", 1pe13.6,
     &          " SCANS=", i5)
      ENDIF

C ..................................................... begin .. ach .. 9/2/95

c...  Display the following statistics every energy timestep:
c...
c...    1. The mean temperature of the ocean
c...    2. The mean salinity of the ocean
c...    3. The global-mean surface heat flux
c...    4. The global-mean surface salinity tendency
c...    5. The mean temperature of each model level
c...    6. The mean salinity of each model level
c...    7. The RMS horizontal velocity of each model level
c...    8. The RMS vertical velocity of each model level
c...    9. The minimum temperature of each model level

      if (nergy .eq. 1) then

      print 914,(tglobe/vglobe),(sglobe*1000/vglobe+35.),
     &          (fhglobe/alevel(1)),(fsglobe/alevel(1))
  914 format(/,' ocean average temp=',f10.8,'  sal=',f11.8,'  fh =',
     &   1pe10.3,'  fs =',1pe10.3)

      do 916 k = 1,km
        tlevel(k,1) = tlevel(k,1)/alevel(k)
        tlevel(k,2) = tlevel(k,2)*1000/alevel(k) + 35.
        tlevel(k,3) = sqrt(tlevel(k,3)/alevel(k))
        tlevel(k,4) = sqrt(tlevel(k,4)/alevel(k))*1000
  916 continue

      print 917,(k,(zdzz(k)/100),alevel(k),(tlevel(k,m),m=1,4),k=1,km)
  917 format(/,'  k   depth     area',8x,'temp',11x,'sal',9x,'rms v',
     &   7x,'rms w',21(/,i3,f8.1,f11.2,f12.4,f14.7,2f12.4),/)
      do 1918 k = 1,km
       print 918,tmin(k),itmin(k),jtmin(k),k
  918  format(' min temperature =',f7.2,'    at point (',2i4,i3,')')
 1918 continue
      print *, " "

      end if

C ........................................................................ end
C
C---------------------------------------------------------------------
C  COMPLETE AND PRINT THE ON-LINE INTEGRALS ON ENERGY TIMESTEPS
C---------------------------------------------------------------------
C
      IF(NERGY.EQ.0) GO TO 390
C
C  1ST, NORMALIZE PREVIOUSLY COMPUTED INTEGRALS BY VOLUME
C

      DO 382 LL=1,8
        ENGINT(LL)=ENGINT(LL)/VOLUME
        ENGEXT(LL)=ENGEXT(LL)/VOLUME
 382  CONTINUE

      DO 383 M=1,NT
        TVAR(M)=TVAR(M)/VOLUME
      DO 383 LL=1,6
        TTDTOT(LL,M)=TTDTOT(LL,M)/VOLUME
 383  CONTINUE

      BUOY=BUOY/VOLUME
C
C  2ND, COMPUTE RESIDUAL TERMS
C
      PLICIN=ENGINT(1)-ENGINT(2)-ENGINT(3)-ENGINT(4)
     &      -ENGINT(5)-ENGINT(6)-ENGINT(7)-ENGINT(8)
      PLICEX=ENGEXT(1)-ENGEXT(2)-ENGEXT(3)-ENGEXT(4)
     &      -ENGEXT(5)-ENGEXT(6)-ENGEXT(7)-ENGEXT(8)

      DO 384 M=1,NT
        TTDTOT(6,M)=TTDTOT(1,M)-TTDTOT(2,M)-TTDTOT(3,M)
     &             -TTDTOT(4,M)-TTDTOT(5,M)
 384  CONTINUE

C
C  3RD, PRINT THE INTEGRALS
C
      PRINT 9100
      PRINT 9101,ENGINT(1),ENGEXT(1),TTDTOT(1,1),TTDTOT(1,2)
      PRINT 9102,ENGINT(2),ENGEXT(2),TTDTOT(2,1),TTDTOT(2,2)
      PRINT 9103,ENGINT(3),ENGEXT(3),TTDTOT(3,1),TTDTOT(3,2)
      PRINT 9104,ENGINT(4),ENGEXT(4),TTDTOT(4,1),TTDTOT(4,2)
      PRINT 9105,ENGINT(5),ENGEXT(5),TTDTOT(5,1),TTDTOT(5,2)
      PRINT 9106,ENGINT(6),ENGEXT(6),TTDTOT(6,1),TTDTOT(6,2)
      PRINT 9109,PLICIN,PLICEX,TVAR(1),TVAR(2)
      PRINT 9107,ENGINT(7),ENGEXT(7)
      PRINT 9108,ENGINT(8),ENGEXT(8)
 9100 FORMAT( 1X,50HWORK BY:              INTERNAL MODE  EXTERNAL MODE,
     &       10X,50H                       TEMPERATURE     SALINITY   )
 9101 FORMAT( 1X,20HTIME RATE OF CHANGE ,2(1PE15.6),
     &       10X,20HTIME RATE OF CHANGE ,2(1PE15.6))
 9102 FORMAT( 1X,20HHORIZONTAL ADVECTION,2(1PE15.6),
     &       10X,20HHORIZONTAL ADVECTION,2(1PE15.6))
 9103 FORMAT( 1X,20HVERTICAL ADVECTION  ,2(1PE15.6),
     &       10X,20HVERTICAL ADVECTION  ,2(1PE15.6))
 9104 FORMAT( 1X,20HHORIZONTAL FRICTION ,2(1PE15.6),
     &       10X,20HHORIZONTAL DIFFUSION,2(1PE15.6))
 9105 FORMAT( 1X,20HVERTICAL FRICTION   ,2(1PE15.6),
     &       10X,20HSURFACE FLUX        ,2(1PE15.6))
c SURFACE FLUX = VERTICAL DIFFUSION
 9106 FORMAT( 1X,20HPRESSURE FORCES     ,2(1PE15.6),
     &       10X,20HTRUNCATION ERROR    ,2(1PE15.6))
 9107 FORMAT( 1X,20HWORK BY WIND        ,2(1PE15.6))
 9108 FORMAT( 1X,20HBOTTOM DRAG         ,2(1PE15.6))
 9109 FORMAT( 1X,20HIMPLICIT EFFECTS    ,2(1PE15.6),
     &       10X,20HCHANGE OF VARIANCE  ,2(1PE15.6))
      TVAR(1)=BUOY-ENGINT(6)-ENGEXT(6)
      DTABS(1)=ENGINT(2)+ENGINT(3)+ENGEXT(2)+ENGEXT(3)
      PRINT 9110,BUOY,TVAR(1),DTABS(1)
 9110 FORMAT(1X,25HWORK BY BUOYANCY FORCES  ,1PE15.6,5X,25HENERGY CONVER
     &SION ERROR  ,1PE15.6,5X,25HNONLINEAR EXCHANGE ERROR ,1PE15.6)
C
C---------------------------------------------------------------------
C  PRINT THE NORTHWARD TRANSPORT OF HEAT AND SALT
C---------------------------------------------------------------------
C
      PRINT 8195
 8195 FORMAT(/,' NORTHWARD TRANSPORT OF HEAT (X10**15 WATTS)',24X,'NORTH
     &WARD TRANSPORT OF SALT (X10**10 CM**3/SEC)',/,6X,'X MEAN  X EDDY
     &Z MEAN  Z EDDY   EKMAN TOT ADV  DIFFUS   TOTAL   X MEAN  X EDDY  Z
     & MEAN  Z EDDY   EKMAN TOT ADV  DIFFUS   TOTAL')
C
C  CONVERT HEAT TRANSPORT TO PETAWATTS, SALT TRNSPT TO 10**10 CM**3/SEC
C

      DO 8198 J=1,JMT
      DO 8198 LL=1,8
        TTN(LL,J,1)=TTN(LL,J,1)*4.186E-15
        TTN(LL,J,2)=TTN(LL,J,2)*1.E-10
 8198 CONTINUE

      DO 8197 JJ=2,JMTM2
        J=JMT-JJ
        PRINT 8196,J,(TTN(I,J,1),I=1,8),(TTN(I,J,2),I=1,8)
 8196   FORMAT(I4,8F8.3,1X,8F8.3)
 8197 CONTINUE

      PRINT 8194
 8194 FORMAT(/,' MERIDIONAL MASS TRANSPORT')
      SCL=1.E12
      CALL MATRIX(TMT,JMT,1,JMT,0,KM,SCL)

      print 475,'GLOBAL MERIDIONAL EDDY TRANSPORT     (SV)'
  475 format(/1x,a70)
      call matrix(tmtedd,jmt,1,jmtm2,0,km,scl)

 390  CONTINUE
C
C---------------------------------------------------------------------
C  INITIATE WRITEOUT OF NEWLY COMPUTED DATA FROM THE FINAL ROW
C---------------------------------------------------------------------
C

      do kk = 1, nt
        do jj = 1, km
          do ii = 1, imt
            odam_t(ii, jmt-1, jj, kk, ndiska) = ta(ii, jj, kk)
          end do
        end do
      end do

      do jj = 1, km
        do ii = 1, imt
          odam_u(ii, jmt-1, jj, ndiska) = ua(ii, jj)
          odam_v(ii, jmt-1, jj, ndiska) = va(ii, jj)
        end do
      end do

c...  Check for conservation of tracers, if required
      if (check_conservation) then

c.....  Calculate mean temperature and salinity of ocean at t+1 time level
        tint = 0.0
        sint = 0.0
        vint = 0.0
        do j = 2, jmtm1
          do i = 2, imtm1
            do k = 1, kmt(i, j)
              dvol = cst(j) * dxt(i) * dyt(j) * dz(k)
              tint = tint + dvol * odam_t(i, j, k, 1, ndiska)
              sint = sint + dvol * odam_t(i, j, k, 2, ndiska)
              vint = vint + dvol
            end do
          end do
        end do
        tbar2 = tint / vint
        sbar2 = 35.0 + 1000.0 * (sint / vint)

c.....  Integrate surface fluxes into ocean for current timestep, and calculate
c.....  equivalent changes in mean temperature and salinity of ocean.
        hint = 0.0
        sint = 0.0
        do j = 2, jmtm1
          do i = 2, imtm1
            if (kmt(i, j) .gt. 0) then
              darea = cst(j) * dxt(i) * dyt(j)
              hint = hint + darea * flux(i, j, 1)
              sint = sint + darea * flux(i, j, 2)
            end if
          end do
        end do
        dtf = 100.0 * c2dtts * hint / (rhocw * volume)
        dsf = 1000.0 * c2dtts * sint * dz(1) / volume

c.....  Display conservation statistics
        write (*, *)
        write (*, *) "ITT = ", itt
        write (*, *)
        write (*, *) "Temperature"
        write (*, *) "-----------"
        write (*, *) "Global mean at t-1 = ", tbar1, " degC"
        write (*, *) "Global mean at t+1 = ", tbar2, " degC"
        write (*, *)
        write (*, *) "Actual change      = ", tbar2-tbar1, " degC"
        write (*, *) "Expected change    = ", dtf, " degC"
        write (*, *)
        write (*, *) "Conservation error = ", tbar2-tbar1-dtf, " degC"
        write (*, *)
        write (*, *) "Salinity"
        write (*, *) "--------"
        write (*, *) "Global mean at t-1 = ", sbar1, " psu"
        write (*, *) "Global mean at t+1 = ", sbar2, " psu"
        write (*, *)
        write (*, *) "Actual change      = ", sbar2-sbar1, " psu"
        write (*, *) "Expected change    = ", dsf, " psu"
        write (*, *)
        write (*, *) "Conservation error = ", sbar2-sbar1-dsf, " psu"
        write (*, *)

      end if

C
C---------------------------------------------------------------------
C  SOLVE FOR THE NEW STREAM FUNCTION
C---------------------------------------------------------------------
C
      CALL RELAX
C
C---------------------------------------------------------------------
C  IF THIS IS THE END OF THE 1ST PASS OF AN EULER BACKWARD TIMESTEP,
C  SET THE INPUT DISC UNITS SO THAT THE PROPER LEVELS ARE FETCHED ON
C  THE NEXT PASS.  THE OUTPUT FOR THE 2ND PASS WILL BE PLACED ON THE
C  "NDISKA" UNIT.  RETURN TO THE TOP OF "STEP" TO DO THE 2ND PASS.
C---------------------------------------------------------------------
C
      IF(MIX.EQ.1 .AND. EB) THEN
        MIX=0
        MXP=1
        NDISKX=NDISKB
        NDISKB=NDISK
        NDISK=NDISKA
        NDISKA=NDISKX
        GO TO 182
      ENDIF
C
C---------------------------------------------------------------------
C  PRINT THE STREAM FUNCTION ON AN ENERGY TIMESTEP
C---------------------------------------------------------------------
C
C
      IF(NERGY.EQ.1) THEN
        PRINT 8000,ITT
 8000   FORMAT(' STREAM FUNCTION IN SVERDRUPS, TS=',I6)
        SCL=1.E12
        CALL MATRIX(P,IMT,2,IMTM1,JMT,0,SCL)
      ENDIF

c...  Add the barotropic streamfunction for this timestep to the output array
      if (save_res) then
        do j = 2, jmtm1
          do i = 2, imtm1
            ohist_res(i, j) = ohist_res(i, j) + p(i, j)
          end do
        end do
      end if

      RETURN
      END
