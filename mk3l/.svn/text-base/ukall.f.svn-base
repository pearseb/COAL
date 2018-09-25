c Add calls to INT to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/15
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Removed several unnecessary "include 'MACHINE.f'"
c SJP 2001/11/22
c
c $Log: ukall.f,v $
c Revision 1.18  2001/10/12 02:06:57  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.17  2001/06/04 02:26:58  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.16  2001/03/07 04:29:01  rot032
c Changes from LDR to bring sulfur cycle to run N63
c
c Revision 1.15  2001/02/28 04:40:57  rot032
c Bring sulfur cycle to run N58.
c
c Revision 1.14  2001/02/22 05:56:38  rot032
c Changes from LDR to sulfur cycle, to bring it to run N52
c
c Revision 1.13  2000/12/08 03:58:53  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.12  2000/11/14 06:50:42  rot032
c Merge of HBG and LDR changes.
c
c Revision 1.11  2000/11/14 03:11:38  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.10.1.1  2000/11/14 06:47:45  rot032
c Changes from LDR for aerosol model
c
c Revision 1.10  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.9  2000/06/21 06:57:51  rot032
c Change order of declarations for NEC.
c
c Revision 1.8  2000/06/20 02:27:10  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.7  1999/05/20 06:23:54  rot032
c HBG changes to V5-2
c
c Revision 1.6  1998/12/10  00:55:46  ldr
c HBG changes to V5-1-21
c
c Revision 1.5  1998/05/26  05:29:37  ldr
c Final changes (inc. ACH salinity stuff) for 10k run.
c
c Revision 1.4  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.3  1997/07/24  05:42:51  ldr
c Tuning of ukconv version and T63 fixes from HBG
c
c Revision 1.2  1997/06/24  06:07:21  ldr
c Changes from HBG to ukconv version to improve SWCF, LWCF in tropics.

C*****************************************************************
CNOTE - To save on the number of subroutines, all of the routines below 
C     - are now kept in one file called "ukall.f"
C*****************************************************************
C     |
C     |
C  CONVECU -(convec2a.f)
C          |
C         QSATU    (qsat1a.f)
C         FLAG_WET (flagw1a.f)
C         LAYER_CN (laycn1a.f)
C         DQS_DTH  (dqsdth1a.f)
C         LIFT_PAR (lifpar1a.f)
C          |    |
C          |   QSATU    (qsat1a.f)
C          |   LATENT_H (latent1a.f)
C          |
C          |
C         CONVEC2  (conv21a.f)
C          |    |
C          |   PARCEL   (parcel1a.f)
C          |    |   |
C          |    |  DETRAIN  (detrai1a.f)
C          |    |   |   |
C          |    |   |  THP_DET  (thpdet1a.f)
C          |    |   |  QSATU    (qsat1a.f)
C          |    |   |  THETAR   (thetar1a.f)
C          |    |   |   |   |
C          |    |   |   |  QSATU    (qsat1a.f)
C          |    |   |   |
C          |    |   |  DET_RATE (detrat1a.f)
C          |    |   |
C          |    |  TERM_CON (termco1a.f)
C          |    |  CLOUD_W  (cloudw1a.f)
C          |    |       |
C          |    |      CON_RAD  (conrad1a.f)
C          |    |
C          |   ENVIRON  (enviro1a.f)
C          |
C          |
C         DD_CALL  (ddcall2a.f)
C          |    |
C          |   FLX_INIT (flxini2a.f)
C          |   LAYER_DD (layerd2a.f)
C          |   DD_INIT  (ddinit2a.f)
C          |    |   |
C          |    |   |
C          |    |  SATCAL   (satcal2a.f)
C          |    |       |
C          |    |      QSATU    (qsat1a.f)
C          |    |      DQS_DTH  (dqsdth1a.f)
C          |    |
C          |   DOWND    (downd2a.f)
C          |        |
C          |       DDRAUGHT (ddraug2a.f) 
C          |        |   |
C          |        |  SATCAL   (satcal2a.f)
C          |        |   |   |
C          |        |   |  QSATU    (qsat1a.f)
C          |        |   |  DQS_DTH  (dqsdth1a.f)
C          |        |   |
C          |        |  CRS_FRZL (crsfrz2a.f)
C          |        |  QSATU    (qsat1a.f)
C          |        |  DEVAP    (devap2a.f)
C          |        |   |   |
C          |        |   |  EVP      (evp2a.f)
C          |        |   |
C          |        |  TERMDD   (termdd2a.f)
C          |        |  DD_ENV   (ddenv2a.f)
C          |        |
C          |       CHG_PHSE (chgphs2a.f)
C          |       PEVP_BCB (pevbcb2a.f)
C          |            |
C          |           QSATU    (qsat1a.f)
C          |           EVP      (evp2a.f)
C          |           SATCAL   (satcal2a.f)
C          |                |
C          |               QSATU    (qsat1a.f)
C          |               DQS_DTH  (dqsdth1a.f)
C          |
C          |
C         COR_ENGY (coreng1a.f)
C
C---------------------------------------------------------------------
C





CLL  SUBROUTINE CHG_PHSE----------------------------------------------
CLL
CLL  PURPOSE : CHANGE OF PHASE ROUTINE FOR POINTS WHERE NO
CLL            DOWNDRAUGHT OCCURING
CLL
CLL            UPDATES POTENTIAL TEMPERATURE OF LAYER K
CLL            AS PRECIPITATION CHANGES PHASE IN SITU
CLL
CLL            ADD LATENT HEATING WHERE PRECIPITATION CROSSES A
CLL            MELTING OR FREEZING LEVEL
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  PROJECT TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CHG_PHSE (NPNTS,RAIN,SNOW,DTHBYDT_KM1,
     &                     EXK,EXKM1,DELPKM1,THE_K,THE_KM1)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include 'C_LHEAT.f'
CHBG  include 'C_R_CP.f'
CHBG  include 'C_G.f'
CHBG  include 'C_0_DG_C.f'
      include 'UKPARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS                ! IN VECTOR LENGTH
C
      INTEGER I                    ! LOOP COUNTER
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL EXK(ln2)                ! IN EXNER RATIO FOR LAYER K
C
      REAL EXKM1(ln2)              ! IN EXNER RATIO FOR LAYER K-1
C
      REAL DELPKM1(ln2)            ! IN PRESSURE DIFFERENCE ACROSS
                                   !    LAYER K-1 (PA)
C
      REAL THE_K(ln2)              ! IN POTENTIAL TEMPERATURE OF
                                   !    ENVIRONMENT IN LAYER K
C
      REAL THE_KM1(ln2)            ! IN POTENTIAL TEMPERATURE OF
                                   !    ENVIRONMENT IN LAYER K-1
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL RAIN(ln2)               ! INOUT
                                   ! IN  AMOUNT OF FALLING RAIN
                                   !     (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     RAIN (KG/M**2/S)
C
      REAL SNOW(ln2)               ! INOUT
                                   ! IN  AMOUNT OF FALLING SNOW
                                   !     (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     SNOW (KG/M**2/S)
C
      REAL DTHBYDT_KM1(ln2)        ! INOUT
                                   ! IN  INCREMENT TO MODEL POTENTIAL
                                   !     TEMPERATURE IN LAYER K-1
                                   ! OUT UPDATED INCREMENT TO MODEL
                                   !     POTENTIAL TEMPERATURE IN LAYER
                                   !     K-1 DUE TO CHANGE OF PHASE
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C---------------------------------------------------------------------
C
      REAL FACTOR                  ! USED IN THE CALCULATION OF
                                   ! CHANGE OF PHASE OF FALLING
                                   ! PRECIPITATION
C
      LOGICAL BPPNWT_K             ! MASK WHERE PRECIPITATION IS LIQUID
                                   ! IN LAYER K
C
      LOGICAL BPPNWT_KM1           ! MASK WHERE PRECIPITATION IS LIQUID
                                   ! IN LAYER K-1
C
CL
CL----------------------------------------------------------------------
CL  ADD LATENT HEATING WHERE PRECIP CROSSES A MELTING OR FREEZING LEVEL
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (11), EQUATION (42)
CL----------------------------------------------------------------------
CL
      DO I=1,NPNTS
        BPPNWT_K = THE_K(I).GT.TM/EXK(I)
        BPPNWT_KM1 = THE_KM1(I).GT.TM/EXKM1(I)
        FACTOR = LF*G/(EXKM1(I)*CP*DELPKM1(I))
C FREEZE
        IF (.NOT.BPPNWT_KM1.AND.(BPPNWT_K.OR.RAIN(I).GT.0.0)) THEN
           DTHBYDT_KM1(I) = DTHBYDT_KM1(I)+RAIN(I)*FACTOR
           SNOW(I) = SNOW(I)+RAIN(I)
           RAIN(I) = 0.0
        END IF
C MELT
        IF (BPPNWT_KM1.AND.(.NOT.BPPNWT_K.OR.SNOW(I).GT.0.0)) THEN
           DTHBYDT_KM1(I) = DTHBYDT_KM1(I)-SNOW(I)*FACTOR
           RAIN(I) = RAIN(I)+SNOW(I)
           SNOW(I) = 0.0
        END IF
      END DO
C
      RETURN
      END
C
CLL  SUBROUTINE CLOUD_W------------------------------------------------
CLL
CLL  PURPOSE : CLOUD MICROPHYSICS ROUTINE
CLL
CLL            CALCULATES PRECIPITATION PRODUCED IN LIFTING PARCEL
CLL            FROM LAYER K TO K+1
CLL
CLL            CALL CON_RAD TO CALCULATE PARAMETERS FOR RADIATION
CLL            CALCULATION
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  3.2  8/07/93 : added convective cloud condensed water diagnostic
CLL               : P Inness
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  PROJECT TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CLOUD_W (K,NPNTS,XPKP1,PREKP1,XSQKP1,BLOWST,
     &                    FLXKP1,XPK,THEKP1,QEKP1,BWKP1,BLAND,
     &                    QSEKP1,BGMKP1,BTERM,CCA,ICCB,ICCT,TCW,DEPTH,
     &                    EKP14,EKP34,DELEXKP1,CCLWP,DELPKP1,
     &                    CCW,CCWI,CCWL,
     &                    qcloud,ipass,thpkp1,exkp1,pkp1,xs,fs,DELPK)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include 'CRITDEP.f'
CHBG  include 'MPARWTR.f'
CHBG  include 'C_G.f'
CHBG  include 'C_R_CP.f'
CHBG  include 'C_EPSLON.f'
      include 'UKPARAMS.f'
      include 'ECPARM.f'

C
C----------------------------------------------------------------------
C VECTOR LENGTH AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER K              ! IN PRESENT MODEL LAYER
C
      INTEGER I              ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL THEKP1(ln2)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEKP1(ln2)        ! IN MIXING RATIO OF CLOUD ENVIRONMENT
                             !    IN LAYER K+1 (KG/KG)
C
      REAL QSEKP1(ln2)       ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL XPK(ln2)          ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
C
      LOGICAL BWKP1(ln2)     ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      LOGICAL BGMKP1(ln2)    ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K+1
C
      LOGICAL BLAND(ln2)     ! IN LAND/SEA MASK
C
      LOGICAL BTERM(ln2)     ! IN MASK FOR PARCELS WHICH TERMINATE IN
                             !    LAYER K+1
C
      LOGICAL BLOWST(ln2)    ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
C
      REAL FLXKP1(ln2)       ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      REAL XSQKP1(ln2)       ! IN EXCESS PARCEL MIXING RATIO IN
                             !    LAYER K+1 (KG/KG)
C
      REAL EKP14(ln2)        ! IN ENTRAINMENT RATE AT LEVEL K+1/4
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
C
      REAL EKP34(ln2)        ! IN ENTRAINEMNT RATE AT LEVEL K+3/4
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
C
      REAL DELEXKP1(ln2)     ! IN DIFFERENCE IN EXNER RATIO ACROSS
                             !    LAYER K+1 (PA)
C
      REAL DELPKP1(ln2)      ! IN PRESSURE DIFFERENCE ACROSS LAYER K+1
C
      real thpkp1(ln2)       ! IN Parcel potential tempeture in layer K+1
C
      real exkp1(ln2)        ! IN Exner function in layer K+1
C
      real pkp1(ln2)         ! IN Air pressure in layer K+1
C
      logical qcloud         ! IN qcloud = T or F
      integer ipass          ! IN ipass  = 1 or 2 (2=leads)

      real xs(ln2)           ! IN sulfate mixing ratio
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL TCW(ln2)          ! INOUT
                             ! IN  TOTAL CONDENSED WATER SUMMED UPTO
                             !     LAYER K (KG/M**2/S)
                             ! OUT TOTAL CONDENSED WATER SUMMED UPTO
                             !     LAYER K+1 (KG/M**2/S)
C
      REAL DEPTH(ln2)        ! INOUT
                             ! IN  DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K (M)
                             ! OUT DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K+1 (M)
C
      REAL CCLWP(ln2)        ! INOUT
                             ! IN  CONDENSED WATER PATH SUMMED UPTO
                             !     LAYER K (KG/M**2)
                             ! OUT CONDENSED WATER PATH SUMMED UPTO
                             !     LAYER K+1 (KG/M**2)

      real fs(ln2,ntrac)     ! INOUT Scavenging Fraction
C
      REAL DELPK(ln2)        ! IN PRESSURE DIFFERENCE ACROSS LAYER K
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL PREKP1(ln2)       ! OUT PRECIPITATION FROM PARCEL AS IT
                             !     RISES FROM LAYER K TO K+1 (KG/M**2/S)
C
      REAL XPKP1(ln2)        ! OUT PARCEL CLOUD WATER IN LAYER K+1
                             !     (KG/KG)
C
      REAL CCA(ln2)          ! OUT CONVECTIVE CLOUD AMOUNT (%)
C
      INTEGER ICCB(ln2)      ! OUT CONVECTIVE CLOUD BASE LEVEL
C
      INTEGER ICCT(ln2)      ! OUT CONVECTIVE CLOUD TOP LEVEL
C
      REAL CCW(ln2)          ! OUT CONVECTIVE CLOUD LIQUID WATER
                             ! (KG/KG) ON MODEL LEVELS
C
      REAL CCWI(ln2)         ! OUT CONVECTIVE CLOUD ICE
                             ! (KG/KG) ON MODEL LEVELS
C
      REAL CCWL(ln2)         ! OUT CONVECTIVE CLOUD LIQUID
                             ! (KG/KG) ON MODEL LEVELS
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE LOCALLY DEFINED
C----------------------------------------------------------------------
C
      REAL DCRIT             ! CRITICAL DEPTH AT WHICH PRECIPITATION
                             ! MAY FORM (M)
C
      REAL XMIN              ! AMOUNT OF CLOUD WATER RETAINED BY THE
                             ! PARCEL ON PRECIPITATION (KG/KG)
C
      REAL EPSS              ! (1.0+EKP14)*(1.0+EKP34)
C
      REAL TEMPRY,FRACOUT

      real tt(ln2)           ! temperature
      real rho(ln2)          ! Air density  (kg/m3)
      real f_so2(ln2)        ! Solubility of SO2

      real xpold(ln2)        ! L.W. mixing ratio of parcel before precip
      integer nt


      REAL ZQTP1,ZE2,ZE3,ZLWC,ZFAC,ZSO4L,ZSO2L
      REAL ZZA,ZZB,ZZP,ZZQ,ZZP2,ZHP,ZQHP,ZHENEFF,P_SO2
C
C----------------------------------------------------------------------
C  EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL CON_RAD
C
C*---------------------------------------------------------------------
CL
CL----------------------------------------------------------------------

CL  CALCULATE CLOUD WATER BEFORE PRECIPITATION
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (2B), EQUATION (13A)
CL----------------------------------------------------------------------
CL
      DO 10 I=1,NPNTS
       EPSS = (1.+EKP14(I))*(1.+EKP34(I))
       XPKP1(I) = (XPK(I)/EPSS) + XSQKP1(I)
       CCWI(I)=0.0
       CCWL(I)=0.0
   10 CONTINUE

c Calculate parcel T, rho

      do i=1,npnts
        tt(i)=thpkp1(i)*exkp1(i)
        rho(i) = pkp1(i) / (r*tt(i))
        f_so2(i)=0.
      enddo

C CALCULATE THE SOLUBILITY OF SO2
C TOTAL SULFATE  IS ONLY USED TO CALCULATE THE PH OF CLOUD WATER
C

      do i=1,npnts
        if(xpkp1(i).gt.0..and.BWKP1(I))then
          ZQTP1=1./tt(i)-1./298.
          ZE2=1.23*EXP(3020.*ZQTP1)
          ZE3=1.2E-02*EXP(2010.*ZQTP1)

          ZLWC=xpkp1(i)        !m.r. of l.w. before precip
          ZFAC=1000./(ZLWC*32.064)
          ZSO4L=xs(i)*ZFAC
          ZSO4L=AMAX1(ZSO4L,0.)
          ZSO2L=xs(i)*ZFAC
          ZSO2L=AMAX1(ZSO2L,0.)
          ZZA=ZE2*8.2E-02*tt(i)*ZLWC*rho(i)*1.E-03
          ZZB=2.5E-06+ZSO4L
          ZZP=(ZZA*ZE3-ZZB-ZZA*ZZB)/(1.+ZZA)
          ZZQ=-ZZA*ZE3*(ZZB+ZSO2L)/(1.+ZZA)
          ZZP=0.5*ZZP
          ZZP2=ZZP*ZZP
          ZHP=-ZZP+SQRT(ZZP2-ZZQ)
          ZQHP=1./ZHP
          ZHENEFF=1.+ZE3*ZQHP
          P_SO2=ZZA*ZHENEFF
          F_SO2(i)=P_SO2/(1.+P_SO2)
          F_SO2(i)=min(max(0.,F_SO2(i)),1.)
c          write(26,'(2f12.3)')tt(i),f_so2(i)
        endif
      enddo

CL
CL----------------------------------------------------------------------
CL STORE CONVECTIVE CLOUD LIQUID WATER BEFORE PRECIPITATION
CLDR Average this later on to include after precipitation value...
CL----------------------------------------------------------------------
CL
      DO I=1,NPNTS
        CCW(I) = XPKP1(I)
      END DO
CL
CL----------------------------------------------------------------------
CL CALCULATE CONVECTIVE CLOUD BASE, CONVECTIVE CLOUD TOP , TOTAL
CL CONDENSED WATER/ICE AND CONVECTIVE CLOUD AMOUNT
CL
CL SUBROUTINE CON_RAD
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (9)
CL----------------------------------------------------------------------
CL
CHBG  CALL CON_RAD(K,XPK,XPKP1,FLXKP1,BTERM,CCA,ICCB,ICCT,
CHBG &            TCW,CCLWP,DELPKP1,BLAND,NPNTS)
CL
CL----------------------------------------------------------------------
CL CALCULATE CLOUD DEPTH AND ASSIGN CRITICAL CLOUD DEPTHS
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (8), EQUATION (34), (35)
CL----------------------------------------------------------------------
CL
      DO 30 I=1,NPNTS
       IF ( BLOWST(I) ) DEPTH(I) = 0.
C
       IF ( BGMKP1(I) )
     &   DEPTH(I) = DEPTH(I) + ( CP * THEKP1(I) *
     &                             (1.0+C_VIRTUAL*QEKP1(I)) *
     &                                          DELEXKP1(I)/G )
C
       IF (.NOT.BWKP1(I)) THEN
          DCRIT = CRITDICE
       ELSE IF (BLAND(I)) THEN
          DCRIT = CRITDLND
       ELSE
          DCRIT = CRITDSEA
       ENDIF
CL
CL----------------------------------------------------------------------
CL CALCULATE PRECIPITATION FROM LAYER K+1 AND ADJUST CLOUD WATER
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (8), EQUATION (36)
CL----------------------------------------------------------------------
CL
       XMIN = MIN (MPARWTR , 0.5*QSEKP1(I))

CHBG   IF (      (DEPTH(I) .GT. DCRIT)
CHBG &     .AND. (XPKP1(I) .GT. XMIN)) THEN
CHBG      PREKP1(I) = (XPKP1(I) - XMIN) * FLXKP1(I) / G
CHBG      XPKP1(I) = XMIN
CHBG   ELSE
CHBG      PREKP1 (I) = 0.
CHBG   ENDIF

       xpold(i)=xpkp1(i)
       IF (      (DEPTH(I) .GT. DCRIT)
     &     .AND. (XPKP1(I) .GT. XMIN)) THEN
          FRACOUT= (XPKP1(I) - XMIN) * FLXKP1(I)
          if(BTERM(I).and.qcloud.and.(ipass.eq.1))then
chbg.. Convection is terminating:
           TEMPRY = FRACOUT/DELPKP1(I)
           if(BWKP1(I))then
chbg.. Send the detrained liquid to LDR cloud scheme
            CCWL(I)= TEMPRY ! (Kgm/Kgm)/sec
           else
chbg.. Send the detrained ice to LDR cloud scheme
            CCWI(I)= TEMPRY ! (Kgm/Kgm)/sec
           endif
          else
chbg.. Not terminating or "qcloud" microphysics not in use
chbg.. (or over leads surface - qcloud i.e. progcloud not used)
chbg..    : Drop all condensate out as precipitation
chbg   (Pa = Kgm/m/sec**2)
           PREKP1(I) = PREKP1(I) + FRACOUT / G ! Kgm/m**2/sec
          endif
          XPKP1(I) = XMIN
       ENDIF
   30  CONTINUE

c Wet deposition scavenging fraction

       do nt=1,ntrac
         if(lwetdep(nt))then
           do i=1,npnts
             if(xpold(i).gt.0..and.BWKP1(I))then
c             if(xpold(i).gt.0.)then
               fs(i,nt)=(xpold(i)-xpkp1(i))/xpold(i)
             endif
           enddo
         endif
       enddo

       do i=1,npnts
         fs(i,itracso2)=fs(i,itracso2)*f_so2(i)
       enddo


CHBG
CHBG--------------------------------------------------------------------
CHBG CALCULATE CONVECTIVE CLOUD BASE, CONVECTIVE CLOUD TOP , TOTAL
CHBG CONDENSED WATER/ICE AND CONVECTIVE CLOUD AMOUNT
CHBG  AFTER PRECIPITATION
CHBG--------------------------------------------------------------------
CHBG
      CALL CON_RAD(K,XPK,XPKP1,FLXKP1,BTERM,CCA,ICCB,ICCT,
     &            TCW,CCLWP,DELPKP1,BLAND,NPNTS)
CHBG
CHBG--------------------------------------------------------------------
CHBG STORE CONVECTIVE CLOUD LIQUID WATER AFTER PRECIPITATION
CHBG--------------------------------------------------------------------
CHBG
      DO I=1,NPNTS
        CCW(I) = 0.5 * ( CCW(I) + XPKP1(I) )
      END DO
C
      RETURN
      END
CLL  SUBROUTINE CON_RAD------------------------------------------------
CLL
CLL  PURPOSE : CALCULATES CONVECTIVE CLOUD TOP, BASE AND
CLL            AMOUNT
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.3   23/12/93 : DG020893 : CHANGE TO CLOUD TOP BECAUSE OF
CLL                               CHANGE TO DETRAINMENT RATE
CLL                               CALCULATION
CLL
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENT NUMBER: P27
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL                  SECTION (9)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CON_RAD (K,XPK,XPKP1,FLXKP1,BTERM,CCA,ICCB,
     &                    ICCT,TCW,CCLWP,DELPKP1,BLAND,NPNTS)
C    &                    ICCT,TCW,CCLWP,DELPKP1,NPNTS)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include 'C_G.f'
      include 'UKPARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTH AND LOOP VARIABLES
C----------------------------------------------------------------------
C
      INTEGER NPNTS        ! IN VECTOR LENGTH
C
      INTEGER K            ! IN PRESENT MODEL LAYER
C
      INTEGER I            ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL XPK(ln2)        ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
C
      REAL XPKP1(ln2)      ! IN PARCEL CLOUD WATER IN LAYER K+1 (KG/KG)
C
      LOGICAL BTERM(ln2)   ! IN MASK FOR POINTS WHERE CONVECTION
                           !    IS ENDING
C
      REAL FLXKP1(ln2)     ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      REAL DELPKP1(ln2)    ! IN PRESSURE DIFFERENCE ACROSS LAYER K+1
C
      LOGICAL BLAND(ln2)     ! IN LAND/SEA MASK
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL TCW(ln2)        ! INOUT
                           ! IN  TOTAL CONDENSED WATER SUMMED TO
                           !     LAYER K (KG/M**2/S)
                           ! OUT TOTAL CONDENSED WATER SUMMED TO
                           !     LAYER K+1 OR IF CONVECTION HAS
                           !     TERMINATED ZEROED (KG/M**2/S)
C
      REAL CCLWP(ln2)      ! INOUT
                           ! IN  TOTAL CLOUD LIQUID WATER PATH
                           !     SUMMED TO LAYER K  (KG/M**2)
                           ! OUT TOTAL CLOUD LIQUID WATER PATH
                           !     SUMMED TO LAYER K+1 (KG/M**2)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE AND OUTPUT
C----------------------------------------------------------------------
C
      REAL CCA(ln2)        ! OUT CONVECTIVE CLOUD AMOUNT (%)
C
      INTEGER ICCB(ln2)     ! OUT CONVECTIVE CLOUD BASE LEVEL
C
      INTEGER ICCT(ln2)     ! OUT CONVECTIVE CLOUD TOP LEVEL
C
C*---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL CALCULATE CLOUD BASE
CL
CL WHEN CLOUD BASE SET ZERO TOTAL CONDENSED WATER
CL---------------------------------------------------------------------
CL
      DO 10 I=1,NPNTS
       IF ( XPK(I) .LE. 0.0 .AND. XPKP1(I) .GT. 0 ) THEN
         ICCB(I)=K+1
         CCLWP(I)=0.0
       END IF
CL
CL---------------------------------------------------------------------
CL CALCULATE CLOUD TOP
CL---------------------------------------------------------------------
CL
       IF( BTERM(I) .AND.
     &      ((XPKP1(I).GT.0.0) .OR.(XPK(I).GT.0.0) ) ) ICCT(I) = K+1
C
       IF ( FLXKP1(I) .GT. 0.0) THEN
CL
CL---------------------------------------------------------------------
CL SUM TOTAL CONDENSED WATER PER SECOND - ASSUMES THAT THE INITIAL
CL CONVECTIVE LAYER IS UNSATURATED
CL---------------------------------------------------------------------
CL
       TCW(I) = TCW(I) + FLXKP1(I) * XPKP1(I) / G
CL
CL---------------------------------------------------------------------
CL SUM CONV CONDENSED WATER PATH - ASSUMES THAT THE INITIAL
CL CONVECTIVE LAYER IS UNSATURATED
CL---------------------------------------------------------------------
CL
       CCLWP(I) = CCLWP(I) + XPKP1(I) * DELPKP1(I) / G
C
       END IF
CL
CL---------------------------------------------------------------------
CL CALCULATE CONVECTIVE CLOUD AMOUNT IF CONVECTION TERMINATES IN
CL LAYER K AND TOTAL CONDENSED WATER PATH OVER A TIME STEP
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (9), EQUATION (37)
CL---------------------------------------------------------------------
CL
       IF( BTERM(I) .AND. TCW(I).GT.0.0 ) THEN
C
CHBG    IF ( TCW(I) .LT. 2.002E-6 ) TCW(I) = 2.002E-6
        IF ( TCW(I) .LT. 2.793E-6 ) TCW(I) = 2.793E-6 ! min 2% cloud
        CCA(I) = 0.7873 + 0.06 * LOG(TCW(I))
chbg    IF ( TCW(I) .LT. 2.365E-6 ) TCW(I) = 2.365E-6 ! min 2% cloud
chbg    CCA(I) = 1.5746 + 0.12 * LOG(TCW(I))
        IF(.NOT.BLAND(I))CCA(I) = CCA(I) * 2.0
        IF (CCA(I) .GT. 0.8) CCA(I) = 0.8
C
        TCW(I) = 0.0
C
       END IF
  10  CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE CONVEC2------------------------------------------------
CLL
CLL  PURPOSE : COMPLETES LIFTING OF THE PARCEL FROM LAYER K TO K+1
CLL
CLL            CALL SUBROUTINE PARCEL AND ENVIRON
CLL
CLL            SUBROUTINE PARCEL CALCULATES AN INITIAL MASS FLUX,
CLL            CARRIES OUT THE DETRAINMENT CALCULATION, TESTS
CLL            TO SEE IF CONVECTION IS TERMINATING AND CALCULATES THE
CLL            PRECIPITATION RATE FROM LAYER K+1
CLL
CLL            SUBROUTINE ENVIRON CALCULATES THE EFFECT OF CONVECTION
CLL            UPON THE LARGE-SCALE ATMOSPHERE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  3.2  8/07/93 : added convective cloud condensed water diagnostic
CLL               : P Inness
CLL
CLL   3.3   23/12/93 : DG060893 : CORRECTION TO REDUCE OVER PREDICTION
CLL                               OF CONVECTIVE SNOW; ALLOW MASS FLUX
CLL                               IN LAYER K-1 WHEN CLOUD UNDERGOING
CLL                               TERMINAL DETRAINMENT TO BE PASSED
CLL                               BACK TO CONVECT CORRECTLY.
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENT NUMBER: P27
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CONVEC2 (NPNTS,NLEV,K,THEK,THEKP1,QEK,QEKP1,QSEKP1,
     &                   DQSKP1,PSTAR,THPK,QPK,THPKP1,QPKP1,XSQKP1,
     &                   RBUOY,QSEK,DQSK,THPI,QPI,XPK,FLXK,BWKP1,BGMKP1,
     &                   BGMK,BLOWST,BLAND,BTERM,DEPTH,PREKP1,DTHEK,
     &                   DQEK,DTHEKP1,DQEKP1,BINIT,CCA,ICCB,ICCT,
     &                   TCW,EKP14,EKP34,AMDETK,PK,PKP1,
     &                   EXK,EXKP1,DELEXKP1,DELPK,DELPKP1,
     &                   CCLWP,CCW,CCWI,CCWL,qcloud,ipass,xs,fs)
C
      IMPLICIT NONE
      include 'PARAMS.f'
      integer ntrac
      parameter (ntrac=3)
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER NLEV           ! IN NUMBER OF MODEL LAYERS
C
      INTEGER I              ! LOOP COUNTER
C
      INTEGER K              ! PRESENT MODEL LAYER
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL THEK(ln2)         ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K (K)
C
      REAL THEKP1(ln2)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEK(ln2)          ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL QEKP1(ln2)        ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL QSEKP1(ln2)       ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL DQSKP1(ln2)       ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT IN LAYER K+1
                             !    (KG/KG)
C
      REAL PSTAR(ln2)        ! IN SURFACE PRESSURE (PA)
C
      REAL THPKP1(ln2)       ! IN PARCEL POTENTIAL TEMPERATURE
                             !    IN LAYER K (K)
C
      REAL QPKP1(ln2)        ! IN PARCEL MIXING RATIO IN LAYER K+1
                             !    (KG/KG)
C
      REAL XSQKP1(ln2)       ! IN EXCESS WATER IN PARCEL AFTER LIFTING
                             !    LAYER K TO K+1 (KG/KG)
C
      REAL RBUOY(ln2)        ! IN PARCEL BUOYANCY IN LAYER K+1 (K)
C
      REAL QSEK(ln2)         ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL DQSK(ln2)         ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT OF LAYER K
                             !    (KG/KG/K)
C
      REAL THPI(ln2)         ! IN INITIAL PARCEL POTENTIAL TEMPERATURE
                             !    (K)
C
      REAL QPI(ln2)          ! IN INITIAL PARCEL MIXING RATIO
                             !    (KG/KG)
C
      LOGICAL BWKP1(ln2)     ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      LOGICAL BGMKP1(ln2)    ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K+1
C
      LOGICAL BLAND(ln2)     ! IN LAND/SEA MASK
C
      LOGICAL BINIT(ln2)     ! IN MASK FOR THOSE POINTS AT WHICH
                             !    CONVECTION IS OCCURING
C
      LOGICAL BLOWST(ln2)    ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
C
      REAL EKP14(ln2)        ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP34(ln2)        ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL AMDETK(ln2)       ! IN MIXING DETRAINMENT COEFFICIENT
                             !    AT LEVEL K MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL PK(ln2)           ! IN PRESSURE AT MID-POINT OF LAYER K
                             !    (PA)
C
      REAL PKP1(ln2)         ! IN PRESSURE AT MID-POINT OF LAYER K+1
                             !    (PA)
C
      REAL EXK(ln2)          ! IN EXNER RATIO AT MID-POINT OF LAYER K
C
      REAL EXKP1(ln2)        ! IN EXNER RATIO AT MID-POINT OF LAYER K+1
C
      REAL DELEXKP1(ln2)     ! IN DIFFERENCE IN EXNER RATIO BETWEEN
                             !    MID-POINTS OF LAYERS K AND K+1
C
      REAL DELPK(ln2)        ! IN DIFFERENCE IN PRESSURE ACROSS LAYER K
                             !    (PA)
C
      REAL DELPKP1(ln2)      ! IN DIFFERENCE IN PRESSURE ACROSS
                             !    LAYER K+1 (PA)
C
      logical qcloud         ! IN qcloud = T or F
      integer ipass          ! IN ipass  = 1 or 2 (2=leads)

      real xs(ln2)           ! IN sulfate mixing ratio
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT BUT WHICH ARE ALSO UPDATED IN THIS ROUTINE
C----------------------------------------------------------------------
C
      REAL THPK(ln2)         ! INOUT
                             ! IN  PARCEL POTENTIAL TEMPERATURE
                             !     IN LAYER K (K)
                             ! OUT PARCEL POTENTIAL TEMPERATURE
                             !     IN LAYER K+1 (K)
C
      REAL QPK(ln2)          ! INOUT
                             ! IN  PARCEL MIXING RATIO IN LAYER K
                             !     (KG/KG)
                             ! OUT PARCEL MIXING RATIO IN LAYER K+1
                             !     (KG/KG)
C
      REAL XPK(ln2)          ! INOUT
                             ! IN  PARCEL CLOUD WATER IN LAYER K
                             !     (KG/KG)
                             ! OUT PARCEL CLOUD WATER IN LAYER K+1
                             !     (KG/KG)
C
      REAL FLXK(ln2)         ! INOUT
                             ! IN  PARCEL MASSFLUX IN LAYER K (PA/S)
                             ! OUT PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      LOGICAL BGMK(ln2)      ! INOUT
                             ! IN  MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K
                             ! OUT MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1
C
      REAL DTHEK(ln2)        ! INOUT
                             ! IN  INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (MAY BE NONE ZERO
                             !     DUE TO A PREVIOUS SPLIT FINAL
                             !     DETRAINEMNT CALCULATION) (K/S)
                             ! OUT UPDATED INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (K/S)
C
      REAL DQEK(ln2)         ! INOUT
                             ! IN  INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K DUE TO CONVECTION
                             !     (MAY BE NONE ZERO DUE TO A
                             !     PREVIOUS SPLIT FINAL DETRAINEMNT
                             !     CALCULATION) (KG/KG/S)
                             ! OUT UPDATED INCREMENT TO MODEL MIXING
                             !     RATIO IN LAYER K DUE TO CONVECTION
                             !     (KG/KG/S)
C
      REAL TCW(ln2)          ! INOUT
                             ! IN  TOTAL CONDENSED WATER SUMMED TO
                             !     LAYER K (KG/M**2/S)
                             ! OUT UPDATED TOTAL CONDENSED WATER
                             !     SUMMED TO LAYER K+1 (KG/M**2/S)
C
      REAL DEPTH(ln2)        ! INOUT
                             ! IN  DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K (M)
                             ! OUT UPDATED DEPTH OF CONVECTIVE CLOUD
                             !     TO LAYER K+1 (M)
C
      REAL CCLWP(ln2)        ! INOUT
                             ! IN  CONDENSED WATER PATH SUMMED TO
                             !     LAYER K (KG/M**2)
                             ! OUT UPDATED CONDENSED WATER PATH
                             !     SUMMED TO LAYER K+1 (KG/M**2)
      
      real fs(ln2,ntrac)     ! INOUT Scavenging Fraction
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      LOGICAL BTERM(ln2)     ! OUT MASK FOR PARCELS WHICH TERMINATE IN
                             !     LAYER K+1
C
      REAL PREKP1(ln2)       ! OUT PRECIPITATION FROM PARCEL AS IT
                             !     RISES FROM LAYER K TO K+1 (KG/M**2/S)
C
      REAL DTHEKP1(ln2)      ! OUT INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K+1 DUE TO
                             !     CONVECTION (K/S)
C
      REAL DQEKP1(ln2)       ! OUT INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K+1 DUE TO CONVECTION
                             !     (KG/KG/S)
C
      REAL CCA(ln2)          ! OUT CONVECTIVE CLOUD AMOUNT (%)
C
      INTEGER ICCB(ln2)       ! OUT CONVECTIVE CLOUD BASE LEVEL
C
      INTEGER ICCT(ln2)       ! OUT CONVECTIVE CLOUD TOP LEVEL
C
      REAL CCW(ln2)          ! OUT CONVECTIVE CLOUD WATER(KG/KG) ON
                             ! MODEL LEVELS
C
      REAL CCWI(ln2)         ! OUT CONVECTIVE CLOUD ICE
                             ! (KG/KG) ON MODEL LEVELS
C
      REAL CCWL(ln2)         ! OUT CONVECTIVE CLOUD LIQUID
                             ! (KG/KG) ON MODEL LEVELS
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C
C ON THE IBM ARRAYS ARE ALLOCATED USING A PARAMETER STATEMENT
C
C ON THE CRAY ARRAYS ARE DYNAMICALLY ALLOCATED
C----------------------------------------------------------------------
C
      REAL THRK(ln2)         ! PARCEL DETRAINMENT POTENTIAL
                             ! TEMPERATURE IN LAYER K (K)
C
      REAL QRK(ln2)          ! PARCEL DETRAINMENT MIXING RATIO
                             ! IN LAYER K (KG/KG)
C
      REAL XPKP1(ln2)        ! PARCEL CLOUD WATER IN LAYER K+1 (KG/KG)
C
      REAL FLXKP1(ln2)       ! PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      REAL DELTAK(ln2)       ! PARCEL FORCED DETRAINMENT RATE
                             ! IN LAYER K MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
C
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL PARCEL,ENVIRON
C
C*---------------------------------------------------------------------
CL
CL----------------------------------------------------------------------
CL COMPLETE LIFTING PARCELS TO LAYER K+1
CL
CL SUBROUTINE PARCEL
CL
CL UM DOCUMENTATION PAPER P27
CL SECTIONS (5),(6),(7),(8),(9)
CL----------------------------------------------------------------------
CL
       CALL PARCEL (K,NPNTS,NLEV,PSTAR,THEKP1,THEK,QEKP1,QEK,
     &              QSEK,QSEKP1,DQSK,DQSKP1,BLAND,BWKP1,
     &              DELTAK,FLXK,THPK,QPK,THRK,QRK,
     &              BTERM,THPKP1,QPKP1,PREKP1,XPK,XPKP1,FLXKP1,
     &              XSQKP1,THPI,QPI,BGMK,BGMKP1,BLOWST,RBUOY,
     &              CCA,ICCB,ICCT,TCW,DEPTH,
     &              EKP14,EKP34,AMDETK,DELPKP1,PK,PKP1,
     &              EXK,EXKP1,DELEXKP1,CCLWP,CCW,CCWI,CCWL,
     &              qcloud,ipass,xs,fs,DELPK)
CL
CL----------------------------------------------------------------------
CL CALCULATE THE EFFECT ON THE ENVIRONMENT (EXCEPT FOR THE
CL THE EVAPORATION OF PRECIPITATION AND CHANGE OF PHASE)
CL
CL SUBROUTINE ENVIRON
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (10)
CL----------------------------------------------------------------------
CL
       CALL ENVIRON (bland,NPNTS,DTHEK,DQEK,DTHEKP1,DQEKP1,
     &               THEK,QEK,DELTAK,FLXK,
     &               THPK,QPK,THRK,QRK,THEKP1,QEKP1,
     &               BTERM,THPKP1,QPKP1,XPK,XPKP1,BWKP1,FLXKP1,
     &               BLOWST,EKP14,EXK,EXKP1,DELPK,DELPKP1,
     &               AMDETK,PK)
CHBG &               AMDETK)
CL
CL----------------------------------------------------------------------
C
      DO 10 I=1,NPNTS
C
C-----------------------------------------------------------------------
C RESET BINIT WHERE CONVECTION HAS TERMINATED
C-----------------------------------------------------------------------
C
        BINIT(I) = .NOT.BTERM(I)
   10 CONTINUE
CL
CL---------------------------------------------------------------------
CL SWAP PARCEL VALUES READY FOR THE NEXT PART OF ASCENT
CL FROM LAYER K+1 TO K+2
CL---------------------------------------------------------------------
CL
        DO 30 I=1,NPNTS
            THPK(I) = THPKP1(I)
            QPK(I) = QPKP1(I)
            XPK(I) = XPKP1(I)
            FLXK(I) = FLXKP1(I)
            BGMK(I) = BGMKP1(I)
 30     CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE CONVECU------------------------------------------------
CLL
CLL  PURPOSE : TOP LEVEL OF THE MASS FLUX CONVECTION SCHEME.
CLL            LOOPS ROUND MODEL LEVELS FORM SURFACE UPWARDS
CLL            A STABILITY TEST IS CARRIED OUT TO DETERMINE WHICH
CLL            POINTS ARE TOO STABLE FOR CONVECTION TO OCCUR
CLL            SUBROUTINE LIFT_PAR AND CONVEC2 ARE CALLED TO CALCULATE
CLL            THE PARCEL ASCENT
CLL            SUBROUTINE POUR IS CALLED TO CALCULATE THE EVAPORATION
CLL            OF FALLING PRECIPITATION
CLL            SUBROUTINE DD_CALL CALLS THE DOWNDRAUGHT CODE
CLL            SUBROUTINE COR_ENGY IS CALLED TO CONSERVE MOIST STATIC
CLL            ENERGY ONCE OTHER CALCULATIONS ARE COMPLETE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  3.2  8/07/93 : added convective cloud condensed water diagnostic
CLL               : P Inness
CLL
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  Dated 05/02/92
CLL
CLL LOGICAL COMPONENTS INCLUDED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CONVECU(qcloud,ipass,
     &                    ttg,xso4,fscav,TH,Q,PSTAR,BLAND,DTHBYDT,
     &                    DQBYDT,RAIN,SNOW,CCA,ICCB,ICCT,CCLWP,EXNER,
     &                    AK,BK,AKM12,BKM12,DELAK,DELBK,TIMESTEP,CCW,
     &                    CCWI,CCWL,FLX,DDFLX,PRECIPX)
C
      IMPLICIT NONE
C
C
C--------------------------------------------------------------------
C MODEL CONSTANTS
C--------------------------------------------------------------------
C
      include 'PARAMS.f'

CHBG  include 'PARXS.f'
CHBG  include 'C_EPSLON.f'
CHBG  include 'C_R_CP.f'
CHBG  include 'XSBMIN.f'
CHBG  include 'MPARB.f'
CHBG  include 'DELTHST.f'
CHBG  include 'C_LHEAT.f'
CHBG  include 'MASSFC.f'
      include 'UKPARAMS.f'

      integer ntrac
      parameter (ntrac=3)
C
C---------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C---------------------------------------------------------------------
C
      INTEGER NP_FIELD            ! LENGTH OF DATA (ALSO USED TO
                                  ! SPECIFY STARTING POINT OF
                                  ! DATA PASSED IN)
      parameter (NP_FIELD=ln2)
C
      INTEGER NPNTS               ! IN FULL VECTOR LENGTH
      parameter (NPNTS=ln2)
C
      INTEGER NLEV                ! IN NUMBER OF MODEL LAYERS
      parameter (NLEV=nl)
C
      INTEGER NCONV               ! NUMBER OF POINTS WHICH PASS
                                  ! INITIAL STABILITY TEST IN LAYER K
C
      INTEGER NINIT               ! NUMBER OF POINTS AT WHICH
                                  ! CONVECTION OCCURS IN LAYER K
C
      INTEGER NTERM               ! NUMBER OF CONVECTING POINTS IN
                                  ! LAYER K AT WHICH CONVECTION IS
                                  ! TERMINATING
C
      INTEGER NCNLV               ! NUMBER OF POINTS AT WHICH CONVECTION
                                  ! OCCURS AT SOME LAYER OF THE DOMAIN
C
      INTEGER I,K                 ! LOOP COUNTERS
C
C
C---------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C---------------------------------------------------------------------
C
      LOGICAL BLAND(NP_FIELD)     ! IN LAND/SEA MASK
C
      REAL PSTAR(NP_FIELD)        ! IN SURFACE PRESSURE (PA)
C
      REAL EXNER(NP_FIELD,NLEV+1) ! IN EXNER RATIO
C
      REAL AK(NLEV),              ! IN HYBRID CO-ORDINATE COEFFICIENTS
     &     BK(NLEV)               !    DEFINE PRESSURE AT MID-POINT
                                  !    OF LAYER K
C
      REAL AKM12(NLEV+1),         ! IN HYBRID CO-ORDINATE COEFFICIENTS
     &     BKM12(NLEV+1)          !    TO DEFINE PRESSURE AT
                                  !    LEVEL K-1/2
C
      REAL DELAK(NLEV),           ! IN DIFFERENCE IN HYBRID CO-ORDINATE
     &     DELBK(NLEV)            !    COEFFICIENTS ACROSS LAYER K
C
      REAL TIMESTEP               ! IN MODEL TIMESTEP (SECS)
C
      logical qcloud              ! IN qcloud = T or F
      integer ipass               ! IN ipass  = 1 or 2 (2=leads)
C
      REAL ttg(NP_FIELD,NLEV)     ! IN MODEL TEMPERATURE (K)

      real xso4(np_field,nlev)    ! IN sulfate mixing ratio
C
C
C---------------------------------------------------------------------
C  VARIABLES WHICH ARE INPUT AND OUTPUT
C---------------------------------------------------------------------
C
      real fscav(np_field,nlev,ntrac) !INOUT tracer scavenging fraction

      REAL TH(NP_FIELD,NLEV)      ! INOUT
                                  ! IN MODEL POTENTIAL TEMPERATURE (K)
                                  ! OUT MODEL POTENTIAL TEMPERATURE
                                  !     AFTER CONVECTION (K)
C
      REAL Q(NP_FIELD,NLEV)       ! INOUT
                                  ! IN MODEL MIXING RATIO (KG/KG)
                                  ! OUT MODEL MIXING RATIO AFTER
                                  !     AFTER CONVECTION (KG/KG)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL DTHBYDT(NP_FIELD,NLEV) ! OUT INCREMENTS TO POTENTIAL
                                  !     TEMPERATURE DUE TO CONVECTION
                                  !     (K/S)
C
      REAL DQBYDT(NP_FIELD,NLEV)  ! OUT INCREMENTS TO MIXING RATIO
                                  !     DUE TO CONVECTION (KG/KG/S)
C
      REAL RAIN(NP_FIELD)         ! OUT SURFACE CONVECTIVE RAINFALL
                                  !     (KG/M**2/S)
C
      REAL SNOW(NP_FIELD)         ! OUT SURFACE CONVECTIVE SNOWFALL
                                  !     (KG/M**2/S)
C
      REAL CCA(NP_FIELD)          ! OUT CONVECTIVE CLOUD AMOUNT (%)
      REAL KCCA(NP_FIELD)         ! 1st CONVECTIVE CLOUD AMOUNT (%)
C
      INTEGER ICCB(NP_FIELD)      ! OUT CONVECTIVE CLOUD BASE LEVEL
      INTEGER KEEB(NP_FIELD)      ! 1st CONVECTIVE CLOUD BASE LEVEL
C
      INTEGER ICCT(NP_FIELD)      ! OUT CONVECTIVE CLOUD TOP LEVEL
      INTEGER KEET(NP_FIELD)      ! 1st CONVECTIVE CLOUD TOP LEVEL
C
      REAL CCLWP(NP_FIELD)        ! OUT CONDENSED WATER PATH (KG/M**2)
C
      REAL CCW(NP_FIELD,NLEV)     ! OUT CONVECTIVE CLOUD LIQUID WATER
                                  ! (KG/KG) ON MODEL LEVELS
C
      REAL CCWI(NP_FIELD,NLEV)    ! OUT CONVECTIVE CLOUD ICE
                                  ! (KG/KG) ON MODEL LEVELS
C
      REAL CCWL(NP_FIELD,NLEV)    ! OUT CONVECTIVE CLOUD LIQUID
                                  ! (KG/KG) ON MODEL LEVELS
C
C----------------------------------------------------------------------
C VARIABLES DEFINED LOCALLY
C
C ON THE IBM ARRAYS ARE ALLOCATED USING A PARAMETER STATEMENT
C
C ON THE CRAY ARRAYS ARE DYNAMICALLY ALLOCATED
C----------------------------------------------------------------------
C
C
C----------------------------------------------------------------------
CHBG NOTE THAT THE BELOW WORK SPACE ALLOCATION (WORK and WORK2) HAVE
CHBG BEEN EQUIVALENCED TO ARRAYS WITH 2ND DIMENSIONS OF 24 and 19 
CHBG RESPECTIVELY. THEY REQUIRE THESE DIMENSIONS AT LEAST REGARDLESS 
CHBG OF THE VALUE OF NLEV*2 AND NLEV+3.
C----------------------------------------------------------------------
C
      REAL WORK(NPNTS,NLEV*2),   ! WORK SPACE
     &     WORK2(NPNTS,NLEV+3)
CHBG
      REAL WORKX(NPNTS,24),WORK2X(NPNTS,19)
      EQUIVALENCE (WORK,WORKX),(WORK2,WORK2X)
CHBG
C
      REAL WORKI(NPNTS,NLEV)     ! FOR SUMMING CCWI IN COR_ENGY
C
      LOGICAL BWORK(NPNTS,4),    ! WORK SPACE FOR 'BIT' MASKS
     &        BWORK2(NPNTS,4)
C
      LOGICAL BCONV(NPNTS)       ! MASK FOR POINTS WHERE STABILITY
                                  ! LOW ENOUGH FOR CONVECTION
                                  ! TO OCCUR
C
      REAL QSE(NPNTS,NLEV)       ! SATURATION MIXING RATIO OF CLOUD
                                  ! ENVIRONMENT (KG/KG)
C
      REAL TT(NPNTS)             ! TEMPORARY STORE FOR TEMPERATURE
                                  ! IN CALCULATION OF SATURATION
                                  ! MIXING RATIO (K)
C
      REAL PT(NPNTS)             ! TEMPORARY STORE FOR PRESSURE
                                  ! IN CALCULATION OF SATURATION
                                  ! MIXING RATIO (PA)
C
      REAL CCAC(NPNTS)            ! COMPRESSED VALUES OF CCA
C
      INTEGER ICCBC(NPNTS)        ! COMPRESSED VALUES OF CCB
C
      INTEGER ICCTC(NPNTS)        ! COMPRESSED VALUES OF CCT
C
      REAL TCW(NPNTS)             ! TOTAL CONDENSED WATER (KG/M**2/S)
C
      REAL TCWC(NPNTS)            ! COMPRESSED VALUES OF TCW
C
      REAL CCLWPC(NPNTS)          ! COMPRESSED VALUE OF CCLWP
C
      REAL DQSTHK(NPNTS)          ! GRADIENT OF SATURATION MIXING
                                  ! RATIO OF CLOUD ENVIRONMENT WITH
                                  ! POTENTIAL TEMPERATURE IN LAYER K
                                  ! (KG/KG/K)
C
      REAL DQSTHKP1(NPNTS)        ! GRADIENT OF SATURATION MIXING
                                  ! RATIO OF CLOUD ENVIRONMENT WITH
                                  ! POTENTIAL TEMPERATURE IN LAYER K+1
                                  ! (KG/KG/K)
C
      REAL PRECIP(NPNTS,NLEV)     ! AMOUNT OF PRECIPITATION
                                  ! FROM EACH LAYER (KG/M**2/S)
                                  ! BEFORE DOWNDRAUGHT CALCULATION 
                                  ! (ZERO AFTER THAT)
C
      REAL PRECIPX(NPNTS,NLEV)    ! FINAL AMOUNT OF PRECIPITATION
                                  ! FROM EACH LAYER (KG/M**2/S)
C
      REAL THPI(NPNTS)            ! INITIAL PARCEL POTENTIAL TEMPERATURE
                                  ! (K)
C
      REAL QPI(NPNTS)             ! INITIAL PARCEL MIXING RATIO
                                  ! (KG/KG)
C
      REAL THP(NPNTS,NLEV)        ! PARCEL POTENTIAL TEMPERATURE
                                  ! IN LAYER K (K)
C
      REAL QP(NPNTS,NLEV)         ! PARCEL MIXING RATIO IN LAYER K
                                  ! (KG/KG)
C
      REAL XPK(NPNTS)             ! PARCEL CLOUD WATER IN LAYER K
                                  ! (KG/KG)
C
      REAL FLX(NPNTS,NLEV)        ! PARCEL MASSFLUX IN LAYER K (PA/S)
C
      REAL DDFLX(NPNTS,NLEV)      ! DOWNDRAUGHT MASSFLUX IN LAYER K (PA/S)
C
      LOGICAL BINIT(NPNTS)        ! MASK FOR POINTS WHERE CONVECTION
                                  ! IS OCCURING
C
      LOGICAL BTERM(NPNTS)        ! MASK FOR POINTS WHERE CONVECTION
                                  ! TERMINATES IN LAYER K+1
C
      LOGICAL BWATER(NPNTS,2:NLEV) ! MASK FOR POINTS AT WHICH
                                      ! PRECIPITATION IS LIQUID
C
      LOGICAL BGMK(NPNTS)         ! MASK FOR POINTS WHERE PARCEL IN
                                  ! LAYER K IS SATURATED
C
      LOGICAL BCNLV(NPNTS)        ! MASK FOR THOSE POINTS AT WHICH
                                  ! CONVECTION HAS OCCURED AT SOME
                                  ! LEVEL OF THE MODEL
C
      REAL DEPTH(NPNTS)           ! DEPTH OF CONVECTIVE CLOUD (M)
C
      REAL FLXMAXK(NPNTS)         ! MAXIMUM INITIL CONVECTIVE MASSFLUX
                                  ! (PA/S)
C
      REAL FLXMAX2(NPNTS)         ! MAXIMUM INITIL CONVECTIVE MASSFLUX
                                  ! (PA/S)
C
      REAL PK(NPNTS)              ! PRESSURE AT MID-POINT OF LAYER K
                                  ! (PA)
C
      REAL PKP1(NPNTS)            ! PRESSURE AT MID-POINT OF LAYER K+1
                                  ! (PA)
C
      REAL DELPK(NPNTS)           ! PRESSURE DIFFERENCE ACROSS LAYER K
                                  ! (PA)
C
      REAL DELPKP1(NPNTS)         ! PRESSURE DIFFERENCE ACROSS LAYER K+1
                                  ! (PA)
C
      REAL DELPKP12(NPNTS)        ! PRESSURE DIFFERENCE BETWEEN
                                  ! LEVELS K AND K+1 (PA)
C
      REAL EKP14(NPNTS),          ! ENTRAINMENT COEFFICIENTS AT LEVELS
     &     EKP34(NPNTS)           ! K+1/4 AND K+3/4 MULTIPLIED BY
                                  ! APPROPRIATE LAYER THICKNESS
C
      REAL AMDETK(NPNTS)          ! MIXING DETRAINMENT COEFFICIENT AT
                                  ! LEVEL K MULTIPLIED BY APPROPRIATE
                                  ! LAYER THICKNESS
C
      REAL EXK(NPNTS)             ! EXNER RATIO AT LEVEL K
C
      REAL EXKP1(NPNTS)           ! EXNER RATIO AT LEVEL K+1
C
      REAL DELEXKP1(NPNTS)        ! DIFFERENCE IN EXNER RATIO
                                  ! ACROSS LAYER K+1
C
      REAL EMINDS(NPNTS)          !
C
      INTEGER INDEX1(NPNTS),      ! INDEX FOR COMPRESS AND
     &        INDEX2(NPNTS),      ! EXPAND
     &        INDEX3(NPNTS),
     &        INDEX4(NPNTS)
C
      REAL EXFULL(NPNTS,NLEV)     ! FULL LEVEL EXNER
C
      real xs(npnts)              ! sulfate mixing ratio
      real fs(npnts,ntrac)        ! Scavenging Fraction
C
C
      REAL FLX2                   ! TEMPORARY STORE FOR MASS FLUX
      REAL SUMPR                  ! RAINFALL CHECKING
CHBG
      LOGICAL FIRSTCALL
      DATA FIRSTCALL/.TRUE./
      SAVE FIRSTCALL
      REAL CXX,CXXL
      SAVE CXX
CHBG
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL QSATU,FLAG_WET,LIFT_PAR,CONVEC2,LAYER_CN,
     &         DQS_DTH,COR_ENGY,DD_CALL
C
      REAL
     &    PU,PL
      include 'P_EXNERC.f'

CHBG
CHBG Check which resolution, and set inital mass flux parameter
CHBG
      IF(FIRSTCALL)THEN
       IF((TRUNC.EQ.'R').AND.(MW.EQ.22))THEN
         CXX=C_R21
       ELSEIF((TRUNC.EQ.'T').AND.(MW.EQ.64))THEN
         CXX=C_T63
       ELSE
         print *,'Initial massflux parameter C=CXX in ukall.f'
         print *,'not set for the resolution=',trunc,mw
         stop
       ENDIF
       FIRSTCALL=.FALSE.
      ENDIF
C*---------------------------------------------------------------------
C
CL
CL---------------------------------------------------------------------
CL CALCULATE AN ARRAY OF SATURATION MIXING RATIOS
CL FIRST CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
CL PRESSURE OF LAYER K
CL
CL SUBROUTINE QSATU
CL UM DOCUMENTATION PAPER P282
CL---------------------------------------------------------------------
CL
      DO 20 K=1,NLEV
       DO 25 I = 1,NPNTS
        PU=PSTAR(I)*BKM12(K+1) + AKM12(K+1)
        PL=PSTAR(I)*BKM12(K) + AKM12(K)
CHBG    TT(I) = TH(I,K)* P_EXNER_C(EXNER(I,K+1),EXNER(I,K),PU,PL,KAPPA)
        TT(I) = ttg(I,K)
        EXFULL(I,K) = P_EXNER_C(EXNER(I,K+1),EXNER(I,K),PU,PL,KAPPA)
        TH(I,K) = ttg(I,K) / EXFULL(I,K)
        PT(I) = AK(K)+BK(K)*PSTAR(I)
   25  CONTINUE
C
       CALL QSATU (QSE(1,K),TT,PT,NPNTS)
C
  20  CONTINUE
CL
CL---------------------------------------------------------------------
CL CALCULATE BIT VECTOR WHERE WATER WILL CONDENSE RATHER THAN ICE
CL SUBROUTINE FLAG_WET
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (2B)
CL---------------------------------------------------------------------
CL
      CALL FLAG_WET(BWATER,TH,EXNER,PSTAR,AKM12,BKM12,
     &                    NPNTS,NLEV)
C
C----------------------------------------------------------------------
C INITIALISE PRECIPITATION, DTH/DT, DQ/DT, CCW ARRAYS
C----------------------------------------------------------------------
C
      DO 40 K=1,NLEV
       DO 40 I=1,NPNTS
        PRECIP(I,K) = 0.0
        PRECIPX(I,K) = 0.0
        CCW(I,K) = 0.0
        CCWI(I,K) = 0.0
        CCWL(I,K) = 0.0
        DTHBYDT(I,K) = 0.0
   40   DQBYDT(I,K) = 0.0
      DO 50 I=1,NPNTS
C
C----------------------------------------------------------------------
C INITIALISE BIT VECTORS FOR POINTS WHICH ARE ALREADY CONVECTING
C AND FOR POINTS AT WHICH CONVECTION OCCURS AT SOME LEVEL OF
C THE ATMOSPHERE
C----------------------------------------------------------------------
C
        BINIT(I) = .FALSE.
        BCNLV(I) = .FALSE.
C
C----------------------------------------------------------------------
C INITIALISE RADIATION DIAGNOSTICS
C----------------------------------------------------------------------
C
       CCA(I) = 0.0
       ICCB(I) = 0
       ICCT(I) = 0
       TCW(I) = 0.0
CHBG added next lines
       CCLWP(I)=0.0
       KEEB(I) = 0
       KEET(I) = 0
       KCCA(I) = 0.0
CHBG
C
C---------------------------------------------------------------------
C INITIALISE SURFACE PRECIPITATION ARRAYS
C---------------------------------------------------------------------
C
       RAIN(I) = 0.0
  50   SNOW(I) = 0.0
CL
CL=====================================================================
CL MAIN LOOP OVER LEVELS - FROM SURFACE TO TOP
CL=====================================================================
CL
      DO 60 K=1,NLEV-1
CL
CL---------------------------------------------------------------------
CL CALCULATE LEVEL PRESSURES, EXNER RATIO FOR MID POINTS, ENTRAINMENT
CL RATES, DETRAINMENTS RATES AND PRESSURE DIFFERENCE ACROS  LAYERS AS
CL A FUNCTION OF GRID-POINT
CL
CL SUBROUTINE LAYER_CN
CL---------------------------------------------------------------------
CL
      CALL LAYER_CN(K,NPNTS,NLEV,EXNER,AK,BK,AKM12,BKM12,
     &              DELAK,DELBK,PSTAR,PK,PKP1,DELPK,DELPKP1,
     &              DELPKP12,EKP14,EKP34,AMDETK,EXK,EXKP1,
     &              DELEXKP1)
CL
CL---------------------------------------------------------------------
CL CALCULATE DQS/DTH FOR LAYERS K AND K+1
CL
CL SUBROUTINE DQS_DTH
CL---------------------------------------------------------------------
CL
      IF (K.EQ.1) THEN
       CALL DQS_DTH(DQSTHK,TH(1,K),QSE(1,K),EXK,NPNTS)
      ELSE
       DO 65 I=1,NPNTS
        DQSTHK(I) = DQSTHKP1(I)
  65   CONTINUE
      END IF
C
       CALL DQS_DTH(DQSTHKP1,TH(1,K+1),QSE(1,K+1),EXKP1,NPNTS)
C
      DO 70 I=1,NPNTS
      BTERM(I) = .FALSE.
C
C---------------------------------------------------------------------
C SET OTHER GIRD-POINT DEPENDENT CONSTANTS
C---------------------------------------------------------------------
C
C---------------------------------------------------------------------
C MAXIMUM INITIAL CONVECTIVE MASSFLUX
C---------------------------------------------------------------------
C
       FLXMAXK(I) = DELPK(I)/((1.0 + EKP14(I)) * TIMESTEP)
C
C---------------------------------------------------------------------
C MAXIMUM CONVECTIVE MASSFLUX AT MID-POINT OF LAYER 2
C---------------------------------------------------------------------
C
      IF (K.EQ.1) FLXMAX2(I) = (PSTAR(I)-PKP1(I)) / TIMESTEP
C
C---------------------------------------------------------------------
C MINIMUM BUOYANCY FOR CONVECTION TO START FROM LAYER K
C---------------------------------------------------------------------
C
       EMINDS(I) = MPARB*DELPKP12(I)/PSTAR(I)
C
C----------------------------------------------------------------------
C SET BIT VECTOR FOR POINTS WHERE CONVECTION HAS OCCURRED AT SOME
C LEVEL OF THE ATMOSPHERE
C-----------------------------------------------------------------------
C
       BCNLV(I) =  BCNLV(I) .OR. BINIT(I)
CL
CL---------------------------------------------------------------------
CL SET INITIAL VALUES FOR POINTS NOT ALREADY INITIATED
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (3), EQUATION(17)
CL---------------------------------------------------------------------
CL
       IF (.NOT.BINIT(I)) THEN
         THPI(I) = TH(I,K) + THPIXS
c         if(bland(i))thpi(i)=th(i,k)+2.*thpixs !LDR test
         QPI(I) = Q(I,K) + QPIXS
         THP(I,K) = TH(I,K) + THPIXS
c         if(bland(i))thp(i,k)=th(i,k)+2.*thpixs !LDR test
         QP(I,K) = Q(I,K) + QPIXS
         XPK(I) = 0.0
         FLX(I,K) = 0.0
         BGMK(I) = .FALSE.
         DEPTH(I) = 0.0
       END IF
CL
CL----------------------------------------------------------------------
CL FORM A BIT VECTOR OF POINTS FOR WHICH CONVECTION MAY BE POSSIBLE
CL FROM LAYER K TO K-1 EITHER BECAUSE STABILITY IS LOW ENOUGH
CL OR BECAUSE CONVECTION OCCURRING FROM LAYER K+1 TO K
CL THIS BIT VECTOR IS USED IN THE FIRST COMPREE OF THE DATA
CL TO CALCULATE PARCEL BUOYANCY IN LAYER K-1
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION(3), EQUATION(16)
CL----------------------------------------------------------------------
CL
        BCONV(I) = BINIT(I) .OR.
     &           ((TH(I,K) - TH(I,K+1) + DELTHST
     &           + MAX(0.0,(Q(I,K)-QSE(I,K+1)))*(LC/(CP*EXKP1(I))))
     &           .GT. 0.)
  70  CONTINUE
CL
CL----------------------------------------------------------------------
CL COMPRESS DOWN POINTS ON THE BASIS OF BIT VECTOR BCONV
CL----------------------------------------------------------------------
CL
      NCONV = 0
      DO 75 I=1,NPNTS
        IF(BCONV(I))THEN
          NCONV = NCONV + 1
          INDEX1(NCONV) = I
        END IF
  75  CONTINUE
C
C----------------------------------------------------------------------
C  WORK SPACE USAGE FOR FIRST COMPRESS ON BASIS OF SIMPLE
C  STABILITY TEST (SECTION (3), EQN(16))
C
C  REFERENCES TO WORK AND BWORK REFER TO STARTING ADDRESS
C
C  LENGTH OF COMPRESSES DATA = NCONV
C
C  WORK(1,1)  = TH(#,K)
C  WORK(1,2)  = TH(#,K+1)
C  WORK(1,3)  = Q(#,K)
C  WORK(1,4)  = Q(#,K+1)
C  WORK(1,5)  = QSE(#,K+1)
C  WORK(1,6)  = DQSTHKP1(#)
C  WORK(1,7)  = THP(#,K)
C  WORK(1,8)  = QP(#,K)
C  WORK(1,9)  = PKP1(#)
C  WORK(1,10) = EXKP1(#)
C  WORK(1,11) = EKP14(#)
C  WORK(1,12) = EKP34(#)
C  WORK(1,13) = PARCEL POT. TEMPERATURE IN LAYER K+1
C  WORK(1,14) = PARCEL MIXING RATIO IN LAYER K+1
C  WORK(1,15) = EXCESS WATER VAPOUR IN PARCEL ABOVE
C               SATURATION AFTER DRY ASCENT
C  WORK(1,16) = PARCEL BUOYANCY IN LAYER K+1
C  WORK(1,17) = DELPKP12(#)
C  WORK(1,18) = PSTAR(#)
C  WORK(1,19) = FLX(#,K)
C  WORK(1,20) = EMINDS(#)
C  WORK(1,21) = FLXMAXK(#)
C  WORK(1,22) = FLXMAX2(#)
C
C  BWORK(1,1) = BWATER(INDEX1(I),K+1)
C  BWORK(1,2) = .TRUE. IF PARCEL SATURATED IN LAYER K+1
C  BWORK(1,3) = .TRUE. IF CONVECTION INITIATE FROM LAYER K+1
C  BWORK(1,4) = BINIT(INDEX1(I))
C----------------------------------------------------------------------
C
      IF (NCONV .NE. 0) THEN
*vdir nodep
        DO 80 I=1,NCONV
          WORK(I,1)  = TH(INDEX1(I),K)
          WORK(I,2)  = TH(INDEX1(I),K+1)
          WORK(I,3)  = Q(INDEX1(I),K)
          WORK(I,4)  = Q(INDEX1(I),K+1)
          WORK(I,5)  = QSE(INDEX1(I),K+1)
          WORK(I,6)  = DQSTHKP1(INDEX1(I))
          WORK(I,7)  = THP(INDEX1(I),K)
          WORK(I,8)  = QP(INDEX1(I),K)
          WORK(I,9)  = PKP1(INDEX1(I))
          WORK(I,10) = EXKP1(INDEX1(I))
          WORK(I,11) = EKP14(INDEX1(I))
          WORK(I,12) = EKP34(INDEX1(I))
          WORK(I,17) = DELPKP12(INDEX1(I))
          WORK(I,18) = PSTAR(INDEX1(I))
          WORK(I,19) = FLX(INDEX1(I),K)
          WORK(I,20) = EMINDS(INDEX1(I))
          WORK(I,21) = FLXMAXK(INDEX1(I))
          WORK(I,22) = FLXMAX2(INDEX1(I))
          BWORK(I,1) = BWATER(INDEX1(I),K+1)
          BWORK(I,4) = BINIT(INDEX1(I))
C
CHBG
          BWORK2(I,1)= BLAND(INDEX1(I))
C
  80    CONTINUE
CL
CL---------------------------------------------------------------------
CL LIFT PARCEL FROM LAYER K TO K-1
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (3) AND (4)
CL---------------------------------------------------------------------
CL
      CALL LIFT_PAR (NCONV,WORK(1,13),WORK(1,14),WORK(1,15),
     &               BWORK(1,2),BWORK(1,1),WORK(1,7),WORK(1,8),
     &               WORK(1,2),WORK(1,4),WORK(1,1),WORK(1,3),
     &               WORK(1,5),WORK(1,6),WORK(1,9),
     &               WORK(1,10),WORK(1,11),WORK(1,12))
C
      DO 110 I=1,NCONV
CL
CL---------------------------------------------------------------------
CL CALCULATE BUOYANCY OF PARCEL IN LAYER K-1
CL---------------------------------------------------------------------
CL
        WORK(I,16) = WORK(I,13)*(1.0 +
     &                            C_VIRTUAL * WORK(I,14))
     &               - WORK(I,2)*(1.0 +
     &                            C_VIRTUAL * WORK(I,4))
C
C----------------------------------------------------------------------
C INITIATE CONVECTION WHERE BUOYANCY IS LARGE ENOUGH
C----------------------------------------------------------------------
C
        BWORK(I,3) = .NOT.BWORK(I,4) .AND. WORK(I,16) .GT.
     &      (WORK(I,20)+ XSBMIN)
C
C----------------------------------------------------------------------
C CALCULATE INITIAL MASSFLUX FROM LAYER K
C----------------------------------------------------------------------
C
CHBG  IF (BWORK(I,3))
CHBG &     WORK(I,19) = 1.0E-3 * WORK(I,18) * (D + C * WORK(I,18) *
CHBG &                    ((WORK(I,16) - XSBMIN) / WORK(I,17)))
      CXXL=CXX
CHBG- Make initial mass flux stronger for land when from BL
      IF(BWORK2(I,1).AND.(K.LE.4))CXXL=CXX*2.0
      IF (BWORK(I,3))
     &     WORK(I,19) = 1.0E-3 * WORK(I,18) * (D + CXXL* WORK(I,18) *
     &                    ((WORK(I,16) - XSBMIN) / WORK(I,17)))
CHBG
  110 CONTINUE
C
C----------------------------------------------------------------------
C LIMIT MASSFLUX IN LOWEST CONVECTING LAYER TO BE <= MASS OF LAYER
C OR
C IF K=1 ADJUST ENTRAINMENT RATE IN BOTTOM HALF OF LAYER 2 SO
C NOT TO AFFECT THE MASS FLUX AT MID-POINT OF LAYER 2
C----------------------------------------------------------------------
C
      IF ( K .EQ. 1 ) THEN
C
       DO I=1,NCONV
C
C--------------------------------------------------------------------
C CARRY OUT CALCULATION IF CONVECTION WAS INITIATED FROM LAYER 1
C--------------------------------------------------------------------
C
        IF ( BWORK(I,3) ) THEN
C
C--------------------------------------------------------------------
C CALCULATE MASS FLUX AT MID-POINT OF LAYER 2 USING STANDARD
C ENTRAINMENT RATES
C--------------------------------------------------------------------
C
         FLX2 = WORK(I,19) * (1.0 + WORK(I,11)) * (1.0 + WORK(I,12))
C
C--------------------------------------------------------------------
C IF MASS FLUX IN LAYER 2 EXCEEDS MASS OF LAYER THEN LIMIT MASS FLUX
C OVER A TIMESTEP TO MASS OF LAYER
C--------------------------------------------------------------------
C
        IF (WORK(I,19) .GT. WORK(I,21)) THEN
C
         WORK(I,19) = WORK(I,21)
C
C--------------------------------------------------------------------
C IF MASS FLUX AT MID-POINT OF LAYER 2 EXCEEDS THE MASS OF THE COLUMN
C DOWN TO THE SURFACE OVER THE TIMESTEP THEN LIMIT MASS FLUX
C--------------------------------------------------------------------
C
        IF ( FLX2 .GT. WORK(I,22)) FLX2 = WORK(I,22)
C
C--------------------------------------------------------------------
C ADJUST ENTRAINMENT RATE IN BOTTOM HALF OF LAYER 2
C--------------------------------------------------------------------
C
       WORK(I,12) = (FLX2/(WORK(I,19) * (1.0 + WORK(I,11)))) - 1.0
       END IF
C
       END IF
      END DO
C
C---------------------------------------------------------------------
C RECALCULATE ASCENT FROM LAYER 1 TO 2 USING ADJUSTED ENTRAINMENT RATE
C---------------------------------------------------------------------
C
      CALL LIFT_PAR (NCONV,WORK(1,13),WORK(1,14),WORK(1,15),
     &               BWORK(1,2),BWORK(1,1),WORK(1,7),WORK(1,8),
     &               WORK(1,2),WORK(1,4),WORK(1,1),WORK(1,3),
     &               WORK(1,5),WORK(1,6),WORK(1,9),
     &               WORK(1,10),WORK(1,11),WORK(1,12))
C
*vdir nodep
       DO I=1,NCONV
C
        IF ( BWORK(I,3) ) THEN
CL
CL---------------------------------------------------------------------
CL RECALCULATE BUOYANCY OF PARCEL IN LAYER K-1
CL---------------------------------------------------------------------
CL
        WORK(I,16) = WORK(I,13)*(1.0 +
     &                            C_VIRTUAL * WORK(I,14))
     &               - WORK(I,2)*(1.0 +
     &                            C_VIRTUAL * WORK(I,4))
C
C----------------------------------------------------------------------
C RESET MASK TO INITIATE CONVECTION WHERE BUOYANCY IS LARGE ENOUGH
C----------------------------------------------------------------------
C
        BWORK(I,3) = .NOT.BWORK(I,4) .AND. WORK(I,16) .GT.
     &      (WORK(I,20)+ XSBMIN)
C
        BWORK(I,4) = BWORK(I,4) .OR. BWORK(I,3)
C
       END IF
C
       FLX(INDEX1(I),K) = WORK(I,19)
C
      END DO
C
C----------------------------------------------------------------------
C END OF CALCULATION FOR LAYER 1
C----------------------------------------------------------------------
C
      ELSE
C
*vdir nodep
       DO I=1,NCONV
C
C----------------------------------------------------------------------
C IF MASS FLUX OUT OF THE INITIAL LAYER IS GREATER THAN THE MASS OF
C THE LAYER OVER THE TIMESTEP THEN LIMIT MASS FLUX TO MASSS OF LAYER
C----------------------------------------------------------------------
C
        IF (BWORK(I,3) .AND. WORK(I,19).GT.WORK(I,21))
     &                 WORK(I,19) = WORK(I,21)
C
        BWORK(I,4) = BWORK(I,4) .OR. BWORK(I,3)
C
        FLX(INDEX1(I),K) = WORK(I,19)
C
       END DO
C
      END IF
C
CL
CL--------------------------------------------------------------------
CL ZERO MIXING DETRAINMENT RATE WHEN CONVECTION STARTS FROM LAYER K
CL--------------------------------------------------------------------
CL
*vdir nodep
      DO I=1,NCONV
       IF ( BWORK(I,3) ) AMDETK(INDEX1(I))=0.0
      END DO
CL
CL--------------------------------------------------------------------
CL COMPRESS DOWN THOSE POINTS WHICH ARE NOT BUOYANT IN LAYER K-1
CL--------------------------------------------------------------------
CL
      NINIT = 0
      DO 115 I=1,NCONV
        IF(BWORK(I,4))THEN
          NINIT = NINIT + 1
          INDEX2(NINIT) = I
      END IF
  115 CONTINUE
C
C
C----------------------------------------------------------------------
C  WORK SPACE USAGE FOR SECOND COMPRESS ON BASIS OF WHETHER
C  PARCEL A PARCEL STARTING FROM LAYER K IS BUOYANT IN LAYER
C  K+1 OR IF CONVECTION ALREADY EXISTS IN LAYER K
C
C  REFERENCES TO WORK, WORK2, BWORK AND BWORK2
C  REFER TO STARTING ADDRESS
C
C  LENGTH OF COMPRESSES DATA = NINIT
C
C  WORK2 AND BWORK2 ARE COMPRESSED DOWN FROM COMPRESSED
C  ARRAYS STORED IN WORK AND BWORK AFTER FIST COMPRESS
C
C  WORK2(1,1)  = TH(#,K)
C  WORK2(1,2)  = TH(#,K+1)
C  WORK2(1,3)  = Q(#,K)
C  WORK2(1,4)  = Q(#,K+1)
C  WORK2(1,5)  = QSE(#,K+1)
C  WORK2(1,6)  = DQSTHKP1(#)
C  WORK2(1,7)  = THP(#,K)
C  WORK2(1,8)  = QP(#,K)
C  WORK2(1,9)  = PKP1(#)
C  WORK2(1,10) = EXKP1(#)
C  WORK2(1,11) = EKP14(#)
C  WORK2(1,12) = EKP34(#)
C  WORK2(1,13) = PARCEL POT. TEMPERATURE IN LAYER K+1
C  WORK2(1,14) = PARCEL MIXING RATIO IN LAYER K+1
C  WORK2(1,15) = EXCESS WATER VAPOUR IN PARCEL ABOVE
C               SATURATION AFTER DRY ASCENT
C  WORK2(1,16) = PARCEL BUOYANCY IN LAYER K+1
C  WORK2(1,17) = NOT USED IN THIS SECTION
C  WORK2(1,18) = PSTAR(#)
C  WORK2(1,19) = FLX(#,K)
C
C  BWORK2(1,1) = BWATER(INDEX1(I),K+1)
C  BWORK2(1,2) = .TRUE. IF PARCEL SATURATED IN LAYER K+1
C  BWORK2(1,3) = .TRUE. IF CONVECTION INITIATE FROM LAYER K+1
C
C  WORK AND BWORK NOW CONTAIN DATA COMPRESSED DOWN
C  FROM FULL LENGTH VECTORS
C
C  WORK(1,1) = not used in this section
C  WORK(1,2) = QSE(#,K)
C  WORK(1,3) = DQSTHK(#)
C  WORK(1,4) = THPI(#)
C  WORK(1,5) = QPI(#)
C  WORK(1,6) = XPK(#)
C  WORK(1,7) = not used in this section
C  WORK(1,8) = DEPTH(#)
C  WORK(1,9) = PRECIP(#,K+1)
C  WORK(1,10) = DTHBYDT(#,K)
C  WORK(1,11) = DQBYDT(#,K)
C  WORK(1,12) = DTHBYDT(#,K+1)
C  WORK(1,13) = DQBYDT(#,K+1)
C  WORK(1,14) = AMDETK(#)
C  WORK(1,15) = NOY USED IN THIS SECTION
C  WORK(1,16) = PK(#)
C  WORK(1,17) = EXK(#)
C  WORK(1,18) = DELEXKP1(#)
C  WORK(1,19) = DELPK(#)
C  WORK(1,20) = DELPKP1(#)
C  WORK(1,21) = CCW(#,K+1)
C  WORK(1,23) = CCWI(#,K+1)
C  WORK(1,24) = CCWL(#,K+1)
C
C  BWORK(1,1) = BGMK(#)
C  BWORK(1,2) = BLAND(#)
C  BWORK(1,3) = BTERM(#)
C  BWORK(1,2) = BLAND(#)
C----------------------------------------------------------------------
C
      IF (NINIT .NE. 0) THEN
C
C-----------------------------------------------------------------------
C FIRST COMPRESS DOWN QUANTITIES FROM PREVIOUSLY COMPRESSED ARRAY
C-----------------------------------------------------------------------
C
*vdir nodep
        DO 120 I=1,NINIT
          WORK2(I,1)  = WORK(INDEX2(I),1)
          WORK2(I,2)  = WORK(INDEX2(I),2)
          WORK2(I,3)  = WORK(INDEX2(I),3)
          WORK2(I,4)  = WORK(INDEX2(I),4)
          WORK2(I,5)  = WORK(INDEX2(I),5)
          WORK2(I,6)  = WORK(INDEX2(I),6)
          WORK2(I,7)  = WORK(INDEX2(I),7)
          WORK2(I,8)  = WORK(INDEX2(I),8)
          WORK2(I,9)  = WORK(INDEX2(I),9)
          WORK2(I,10) = WORK(INDEX2(I),10)
          WORK2(I,11) = WORK(INDEX2(I),11)
          WORK2(I,12) = WORK(INDEX2(I),12)
          WORK2(I,13) = WORK(INDEX2(I),13)
          WORK2(I,14) = WORK(INDEX2(I),14)
          WORK2(I,15) = WORK(INDEX2(I),15)
          WORK2(I,16) = WORK(INDEX2(I),16)
          WORK2(I,17) = WORK(INDEX2(I),17)
          WORK2(I,18) = WORK(INDEX2(I),18)
          WORK2(I,19) = WORK(INDEX2(I),19)
          BWORK2(I,1) = BWORK(INDEX2(I),1)
          BWORK2(I,2) = BWORK(INDEX2(I),2)
          BWORK2(I,3) = BWORK(INDEX2(I),3)
  120   CONTINUE
C
C----------------------------------------------------------------------
C COMPRESS DOWN REST OF DATA FROM FULL ARRAYS
C
C FIRST EXPAND BACK BWORK(1,2) (=BINIT) BACK TO FULL VECTORS
C----------------------------------------------------------------------
C
CDIR$ IVDEP
*vdir nodep
      DO 130 I=1,NCONV
        BINIT(INDEX1(I)) = BWORK(I,4)
  130 CONTINUE
C
      NINIT = 0
      DO 135 I=1,NPNTS
        IF(BINIT(I))THEN
          NINIT = NINIT + 1
          INDEX3(NINIT) = I
        END IF
  135 CONTINUE
C
*vdir nodep
      DO 140 I=1,NINIT
        xs(i)     = xso4(index3(i),k+1)
        fs(i,:)     = fscav(index3(i),k+1,:)
        WORK(I,2) = QSE(INDEX3(I),K)
        WORK(I,3) = DQSTHK(INDEX3(I))
        WORK(I,4) = THPI(INDEX3(I))
        WORK(I,5) = QPI(INDEX3(I))
        WORK(I,6) = XPK(INDEX3(I))
        WORK(I,8) = DEPTH(INDEX3(I))
        CCAC(I)    = CCA(INDEX3(I))
        ICCBC(I)   = ICCB(INDEX3(I))
        ICCTC(I)   = ICCT(INDEX3(I))
        TCWC(I)    = TCW(INDEX3(I))
        CCLWPC(I)  = CCLWP(INDEX3(I))
        BWORK(I,1) = BGMK(INDEX3(I))
        BWORK(I,2) = BLAND(INDEX3(I))
        WORK(I,10) = DTHBYDT(INDEX3(I),K)
        WORK(I,11) = DQBYDT(INDEX3(I),K)
        WORK(I,12) = DTHBYDT(INDEX3(I),K+1)
        WORK(I,13) = DQBYDT(INDEX3(I),K+1)
        WORK(I,14) = AMDETK(INDEX3(I))
        WORK(I,16) = PK(INDEX3(I))
        WORK(I,17) = EXK(INDEX3(I))
        WORK(I,18) = DELEXKP1(INDEX3(I))
        WORK(I,19) = DELPK(INDEX3(I))
        WORK(I,20) = DELPKP1(INDEX3(I))
C
        BWORK(I,4) = .TRUE.
  140 CONTINUE
CL
CL----------------------------------------------------------------------
CL CALCULATE REST OF PARCEL ASCENT AND EFFECT OF CONVECTION
CL UPON THE LARGE-SCALE ATMOSPHERE
CL
CL SUBROUTINE CONVEC2
CL
CL UM DOCUMENTATION PAPER P27
CL SECTIONS (5),(6),(7),(8),(9),(10)
CL----------------------------------------------------------------------
CL
         CALL CONVEC2 (NINIT,NLEV,K,WORK2(1,1),WORK2(1,2),WORK2(1,3),
     &                WORK2(1,4),WORK2(1,5),WORK2(1,6),WORK2(1,18),
     &                WORK2(1,7),WORK2(1,8),WORK2(1,13),WORK2(1,14),
     &                WORK2(1,15),WORK2(1,16),WORK(1,2),WORK(1,3),
     &                WORK(1,4),WORK(1,5),WORK(1,6),WORK2(1,19),
     &                BWORK2(1,1),BWORK2(1,2),BWORK(1,1),BWORK2(1,3),
     &                BWORK(1,2),BWORK(1,3),WORK(1,8),WORK(1,9),
     &                WORK(1,10),WORK(1,11),WORK(1,12),WORK(1,13),
     &                BWORK(1,4),CCAC,ICCBC,ICCTC,TCWC,
     &                WORK2(1,11),WORK2(1,12),WORK(1,14),
     &                WORK(1,16),WORK2(1,9),WORK(1,17),WORK2(1,10),
     &                WORK(1,18),WORK(1,19),WORK(1,20),
     &                CCLWPC,WORK(1,21),WORK(1,23),WORK(1,24),
     &                qcloud,ipass,xs,fs)
CL
CL---------------------------------------------------------------------
CL EXPAND REQUIRED VECTORS BACK TO FULL FIELDS
CL----------------------------------------------------------------------
CL
      DO 145 I=1,NPNTS
        THP(I,K+1) = 0.0
        QP(I,K+1) = 0.0
        XPK(I) = 0.0
        FLX(I,K+1)= 0.0
        DEPTH(I) = 0.0
        PRECIP(I,K+1) = 0.0
        BGMK(I) = .FALSE.
        BTERM(I) = .FALSE.
        BINIT(I) = .FALSE.
  145 CONTINUE
C
CDIR$ IVDEP
*vdir nodep
      DO 150 I=1,NINIT
        xso4(index3(i),k+1) = xs(i)
        fscav(index3(i),k+1,:) = fs(i,:)
        THP(INDEX3(I),K+1) = WORK2(I,7)
        QP(INDEX3(I),K+1) = WORK2(I,8)
        XPK(INDEX3(I)) = WORK(I,6)
        FLX(INDEX3(I),K+1) = WORK2(I,19)
        DEPTH(INDEX3(I)) = WORK(I,8)
        PRECIP(INDEX3(I),K+1) = WORK(I,9)
        DTHBYDT(INDEX3(I),K) = WORK(I,10)
        DQBYDT(INDEX3(I),K) = WORK(I,11)
        DTHBYDT(INDEX3(I),K+1) = WORK(I,12)
        DQBYDT(INDEX3(I),K+1) = WORK(I,13)
        CCA(INDEX3(I)) = CCAC(I)
        ICCB(INDEX3(I)) = ICCBC(I)
        ICCT(INDEX3(I)) = ICCTC(I)
        TCW(INDEX3(I)) = TCWC(I)
        CCLWP(INDEX3(I)) = CCLWPC(I)
        CCW(INDEX3(I),K+1) = WORK(I,21)
        CCWI(INDEX3(I),K+1) = WORK(I,23)
        CCWL(INDEX3(I),K+1) = WORK(I,24)
C
        BGMK(INDEX3(I)) = BWORK(I,1)
        BTERM(INDEX3(I)) = BWORK(I,3)
        BINIT(INDEX3(I)) = BWORK(I,4)
  150 CONTINUE
C
       END IF
C
      END IF
CL
CL---------------------------------------------------------------------
CL DOWNDRAUGHT CALCULATION
CL
CL CARRIED OUT FOR THOSE CLOUD WHICH ARE TERMINATING
CL
CL SUBROUTINE DD_CALL
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (11)
CL---------------------------------------------------------------------
CL
C
      NTERM = 0
      DO 160 I=1,NPNTS
        IF (BTERM(I)) THEN
         NTERM = NTERM + 1
        END IF
  160 CONTINUE
C
      IF (NTERM .NE. 0) THEN
C
         CALL DD_CALL (NPNTS,K,THP(1,1),QP(1,1),TH(1,1),Q(1,1),
     &                 DTHBYDT(1,1),DQBYDT(1,1),FLX(1,1),PSTAR,
     &                 AK,BK,AKM12,BKM12,DELAK,DELBK,EXNER(1,1),
     &                 PRECIP(1,1),RAIN,SNOW,ICCB,ICCT,BWATER(1,2),
     &                 BTERM,BGMK,TIMESTEP,CCA,
     &                 DDFLX,PRECIPX(1,1))
C
C---------------------------------------------------------------------
C ADJUSTMENT TO CLOUD BASE, TOP AND AMOUNT
C
C IF CLOUD BASE AND TOP ARE EQUAL THEN ERRORS OCCUR IN RADIATION SCHEME
C
C ONLY OCCURS IF CONVECTION SATURATES UPON FORCED DETRAINMENT
C
C WHEN OCCURS ZERO CLOUD BASE, TOP AND AMOUNT
C
C---------------------------------------------------------------------
C
      DO I=1,NPNTS
        IF (BTERM(I) .AND. ICCB(I) .EQ. ICCT(I)) THEN
          ICCB(I) = 0.0
          ICCT(I) = 0.0
          CCA(I) = 0.0
          TCW(I) = 0.0
          CCLWP(I) = 0.0
        END IF
      END DO
C
C---------------------------------------------------------------------
C RESET BTERM TO FALSE
C---------------------------------------------------------------------
C
      DO 200 I=1,NPNTS
  200  BTERM(I) = .FALSE.
C
      END IF
CHBG
CHBG THERE MAY BE TWO LOTS OF CONVECTION IN A COLUMN. TRAP THESE
CHBG BY USING KEEB,KEET,KCCA TO STORE THE FIRST OCCURRENCE
CHBG
      DO I=1,NPNTS
        IF((KEEB(I).EQ.0).AND.(ICCB(I).LT.ICCT(I)))THEN
          KEEB(I)=ICCB(I)
          KEET(I)=ICCT(I)
          KCCA(I)=CCA(I)
        ENDIF
      ENDDO
CL
CL=====================================================================
CL END OF MAIN LOOP
CL=====================================================================
CL

  60  CONTINUE
CHBG
CHBG CHECK FOR TWO SETS OF CONVECTION PER GRID POINT :-
CHBG IF TWO CONVECTIVE ELEMENTS:    
CHBG   JOIN THEN IF THE TWO CONVECTIVE ELEMENTS ARE NO MORE 
CHBG   THAN 1 MODEL LEVEL APART; ALSO OVERLAP THE CONVECTIVE CLOUD.
CHBG   IF THEY ARE MORE THAN 1 MODEL LEVEL APART
CHBG   USE THE LOWER CONVECTIVE ELEMENT & CLOUD.
CHBG IF ONLY ONE CONVECTIVE ELEMENT, THEN ICCB,ICCT,CCA
CHBG WILL CONTAIN THE RELEVANT INFORMATION FOR 1 ELEMENT (& KEEB=ICCB).
CHBG
      DO I=1,NPNTS
CHBG TWO ELEMENTS SHOWN BY KEEB < ICCB
        IF(KEEB(I).LT.ICCB(I))THEN
          IF((ICCB(I)-KEET(I)).LE.1)THEN
            ICCB(I)=KEEB(I)
            CCA(I)=CCA(I)+KCCA(I)-CCA(I)*KCCA(I)
          ELSE
            ICCB(I)=KEEB(I)
            ICCT(I)=KEET(I)
            CCA(I)=KCCA(I)
          ENDIF
        ENDIF
      ENDDO
CHBG
CHBG MAKE SURE ALL CONVECTIVE LEVELS HAVE A MINIMUM OF CCW
CHBG
c     DO I=1,NPNTS
c     IF(ICCB(I).NE.0)THEN
c       DO K=ICCB(I),ICCT(I)
c         IF(CCW(I,K).LT.1.0E-06)CCW(I,K)=1.0E-06
c       ENDDO
c     ENDIF
c     ENDDO
      DO K=1,NLEV
      DO I=1,NPNTS
        if((k.ge.ICCB(I)).and.(k.le.ICCT(I)).and.
     &   (CCW(I,K).LT.1.0E-06))CCW(I,K)=1.0E-06
      ENDDO
      ENDDO
CHBG
CHBG RETAIN THE FLUX OF CONVECTIVE RAIN PER MODEL LEVEL
CHBG FOR THE QCLOUD SCHEME (PRECIPX). THERE IS NO PRECIP
CHBG FROM LEVEL K=1 SO PUT FLUX THROUGH LEVEL K=1 EQUAL TO
CHBG FLUX THROUGH LEVEL K=2
CHBG
      DO I=1,NPNTS
       SUMPR=RAIN(I)+SNOW(I)
       IF(SUMPR.GT.0.0)THEN
        PRECIPX(I,1)=PRECIPX(I,2)
c       write(6,1617)I,SUMPR,ICCB(I),ICCT(I)
c1617 format(i3,2x,e13.6,2(2x,i2))
c       write(6,1618)(PRECIPX(I,K),K=1,nl)
c1618 format(6e12.5)
       ENDIF
      ENDDO
CL
CL---------------------------------------------------------------------
CL BALANCE ENERGY BUDGET BY APPLYING CORRECTION TO THE TEMPERATURES
CL
CL SUBROUTINE COR_ENGY
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (12)
CL---------------------------------------------------------------------
CL
      NCNLV = 0
      DO 210 I=1,NPNTS
        IF(BCNLV(I))THEN
          NCNLV = NCNLV + 1
          INDEX4(NCNLV) = I
        END IF
  210 CONTINUE
C
C
C----------------------------------------------------------------------
C WORK SPACE USAGE FOR ENERGY CORRECTION CALCULATION
C
C  REFERENCES TO WORK AND WORK2
C  REFER TO STARTING ADDRESS
C
C  LENGTH OF COMPRESSES DATA = NCNLV
C
C  WORK(1,1 TO NLEV)        = DTHBYDT(#,1 TO NLEV)
C  WORK(1,NLEV+1 TO 2*NLEV) = DQBYDT(#,1 TO NLEV)
C  WORK2(1,1 TO NLEV+1)     = EXNER(#,1 TO NLEV+1)
C  WORK2(1,NLEV+2)          = TH(#,1)
C  WORK2(1,NLEV+3)          = PSTAR(#)
C----------------------------------------------------------------------
C
      IF (NCNLV .NE. 0)THEN
        DO 220 K=1,NLEV
*vdir nodep
          DO 220 I=1,NCNLV
            WORK(I,K)             = DTHBYDT(INDEX4(I),K)
            WORK(I,NLEV+K)        = DQBYDT(INDEX4(I),K)
            WORK2(I,K)            = EXNER(INDEX4(I),K)
            WORKI(I,K)            = CCWI(INDEX4(I),K)
  220   CONTINUE
C
*vdir nodep
      DO 225 I=1,NCNLV
        WORK2(I,NLEV+1) = EXNER(INDEX4(I),NLEV+1)
        WORK2(I,NLEV+2) = SNOW(INDEX4(I))
        WORK2(I,NLEV+3) = PSTAR(INDEX4(I))
  225 CONTINUE
C
C
        CALL COR_ENGY (NCNLV,NLEV,WORK(1,1),WORK(1,NLEV+1),
     &                 WORK2(1,NLEV+2),WORK2(1,1),
     &                 WORK2(1,NLEV+3),DELAK,DELBK,AKM12,BKM12,WORKI)
C
        DO 230 K=1,NLEV
CDIR$ IVDEP
*vdir nodep
          DO 230 I=1,NCNLV
            DTHBYDT(INDEX4(I),K) = WORK(I,K)
  230   CONTINUE
CL
CL---------------------------------------------------------------------
CL  UPDATE MODEL POTENTIAL TEMPERATURE AND MIXING RATIO
CL  WITH INCREMENTS DUE TO CONVECTION
CL---------------------------------------------------------------------
CL
        DO 250 K=1,NLEV
         DO 250 I=1,NPNTS
           TH(I,K) = TH(I,K) + DTHBYDT(I,K) * TIMESTEP
           ttg(I,K) = TH(I,K) * EXFULL(I,K)
           Q(I,K) = Q(I,K) + DQBYDT(I,K) * TIMESTEP
           CCWI(I,K) = CCWI(I,K) * TIMESTEP
           CCWL(I,K) = CCWL(I,K) * TIMESTEP
  250   CONTINUE
C
      END IF
C
      RETURN
      END
C
CLL  SUBROUTINE COR_ENGY-----------------------------------------------
CLL
CLL  PURPOSE : TO ADJUST THE POTENTIAL TEMPERATURE INCREMENTS
CLL            TO ENSURE THE CONSERVATION OF MOIST STATIC ENERGY
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL                  SECTION (12)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE COR_ENGY(NCORE,NLEV,DTHBYDT,DQBYDT,SNOW,
     &                   EXNER,PSTAR,DELAK,DELBK,AKH,BKH,CCWI)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include 'C_R_CP.f'
CHBG  include 'C_G.f'
CHBG  include 'C_LHEAT.f'
      include 'UKPARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTH AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NCORE               ! IN VECTOR LENGTHS
C
      INTEGER NLEV                ! IN NUMBER OF MODEL LAYERS
C
      INTEGER I,K                 ! LOOP COUNTERS
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL DQBYDT(ln2,NLEV)       ! IN INCREMENT TO MODEL MIXING
                                  !    RATIO DUE TO CONVECTION
                                  !    (KG/KG/S)
C
      REAL SNOW(ln2)              ! IN SNOW AT SURFACE (KG/M**2/S)
C
      REAL EXNER(ln2,NLEV+1)      ! IN EXNER RATIO
C
      REAL PSTAR(ln2)             ! IN SURFACE PRESSURE (PA)
C
      REAL DELAK(NLEV),           ! IN DIFFERENCE IN HYBRID CO-ORDINATE
     &     DELBK(NLEV)            !    COEFFICIENTS A AND B
                                  !    ACROSS LAYER K
C
      REAL AKH(NLEV+1)              ! IN Hybrid coordinate A at
                                    !    layer boundary
      REAL BKH(NLEV+1)              ! IN Hybrid coordinate B at
                                    !    layer boundary
C
      REAL CCWI(ln2,NLEV)         ! IN Increment to model cloud ice mixing
                                  !    ratio due to convection
                                  !    (kg/kg/s)
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL DTHBYDT(ln2,NLEV)      ! INOUT
                                  ! IN  INCREMENT TO MODEL POTENTIAL
                                  !     TEMPERATURE DUE TO CONVECTION
                                  !     (K/S)
                                  ! OUT CORRECTED INCREMENT TO MODEL
                                  !     POTENTIAL TEMPERATURE DUE TO
                                  !     CONVECTION (K/S)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE LOCALLY DEFINED
C
C ON THE IBM ARRAYS ARE ALLOCATED USING A PARAMETER STATEMENT
C
C ON THE CRAY ARRAYS ARE DYNAMICALLY ALLOCATED
C----------------------------------------------------------------------
C
      REAL QSUM(ln2)              ! SUMMATION OF INCREMENTS TO MODEL
                                  ! MIXING RATIO DUE TO CONVECTION
                                  ! IN THE VERTICAL, WEIGHTED
                                  ! ACCORDING TO THE MASS OF THE
                                  ! LAYER (KG/M**2/S)
C
      REAL QISUM(ln2)             ! Summation of increments to model
                                  ! cloud ice mixing ratio due to convection
                                  ! in the vertical, weighted
                                  ! according to the mass of the
                                  ! layer (kg/m**2/s)
C
      REAL TSPOS(ln2)             ! SUMMATION OF POSITIVE INCREMENTS
                                  ! TO MODEL POTENTIAL TEMPERATURE
                                  ! DUE TO CONVECTION WITH HEIGHT,
                                  ! WEIGHTED ACCORDING TO THE MASS
                                  ! OF THE LAYER (K/M**2/S)
C
      REAL TSNEG(ln2)             ! SUMMATION OF NEGATIVE INCREMENTS
                                  ! TO MODEL POTENTIAL TEMPERATURE
                                  ! DUE TO CONVECTION WITH HEIGHT,
                                  ! WEIGHTED ACCORDING TO THE MASS
                                  ! OF THE LAYER (K/M**2/S)
C
      REAL TERR(ln2)              ! SUMMATION OF ALL INCREMENTS TO
                                  ! MODEL POTENTIAL TEMPERATURE
                                  ! DUE TO CONVECTION WITH HEIGHT,
                                  ! WEIGHTED ACCORDING TO THE MASS
                                  ! OF THE LAYER (K/M**2/S)
C
      LOGICAL BPOSER(ln2)         ! MASK FOR POINTS IN LAYER K AT WHICH
                                  ! INCREMENTS TO MODEL POTENTIAL
                                  ! TEMPERATURE DUE TO CONVECTION ARE
                                  ! POSITIVE
C
      LOGICAL BCORR(ln2)          ! MASK FOR POINTS AT WHICH ENTHALPY
                                  ! CORRECTION IS NECESSARY
C
      REAL DELPK                  ! DIFFERENCE IN PRESSURE ACROSS A
                                  ! LAYER (PA)
C
      REAL EXTEMPK                ! EXNER RATIO AT THE MID-POINT OF
                                  ! LAYER K
C

      REAL
     &    PU,PL
      include 'P_EXNERC.f'

C*----------------------------------------------------------------------
CL
CL----------------------------------------------------------------------
CL  SUM UP MIXING RATIO AND +VE AND -VE TEMPERATURE INCREMENTS
CL----------------------------------------------------------------------
CL
      DO 20 I=1,NCORE
       QSUM (I) = 0.0
       QISUM(I) = 0.0
       TSPOS(I) = 0.0
       TSNEG(I) = 0.0
   20  CONTINUE
C
      DO 40 K=1,NLEV
       DO 30 I=1,NCORE
C
        DELPK = -DELAK(K) - DELBK(K)*PSTAR(I)
C
        PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
        PL=PSTAR(I)*BKH(K) + AKH(K)
        EXTEMPK  = P_EXNER_C(EXNER(I,K+1),EXNER(I,K),PU,PL,KAPPA)
C
        QSUM(I) = QSUM(I) + DQBYDT(I,K)*DELPK
        QISUM(I) = QISUM(I) + CCWI(I,K)*DELPK
C
        IF (DTHBYDT(I,K) .GT. 0.0) THEN
           TSPOS(I) = TSPOS(I) + DTHBYDT(I,K)*(CP*DELPK*EXTEMPK)
        ELSE
           TSNEG(I) = TSNEG(I) + DTHBYDT(I,K)*(CP*DELPK*EXTEMPK)
        ENDIF
   30  CONTINUE
   40 CONTINUE
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE ERROR AND APPLY THE NECESSARY CORRECTION
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (12), EQUATION (48), (49)
CL----------------------------------------------------------------------
CL
      DO 50 I=1,NCORE
C
       TERR(I) = LC*QSUM(I) - LF*G*SNOW(I) + TSPOS(I) + TSNEG(I)
     &           -LF*QISUM(I)
C
       BPOSER(I) = TERR(I) .GT. 0.0
C
       IF (BPOSER(I) .AND. (TSPOS(I) .EQ. 0.0)) THEN
          BPOSER(I) = .FALSE.
       ELSE IF (.NOT.BPOSER(I) .AND. (TSNEG(I) .EQ. 0.0)) THEN
          BPOSER(I) = .TRUE.
       ENDIF
C
       BCORR(I) = (TSPOS(I) .NE. 0.0) .OR. (TSNEG(I) .NE. 0.0)
C
       IF (BPOSER(I) .AND. BCORR(I)) THEN
          TERR(I) = 1. - TERR(I)/TSPOS(I)
       ELSE IF (.NOT.BPOSER(I) .AND. BCORR(I)) THEN
          TERR(I) = 1. - TERR(I)/TSNEG(I)
       ENDIF
C
  50  CONTINUE
C
      DO 100 K=1,NLEV
CDIR$ IVDEP
*vdir nodep
       DO 100 I=1,NCORE
        IF (BCORR(I) .AND. (( BPOSER(I) .AND.
     &   (DTHBYDT(I,K) .GT. 0.0)) .OR. ( .NOT.BPOSER(I)
     &   .AND. (DTHBYDT(I,K) .LT. 0.0))))
     &       DTHBYDT(I,K) = DTHBYDT(I,K)*TERR(I)
  100      CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE CRS_FRZL------------------------------------------
CLL
CLL  PURPOSE : CHANGE OF PHASE ROUTINE WHERE PRECIPITATION
CLL            CROSSES A MELTING OR FREEZING LEVEL
CLL
CLL            ADD LATENT HEATING
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE CRS_FRZL (NPNTS,RAIN,SNOW,THDD_KM1,EXKM1,
     &                     FLX_DD_KM1,BDDWT_KM1)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include'C_LHEAT.f'
CHBG  include'C_R_CP.f'
CHBG  include'C_G.f'
CHBG  include'C_0_DG_C.f'
      include 'UKPARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS                ! IN VECTOR LENGTH
C
      INTEGER I                    ! LOOP COUNTER
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      LOGICAL BDDWT_KM1(ln2)       ! IN MASK FOR POINTS WHERE
                                   !    PRECIPITATION IS LIQUID
                                   !    IN LAYER K+1
C
      REAL EXKM1(ln2)              ! IN EXNER RATIO FOR LAYER K-1
C
      REAL FLX_DD_KM1(ln2)         ! IN MASS FLUX OF LAYER K-1 (PA/S)
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL RAIN(ln2)               ! INOUT
                                   ! IN  AMOUNT OF FALLING RAIN
                                   !     DESCENDING FROM LAYER
                                   !     K-1 TO K-2 (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     RAIN (KG/M**2/S)
C
      REAL SNOW(ln2)               ! INOUT
                                   ! IN  AMOUNT OF FALLING SNOW
                                   !     DESCENDING FROM LAYER
                                   !     K-1 TO K-2 (KG/M**2/S)
                                   ! OUT UPDATED AMOUNT OF FALLING
                                   !     SNOW (KG/M**2/S)
C
      REAL THDD_KM1(ln2)           ! INOUT
                                   ! IN  DOWNDRAUGHT POTENTIAL
                                   !     TEMPERATURE IN LAYER K-1 (K)
                                   ! OUT UPDATED DOWNDRAUGHT POTENTIAL
                                   !     TEMPERATURE IN LAYER K-1
                                   !     DUE TO CHANGE OF PHASE (K)
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C---------------------------------------------------------------------
C
      REAL FACTOR                  ! USED IN THE CALCULATION OF
                                   ! CHANGE OF PHASE OF FALLING
                                   ! PRECIPITATION
C
      REAL PRECIP_FRE              ! FREEZING PRECIPITATION
C
      REAL PRECIP_MELT             ! MELTING PRECIPITATION
C
CL
CL----------------------------------------------------------------------
CL  ADD LATENT HEATING WHERE PRECIP CROSSES A MELTING OR FREEZING LEVEL
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (11), EQUATION (42)
CL----------------------------------------------------------------------
CL
       DO I=1,NPNTS
        FACTOR = LF*G/(EXKM1(I)*CP*FLX_DD_KM1(I))
C
        IF (.NOT.BDDWT_KM1(I).AND.RAIN(I).GT.0.0.AND.THDD_KM1(I).LT.
     &       TM/EXKM1(I)) THEN
C FREEZE
          PRECIP_FRE = (TM/EXKM1(I)-THDD_KM1(I))/ FACTOR
          IF (PRECIP_FRE.GT.RAIN(I)) PRECIP_FRE = RAIN(I)
          THDD_KM1(I) = THDD_KM1(I)+PRECIP_FRE*FACTOR
          RAIN(I) = RAIN(I)-PRECIP_FRE
          SNOW(I) = SNOW(I)+PRECIP_FRE
        ENDIF
C
        IF (BDDWT_KM1(I).AND.SNOW(I).GT.0.0) THEN
C MELT
          PRECIP_MELT = (THDD_KM1(I)-TM/EXKM1(I))/FACTOR
          IF (PRECIP_MELT.GT.SNOW(I)) PRECIP_MELT = SNOW(I)
          THDD_KM1(I) = THDD_KM1(I)-PRECIP_MELT*FACTOR
          RAIN(I) = RAIN(I)+PRECIP_MELT
          SNOW(I) = SNOW(I)-PRECIP_MELT
        END IF
      END DO
C
      RETURN
      END
C
CLL  SUBROUTINE DD_CALL------------------------------------------------
CLL
CLL  PURPOSE : CALCULATE INITIAL DOWNDRAUGHT MASSFLUX
CLL
CLL            RESET EN/DETRAINMENT RATES FOR DOWNDRAUGHT
CLL
CLL            COMPRESS/EXPAND VARIABLES
CLL
CLL            INITIALISE DOWNDRAUGHT
CLL
CLL            CALL DOWNDRAUGHT ROUTINE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.3   23/12/93 : DG060893 : CORRECTION TO REDUCE OVER PREDICTION
CLL                               OF CONVECTIVE SNOW; TO PASS ADDITIONAL
CLL                               DATA DOWN TO DOWN2A AND PREVENT DD
CLL                               FORMING BELOW UPDRAUGHT BASE
CLL
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DD_CALL (NPNTS,KCT,THP,QP,THE,QE,DTHBYDT,
     &                    DQBYDT,FLX,PSTAR,AK,BK,AKM12,BKM12,DELAK,
     &                    DELBK,EXNER,PRECIP,RAIN,SNOW,ICCB,ICCT,
     &                    BWATER,BTERM,BGMK,TIMESTEP,CCA,
     &                    DDFLX,PRECIPX)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                  ! LOOP COUNTER
C
      INTEGER K                  ! PRESENT MODEL LAYER
C
      INTEGER NPNTS              ! IN NUMBER OF POINTS
C
      INTEGER NDD                ! COMPRESSED VECTOR LENGTH FOR
                                 ! DOWNDRAUGHT CALCULATION
C
      INTEGER NDDON_TMP          ! NUMBER OF POINTS WITH ACTIVE
                                 ! DOWNDRAUGHT
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      INTEGER KCT                ! IN CONVECTIVE CLOUD TOP LAYER
C
      REAL AK(KCT+1)             ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BK(KCT+1)             ! IN ) MID-LAYER OF LAYER K
C
      REAL AKM12(KCT+2)          ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BKM12(KCT+2)          ! IN ) LOWER LAYER BOUNDARY OF LAYER K
C
      REAL DELAK(KCT+1)          ! IN ) HYBRID CO-ORDINATE VALUES FOR
      REAL DELBK(KCT+1)          ! IN ) THICKNESS OF LAYER K
C
      REAL EXNER(ln2,KCT+2)      ! IN EXNER FUNCTION AT LAYER BOUNDARIES
                                 !    STARTING AT LEVEL K-1/2
C
      REAL THP(ln2,KCT+1)        ! IN POTENTIAL TEMPERATURE OF
                                 !    PARCEL (K)
C
      REAL QP(ln2,KCT+1)         ! IN MODEL MIXING RATIO (KG/KG)
C
      REAL THE(ln2,KCT+1)        ! IN MODEL ENVIRONMENTAL POTENTIAL
                                 !    TEMPERATURE (K)
C
      REAL QE(ln2,KCT+1)         ! IN ENVIRONMENT MIXING RATIO
                                 !    (KG/KG)
C
      REAL FLX(ln2,KCT+1)        ! IN CONVECTIVE MASSFLUX (PA/S)
C
      REAL PSTAR(ln2)            ! IN SURFACE PRESSURE (PA)
C
      REAL PRECIP(ln2,KCT+1)     ! IN PRECIPITATION ADDED WHEN
                                 !    DESCENDING FROM LAYER K TO K-1
                                 !    (KG/M**2/S)
C
      REAL PRECIPX(ln2,KCT+1)    ! FINAL AMOUNT OF PRECIPITATION
                                 ! FROM EACH LAYER (KG/M**2/S)
C
      INTEGER ICCB(ln2)          ! IN CLOUD BASE LEVEL
C
      INTEGER ICCT(ln2)          ! IN CLOUD TOP LEVEL
C
      REAL CCA(ln2)              ! IN CONVECTIVE CLOUD AMOUNT
C
      LOGICAL BWATER(ln2,2:KCT+1)!IN  MASK FOR THOSE POINTS AT WHICH
                                 !     CONDENSATE IS WATER IN LAYER K
C
      LOGICAL BTERM(ln2)         ! IN MASK FOR THOSE POINTS WHERE
                                 !    UPDRAUGHT IS TERMINATING
C
      LOGICAL BGMK(ln2)          ! IN MASK FOR POINTS WHERE PARCEL IN
                                 !    LAYER K IS SATURATED
C
      REAL TIMESTEP
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------

      REAL DTHBYDT(ln2,KCT+1)    ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE (K/S)
C
      REAL DQBYDT(ln2,KCT+1)     ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING
                                 !     RATIO (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     MIXING RATIO (KG/KG/S)
C
      REAL DDFLX(ln2,KCT+1)      ! DOWNDRAUGHT MASSFLUX IN LAYER K (PA/S)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL RAIN(ln2)        ! OUT RAINFALL AT SURFACE (KG/M**2/S)
C
      REAL SNOW(ln2)        ! OUT SNOWFALL AT SURFACE (KG/M**2/S)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL EXNER_KM12_C(ln2)     ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K
C
      REAL EXNER_KP12_C(ln2)     ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K+1
C
      REAL EXNER_KM32_C(ln2)     ! COMPRESSED EXNER FUNCTION AT
                                 ! LAYER K-1
C
      REAL PK(ln2)               ! PRESSURE OF LAYER K (PA)
C
      REAL P_KM1(ln2)            ! PRESSURE OF LAYER K-1 (PA)
C
      REAL EXK(ln2)              ! EXNER RATIO FOR LAYER K
C
      REAL EXKM1(ln2)            ! EXNER RATIO FOR LAYER K-1
C
      REAL DELPK(ln2)            ! PRESSURE DIFFERENCE ACROSS LAYER K
                                 ! (PA)
C
      REAL DELPKM1(ln2)          ! PRESSURE DIFFERENCE ACROSS
                                 ! LAYER K-1 (PA)
C
      REAL AMDETK(ln2)           ! MIXING DETRAINMENT AT LEVEL K
                                 ! MULTIPLIED BY APPROPRIATE LAYER
                                 ! THICKNESS
C
      REAL EKM14(ln2)            ! EXNER RATIO AT LAYER K-1/4
C
      REAL EKM34(ln2)            ! EXNER RATIO AT LAYER K-3/4
C
      LOGICAL BWATER_K_C(ln2)    ! COMPRESSED MASK FOR THOSE
                                 ! POINTS AT WHICH CONDENSATE
                                 ! IS WATER IN LAYER K
C
      REAL PRECIP_K_C(ln2)       ! COMPRESSED PRECIPITATION
                                 ! ADDED WHEN DESCENDING FROM
                                 ! LAYER K TO K-1 (KG/M**2/S)
C
      REAL Q_K_C(ln2)            ! COMPRESSED PARCEL MIXING RATIO
                                 ! OF LAYER K (KG/KG)
C
      REAL TH_K_C(ln2)           ! COMPRESSED PARCEL POTENTIAL
                                 ! TEMPERATURE OF LAYER K (K)
C
      REAL PSTAR_C(ln2)          ! COMPRESSED SURFACE PRESSURE (PA)
C
      INTEGER ICCB_C(ln2)        ! COMPRESSED CLOUD BASE LEVEL
C
      REAL DTHBYDT_K_C(ln2)      ! COMPRESSED INCREMENT TO MODEL
                                 ! POTENTIAL TEMPERATURE OF LAYER K
                                 ! (K/S)
C
      REAL DTHBYDT_KM1_C(ln2)    ! COMPRESSED INCREMENT TO MODEL
                                 ! POTENTIAL TEMPERATURE OF LAYER K-1
                                 ! (K/S)
C
      REAL DQBYDT_K_C(ln2)       ! COMPRESSED INCREMENT TO MODEL
                                 ! MIXING RATIO OF LAYER K (KG/KG/S)
C
      REAL DQBYDT_KM1_C(ln2)     ! COMPRESSED INCREMENT TO MODEL
                                 ! MIXING RATIO OF LAYER K-1 (KG/KG/S)
C
      REAL DELTD(ln2)            ! COOLING NECESSARY TO
                                 ! ACHIEVE SATURATION (K)
C
      REAL DELQD(ln2)            ! MOISTENING NECESSARY TO
                                 ! ACHIEVE SATURATION (KG/KG)
C
      REAL QDD_K(ln2)            ! MIXING RATIO OF DOWNDRAUGHT IN
                                 ! LAYER K (KG/KG)
C
      REAL THDD_K(ln2)           ! MODEL POTENTIAL TEMPERATURE
                                 ! OF DOWNDRAUGHT IN LAYER K (K)
C
      REAL FLX_DD_K(ln2)         ! DOWNDRAUGHT INITIAL MASS FLUX
                                 ! (PA/S)
C
      REAL FLX_DD_K_C(ln2)       ! COMPRESSED DOWNDRAUGHT INITIAL
                                 ! MASS FLUX (PA/S)
C
      LOGICAL BDDI(ln2)          ! MASK FOR POINTS WHERE DOWNDRAUGHT
                                 ! MIGHT OCCUR
C
      LOGICAL BDDI_C(ln2)        ! COMPRESSED MASK FOR POINTS WHERE
                                 ! DOWNDRAUGHT MAY INITIATE
C
      INTEGER INDEX1(ln2)        ! INDEX FOR COMPRESS AND EXPAND
C
      REAL QE_K_C(ln2)           ! COMPRESSED ENVIRONMENT MIXING
                                 ! RATIO OF LAYER K (KG/KG)
C
      REAL QE_KM1_C(ln2)         ! COMPRESSED ENVIRONMENT MIXING
                                 ! RATIO OF LAYER K-1 (KG/KG)
C
      REAL THE_K_C(ln2)          ! COMPRESSED POTENTIAL TEMPERATURE
                                 ! OF ENVIRONMENT IN LAYER K (K)
C
      REAL THE_KM1_C(ln2)        ! COMPRESSED POTENTIAL TEMPERATURE
                                 ! OF ENVIRONMENT IN LAYER K-1 (K)
C
      REAL RAIN_C(ln2)           ! COMPRESSED SURFACE RAINFALL
                                 ! (KG/M**2/S)
C
      REAL SNOW_C(ln2)           ! COMPRESSED SURFACE SNOWFALL
                                 ! (KG/M**2/S)
C
      REAL FLX_UD_K_C(ln2)       ! UPDRAUGHT MASS FLUX AT LAYER K
C
      REAL RAIN_ENV(ln2)         ! AMOUNT OF RAINFALL PASSING THROUGH
                                 ! ENVIRONMENT (KG/M**2/S)
C
      REAL SNOW_ENV(ln2)         ! AMOUNT OF SNOWFALL PASSING THROUGH
                                 ! ENVIRONMENT (KG/M**2/S)
C
      REAL RAIN_DD(ln2)          ! AMOUNT OF RAINFALL PASSING THROUGH
                                 ! DOWNDRAUGHT (KG/M**2/S)
C
      REAL SNOW_DD(ln2)          ! AMOUNT OF SNOWFALL PASSING THROUGH
                                 ! DOWNDRAUGHT (KG/M**2/S)
C
      LOGICAL BDD_START(ln2)     ! MASK FOR THOSE POINT WHERE
                                 ! DOWNDRAUGHT IS ABLE TO START
                                 ! FROM LEVEL K
C
      LOGICAL BDD_START_C(ln2)   ! COMPRESSED MASK FOR THOSE POINT
                                 ! WHERE DOWNDRAUGHT IS ABLE TO START
                                 ! FROM LEVEL K
C
      LOGICAL BDDWT_K(ln2)       ! MASK FOR POINTS IN DOWNDRAUGHT
                                 ! WHERE PPT IN LAYER K IS LIQUID
C
      LOGICAL BDDWT_K_C(ln2)     ! COMPRESSED MASK FOR POINTS IN DD
                                 ! WHERE PPT IN LAYER K IS LIQUID
C
      LOGICAL BDDWT_KM1(ln2)     ! MASK FOR POINTS IN DOWNDRAUGHT
                                 ! WHERE PPT IN LAYER K-1 IS LIQUID
C
      LOGICAL BDDWT_KM1_C(ln2)   ! COMPRESSED MASK FOR POINTS IN DD
                                 ! WHERE PPT IN LAYER K-1 IS LIQUID
C
      LOGICAL BDD_ON(ln2)        ! MASK FOR THOSE POINTS WHERE DD
                                 ! CONTINUES FROM LAYER K+1
C
      LOGICAL BDD_ON_C(ln2)      ! COMPRESSED MASK FOR POINTS WHERE DD
                                 ! CONTINUES FROM LAYER K+1
C
      INTEGER KMIN(ln2)          ! FREEZING LEVEL WHERE ENTRAINMENT
                                 ! RATES ARE INCREASED
C
      REAL FLX_STRT(ln2)         ! MASSFLUX AT LEVEL WHERE DOWNDRAUGHT
                                 ! STARTS (PA/S)
C
      REAL FLX_STRT_C(ln2)       ! COMPRESSED VALUE OF FLX_STRT
C
      REAL CCA_C(ln2)            ! COMPRESSED CONVECTIVE CLOUD AMOUNT
C
      REAL LR_UD_REF(ln2)        ! PRECIPITATION MIXING RATIO AT LOWEST
                                 ! PRECIPITATING LEVEL OF UD
C
C
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL FLX_INIT, LAYER_DD, DD_INIT, DOWND
C
C-----------------------------------------------------------------------
C CALCULATE INDEX FOR COMPRESS ON BASIS OF BTERM
C-----------------------------------------------------------------------
C
      NDD = 0
      DO I=1,NPNTS
       IF (BTERM(I)) THEN
          NDD = NDD+1
          INDEX1(NDD) = I
       END IF
      END DO
C
C----------------------------------------------------------------------
C INITIALISE LOGICAL ARRAYS AS FALSE
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       FLX_DD_K(I)=0.0
       BDDI(I) = .FALSE.
       BDD_START(I) = .FALSE.
       BDDWT_K(I) = .FALSE.
       BDDWT_KM1(I) = .FALSE.
       BDD_ON(I) = .FALSE.
C
C-----------------------------------------------------------------------
C CALCULATE MASK FOR THOSE POINT WHERE DOWNDRAUGHT MIGHT OCCUR
C AND LEVEL AT WHICH IT MIGHT INITIATE
C-----------------------------------------------------------------------
C
        IF (KCT .GE. 4 .AND. BTERM(I) .AND. BGMK(I) .AND. (KCT-ICCB(I))
     &       .GT. 2)  BDDI(I) = .TRUE.
      END DO
C
C----------------------------------------------------------------------
C CALCULATE INITIAL DOWNDRAUGHT MASS FLUX
C-----------------------------------------------------------------------
C
      IF (KCT .GE. 4)
     &  CALL FLX_INIT (NPNTS,KCT,ICCB,ICCT,FLX,FLX_DD_K,BDDI,FLX_STRT)
C
C-----------------------------------------------------------------------
C COMPRESS ALL INPUT ARRAYS FOR THE DOWNDRAUGHT CALCULATION
C-----------------------------------------------------------------------
      DO 10 K = KCT+1,2,-1
C
*vdir nodep
         DO I=1,NDD
            TH_K_C(I) = THP(INDEX1(I),K)
            Q_K_C(I) = QP(INDEX1(I),K)
            THE_K_C(I) = THE(INDEX1(I),K)
            THE_KM1_C(I) = THE(INDEX1(I),K-1)
            QE_K_C(I) = QE(INDEX1(I),K)
            QE_KM1_C(I) = QE(INDEX1(I),K-1)
            DTHBYDT_K_C(I) = DTHBYDT(INDEX1(I),K)
            DTHBYDT_KM1_C(I) = DTHBYDT(INDEX1(I),K-1)
            DQBYDT_K_C(I) = DQBYDT(INDEX1(I),K)
            DQBYDT_KM1_C(I) = DQBYDT(INDEX1(I),K-1)
            EXNER_KM12_C(I) = EXNER(INDEX1(I),K)
            EXNER_KP12_C(I) = EXNER(INDEX1(I),K+1)
            EXNER_KM32_C(I) = EXNER(INDEX1(I),K-1)
            PRECIP_K_C(I) = PRECIP(INDEX1(I),K)
            FLX_UD_K_C(I) = FLX(INDEX1(I),K)
            BWATER_K_C(I) = BWATER(INDEX1(I),K)
         END DO
         IF (K.EQ.KCT+1) THEN
*vdir nodep
          DO I=1,NDD
            FLX_DD_K_C(I) = FLX_DD_K(INDEX1(I))
            FLX_STRT_C(I) = FLX_STRT(INDEX1(I))
            PSTAR_C(I) = PSTAR(INDEX1(I))
            ICCB_C(I) = ICCB(INDEX1(I))
            BDDI_C(I) = BDDI(INDEX1(I))
            BDD_START_C(I) = BDD_START(INDEX1(I))
            RAIN_C(I) = RAIN(INDEX1(I))
            SNOW_C(I) = SNOW(INDEX1(I))
            BDDWT_K_C(I) = BDDWT_K(INDEX1(I))
            BDDWT_KM1_C(I) = BDDWT_KM1(INDEX1(I))
            BDD_ON_C(I) = BDD_ON(INDEX1(I))
            CCA_C(I) = CCA(INDEX1(I))
            LR_UD_REF(I) = 0.0
          END DO
         END IF
C
C----------------------------------------------------------------------
C IF BELOW CONVECTIVE CLOUD BASE DOWNDRAUGHT NOT ALLOWED TO FORM
C----------------------------------------------------------------------
C
      DO I=1,NDD
       IF (K.LT.ICCB_C(I)) BDDI_C(I)=.FALSE.
      END DO
C
C-----------------------------------------------------------------------
C RESET EN/DETRAINMENT RATES FOR DOWNDRAUGHT
C-----------------------------------------------------------------------
C
      CALL LAYER_DD (NDD,K,KCT,THE_K_C,THE_KM1_C,FLX_STRT_C,AK,BK,
     &               AKM12,BKM12,DELAK,DELBK,EXNER_KM12_C,EXNER_KP12_C,
     &               EXNER_KM32_C,PSTAR_C,PK,P_KM1,DELPK,DELPKM1,EXK,
     &               EXKM1,AMDETK,EKM14,EKM34,KMIN,BDDI_C)
C
C-----------------------------------------------------------------------
C INITIALISE DOWNDRAUGHT
C DOWNDRAUGHT NOT ALLOWED TO FORM FROM CLOUD TOP LAYER (KCT+1)
C OR FROM BELOW CLOUD BASE
C-----------------------------------------------------------------------
C
      IF (KCT .GE. 4 .AND. K.LT.KCT+1)
     & CALL DD_INIT(NDD,TH_K_C,Q_K_C,THE_K_C,QE_K_C,PK,EXK,THDD_K,
     &              QDD_K,DELTD,DELQD,BDD_START_C,K,BDDI_C,BDD_ON_C)
C
C-----------------------------------------------------------------------
C UPDATE MASK FOR WHERE DOWNDRAUGHT OCCURS
C-----------------------------------------------------------------------
C
      DO I=1,NDD
c       IF (BDD_START_C(I).OR.BDD_ON_C(I)) BDD_ON_C(I)=.TRUE.
        BDD_ON_C(I)=BDD_ON_C(I).OR.BDD_START_C(I)
      END DO
C
      NDDON_TMP = 0
      DO I=1,NDD
        IF (BDD_ON_C(I)) THEN
          NDDON_TMP = NDDON_TMP+1
        END IF
      END DO
C
C-----------------------------------------------------------------------
C COLLECT MASSFLUX WHERE DOWNDRAUGHT OCCURS
C-----------------------------------------------------------------------
C
c     IF (K.EQ.KCT+1) THEN
c     DO I=1,NDD
c       DDFLX(INDEX1(I),K) = DDFLX(INDEX1(I),K) + FLX_DD_K_C(I)
c     END DO
c     ENDIF
C
C-----------------------------------------------------------------------
C CALL DOWNDRAUGHT ROUTINE
C-----------------------------------------------------------------------
C

      CALL DOWND(NDD,K,KCT,THDD_K,QDD_K,THE_K_C,THE_KM1_C,QE_K_C,
     &           QE_KM1_C,DTHBYDT_K_C,DTHBYDT_KM1_C,DQBYDT_K_C,
     &           DQBYDT_KM1_C,FLX_DD_K_C,P_KM1,DELPK,DELPKM1,EXK,
     &           EXKM1,DELTD,DELQD,AMDETK,EKM14,EKM34,PRECIP_K_C,
     &           RAIN_C,SNOW_C,ICCB_C,BWATER_K_C,BDD_START_C,
     &           BDDWT_K_C,BDDWT_KM1_C,BDD_ON_C,RAIN_ENV,SNOW_ENV,
     &           RAIN_DD,SNOW_DD,FLX_UD_K_C,TIMESTEP,CCA_C,NDDON_TMP,
     &           LR_UD_REF)
C
C-----------------------------------------------------------------------
C DECOMPRESS/EXPAND THOSE VARIABLES WHICH ARE TO BE OUTPUT
C-----------------------------------------------------------------------
C
C COLLECT MASSFLUX WHERE DOWNDRAUGHT OCCURS
*vdir nodep
      DO I=1,NDD
        DDFLX(INDEX1(I),K-1) = DDFLX(INDEX1(I),K-1) + FLX_DD_K_C(I)
      END DO
CDIR$ IVDEP
*vdir nodep
        DO I=1,NDD
         DTHBYDT(INDEX1(I),K) = DTHBYDT_K_C(I)
         DTHBYDT(INDEX1(I),K-1) = DTHBYDT_KM1_C(I)
         DQBYDT(INDEX1(I),K) = DQBYDT_K_C(I)
         DQBYDT(INDEX1(I),K-1) = DQBYDT_KM1_C(I)
         PRECIPX(INDEX1(I),K)=PRECIPX(INDEX1(I),K)+
     & RAIN_DD(I)+RAIN_ENV(I)+SNOW_DD(I)+SNOW_ENV(I)
         IF (K.EQ.2) THEN
          RAIN(INDEX1(I)) = RAIN_C(I)
          SNOW(INDEX1(I)) = SNOW_C(I)
         END IF
         PRECIP(INDEX1(I),K) = PRECIP_K_C(I)
        END DO
C
C----------------------------------------------------------------------
C   END OF MAIN K LOOP
C----------------------------------------------------------------------
C
 10   CONTINUE
C
      RETURN
      END
C
CLL  SUBROUTINE DD_ENV-------------------------------------------------
CLL
CLL  PURPOSE : CALCULATE THE EFFECT OF THE DOWNDRAUGHT
CLL            ON THE LARGE_SCALE ATMOSPHERE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DD_ENV (NPNTS,THDD_K,THDD_KM1,QDD_K,QDD_KM1,THE_K,
     &                   THE_KM1,QE_K,QE_KM1,DTHBYDT_K,DTHBYDT_KM1,
     &                   DQBYDT_K,DQBYDT_KM1,FLX_DD_K,FLX_DD_KM1,DELPK,
     &                   DELPKM1,DELTD,DELQD,AMDETK,EKM14,
     &                   B_DD_END,BDD_START,BDD_ON)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS                 ! IN VECTOR LENGTH
C
      INTEGER I                     ! LOOP COUNTER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THDD_K(ln2)              ! IN DOWNDRAUGHT POTENTIAL
                                    !    TEMPERATURE IN LAYER K (K)
C
      REAL THDD_KM1(ln2)            ! IN DOWNDRAUGHT POTENTIAL
                                    !    TEMPERATURE IN LAYER K-1 (K)
C
      REAL QDD_K(ln2)               ! IN DOWNDRAUGHT MIXING RATIO
                                    !    AT LAYER K (KG/KG)
C
      REAL QDD_KM1(ln2)             ! IN DOWNDRAUGHT MIXING RATIO
                                    !    AT LAYER K-1 (KG/KG)
C
      REAL THE_K(ln2)               ! IN POTENTIAL TEMPERATURE OF
                                    !    ENVIRONMENT IN LAYER K (K)
C
      REAL THE_KM1(ln2)             ! IN POTENTIAL TEMPERATURE OF
                                    !    ENVIRONMENT IN LAYER K-1 (K)
C
      REAL QE_K(ln2)                ! IN MIXING RATIO AT LAYER K (KG/KG)
C
      REAL QE_KM1(ln2)              ! IN MIXING RATIO AT LAYER K-1
                                    !    (KG/KG)
C
      REAL FLX_DD_K(ln2)            ! IN MASS FLUX IN LAYER K (PA/S)
C
      REAL FLX_DD_KM1(ln2)          ! IN MASS FLUX IN LAYER K-1 (PA/S)
C
      REAL DELPK(ln2)               ! IN DIFFERENCE IN PRESSURE ACROSS
                                    !    LAYER K (PA)
C
      REAL DELPKM1(ln2)             ! IN DIFFERENCE IN PRESSURE ACROSS
                                    !    LAYER K-1 (PA)
C
      REAL DELTD(ln2)               ! IN COOLING NECESSARY TO ACHIEVE
                                    !    SATURATION (K)
C
      REAL DELQD(ln2)               ! IN MOISTENING NECESSARY TO ACHIEVE
                                    !    SATURATION (KG/KG)
C
      REAL AMDETK(ln2)              ! IN MIXING DETRAINMENT AT LEVEL K
                                    !    MULTIPLIED BY APPROPRIATE LAYER
                                    !    THICKNESS
C
      REAL EKM14(ln2)               ! IN EXNER RATIO AT LAYER K-1/4
C
      LOGICAL B_DD_END(ln2)         ! IN MASK FOR THOSE POINTS WHERE
                                    !    DOWNDRAUGHT IS TERMINATING
C
      LOGICAL BDD_START(ln2)        ! IN MASK FOR THOSE POINTS WHERE
                                    !    DOWNDRAUGHT IS STARTING
C
      LOGICAL BDD_ON(ln2)           ! IN MASK FOR THOSE POINTS WHERE
                                    !    DOWNDRAUGHT IS ON
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL DTHBYDT_K(ln2)           ! INOUT
                                    ! IN  INCREMENT TO POTENTIAL
                                    !     TEMPERATURE IN LAYER K (K/S)
                                    ! OUT UPDATED INCREMENT TO POTENTIAL
                                    !     TEMPERATURE IN LAYER K (K/S)
C
      REAL DTHBYDT_KM1(ln2)         ! INOUT
                                    ! IN  INCREMENT TO POTENTIAL
                                    !     TEMPERATURE AT LAYER K-1 (K/S)
                                    ! OUT UPDATED INCREMENT TO POTENTIAL
                                    !     TEMPERATURE AT LAYER K-1 (K/S)
C
      REAL DQBYDT_KM1(ln2)          ! INOUT
                                    ! IN  INCREMENT TO MIXING RATIO AT
                                    !     LAYER K-1 (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MIXING
                                    !     RATIO AT LAYER K-1 (KG/KG/S)
C
      REAL DQBYDT_K(ln2)            ! INOUT
                                    ! IN  INCREMENT TO MIXING RATIO
                                    !     AT LAYER K (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MIXING
                                    !     RATIO AT LAYER K (KG/KG/S)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL TEMPRY                   ! USED IN CALCULATIONS OF THE
                                    ! EFFECT ON THE ENVIRONMENT
C
C-----------------------------------------------------------------------
C CALCULATE THE EFFECT ON THE ENVIRONMENT IN LAYER K
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       IF (BDD_ON(I)) THEN
C
C-----------------------------------------------------------------------
C SUBTRACT ENERGY USED TO FORM DOWNDRAUGHT
C-----------------------------------------------------------------------
C
       TEMPRY = FLX_DD_K(I)/DELPK(I)
       IF (BDD_START(I)) THEN
         DTHBYDT_K(I) = DTHBYDT_K(I)-TEMPRY*DELTD(I)
         DQBYDT_K(I) = DQBYDT_K(I)-TEMPRY*DELQD(I)
       END IF
C
C-----------------------------------------------------------------------
C EFFECT OF CONVECTION AND DOWNDRAUGHT UPON POTENTIAL TEMPERATURE OF
C LAYER K
C-----------------------------------------------------------------------
C
       DTHBYDT_K(I) = DTHBYDT_K(I) + TEMPRY * (
     &
     &          (1.0+EKM14(I)) * (1.0-AMDETK(I)) *      ! COMPENSATING
     &           (THE_KM1(I)-THE_K(I))                  ! SUBSIDENCE
     &        +
     &          AMDETK(I)* (THDD_K(I)-THE_K(I))         ! MIXING
     &        )                                         ! DETRAINMENT
C
C-----------------------------------------------------------------------
C EFFECT OF CONVECTION AND DOWNDRAUGHT UPON MIXING RATIO OF
C LAYER K
C-----------------------------------------------------------------------
C
       DQBYDT_K(I) = DQBYDT_K(I) + TEMPRY * (
     &
     &      (1.0+EKM14(I)) * (1.0-AMDETK(I)) *       ! COMPENSATING
     &      (QE_KM1(I)-QE_K(I))                      ! SUBSIDENCE
     &    +
     &      AMDETK(I)* (QDD_K(I)-QE_K(I))            ! MIXING
     &    )                                          ! DETRAINMENT
C
C-----------------------------------------------------------------------
C TERMINAL DETRAINMENT AND SUBSIDENCE IN TERMINAL LAYER
C-----------------------------------------------------------------------
C
         IF (B_DD_END(I)) THEN
           TEMPRY = FLX_DD_KM1(I)/DELPKM1(I)
           DTHBYDT_KM1(I) = DTHBYDT_KM1(I)+TEMPRY*
     &                      (THDD_KM1(I)-THE_KM1(I))
           DQBYDT_KM1(I) = DQBYDT_KM1(I)+TEMPRY*(QDD_KM1(I)-QE_KM1(I))
         END IF
C
       END IF
      END DO
C
      RETURN
      END
C
CLL  SUBROUTINE DD_INIT------------------------------------------------
CLL
CLL  PURPOSE : ROUTINE TO INITIALISE THE DOWNDRAUGHT
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DD_INIT(NPNTS,TH_UD_K,Q_UD_K,THE_K,QE_K,PK,EXK,THDD_K,
     &                   QDD_K,DELTD,DELQD,BDD_START,K,BDDI,BDD_ON)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                 ! LOOP COUNTER
C
      INTEGER NPNTS             ! VECTOR LENGTH
C
      INTEGER K                 ! IN PRESENT MODEL LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THE_K(ln2)           ! IN POTENTIAL TEMPERATURE OF
                                !    ENVIRONMENT IN LAYER K (K)
C
      REAL TH_UD_K(ln2)         ! IN PARCEL POTENTIAL TEMPERATURE OF
                                !    UPDRAUGHT, LAYER K (K)
C
      REAL QE_K(ln2)            ! IN MIXING RATIO OF ENVIRONMENT IN
                                !    LAYER K (KG/KG)
C
      REAL Q_UD_K(ln2)          ! IN PARCEL MIXING RATIO OF UPDRAUGHT,
                                !    LAYER K (KG/KG)
C
      REAL EXK(ln2)             ! IN EXNER RATIO OF LAYER K
C
      REAL PK(ln2)              ! IN PRESSURE OF LAYER K (PA)
C
      LOGICAL BDDI(ln2)         ! IN MASK FOR THOSE POINTS WHERE
                                !    DOWNDRAUGHT MAY INITIATE
C
      LOGICAL BDD_ON(ln2)       ! IN MASK FOR THOSE POINTS WHERE
                                !    DOWNDRAUGHT IS ON
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL BDD_START(ln2)    ! INOUT
                                ! IN  MASK FOR THOSE POINT WHERE
                                !     DOWNDRAUGHT MAY START
                                ! OUT MASK FOR THOSE POINTS WHERE
                                !
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL THDD_K(ln2)          ! OUT DOWNDRAUGHT POTENTIAL TEMPERATURE
                                !     OF LAYER K
C
      REAL QDD_K(ln2)           ! OUT DOWNDRAUGHT MIXING RATIO OF
                                !     LAYER K
C
      REAL DELTD(ln2)           ! OUT COOLING NECESSARY TO ACHIEVE
                                !     SATURATION
C
      REAL DELQD(ln2)           ! OUT MOISTENING NECESSARY TO ACHIEVE
                                !     SATURATION
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL TH_MEAN(ln2)         ! MEAN POTENTIAL TEMPERATURE USED IN
                                ! CALCULATION OF SATURATED DOWNDRAUGHT
                                ! POTENTIAL TEMPERATURE IN LAYER K
C
      REAL Q_MEAN(ln2)          ! MEAN MIXING RATIO USED IN CALCULATION
                                ! OF SATURATED DOWNDRAUGHT
                                ! MIXING RATIO FOR LAYER K
C
      REAL THDDS(ln2)           ! SATURATED DOWNDRAUGHT POTENTIAL
                                ! TEMPERATURE IN LAYER K (K)
C
      REAL QDDS(ln2)            ! SATURATED DOWNDRAUGHT MIXING RATIO
                                ! IN LAYER K (KG/KG)
C
      REAL BUOY                 ! BUOYANCY OF PARCEL IN LAYER K
C
      REAL THDD_V               ! VIRTUAL POTENTIAL TEMPERATURE OF
                                ! PARCEL IN LAYER K
C
      REAL THE_V                ! VIRTUAL POTENTIAL TEMPERATURE OF
                                ! ENVIRONMENT IN LAYER K
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL SATCAL
C
C-----------------------------------------------------------------------
C CALCULATE MEAN TEMPERATURE AND MIXING RATIO
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       TH_MEAN(I) = (THE_K(I)+TH_UD_K(I))*0.5
       Q_MEAN(I) = (QE_K(I)+Q_UD_K(I))*0.5
      END DO
C
C-----------------------------------------------------------------------
C CALCULATE SATURATED DOWNDRAUGHT POTENTIAL TEMPERATURE FOR LAYER K
C-----------------------------------------------------------------------
C
      CALL SATCAL(NPNTS,TH_MEAN,PK,QDDS,THDDS,EXK,Q_MEAN,THE_K)
C
C-----------------------------------------------------------------------
C IS SATURATED PARCEL NEGATIVELY BUOYANT COMPARED TO ENVIRONMENT
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       IF (.NOT. BDD_ON(I) .AND. BDDI(I) .AND. K.GE.4) THEN
          THDD_V = THDDS(I)*(1.0+0.61*QDDS(I))
          THE_V = THE_K(I)*(1.0+0.61*QE_K(I))
          BUOY = THDD_V - THE_V
C
          IF (BUOY .LT. 0.5 ) THEN
C
C-----------------------------------------------------------------------
C INITIATE DOWNDRAUGHT
C-----------------------------------------------------------------------
C
             THDD_K(I) = THDDS(I)
             QDD_K(I) = QDDS(I)
             BDD_START(I) = .TRUE.
C
C-----------------------------------------------------------------------
C CALCULATE COOLING AND MOISTENING TO ACHIEVE SATURATION
C-----------------------------------------------------------------------
C
             DELTD(I) = THDDS(I)-THE_K(I)
             DELQD(I) = QDDS(I)-QE_K(I)
          END IF
       END IF
      END DO
C
      RETURN
      END
C
CLL  SUBROUTINE DDRAUGHT-----------------------------------------------
CLL
CLL  PURPOSE : DOWNDRAUGHT ROUTINE
CLL
CLL            CONVECTIVE DOWNDRAUGHT BASED ON PARCEL THEORY
CLL
CLL            CARRY OUT DRY DESCENT
CLL
CLL            CALCULATE SUBSATURATION
CLL
CLL            CALCULATE EFFECT ON THE ENVIRONMENT
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DDRAUGHT (NPNTS,K,KCT,THDD_K,QDD_K,THE_K,THE_KM1,QE_K,
     &                     QE_KM1,DTHBYDT_K,DTHBYDT_KM1,DQBYDT_K,
     &                     DQBYDT_KM1,FLX_DD_K,P_KM1,DELPK,
     &                     DELPKM1,EXK,EXKM1,DELTD,DELQD,AMDETK,EKM14,
     &                     EKM34,RAIN,SNOW,BDD_START,BDDWT_K,
     &                     BDDWT_KM1,BDD_ON,B_DD_END,CCA)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
CHBG  include 'C_0_DG_C.f'
CHBG  include 'DDKMDET.f'
      include 'UKPARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                     ! LOOP COUNTER
C
      INTEGER NPNTS                 ! IN NUMBER OF POINTS
C
      INTEGER K                     ! IN PRESENT MODEL LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      INTEGER KCT                   ! IN CONVECTIVE CLOUD TOP
C
      REAL THE_KM1(ln2)             ! IN POTENTIAL TEMPERATURE OF
                                    !    ENVIRONMENT IN LAYER K-1 (K)
C
      REAL QE_KM1(ln2)              ! IN MIXING RATIO OF ENVIRONMENT IN
                                    !    LAYER K-1 (KG/KG)
C
      REAL P_KM1(ln2)               ! IN PRESSURE OF LAYER K-1 (PA)
C
      REAL DELPK(ln2)               ! IN CHANGE IN PRESSURE ACROSS
                                    !    LAYER K (PA)
C
      REAL DELPKM1(ln2)             ! IN CHANGE IN PRESSURE ACROSS
                                    !    LAYER K-1 (PA)
C
      REAL EXK(ln2)                 ! IN EXNER RATIO IN LAYER K
C
      REAL EXKM1(ln2)               ! IN EXNER RATIO IN LAYER K-1
C
      REAL AMDETK(ln2)              ! IN MIXING DETRAINMENT RATE
C
      REAL EKM14(ln2)               ! IN EXNER RATIO AT LAYER K-1/4
C
      REAL EKM34(ln2)               ! IN EXNER RATIO AT LAYER K-3/4
C
      REAL DELTD(ln2)               ! IN COOLING NECESSARY TO ACHIEVE
                                    !    SATURATION (K)
C
      REAL DELQD(ln2)               ! IN MOISTENING NECESSARY TO ACHIEVE
                                    !    SATURATION (KG/KG)
C
      LOGICAL BDDWT_K(ln2)          ! IN MASK FOR THOSE POINTS IN
                                    !    DOWNDRAUGHT WHERE PRECIPITATION
                                    !    IS LIQUID IN LAYER K
C
      LOGICAL BDDWT_KM1(ln2)        ! IN MASK FOR THOSE POINTS IN
                                    !    DOWNDRAUGHT WHERE PRECIPITATION
                                    !    IS LIQUID IN LAYER K-1
C
      REAL CCA(ln2)                 ! IN CONVECTIVE CLOUD AMOUNT
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL THDD_K(ln2)              ! INOUT
                                    ! IN  POTENTIAL TEMPERATURE OF
                                    !     DOWNDRAUGHT IN LAYER K (K)
                                    ! OUT POTENTIAL TEMPERATURE RESET
                                    !     FOR NEXT LAYER (K)
C
      REAL QDD_K(ln2)               ! INOUT
                                    ! IN  DOWNDRAUGHT MIXING RATIO OF
                                    !     LAYER K (KG/KG)
                                    ! OUT MIXING RATIO RESET FOR NEXT
                                    !     LAYER (KG/KG)
C
      REAL THE_K(ln2)               ! INOUT
                                    ! IN  POTENTIAL TEMPERATURE OF
                                    !     ENVIRONMENT IN LAYER K (K)
                                    ! OUT ENVIRONMENT POTENTIAL
                                    !     TEMPERATURE RESET FOR NEXT
                                    !     LAYER (K)
C
      REAL QE_K(ln2)                ! INOUT
                                    ! IN  MIXING RATIO OF ENVIRONMENT
                                    !     LAYER K (KG/KG)
                                    ! OUT ENVIRONMENT MIXING RATIO
                                    !     RESET FOR NEXT LAYER (KG/KG)
C
      REAL FLX_DD_K(ln2)            ! INOUT
                                    ! IN  DOWNDRAUGHT MASS FLUX OF
                                    !     LAYER K (PA/S)
                                    ! OUT DOWNDRAUGHT MASS FLUX RESET
                                    !      FOR NEXT LAYER (PA/S)
C
      REAL RAIN(ln2)                ! INOUT
                                    ! IN  AMOUNT OF RAIN (KG/M**2/S)
                                    ! OUT UPDATED RAINFALL (KG/M**2/S)
C
      REAL SNOW(ln2)                ! INOUT
                                    ! IN  AMOUNT OF SNOW(KG/M**2/S)
                                    ! OUT UPDATED SNOWFALL (KG/M**2/S)
C
      REAL DTHBYDT_K(ln2)           ! INOUT
                                    ! IN  INCREMENT TO MODEL POTENTIAL
                                    !     TEMPERATURE OF LAYER K (K/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     POTENTIAL TEMPERATURE IN
                                    !     LAYER K (K/S)
C
      REAL DTHBYDT_KM1(ln2)         ! INOUT
                                    ! IN  INCREMENT TO MODEL POTENTIAL
                                    !     TEMPERATURE IN LAYER K-1 (K/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     POTENTIAL TEMPERATURE IN
                                    !     LAYER K-1 (K/S)
C
      REAL DQBYDT_K(ln2)            ! INOUT
                                    ! IN  INCREMENT TO MODEL MIXING
                                    !     RATIO IN LAYER K (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     MIXING RATIO IN LAYER K
                                    !     (KG/KG/S)
C
      REAL DQBYDT_KM1(ln2)          ! INOUT
                                    ! IN  INCREMENT TO MODEL MIXING
                                    !     RATIO IN LAYER K-1 (KG/KG/S)
                                    ! OUT UPDATED INCREMENT TO MODEL
                                    !     MIXING RATIO IN LAYER K-1
                                    !     (KG/KG/S)
C
      LOGICAL BDD_ON(ln2)           ! INOUT
                                    ! IN  MASK FOR THOSE POINTS WHERE DD
                                    !     HAS CONTINUED FROM LAYER K+1
                                    ! OUT MASK FOR THOSE POINTS WHERE DD
                                    !     CONTINUES TO LAYER K-1
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL BDD_START(ln2)        ! OUT MASK FOR THOSE POINTS WHERE
                                    !     DOWNDRAUGHT MAY START IN
                                    !     LAYER K-1
C
      LOGICAL B_DD_END(ln2)         ! OUT MASK FOR THOSE POINTS WHERE
                                    !     DOWNDRAUGHT IS ENDING IN
                                    !     LAYER K-1
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL THDD_KM1(ln2)            ! POTENTIAL TEMPERATURE OF
                                    ! DOWNDRAUGHT IN LAYER K-1 (K)
C
      REAL QDD_KM1(ln2)             ! DOWNDRAUGHT MIXING RATIO OF
                                    ! LAYER K-1 (KG/KG)
C
      REAL QSATDD(ln2)              ! SATURATED DOWNDRAUGHT MIXING
                                    ! RATIO (KG/KG)
C
      REAL TDD_KM1(ln2)             ! TEMPERATURE OF DOWNDRAUGHT
                                    ! IN LAYER K-1 (K)
C
      REAL THDDS(ln2)               ! POTENTIAL TEMPERATURE OF
                                    ! SATURATED DOWNDRAUGHT (K)
C
      REAL QDDS(ln2)                ! SATURATED DOWNDRAUGHT MIXING
                                    ! RATIO (KG/KG)
C
      REAL FLX_DD_KM1(ln2)          ! DOWNDRAUGHT MASS FLUX IN
                                    ! LAYER K-1 (PA/S)
C
      REAL RAIN_TMP(ln2)            ! LIQUID PRECIPITATION STORE
C
      REAL SNOW_TMP(ln2)            ! SNOW STORE
C
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL SATCAL, CRS_FRZL, QSATU, DEVAP, TERMDD,
     &         DD_ENV, EVP
C
C-----------------------------------------------------------------------
C CALCULATE MASK FOR THOSE POINTS IN DOWNDRAUGHT WHERE PRECIPITATION
C IS LIQUID
C
C STORE PRECIPITATION IN LAYER K IN TEMPORARY VARIABLES
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        IF (K .EQ. KCT+1 .OR. BDD_START(I)) THEN
          BDDWT_K(I) = EXK(I)*THDD_K(I) .GT. TM
        ELSE
          BDDWT_K(I) = BDDWT_KM1(I)
        END IF
          RAIN_TMP(I) = RAIN(I)
          SNOW_TMP(I) = SNOW(I)
C
C-----------------------------------------------------------------------
C DRY DESCENT FROM LAYER K TO K-1
C
C ENTRAINMENT CALCULATION
C-----------------------------------------------------------------------
C
          THDD_KM1(I) = (THDD_K(I)+(EKM14(I)*THE_K(I)) +
     &                  (1.0+EKM14(I))*EKM34(I)*THE_KM1(I)) /
     &                  ((1.0+EKM14(I))*(1.0+EKM34(I)))
          QDD_KM1(I) = (QDD_K(I)+(EKM14(I)*QE_K(I)) +
     &                 (1.0+EKM14(I))*EKM34(I)*QE_KM1(I))/
     &                 ((1.0+EKM14(I))*(1.0+EKM34(I)))
C
C-----------------------------------------------------------------------
C UPDATE MASS FLUX  AND CALCULATE TEMPERATURE OF LAYER K-1
C-----------------------------------------------------------------------
C
          FLX_DD_KM1(I) = FLX_DD_K(I)*(1.0+EKM34(I))*(1.0+EKM14(I))*
     &                (1.0-AMDETK(I))
C
          TDD_KM1(I) = THDD_KM1(I)*EXKM1(I)
      END DO
C
C-----------------------------------------------------------------------
C CALCULATE SUBSATURATION
C CALCULATE TEMPERATURE IF BROUGHT TO SATURATION
C-----------------------------------------------------------------------
C
       CALL SATCAL(NPNTS,THDD_KM1,P_KM1,QDDS,THDDS,
     &             EXKM1,QDD_KM1,THE_KM1)
C
      DO I=1,NPNTS
        BDDWT_KM1(I) = THDDS(I) .GT. TM/EXKM1(I)
      END DO
C
C-----------------------------------------------------------------------
C CALCULATE CHANGE OF PHASE DUE TO DOWNDRAUGHT SATURATION TEMPERATURE
C-----------------------------------------------------------------------
C
       CALL CRS_FRZL (NPNTS,RAIN,SNOW,THDD_KM1,EXKM1,FLX_DD_KM1,
     &                BDDWT_KM1)
C
      DO I=1,NPNTS
        TDD_KM1(I) = THDD_KM1(I)*EXKM1(I)
      END DO
C
C-----------------------------------------------------------------------
C RECALCULATE SUBSATURATION TEMPERATURE
C-----------------------------------------------------------------------
C
       CALL SATCAL(NPNTS,THDD_KM1,P_KM1,QDDS,THDDS,
     &             EXKM1,QDD_KM1,THE_KM1)
C
C-----------------------------------------------------------------------
C CALCULATE MOISTURE SUBSATURATION
C-----------------------------------------------------------------------
C
       CALL QSATU(QSATDD,TDD_KM1,P_KM1,NPNTS)
C
C-----------------------------------------------------------------------
C EVAPORATION CALCULATION AND ADJUSTMENT OF DOWNDRAUGHT TEMPERATURE
C AND MOISTURE
C-----------------------------------------------------------------------
C
       CALL DEVAP (NPNTS,THDD_K,THDD_KM1,QDD_KM1,THDDS,QDDS,
     &             FLX_DD_KM1,EXK,EXKM1,QSATDD,RAIN,SNOW,
     &             DELPKM1,BDDWT_KM1,CCA,P_KM1)
C
C-----------------------------------------------------------------------
C CHECK IF PARCEL STILL NEGATIVELY BUOYANT SUCH THAT DOWNDRAUGHT CAN
C CONTINUE TO K-1
C-----------------------------------------------------------------------
C
       CALL TERMDD (NPNTS,BDD_START,THDD_KM1,QDD_KM1,THE_KM1,
     &              QE_KM1,K,B_DD_END,BDD_ON)
C
C-----------------------------------------------------------------------
C CALCULATE THE EFFECT ON THE ENVIRONMENT IN LAYER K
C-----------------------------------------------------------------------
C
       CALL DD_ENV (NPNTS,THDD_K,THDD_KM1,QDD_K,QDD_KM1,THE_K,THE_KM1,
     &              QE_K,QE_KM1,DTHBYDT_K,DTHBYDT_KM1,DQBYDT_K,
     &              DQBYDT_KM1,FLX_DD_K,FLX_DD_KM1,DELPK,DELPKM1,
     &              DELTD,DELQD,AMDETK,EKM14,B_DD_END,
     &              BDD_START,BDD_ON)
C
C-----------------------------------------------------------------------
C RESET DOWNDRAUGHT BIT VECTORS
C
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       BDD_START(I) = .FALSE.
       IF (.NOT. BDD_ON(I)) THEN
         RAIN(I) = RAIN_TMP(I)
         SNOW(I) = SNOW_TMP(I)
       END IF
       IF (B_DD_END(I)) BDD_ON(I) = .FALSE.
      END DO
C
C-----------------------------------------------------------------------
C SWITCH POTENTIAL TEMPERATURE, MIXING RATIO AND MASS FLUX FOR
C CALCULATION AT NEXT MODEL LAYER
C-----------------------------------------------------------------------
C
      IF (K.GT.2) THEN
        DO I=1,NPNTS
         IF (BDD_ON(I)) THEN
          THDD_K(I) = THDD_KM1(I)
          QDD_K(I) = QDD_KM1(I)
          FLX_DD_K(I) = FLX_DD_KM1(I)
         END IF
        END DO
      END IF

      RETURN
      END
C
CLL  SUBROUTINE DETRAIN------------------------------------------------
CLL
CLL  PURPOSE : FORCED DETRAINMENT CALCULATION
CLL
CLL            SUBROUTINE THP_DET CALCULATES THE POTENTIAL
CLL            TEMPERATURE OF THE PARCEL IN LAYER K+1
CLL            AFTER FORCED DETRAINMENT
CLL
CLL            SUBROUTINE THETAR CALCULATES THE POTENTIAL TEMPERATURE
CLL            OF THE AIR IN LAYER K UNDERGOING FORCED DETRAINMENT
CLL
CLL            SUBROUTINE DET_RATE CALCULATES THE FORCED DETRAINMENT
CLL            RATE OF THE ENSEMBLE IN LAYER K
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.3   23/12/93  : DG020893 : TO MAKE CALCULATIONS OF FORCED
CLL                                DETRAINMENT RATE LESS PRONE TO
CLL                                FAILURE
CLL
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  SYSTEM TASK :
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DETRAIN (NPNTS,THEK,QEK,THPK,QPK,QSEK,DQSK,BGMK,
     &                     THEKP1,QEKP1,THPKP1,QPKP1,QSEKP1,DQSKP1,
     &                     BGMKP1,BWKP1,XSQKP1,
     &                     DELTAK,THRK,QRK,EKP14,EKP34,PK,PKP1,
     &                     EXK,EXKP1)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER I              ! LOOP COUNTER
C
      INTEGER NREDO          ! NUMBER OF POINTS FOR WHICH FORCED
                             ! DETRAINMENT CALCULATION MUST BE
                             ! AS THE PROCESSES EITHER CAUSES THE
                             ! PARCEL TO BECOME SATURATED OR
                             ! SUB-SATURATED
C
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      REAL THEK(ln2)         ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K (K)
C
      REAL THEKP1(ln2)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEK(ln2)          ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL QEKP1(ln2)        ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL QSEKP1(ln2)       ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL DQSKP1(ln2)       ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT IN LAYER K+1
                             !    (KG/KG/K)
C
      REAL THPK(ln2)         ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K (K)
C
      REAL QPK(ln2)          ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
C
      REAL QSEK(ln2)         ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL DQSK(ln2)         ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT OF LAYER K
                             !    (KG/KG/K)
C
      LOGICAL BWKP1(ln2)     ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      LOGICAL BGMK(ln2)      ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K
C
      REAL EKP14(ln2)        ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP34(ln2)        ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+3/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EXKP1(ln2)        ! IN EXNER RATIO AT LEVEL K+1
C
      REAL EXK(ln2)          ! IN EXNER RATIO AT LEVEL K
C
      REAL PKP1(ln2)         ! IN PRESSURE AT LEVEL K+1 (PA)
C
      REAL PK(ln2)           ! IN PRESSURE AT LEVEL K (PA)
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL THPKP1(ln2)       ! INOUT
                             ! IN  PARCEL POTENTIAL TEMPERATURE IN
                             !     LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (K)
                             ! OUT ADJUSTED PARCEL POTENTIAL
                             !     IN LAYER K+1 AFTER FORCED
                             !     DETRAINMENT (K)
C
      REAL QPKP1(ln2)        ! INOUT
                             ! IN  PARCEL MIXING RATIO IN
                             !     LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (KG/KG)
                             ! OUT ADJUSTED PARCEL POTENTIAL
                             !     IN LAYER K+1 AFTER FORCED
                             !     DETRAINMENT (KG/KG)
C
      REAL XSQKP1(ln2)       ! INOUT
                             ! IN  EXCESS WATER IN PARCEL AFTER
                             !     LIFTING FROM LAYER K TO K+1 AFTER
                             !     ENTRAINMENT AND LATENT HEATING
                             !     (KG/KG)
                             ! OUT EXCESS WATER IN PARCEL IN LAYER
                             !     K+1 AFTER FORCED DETRAINMENT
                             !     (KG/KG)
C
      LOGICAL BGMKP1(ln2)    ! INOUT
                             ! IN  MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1 AFTER
                             !     ENTRAINMENT AND LATENT HEATING
                             ! OUT MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1 AFTER
                             !     FORCED DETRAINMENT
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL THRK(ln2)         ! OUT PARCEL DETRAINMENT POTENTIAL
                             !     TEMPERATURE IN LAYER K (K)
C
      REAL QRK(ln2)          ! OUT PARCEL DETRAINMENT MIXING RATIO
                             !     IN LAYER K (KG/KG)
C
      REAL DELTAK(ln2)       ! OUT PARCEL FORCED DETRAINMENT RATE
                             !     IN LAYER K
C
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C
C ON THE IBM ARRAYS ARE ALLOCATED USING A PARAMETER STATEMENT
C
C ON THE CRAY ARRAYS ARE DYNAMICALLY ALLOCATED
C-----------------------------------------------------------------------
C
      LOGICAL BDETK(ln2)     ! MASK FOR PARCELS WHICH ARE
                             ! UNDERGOING FORCED DETRAINMENT
                             ! IN THEIR ASCENT FROM LAYER K
                             ! TO K+1
C
      REAL XSQR(ln2)         ! EXCESS PARCEL WATER VAPOUR
                             ! DURING DETRAINMENT (KG/KG)
C
      REAL THPKP1W(ln2)   ,  ! TEMPORARY STOREAGE FOR PARCEL
     &     QPKP1W(ln2)   ,   ! POTENTIAL TEMPERATURE (K), MIXING
     &     XSQK1W(ln2)       ! RATIO (KG/KG), EXCESS WATER VAPOUR
      LOGICAL BGKP1W(ln2)    ! (KG/KG) AND MASK FOR SATURATION
                             ! IN LAYER K+1
C
      LOGICAL BRECAL(ln2)    ! MASK FOR THOSE POINTS AT WHICH THE
                             ! THE DETRAINMENT CALCULATION NEEDS
                             ! REPEATING
C
      REAL TT(ln2)           ! TEMPORARY STORE FOR TEMPERATURE
                             ! FOR THE CALCULATION OF SATURATED
                             ! MIXING RATIO (K)
C
      REAL EPSS              ! (1+EKP14)*(1+EKP34)
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL THP_DET,QSATU,THETAR,DET_RATE
C
C*---------------------------------------------------------------------
C
      DO 10 I=1,NPNTS
C
C----------------------------------------------------------------------
C AT START OF ROUTINE FORCED DETARINMENT DONE AT ALL POINTS SO
C SET ARRAY BDETK EQUAL TO .TRUE.
C SET FORCED DETRAINMENT RATE EQUAL TO ZERO
C----------------------------------------------------------------------
C
       BDETK(I) = .TRUE.
       DELTAK(I) = 0.0
C
C-----------------------------------------------------------------------
C   SAVE THE CURRENT VALUES OF QPKP1, XSQKP1 AND BGMKP1
C-----------------------------------------------------------------------
C
       THPKP1W(I) = THPKP1(I)
       QPKP1W(I) = QPKP1(I)
       XSQK1W(I) = XSQKP1(I)
       BGKP1W(I) = BGMKP1(I)
C
C-----------------------------------------------------------------------
C   ADD THE EXCESS WATER VAPOUR BACK INTO THE DETRAINING PARCELS
C-----------------------------------------------------------------------
C
       QPKP1(I) = QPKP1(I) + XSQKP1(I)
   10 CONTINUE
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE ENSEMBLE AVERAGE POTENTIAL TEMPERATURE IN LAYER K+1
CL  AT THE POINTS WHERE DETRAINMENT IS TAKING PLACE
CL
CL  SUBROUTINE THP_DET
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (28)
CL----------------------------------------------------------------------
CL
      CALL THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,DQSKP1,
     &              BGMKP1,BDETK)
CL
CL---------------------------------------------------------------------
CL  CHECK TO SEE IF SUFFICIENT EXCESS WATER VAPOUR IN THE
CL  INITIAL DRY ASCENT TO ALLOW PARCEL TO BE SATURATED
CL  IN LAYER K+1 AFTER FORCED DETRAINMENT
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (29)
CL
CL  NOTE : ONLY ALLOW PARCEL TO BE SATURATED IN LAYER K+1 IF
CL         SATURATED INITIALLY.  IT IS POSSIBLE FOR SMALL
CL         SUPERSATURATIONS TO IF SUBROUTINE LATENT_H CAUSES
CL         PARCEL TO BE COME UNSATURATED.  IN THIS CASE TREAT
CL         THE PARCEL AS UNSATURATED IN LAYER K+1
CL---------------------------------------------------------------------
CL
C
C-----------------------------------------------------------------------
C   CALCULATE THE EXCESS WATER VAPOUR IN LAYER K+1 AND RECALCULATE
C   BGMKP1 AND QPKP1.
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
C PRESSURE OF LAYER K FOR CALCULATION OF SATURATED
C MIXING RATIO
C-----------------------------------------------------------------------
C
      DO 25 I = 1,NPNTS
       TT(I) = THPKP1(I)*EXKP1(I)
   25 CONTINUE
      CALL QSATU (XSQKP1,TT,PKP1,NPNTS)
C
      DO 30 I=1,NPNTS
       XSQKP1(I) = QPKP1(I) - XSQKP1(I)
C
       BRECAL(I) = BGMKP1(I)
C
C----------------------------------------------------------------------
C ONLY ALLOW PARCEL TO BE SATURATED IN INITIAL BGMKP1 = .TRUE.
C (STORED IN BRECAL AT THIS POINT)
C----------------------------------------------------------------------
C
       IF ( BGMK(I) .OR.( (XSQKP1(I) .GT. 0.) .AND. BRECAL(I) ) ) THEN
         BGMKP1(I) = .TRUE.
       ELSE
         BGMKP1(I) = .FALSE.
         XSQKP1(I) = 0.0
       END IF
C
       QPKP1(I) = QPKP1(I) - XSQKP1(I)
CL
CL----------------------------------------------------------------------
CL  RECALCULATE THE ENSEMBLE AVERAGE POTENTIAL TEMPERATURE AT POINTS
CL  WHERE THE ENSEMBLE HAS BECOME UNSATURATED.
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (28)
CL----------------------------------------------------------------------
CL
       BRECAL(I) = BDETK(I) .AND. BRECAL(I) .AND. .NOT.BGMKP1(I)
   30 CONTINUE
C
      CALL THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,DQSKP1,
     &             BGMKP1,BRECAL)
CL
CL----------------------------------------------------------------------
CL  BECAUSE OF THE REMOVAL OF LATENT HEATING, THE NEW PARCEL POTENTIAL
CL  TEMPERATURE MAY BE LOWER THAN ITS VALUE BEFORE THE DETRAINMENT
CL  CALCULATION. IN THIS CASE ABANDON THE DETRAINMENT CALCULATION.
CL----------------------------------------------------------------------
CL
      DO 90 I=1,NPNTS
       BDETK(I) = THPKP1(I) .GT. THPKP1W(I)
   90 CONTINUE
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE POTENTIAL TEMPERATURE AND MIXING RATIO  OF DETRAINING
CL  AIR AND THE EXCESS WATER VAPOUR CONDESED FROM DETRAINING AIR
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (26)
CL----------------------------------------------------------------------
CL
      CALL THETAR (NPNTS,THRK,QRK,XSQR,BGMK,THEK,QEK,QPK,QSEK,DQSK,
     &             EXK,PK)
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE DETRAINMENT RATE, DELTAK.
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (31)
CL----------------------------------------------------------------------
CL
      CALL DET_RATE (NPNTS,DELTAK,THRK,XSQR,THPK,THEK,THEKP1,
     &             XSQKP1,THPKP1,BWKP1,BDETK,EKP14,EKP34,EXK,EXKP1)
C
      NREDO = 0
CL
CL----------------------------------------------------------------------
CL  ADD WATER VAPOUR WHICH WAS REMOVED FROM DETRAINING AIR INTO XSQKP1
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION 86), EQUATION (11C)
CL----------------------------------------------------------------------
CL
      DO 120 I=1,NPNTS
C
       EPSS = (1.+EKP14(I))*(1.+EKP34(I))
C
       IF (BDETK(I))
     & XSQKP1(I) = XSQKP1(I) + (DELTAK(I)*XSQR(I)/
     &               (EPSS*(1.-DELTAK(I))))
CL
CL----------------------------------------------------------------------
CL  IF THE EXCESS WATER VAPOUR IN LAYER K+1 IS LESS THAN ZERO
CL  I.E. THE PARCEL HAS BECOME UNSATURATED THROUGH THE FORCED
CL  DETRAINMENT PROCESS THEN ABANDON THE CALCULATION
CL----------------------------------------------------------------------
CL
       BRECAL(I) = BGMKP1(I)
C
       BGMKP1(I) = XSQKP1(I) .GT. 0.
C
       BRECAL(I) = BDETK(I) .AND. BRECAL(I) .AND. .NOT.BGMKP1(I)
C
       IF (BRECAL(I)) THEN
          QPKP1(I)  = QPKP1(I) + XSQKP1(I)
     &               - (DELTAK(I)*XSQR(I)/(EPSS*(1.-DELTAK(I))))
          XSQKP1(I) = 0.
       ENDIF
C
C----------------------------------------------------------------------
C COUNT POINTS AT WHICH DETRAINMENT CALCULATION NEEDS REPEATING
C----------------------------------------------------------------------
C
       IF (BRECAL(I)) NREDO = NREDO + 1
  120 CONTINUE
CL
CL---------------------------------------------------------------------
CL  REPEAT CALCULATION OF PARCEL POTENTIAL TEMPERATURE, DETRAINMENT
CL  RATE AND EXCESS PARCEL WATER IF THE PARCEL BECOMES UNSATURATED
CL  IN LAYER K+1 AFTER FORCED DETARINMENT
CL---------------------------------------------------------------------
CL
      IF (NREDO .GT. 0) THEN
C
C----------------------------------------------------------------------
C  CALCULATE NEW PARCEL POTENTIAL TEMPERATURE IN LAYER K+1
C  AFTER FORCED DETRAINMENT
C----------------------------------------------------------------------
C
        CALL THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,DQSKP1,
     &               BGMKP1,BRECAL)
C
C----------------------------------------------------------------------
C  CHECK IF FORCED DETRAINMENT STILL POSSIBLE AND RESET RECALCUATION
C  MASK TO FALSE IF IT IS NOT
C----------------------------------------------------------------------
C
        DO 130 I=1,NPNTS
          IF (BRECAL(I)) THEN
            BDETK(I) = THPKP1(I) .GT. THPKP1W(I)
            BRECAL(I) = BDETK(I)
          END IF
  130   CONTINUE
C
C----------------------------------------------------------------------
C  RCALCULATE FORCED DETRAINEMNT RATE
C----------------------------------------------------------------------
C
        CALL DET_RATE (NPNTS,DELTAK,THRK,XSQR,THPK,THEK,THEKP1,
     &             XSQKP1,THPKP1,BWKP1,BRECAL,EKP14,EKP34,EXK,EXKP1)
C
C----------------------------------------------------------------------
C  RECALCULATE EXCESS WATER VAPOUR IN LAYER K+1
C  AFTER FORCED DETRAINMENT
C----------------------------------------------------------------------
C
        DO 140 I=1,NPNTS
         IF (BRECAL(I)) THEN
            EPSS = (1.+EKP14(I))*(1.+EKP34(I))
            XSQKP1(I) = XSQKP1(I) + (DELTAK(I)*XSQR(I)/
     &                       (EPSS*(1.-DELTAK(I))))
         END IF
  140   CONTINUE
C
      END IF
CL
CL----------------------------------------------------------------------
CL  MAKE SURE THAT THE DETRAINMENT RATE IS BETWEEN 0 AND 1
CL
CL  IF <0 THEN NO DETRAINMENT OCCURS AND ORIGINAL VALUES ARE
CL  RESTORED
CL
CL  IF >1 THEN SET TO 1 AND THRK = THPK, QRK = QPK AND VALUES
CL  IN LAYER K+1 ARE RESTORED.  ALTHOUGH THESE ARE NOT USED
CL  IN ANY THERMODYNAMIC CALCULATION THEY ARE USED TO SPECIFY
CL CLOUD TOP IN SUBROUTIBE CONRAD
CL----------------------------------------------------------------------
CL
      DO 180 I=1,NPNTS
C
       IF (BDETK(I)) THEN
C
        IF (DELTAK(I).LE.0.0) THEN
           BDETK(I) = .FALSE.
           THPKP1 (I) = THPKP1W(I)
           QPKP1 (I) = QPKP1W(I)
           XSQKP1(I) = XSQK1W(I)
           BGMKP1(I) = BGKP1W(I)
           DELTAK(I) = 0.0
        END IF
        IF (DELTAK(I).GT.1.0) THEN
           DELTAK(I) = 1.0
           THRK(I) = THPK(I)
           QRK(I) = QPK(I)
           THPKP1 (I) = THPKP1W(I)
           QPKP1 (I) = QPKP1W(I)
           XSQKP1(I) = XSQK1W(I)
           BGMKP1(I) = BGKP1W(I)
        END IF
C
       ENDIF
  180  CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE DET_RATE-----------------------------------------------
CLL
CLL  PURPOSE : CALCULATES THE FORCED DETRAINMENT RATE IN LAYER K
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  SYSTEM TASK :
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL                  SECTION (6), EQUATION (31)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DET_RATE (NPNTS,DELTAK,THRK,XSQR,THPK,THEK,THEKP1,
     &                   XSQKP1,THPKP1,BWKP1,BCALC,EKP14,EKP34,
     &                   EXK,EXKP1)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
CHBG  include 'C_R_CP.f'
CHBG  include 'C_LHEAT.f'
      include 'UKPARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS            ! VECTOR LENGTH
C
      INTEGER I                ! LOOP COUNTER
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THRK(ln2)           ! IN PARCEL DETRAINMENT POTENTIAL
                               !    TEMPERATURE IN LAYER K (K)
C
      REAL XSQR(ln2)           ! IN EXCESS WATER VAPOUR OF THE
                               !    DETRAINING AIR IN LAYER K (KG/KG)
C
      REAL THPK(ln2)           ! IN PARCEL POTENTIAL TEMPERATURE
                               !    IN LAYER K (K)
C
      REAL THEK(ln2)           ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                               !    IN LAYER K (K)
C
      REAL THEKP1(ln2)         ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                               !    IN LAYER K+1 (K)
C
      REAL XSQKP1(ln2)         ! IN EXCESS WATER VAPOUR OF THE PARCEL
                               !    IN LAYER K+1 (KG/KG)
C
      REAL THPKP1(ln2)         ! IN PARCEL POTENTIAL TEMPERATURE
                               !    IN LAYER K+1 (K)
C
      LOGICAL BCALC(ln2)       ! IN MASK FOR POINTS AT WHICH
                               !    CALCULATIONS OF THIS ROUTINE
                               !    ARE NEEDED
C
      LOGICAL BWKP1(ln2)       ! IN MASK FOR THOSE POINTS AT WHICH
                               !    CONDENSATE IS LIQUID IN LAYER K+1
C
      REAL EKP14(ln2)          ! IN ENTRAINEMNT RATE FOR LEVEL K+1/4
                               !    MULTIPLIED BY APPROPRIATE LAYER
                               !    THICKNESS
C
      REAL EKP34(ln2)          ! IN ENTRAINEMNT RATE FOR LEVEL K+3/4
                               !    MULTIPLIED BY APPROPRIATE LAYER
                               !    THICKNESS
C
      REAL EXK(ln2)            ! IN EXNER RATIO FOR LEVEL K
C
      REAL EXKP1(ln2)          ! IN EXNER RATIO FOR LEVEL K+1
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL DELTAK(ln2)         ! OUT PARCEL FORCED DETRAINMENT RATE
                               !     IN LAYER K MULTIPLIED BY
                               !     APPROPRIATE LAYER THICKNESS
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL EL                  ! LATENT HEAT OF CONDENSATION OR
                               ! (CONDENSATION + FUSION) (J/KG)
C
      REAL EPSS                ! (1+EKP14)*(1+EKP34)
C
C*---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL  NO SIGNIFICANT STRUCTURE
CL---------------------------------------------------------------------
CL
C
      DO 10 I=1,NPNTS
       EPSS = (1. + EKP14(I)) * (1. + EKP34(I))
C
C-----------------------------------------------------------------------
C   CREATE A VECTOR OF LATENT HEATS
C-----------------------------------------------------------------------
C
       IF (BWKP1(I)) THEN
          EL = LC
       ELSE
          EL = LCPLF
       ENDIF
C
C-----------------------------------------------------------------------
C   CALCULATE DETRAINMENT RATES
C-----------------------------------------------------------------------
C
       IF (BCALC(I)) THEN
          DELTAK(I) = EKP14(I)*THEK(I)
     &        + EKP34(I)*(1.+EKP14(I))*THEKP1(I)
     &        - EPSS*(THPKP1(I) - EL/(EXKP1(I)*CP) * XSQKP1(I))
C
          DELTAK(I) =   (DELTAK(I) + THPK(I))
     &              /(DELTAK(I) + THRK(I) - EL/(EXK(I)*CP) * XSQR(I))
C
C----------------------------------------------------------------------
C  FROM A THEORETICAL VIEW POINT DELTAK CANNOT = 1 . HOWEVER
C  BECAUSE OF APPROXIMATION USED IN THE CALCULATION NUMERICALLY IT
C  MAY BE POSSIBLE.  HENCE IF DELTAK = 1 SET IT TO SLIGHTLY SMALLER
C  THAN 1
C----------------------------------------------------------------------
C
          IF (DELTAK(I).EQ.1.0) DELTAK(I) = 0.99999
C
       ENDIF
   10 CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE DEVAP--------------------------------------------------
CLL
CLL  PURPOSE : EVAPORATION ROUTINE
CLL
CLL            CARRIES OUT EVAPORATION AND UPDATES PRECIPITATION
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1990
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DEVAP(NPNTS,THDD_K,THDD_KM1,QDD_KM1,THDDS,QDDS,
     &                 FLX_DD_KM1,EXK,EXKM1,QSATDD,RAIN,SNOW,
     &                 DELPKM1,BDDWT_KM1,CCA,PKM1)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
CHBG  include 'C_G.f'
CHBG  include 'C_LHEAT.f'
CHBG  include 'C_R_CP.f'
CHBG  include 'DDAREA.f'
      include 'UKPARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I               ! LOOP COUNTER
C
      INTEGER NPNTS           ! IN VECTOR LENGTH
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THDDS(ln2)         ! IN SATURATED POTENTIAL
                              !    TEMPERATURE OF DOWNDRAUGHT
                              !    (K)
C
      REAL QDDS(ln2)          ! IN MIXING RATIO OF SATURATED
                              !    DOWNDRAUGHT (KG/KG)
C
      REAL FLX_DD_KM1(ln2)    ! IN DOWNDRAUGHT MASS FLUX IN
                              !    LAYER K-1 (PA/S)
C
      REAL THDD_K(ln2)        ! IN POTENTIAL TEMPERATURE OF
                              !    DOWNDRAUGHT IN LAYER K (K)
C
      REAL EXK(ln2)           ! IN EXNER RATIO OF LAYER K
C
      REAL EXKM1(ln2)         ! IN EXNER RATIO OF LAYER K-1
C
      REAL QSATDD(ln2)        ! IN SATURATED DOWNDRAUGHT
                              !    MIXING RATIO (KG/KG)
C
      REAL DELPKM1(ln2)       ! IN CHANGE IN PRESSURE ACROSS
                              !    LAYER K-1 (PA)
C
      LOGICAL BDDWT_KM1(ln2)  ! IN MASK WHERE PRECIPITATION IN
                              !    DOWNDRAUGHT IS LIQUID
C
      REAL CCA(ln2)           ! IN CONVECTIVE CLOUD AMOUNT
C
      REAL PKM1(ln2)          ! IN PRESSURE OF LAYER K-1
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL THDD_KM1(ln2)      ! INOUT
                              ! IN  POTENTIAL TEMPERATURE OF
                              !     DOWNDRAUGHT IN LAYER K-1 (K)
                              ! OUT UPDATED POTENTIAL TEMPERATURE
                              !     OF DOWNDRAUGHT IN LAYER K-1
                              !     AFTER EVAPORATION OR
                              !     SATURATION (K)
C
      REAL QDD_KM1(ln2)       ! INOUT
                              ! IN  MODEL MIXING RATIO OF
                              !     DOWNDRAUGHT IN LAYER K-1
                              !     (KG/KG)
                              ! OUT UPDATED MODEL MIXING RATIO
                              !     OF DOWNDRAUGHT IN LAYER K-1
                              !     AFTER EVAPORATION OR
                              !     SATURATION (KG/KG)
C
      REAL RAIN(ln2)          ! INOUT
                              ! IN  AMOUNT OF RAIN (KG/M**2/S)
                              ! OUT UPDATED RAINFALL (KG/M**2/S)
C
      REAL SNOW(ln2)          ! INOUT
                              ! IN  AMOUNT OF SNOW (KG/M**2/S)
                              ! OUT UPDATED SNOWFALL (KG/M**2/S)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE LOCALLY DEFINED
C-----------------------------------------------------------------------
C
C
      REAL TEVP(ln2)          ! TEMPERATURE USED IN EVAPORATION
                              ! CALCULATION (K)
C
      LOGICAL BEVAP(ln2)      ! MASK FOR THOSE POINTS AT WHICH
                              ! EVAPORATION CALCULATION IS TO
                              ! BE CARRIED OUT
C
      LOGICAL BSAT(ln2)       ! MASK FOR THOSE POINTS WHICH
                              ! ARE SUBSATURATED
C
      REAL EVAP_RAIN(ln2)     ! AMOUNT OF EVAPORATION OF RAIN
C
      REAL SUB_SNOW(ln2)      ! AMOUNT OF SNOW SUBLIMATION
C
      REAL DELQ(ln2)          ! DIFFERENCE IN MIXING RATIOS
                              ! (KG/KG)
C
      REAL DELTH(ln2)         ! INCREMENT TO DOWNDRAUGHT POTENTIAL
                              ! TEMPERATURE IN LAYER K-1 DUE TO
                              ! EVAPORATION
C
      REAL DELQE(ln2)         ! INCREMENT TO DOWNDRAUGHT MIXING RATIO
                              ! IN LAYER K-1 DUE TO EVAPORATION
C
      REAL DELTHS(ln2)        ! SATURATED POTENTIAL TEMPERATURE MINUS
                              ! POTENTIAL TEMPERATURE OF DOWNDRAUGHT
C
      REAL FACTOR(ln2)        ! DELTHS / DELTH
C
      REAL PINCR(ln2)         ! INCREASE IN PRECIPITATION IF PARCEL
                              ! SUPERSATURATES
C
      REAL RHO(ln2)           ! DENSITY OF AIR IN PARCEL
C
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL EVP
C
C-----------------------------------------------------------------------
C CHECK IF EVAPORATION POSSIBLE
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       DELQ(I) = QSATDD(I)-QDD_KM1(I)
C
       BEVAP(I) =((RAIN(I).GT.0.0) .OR. (SNOW(I).GT.0.0))
     &             .AND. (DELQ(I).GT.0.0)
       BSAT(I) = DELQ(I) .LT. 0.0
C
C-----------------------------------------------------------------------
C CALCULATE TEMPERATURE USED IN CALCULATION OF EVAPORATION CONSTANTS
C BASED ON TEMPERATURE OF PARCEL AFTER UNSATURATED DESCENT
C-----------------------------------------------------------------------
C
        IF (BEVAP(I)) THEN
          TEVP(I) = ((THDD_K(I)*EXK(I))+(THDD_KM1(I)*EXKM1(I)))*0.5
          RHO(I) = PKM1(I) / (R*TEVP(I))
        END IF
      END DO
C
C-----------------------------------------------------------------------
C EVAPORATION CALCULATION - CALCULATE RATES FOR RAIN AND SNOW
C-----------------------------------------------------------------------
C
      CALL EVP(NPNTS,RAIN,TEVP,CCA,RHO,DELQ,DELPKM1,EVAP_RAIN,
     &         BEVAP,1,DDCLDFRA)
C
      CALL EVP(NPNTS,SNOW,TEVP,CCA,RHO,DELQ,DELPKM1,SUB_SNOW,
     &         BEVAP,2,DDCLDFRA)
C
      DO I=1,NPNTS
       IF (BEVAP(I)) THEN
C
C-----------------------------------------------------------------------
C ADJUST EVAPORATION AND SUBLIMNATION RATES BACK TO GRID BOX MEAN
C-----------------------------------------------------------------------
C
        EVAP_RAIN(I) = EVAP_RAIN(I)*CCA(I)*DDCLDFRA
        SUB_SNOW(I) = SUB_SNOW(I)*CCA(I)*DDCLDFRA
C
C-----------------------------------------------------------------------
C CHECK IF PARCEL SUPERSATURATED
C-----------------------------------------------------------------------
C
        DELTH(I) = -((LC*EVAP_RAIN(I))+((LC+LF)*SUB_SNOW(I)))*G/
     &           (CP*EXKM1(I)*FLX_DD_KM1(I))
        DELQE(I) = (EVAP_RAIN(I)+SUB_SNOW(I))*G/FLX_DD_KM1(I)
C
        DELTHS(I) = THDDS(I)-THDD_KM1(I)
        IF (DELTH(I).LT.DELTHS(I)) THEN
C
C-----------------------------------------------------------------------
C ADJUST EVAP AND SUBLIMATION RATES TO GIVE SATURATION
C-----------------------------------------------------------------------
C
          FACTOR(I) = DELTHS(I)/DELTH(I)
          DELTH(I) = DELTHS(I)
          DELQE(I) = DELQE(I)*FACTOR(I)
          EVAP_RAIN(I) = EVAP_RAIN(I)*FACTOR(I)
          SUB_SNOW(I) = SUB_SNOW(I)*FACTOR(I)
        END IF
C
C-----------------------------------------------------------------------
C UPDATE T,Q AND PRECIPITATION
C-----------------------------------------------------------------------
C
        RAIN(I) = RAIN(I)-EVAP_RAIN(I)
        SNOW(I) = SNOW(I)-SUB_SNOW(I)
        THDD_KM1(I) = THDD_KM1(I)+DELTH(I)
        QDD_KM1(I) = QDD_KM1(I)+DELQE(I)
C
C-----------------------------------------------------------------------
C PARCEL IS SUPERSATURATED BEFORE EVAPORATION OCCURS
C BRING PARCEL TO SATURATION AND PRECIPITATE WATER
C-----------------------------------------------------------------------
C
      ENDIF
      IF (BSAT(I)) THEN
         PINCR(I) = (QDD_KM1(I)-QDDS(I))*FLX_DD_KM1(I)/G
         QDD_KM1(I) = QDDS(I)
         THDD_KM1(I) = THDDS(I)
         IF (BDDWT_KM1(I)) THEN
           RAIN(I) = RAIN(I)+PINCR(I)
         ELSE
           SNOW(I) = SNOW(I)+PINCR(I)
         END IF
      END IF
      END DO
C
      RETURN
      END
C
CLL  SUBROUTINE DOWND--------------------------------------------------
CLL
CLL  PURPOSE : CALL DOWNDRAUGHT CALCULATION
CLL
CLL            CHANGE OF PHASE CALCULATION WHERE NO DOWNDRAUGHT OCCURS
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.3   23/12/93 : DG060893 : CORRECTION TO PREVENT OVER PREDICTION
CLL                               OF SNOW SHOWERS; CHANGE TO CALCULATION
CLL                               OF AMOUNT OF SNOW WHICH FALLS THROUGH
CLL                               THE DOWNDRAUGHT
CLL
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DOWND (NPNTS,K,KCT,THDD_K,QDD_K,THE_K,THE_KM1,QE_K,
     &                  QE_KM1,DTHBYDT_K,DTHBYDT_KM1,DQBYDT_K,
     &                  DQBYDT_KM1,FLX_DD_K,P_KM1,DELPK,DELPKM1,EXK,
     &                  EXKM1,DELTD,DELQD,AMDETK,EKM14,EKM34,PRECIP_K,
     &                  RAIN,SNOW,ICCB,BWATER_K,BDD_START,
     &                  BDDWT_K,BDDWT_KM1,BDD_ON,RAIN_ENV,SNOW_ENV,
     &                  RAIN_DD,SNOW_DD,FLX_UD_K,TIMESTEP,CCA,NDDON_A,
     &                  LR_UD_REF)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C

CHBG  include 'C_0_DG_C.f'
CHBG  include 'C_G.f'
      include 'UKPARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                  ! LOOP COUNTER
C
      INTEGER K                  ! IN PRESENT MODEL LAYER
C
      INTEGER NPNTS              ! IN NUMBER OF POINTS
C
      INTEGER NDDON,NDDON_A      ! NUMBER OF POINTS AT WHICH
                                 ! DOWNDRAUGHT DOES OCCUR
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      INTEGER KCT                ! IN CONVECTIVE CLOUD TOP LAYER
C
      REAL THDD_K(ln2)           ! IN MODEL POTENTIAL TEMPERATURE
                                 !    OF DOWNDRAUGHT IN LAYER K (K)
C
      REAL QDD_K(ln2)            ! IN MIXING RATIO OF DOWNDRAUGHT IN
                                 !    LAYER K (KG/KG)
C
      REAL THE_K(ln2)            ! IN POTENTIAL TEMPERATURE OF
                                 !    ENVIRONMENT IN LAYER K (K)
C
      REAL THE_KM1(ln2)          ! IN POTENTIAL TEMPERATURE OF
                                 !    ENVIRONMENT IN LAYER K-1 (K)
C
      REAL QE_K(ln2)             ! IN MIXING RATIO OF ENVIRONMENT IN
                                 !    LAYER K (KG/KG)
C
      REAL QE_KM1(ln2)           ! IN MIXING RATIO OF ENVIRONMENT IN
                                 !    LAYER K-1 (KG/KG)
C
      REAL FLX_DD_K(ln2)         ! IN DOWNDRAUGHT MASS FLUX OF LAYER K
                                 !    (PA/S)
C
      REAL P_KM1(ln2)            ! IN PRESSURE OF LAYER K-1 (PA)
C
      REAL DELPK(ln2)            ! IN PRESSURE DIFFERENCE ACROSS
                                 !    LAYER K (PA)
C
      REAL DELPKM1(ln2)          ! IN PRESSURE DIFFERENCE ACROSS
                                 !    LAYER K-1 (PA)
C
      REAL EXK(ln2)              ! IN EXNER RATIO FOR LAYER K
C
      REAL EXKM1(ln2)            ! IN EXNER RATIO FOR LAYER K-1
C
      REAL PRECIP_K(ln2)         ! IN PRECIPITATION ADDED WHEN
                                 !    DESCENDING FROM LAYER K TO K-1
                                 !    (KG/M**2/S)
C
      REAL AMDETK(ln2)           ! IN MIXING DETRAINMENT AT LEVEL K
                                 !    MULTIPLIED BY APPROPRIATE LAYER
                                 !    THICKNESS
C
      REAL EKM14(ln2)            ! IN EXNER RATIO AT LAYER K-1/4
C
      REAL EKM34(ln2)            ! IN EXNER RATIO AT LAYER K-3/4
C
      REAL DELTD(ln2)            ! IN COOLING NECESSARY TO
                                 !    ACHIEVE SATURATION (K)
C
      REAL DELQD(ln2)            ! IN MOISTENING NECESSARY TO
                                 !    ACHIEVE SATURATION (KG/KG)
C
CHBG  REAL ICCB(ln2)             ! IN CLOUD BASE LEVEL
      INTEGER ICCB(ln2)          ! IN CONVECTIVE CLOUD BASE
C
      LOGICAL BWATER_K(ln2)      ! IN MASK FOR THOSE POINTS AT WHICH
                                 !    CONDENSATE IS WATER IN LAYER K
C
      LOGICAL BDDWT_K(ln2)       ! IN MASK FOR THOSE POINTS IN
                                 !    DOWNDRAUGHT WHERE PRECIPITATION
                                 !    IS LIQUID IN LAYER K
C
      LOGICAL BDDWT_KM1(ln2)     ! IN MASK FOR THOSE POINTS IN
                                 !    DOWNDRAUGHT WHERE PRECIPITATION
                                 !    IS LIQUID IN LAYER K-1
C
      REAL RAIN_ENV(ln2)         ! IN AMOUNT OF RAIN FALLING THROUGH
                                 !    THE ENVIRONMENT
C
      REAL SNOW_ENV(ln2)         ! IN AMOUNT OF SNOW FALLING THROUGH
                                 !    THE ENVIRONMENT
C
      REAL RAIN_DD(ln2)          ! IN AMOUNT OF RAIN FALLING THROUGH
                                 !    THE DOWNDRAUGHT
C
      REAL SNOW_DD(ln2)          ! IN AMOUNT OF SNOW FALLING THROUGH
                                 !    THE DOWNDRAUGHT
C
      REAL FLX_UD_K(ln2)         ! IN UPDRAUGHT MASSFLUX AT LAYER K
C
      REAL TIMESTEP              ! IN MODEL TIMESTEP (S)
C
      REAL CCA(ln2)              ! IN CONVECTIVE CLOUD AMOUNT
C
      REAL LR_UD_REF(ln2)        ! IN UD PPN MIXING RATION IN LOWEST
                                 !    PRECIPITATING LAYER IN UD
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL BDD_START(ln2)     ! INOUT
                                 ! IN  MASK FOR THOSE POINTS WHERE
                                 !     DOWNDRAUGHT MAY FORM IN LAYER K
                                 ! OUT MASK FOR THOSE POINTS WHERE
                                 !     DOWNDRAUGHT MAY FORM IN LAYER
                                 !     K-1
C
      REAL DTHBYDT_K(ln2)        ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE OF LAYER K (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE OF LAYER K
                                 !     (K/S)
C
      REAL DTHBYDT_KM1(ln2)      ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE OF LAYER K-1 (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE OF
                                 !     LAYER K-1 (K/S)
C
      REAL DQBYDT_K(ln2)         ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING
                                 !     RATIO OF LAYER K (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     MIXING RATIO OF LAYER K (KG/KG/S)
C
      REAL DQBYDT_KM1(ln2)       ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING
                                 !     RATIO OF LAYER K-1 (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE OF
                                 !     LAYER K-1 (KG/KG/S)
C
      REAL RAIN (ln2)            ! INOUT
                                 ! IN  INITIALISED RAINFALL (KG/M**2/S)
                                 ! OUT SURFACE RAINFALL (KG/M**2/S)
C
      REAL SNOW(ln2)             ! INOUT
                                 ! IN  INITIALISED SNOWFALL (KG/M**2/S)
                                 ! OUT SURFACE SNOWFALL (KG/M**2/S)
C
      LOGICAL BDD_ON(ln2)        ! INOUT
                                 ! IN  MASK FOR THOSE POINTS WHERE DD
                                 !     HAS CONTINUED FROM PREVIOUS LAYER
                                 ! OUT MASK FOR THOSE POINTS WHERE DD
                                 !     CONTINUES TO LAYER K-1
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL WORK(ln2,24)          !  WORK SPACE
C
      LOGICAL BWORK(ln2,5)       !  WORK SPACE FOR 'BIT' MASKS
C
      INTEGER INDEX1(ln2)        !  INDEX FOR COMPRESS AND
C
      LOGICAL B_DD_END(ln2)      !  MASK FOR POINTS WHERE DOWNDRAUGHT
                                 ! HAS ENDED
C
      REAL FACTOR                !  PROPORTION OF RAINFALL GOING INTO
                                 !  DOWNDRAUGHT FROM UD
C
      REAL FACTOR_ENV            !  PROPORTION OF RAINFALL GOING INTO
                                 !  DD FROM FALLING PPN
C
      REAL PPN_DD_REF            !  REFERENCE DD PPN MASS
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL CHG_PHSE, PEVP_BCB, DDRAUGHT
C
C-----------------------------------------------------------------------
C START OF MAIN LOOP
C   UPDATE PRECIPITATION AND CALCULATE MASK FOR WHERE PRECIPITATION
C   IS LIQUID
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        B_DD_END(I) = .FALSE.
      END DO
C
      IF (K.EQ.KCT+1) THEN
        DO I=1,NPNTS
         RAIN_DD(I) = 0.0
         RAIN_ENV(I) = 0.0
         SNOW_DD(I) = 0.0
         SNOW_ENV(I) = 0.0
        END DO
      END IF
C
C----------------------------------------------------------------------
C INJECTION OF PRECIPITATION FROM UD AT LEVEL K
C----------------------------------------------------------------------
C
      DO I=1,NPNTS
       FACTOR= 0.0
       IF (BDD_ON(I) .AND. FLX_UD_K(I).GT.0.0) THEN
        FACTOR = G * FLX_DD_K(I)/FLX_UD_K(I)
        FACTOR = MIN(FACTOR,1.0)
       END IF
c
       IF (BWATER_K(I)) THEN
        RAIN_DD(I) = RAIN_DD(I) + PRECIP_K(I)*FACTOR
        RAIN_ENV(I) = RAIN_ENV(I) + PRECIP_K(I)*(1.0-FACTOR)
       ELSE
        SNOW_DD(I) = SNOW_DD(I) + PRECIP_K(I)*FACTOR
        SNOW_ENV(I) = SNOW_ENV(I) + PRECIP_K(I)*(1.0-FACTOR)
       END IF
c
      END DO
C
C----------------------------------------------------------------------
C INTERACTION OF DOWNDRAUGHT WITH RESERVE OF PRECIPITATION OUTSIDE
C DOWNDRAUGHT
C
C BASED UPON CONTINUITY OF PRECIPITATION MIXING RATIO WITHIN
C DOWNDRAUGHT - EITHER AFTER INJECTION OF RAIN FROM UD IN LEVEL
C K OR WITH PPN MIXING RATIO IN LOWEST PRECIPITATING LAYER
C
C IF DOWNDRAUGHT INCREASES IN MASS THEN WATER INJECTED
C IF DOWNDRAUGHT DECREASES IN MASS THEN WATER IS REMOVED
C
C----------------------------------------------------------------------
C
      if(NDDON_A.gt.0)then
      DO I=1,NPNTS
C
       IF (BDD_ON(I)) THEN
C
        FACTOR_ENV = 0.0
        IF (PRECIP_K(I).GT.0.0) THEN
C
C---------------------------------------------------------------------
C CALCULATE NEW REFERENCE PPN MIXING RATIO
C DD PPN MIXING RATIO IN LAYER KM1 BASED ON CONTINUITY
C WITH THAT IN LAYER K
C---------------------------------------------------------------------
C
         LR_UD_REF(I) = G * PRECIP_K(I)/FLX_UD_K(I)
         PPN_DD_REF = RAIN_DD(I)+SNOW_DD(I)
        ELSE
C
C---------------------------------------------------------------------
C DD PPN MIXING RATIO IN LAYER KM1 BASED ON CONTINUITY
C WITH THAT IN LAST PRECIPITATING UD LAYER
C---------------------------------------------------------------------
C
         PPN_DD_REF = LR_UD_REF(I) * FLX_DD_K(I)
        END IF
C
C--------------------------------------------------------------------
C INJECT PPN INTO DD FROM PPN FALLING OUTSIDE OF THE DD
C--------------------------------------------------------------------
C
        IF (BDD_ON(I) .AND. (RAIN_ENV(I).GT.0.0 .OR.
     &                             SNOW_ENV(I).GT.0.0)) THEN
         FACTOR_ENV = ( (PPN_DD_REF * (1.0+EKM14(I))*
     &                    (1.0+EKM34(I))*(1.0-AMDETK(I))) -
     &                         (RAIN_DD(I)+SNOW_DD(I)) ) /
     &                          (RAIN_ENV(I)+SNOW_ENV(I))
         FACTOR_ENV = MIN(FACTOR_ENV,1.0)
         FACTOR_ENV = MAX(FACTOR_ENV,-1.0)
        END IF
C
        IF (FACTOR_ENV.GT.0.0) THEN
         RAIN_DD(I) = RAIN_DD(I) + RAIN_ENV(I)*FACTOR_ENV
         RAIN_ENV(I) = RAIN_ENV(I) * (1.0-FACTOR_ENV)
         SNOW_DD(I) = SNOW_DD(I) + SNOW_ENV(I)*FACTOR_ENV
         SNOW_ENV(I) = SNOW_ENV(I) * (1.0-FACTOR_ENV)
        ELSE
         RAIN_ENV(I) = RAIN_ENV(I) - RAIN_DD(I)*FACTOR_ENV
         RAIN_DD(I) = RAIN_DD(I) * (1.0+FACTOR_ENV)
         SNOW_ENV(I) = SNOW_ENV(I) - SNOW_DD(I)*FACTOR_ENV
         SNOW_DD(I) = SNOW_DD(I) * (1.0+FACTOR_ENV)
        END IF
C
       END IF
      END DO
      endif
C
C--------------------------------------------------------------------
C ZERO PRECIPITATION RATE IN LAYER K
C--------------------------------------------------------------------
C
      DO I=1,NPNTS
       PRECIP_K(I) = 0.0
      END DO
C
      if(NDDON_A.gt.0)then
C
C-----------------------------------------------------------------------
C COMPRESS OUT ON BASIS OF BIT VECTOR BDDON - THOSE POINTS WITH A
C DOWNDRAUGHT
C-----------------------------------------------------------------------
C
      NDDON=0 ! same as NDDON_A
      DO I=1,NPNTS
        IF (BDD_ON(I)) THEN
           NDDON = NDDON+1
           INDEX1(NDDON) = I
        END IF
      END DO
C
c     IF (NDDON .NE. 0) THEN
*vdir nodep
         DO I=1,NDDON
          WORK(I,1) = THDD_K(INDEX1(I))
          WORK(I,2) = QDD_K(INDEX1(I))
          WORK(I,3) = THE_K(INDEX1(I))
          WORK(I,4) = THE_KM1(INDEX1(I))
          WORK(I,5) = QE_K(INDEX1(I))
          WORK(I,6) = QE_KM1(INDEX1(I))
          WORK(I,7) = DTHBYDT_K(INDEX1(I))
          WORK(I,8) = DTHBYDT_KM1(INDEX1(I))
          WORK(I,9) = DQBYDT_K(INDEX1(I))
          WORK(I,10) = DQBYDT_KM1(INDEX1(I))
          WORK(I,11) = FLX_DD_K(INDEX1(I))
          WORK(I,12) = P_KM1(INDEX1(I))
          WORK(I,13) = DELPK(INDEX1(I))
          WORK(I,14) = DELPKM1(INDEX1(I))
          WORK(I,15) = EXK(INDEX1(I))
          WORK(I,16) = EXKM1(INDEX1(I))
          WORK(I,17) = DELTD(INDEX1(I))
          WORK(I,18) = DELQD(INDEX1(I))
          WORK(I,19) = AMDETK(INDEX1(I))
          WORK(I,20) = EKM14(INDEX1(I))
          WORK(I,21) = EKM34(INDEX1(I))
          WORK(I,22) = RAIN_DD(INDEX1(I))
          WORK(I,23) = SNOW_DD(INDEX1(I))
          WORK(I,24) = CCA(INDEX1(I))
          BWORK(I,1) = BDD_START(INDEX1(I))
          BWORK(I,2) = BDDWT_K(INDEX1(I))
          BWORK(I,3) = BDDWT_KM1(INDEX1(I))
          BWORK(I,4) = BDD_ON(INDEX1(I))
          BWORK(I,5) = B_DD_END(INDEX1(I))
      END DO
C
C-----------------------------------------------------------------------
C START DOWNDRAUGHT CALCULATION
C-----------------------------------------------------------------------
C
C
         CALL DDRAUGHT (NDDON,K,KCT,WORK(1,1),WORK(1,2),WORK(1,3),
     &                  WORK(1,4),WORK(1,5),WORK(1,6),WORK(1,7),
     &                  WORK(1,8),WORK(1,9),WORK(1,10),WORK(1,11),
     &                  WORK(1,12),WORK(1,13),WORK(1,14),
     &                  WORK(1,15),WORK(1,16),WORK(1,17),WORK(1,18),
     &                  WORK(1,19),WORK(1,20),WORK(1,21),WORK(1,22),
     &                  WORK(1,23),BWORK(1,1),BWORK(1,2),BWORK(1,3),
     &                  BWORK(1,4),BWORK(1,5),WORK(1,24))
C
C-----------------------------------------------------------------------
C EXPAND REQUIRED VECTORS BACK TO FULL FIELDS
C-----------------------------------------------------------------------
C
*vdir nodep
      DO I=1,NDDON
       THDD_K(INDEX1(I)) = WORK(I,1)
       QDD_K(INDEX1(I)) = WORK(I,2)
       DTHBYDT_K(INDEX1(I)) = WORK(I,7)
       DTHBYDT_KM1(INDEX1(I)) = WORK(I,8)
       DQBYDT_K(INDEX1(I)) = WORK(I,9)
       DQBYDT_KM1(INDEX1(I)) = WORK(I,10)
       FLX_DD_K(INDEX1(I)) = WORK(I,11)
       RAIN_DD(INDEX1(I)) = WORK(I,22)
       SNOW_DD(INDEX1(I)) = WORK(I,23)
       BDD_START(INDEX1(I)) = BWORK(I,1)
       BDDWT_K(INDEX1(I)) = BWORK(I,2)
       BDDWT_KM1(INDEX1(I)) = BWORK(I,3)
       BDD_ON(INDEX1(I)) = BWORK(I,4)
       B_DD_END(INDEX1(I)) = BWORK(I,5)
      END DO
c     END IF
      endif
C
C-----------------------------------------------------------------------
C RESET PRECIPITATION FALLING THROUGH ENVIRONMENT IF DOWNDRAUGHT
C DID NOT FORM
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
c     FACTOR=99.0
        IF (.NOT.BDD_ON(I).AND..NOT.B_DD_END(I)) THEN
          RAIN_ENV(I) = RAIN_ENV(I)+RAIN_DD(I)
          SNOW_ENV(I) = SNOW_ENV(I)+SNOW_DD(I)
          RAIN_DD(I) = 0.0
          SNOW_DD(I) = 0.0
c         FACTOR=-88.0
        END IF
      END DO
C
C-----------------------------------------------------------------------
C CARRY OUT CHANGE OF PHASE CALCULATION FOR PRECIPITATION FALLING
C THROUGH ENVIRONMENT
C-----------------------------------------------------------------------
C
         CALL CHG_PHSE (NPNTS,RAIN_ENV,SNOW_ENV,DTHBYDT_KM1,
     &                  EXK,EXKM1,DELPKM1,THE_K,THE_KM1)
C
C-----------------------------------------------------------------------
C EVAPORATE RAIN FALLING THROUGH ENVIRONMENT IF LAYER K BELOW
C CLOUD BASE
C-----------------------------------------------------------------------
C
         CALL PEVP_BCB (NPNTS,K-1,ICCB,THE_KM1,P_KM1,QE_KM1,DELPKM1,
     &                  RAIN_ENV,SNOW_ENV,DTHBYDT_KM1,DQBYDT_KM1,
     &                  EXKM1,TIMESTEP,CCA)
C
C-----------------------------------------------------------------------
C RESET PRECIPITATION FALLING THROUGH ENVIRONMENT IF DOWNDRAUGHT
C TERMINATES
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        IF (B_DD_END(I)) THEN
          RAIN_ENV(I) = RAIN_ENV(I)+RAIN_DD(I)
          SNOW_ENV(I) = SNOW_ENV(I)+SNOW_DD(I)
          RAIN_DD(I) = 0.0
          SNOW_DD(I) = 0.0
        END IF
      END DO
C
C-----------------------------------------------------------------------
C UPDATE RAIN AND SNOW
C-----------------------------------------------------------------------
C
       IF (K.EQ.2) THEN
         DO I=1,NPNTS
           RAIN(I) = RAIN(I)+RAIN_DD(I)+RAIN_ENV(I)
           SNOW(I) = SNOW(I)+SNOW_DD(I)+SNOW_ENV(I)
         END DO
       END IF
C
      RETURN
      END
C
CLL  SUBROUTINE DQS_DTH-------------------------------------------------
CLL
CLL  PURPOSE : CALCULATES GARDIENT OF SATURATION MIXING RATIO
CLL            WITH POTENTIAL TEMPERATURE FORM THE
CLL            CLAUSIUS-CLAPEYRON EQUATION
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  SYSTEM TASK :
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL                  SECTION(4), EQUATION (20)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE DQS_DTH (DQS,THEK,QSEK,EXK,NPNTS)
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include 'RV.f'
CHBG  include 'C_LHEAT.f'
      include 'UKPARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS        ! IN VECTOR LENGTH
C
      INTEGER I            ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL THEK(ln2)       ! IN POTENTIAL TEMPERATURE FOR LAYER K (K)
C
      REAL QSEK(ln2)       ! IN SATURATION MIXING RATIO FOR LAYER K (K)
C
      REAL EXK(ln2)        ! IN EXNER RATIO FOR LAYER K
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL DQS(ln2)        ! OUT GRADIENT OF SATURATION MIXING RATIO
                           !     WITH POTENTIAL TEMPERATURE FOR LAYER K
                           !     (KG/KG/S)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL EL              ! LATENT HEATING OF CONDENSATION OR
                           ! (CONDENSATION PLUS HEATING) (J/KG)
C
C*---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL NO SIGNIFICANT STRUCTURE
CL---------------------------------------------------------------------
CL
      DO 10 I=1,NPNTS
C
C-----------------------------------------------------------------------
C  CREATE A VECTOR OF LATENT HEATS ACCORDING TO WHETHER QSAT IS WRT
C  ICE OR WATER
C-----------------------------------------------------------------------
C
       IF (EXK(I)*THEK(I) .GT. 273.) THEN
          EL = LC
       ELSE
          EL = LCPLF
       ENDIF
C
C-----------------------------------------------------------------------
C CALCULATE D(QSAT)/D(THETA)
C-----------------------------------------------------------------------
C
       DQS(I) = EL*QSEK(I)/(EXK(I)*RV*THEK(I)*THEK(I))
   10  CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE ENVIRON------------------------------------------------
CLL
CLL  PURPOSE : CALCULATE THE EFFECT OF CONVECTION UPON THE
CLL            LARGE-SCALE ATMOSPHERE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  SYSTEM TASK :
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE ENVIRON (bland,NPNTS,DTHEK,DQEK,DTHEKP1,DQEKP1,
     &                    THEK,QEK,DELTAK,FLXK,THPK,QPK,
     &                    THRK,QRK,THEKP1,QEKP1,BTERM,THPKP1,
     &                    QPKP1,XPK,XPKP1,BWKP1,FLXKP1,BLOWST,
     &                    EKP14,EXK,EXKP1,DELPK,DELPKP1,AMDETK,PK)
CHBG &                    EKP14,EXK,EXKP1,DELPK,DELPKP1,AMDETK)
C
      IMPLICIT NONE
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
      include 'PARAMS.f'
CHBG  include 'C_R_CP.f'
CHBG  include 'C_LHEAT.f'
CHBG  include 'PARXS.f'
      include 'UKPARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER I              ! LOOP COUNTER
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE INPUT
C-----------------------------------------------------------------------
C
      LOGICAL BLAND(ln2)     ! IN LAND/SEA MASK

      REAL THEK(ln2)         ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K (K)
C
      REAL THEKP1(ln2)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEK(ln2)          ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL QEKP1(ln2)        ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL THPK(ln2)         ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K (K)
C
      REAL QPK(ln2)          ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
C
      REAL THPKP1(ln2)       ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K+1 (K)
C
      REAL QPKP1(ln2)        ! IN PARCEL MIXING RATIO IN LAYER K+1
                             !    (KG/KG)
C
      REAL XPK(ln2)          ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
C
      REAL FLXK(ln2)         ! IN PARCEL MASSFLUX IN LAYER K (PA/S)
C
      LOGICAL BWKP1(ln2)     ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      LOGICAL BTERM(ln2)     ! IN MASK FOR PARCELS WHICH TERMINATE IN
                             !    LAYER K+1
C
      LOGICAL BLOWST(ln2)    ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
C
      REAL THRK(ln2)         ! IN PARCEL DETRAINMENT POTENTIAL
                             !    TEMPERATURE IN LAYER K (K)
C
      REAL QRK(ln2)          ! IN PARCEL DETRAINMENT MIXING RATIO
                             !    IN LAYER K (KG/KG)
C
      REAL XPKP1(ln2)        ! IN PARCEL CLOUD WATER IN LAYER K+1
                             !    (KG/KG)
C
      REAL FLXKP1(ln2)       ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      REAL DELTAK(ln2)       ! IN PARCEL FORCED DETRAINMENT RATE
                             !    IN LAYER K MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP14(ln2)        ! IN ENTRAINMENT RATE FOR LEVEL K+1/4
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
C
      REAL EXK(ln2)          ! IN EXNER RATIO FOR MID-POINT OF LAYER K
C
      REAL EXKP1(ln2)        ! IN EXNER RATIO FOR MID-POINT OF
                             !    LAYER K+1
C
      REAL DELPK(ln2)        ! IN PRESSURE DIFFERENCE ACROSS LAYER K
                             !    (PA)
C
      REAL DELPKP1(ln2)      ! IN PRESSURE DIFFERENCE ACROSS LAYER K+1
                             !    (PA)
C
      REAL AMDETK(ln2)       ! IN MIXING DETRIANMENT AT LEVEL K
                             !    MULTIPLIED BY APPROPRIATE LAYER
                             !    THICKNESS
CHBG
      REAL PK(ln2)           ! IN PRESSURE AT MID-POINT OF LAYER K
                             !    (PA)
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL DTHEK(ln2)        ! INOUT
                             ! IN  INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (MAY BE NONE ZERO
                             !     DUE TO A PREVIOUS SPLIT FINAL
                             !     DETRAINMENT CALCULATION) (K/S)
                             ! OUT UPDATED INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K DUE TO
                             !     CONVECTION (K/S)
C
      REAL DQEK(ln2)         ! INOUT
                             ! IN  INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K DUE TO CONVECTION
                             !     (MAY BE NONE ZERO
                             !     DUE TO A PREVIOUS SPLIT FINAL
                             !     DETRAINMENT CALCULATION) (KG/KG/S)
                             ! OUT UPDATED INCREMENT TO MODEL MIXING
                             !     RATIO IN LAYER K DUE TO
                             !     CONVECTION (KG/KG/S)
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL DTHEKP1(ln2)      ! OUT INCREMENT TO MODEL POTENTIAL
                             !     TEMPERATURE IN LAYER K+1 DUE TO
                             !     CONVECTION (K/S)
C
      REAL DQEKP1(ln2)       ! OUT INCREMENT TO MODEL MIXING RATIO
                             !     IN LAYER K+1 DUE TO CONVECTION
                             !     (KG/KG/S)
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL EL                ! LATENT HEAT OF CONDENSATION OR
                             ! (CONDENSATION + FUSION) (J/KG)
C
      REAL TEMPRY            ! TEMPORARY ARRAY
CHBG
CSJP      REAL XXDTHE,XXDQE,PKP1FAC,PKFAC
      REAL FRACENV  ! FRACTION OF DETRAINED LIQUID WATER INTO ENVIRON
CHBG
      FRACENV=1.00
C
C*---------------------------------------------------------------------
C
      DO 10 I=1,NPNTS
C
C-----------------------------------------------------------------------
C   CREATE A VECTOR OF LATENT HEATS
C-----------------------------------------------------------------------
C
       IF (BWKP1(I)) THEN
          EL = LC
       ELSE
          EL = LCPLF
       ENDIF
C
C----------------------------------------------------------------------
C CALCULATE PARCEL MASSFLUX DIVIDED BY THE THICKNESS OF LAYER K
C THIS VALUE IS USED IN SEVERAL PLACES IN THE SUBROUTINE
C----------------------------------------------------------------------
C
       TEMPRY = FLXK(I)/DELPK(I)
C
       IF (BLOWST(I)) THEN
CL
CL----------------------------------------------------------------------
CL AT THE LOWEST CONVECTIVE LAYER, THE PARCEL MASS FLUX IS A FLUX FROM
CL THE ENVIRONMENT. IE. THE INITIAL MASS FLUX IS ENTRAINED WITH EXCESS
CL POTENTIAL TEMPERATURE AND MIXING RATIO TPIXS, QPIXS
CL
CL UM DOCUMENTATIO PAPER P27
CL SECTION (10), EQUATION (39)
CL----------------------------------------------------------------------
CL
         DTHEK(I) = DTHEK(I) - TEMPRY*THPIXS
c         if(bland(i))dthek(i)=dthek(i) - tempry*thpixs !Doubling for land
         DQEK(I) = DQEK(I) - TEMPRY*QPIXS
       ENDIF
CL
CL---------------------------------------------------------------------
CL EFFECT OF CONVECTION UPON POTENTIAL TEMPERATURE OF LAYER K
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (10), EQUATION (38A)
CL--------------------------------------------------------------------
CL
CH     DTHEK(I) = DTHEK(I) + TEMPRY * (
CH   &
CH   &           (1+EKP14(I)) * (1.0-DELTAK(I)) *        ! COMPENSATING
CH   &           (1-AMDETK(I)) * (THEKP1(I)-THEK(I))     ! SUBSIDENCE
CH   &         +
CH   &           DELTAK(I) * (1.0-AMDETK(I)) *           ! FORCED
CH   &           (THRK(I)-THEK(I)-                       ! DETRAINMENT
CH   &                    ((EL/CP)*XPK(I)/EXK(I)))
CH   &         +
CH   &           AMDETK(I) * (THPK(I)-THEK(I)-           ! MIXING
CH   &                    ((EL/CP)*XPK(I)/EXK(I)))       ! DETRAINMENT
CH   &         )
       DTHEK(I) = DTHEK(I) + TEMPRY * (
     &
     &           (1+EKP14(I)) * (1.0-DELTAK(I)) *        ! COMPENSATING
     &           (1-AMDETK(I)) * (THEKP1(I)-THEK(I))     ! SUBSIDENCE
     &         +
     &           DELTAK(I) * (1.0-AMDETK(I)) *           ! FORCED
     &           (THRK(I)-THEK(I)-                       ! DETRAINMENT
     &            FRACENV*((EL/CP)*XPK(I)/EXK(I)))
     &         +
     &           AMDETK(I) * (THPK(I)-THEK(I)-           ! MIXING
     &            FRACENV*((EL/CP)*XPK(I)/EXK(I)))       ! DETRAINMENT
     &         )
CL
CL---------------------------------------------------------------------
CL EFFECT OF CONVECTION UPON MIXING RATIO OF LAYER K
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (10), EQUATION (38B)
CL--------------------------------------------------------------------
CL
CH     DQEK(I) = DQEK(I) + TEMPRY * (
CH   &
CH   &           (1+EKP14(I)) * (1.0-DELTAK(I)) *        ! COMPENSATING
CH   &           (1-AMDETK(I)) * (QEKP1(I)-QEK(I))       ! SUBSIDENCE
CH   &         +
CH   &           DELTAK(I) * (1.0-AMDETK(I)) *           ! FORCED
CH   &           (QRK(I)-QEK(I)+XPK(I))                  ! DETRAINMENT
CH   &         +
CH   &           AMDETK(I) * (QPK(I)-QEK(I)+             ! MIXING
CH   &                                XPK(I))            ! DETRAINMENT
CH   &         )
       DQEK(I) = DQEK(I) + TEMPRY * (
     &
     &           (1+EKP14(I)) * (1.0-DELTAK(I)) *        ! COMPENSATING
     &           (1-AMDETK(I)) * (QEKP1(I)-QEK(I))       ! SUBSIDENCE
     &         +
     &           DELTAK(I) * (1.0-AMDETK(I)) *           ! FORCED
     &           (QRK(I)-QEK(I)+FRACENV*XPK(I))          ! DETRAINMENT
     &         +
     &           AMDETK(I) * (QPK(I)-QEK(I)+             ! MIXING
     &                        FRACENV*XPK(I))            ! DETRAINMENT
     &         )
CL
CL----------------------------------------------------------------------
CL TERMINAL DETRAINMENT AND SUBSIDENCE IN TERMINAL LAYER
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (10), EQUATION (40)
CL--------------------------------------------------------------------
CL
       IF ( BTERM(I) ) THEN
          TEMPRY = FLXKP1(I)/DELPKP1(I)
CHBG  USE UKMO SPLIT FINAL DETRAINMENT
          DTHEKP1(I) = DTHEKP1(I) + TEMPRY*((THPKP1(I)-THEKP1(I))
     &                                   - EL*XPKP1(I)/(EXKP1(I)*CP))
          DQEKP1(I)  = DQEKP1(I) + TEMPRY*(QPKP1(I)-QEKP1(I)
     &                                                + XPKP1(I))
CHBG OR USE SPLIT FINAL DETRAINMENT ONLY AT LOW LEVELS (FADE OUT OVER 
CHBG   700-500MB). PK IS LEVEL PRESSURE IN Pa = 100xMbs.

CHBG      XXDTHE = TEMPRY*((THPKP1(I)-THEKP1(I))
CHBG &                                   - EL*XPKP1(I)/(EXKP1(I)*CP))
CHBG      XXDQE  = TEMPRY*(QPKP1(I)-QEKP1(I) + XPKP1(I))
CHBG      PKP1FAC=MIN(1.0,MAX(0.0,(PK(I)-50000.)/(70000.-50000.)))
CHBG      DTHEKP1(I) = DTHEKP1(I) + PKP1FAC * XXDTHE
CHBG      DQEKP1(I)  = DQEKP1(I)  + PKP1FAC * XXDQE
CHBG      PKFAC = (1.0-PKP1FAC) * DELPKP1(I)/DELPK(I)
CHBG      DTHEK(I)   = DTHEK(I)   + PKFAC * XXDTHE
CHBG      DQEK(I)    = DQEK(I)    + PKFAC * XXDQE
       ENDIF

  10   CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE EVP----------------------------------------------------
CLL
CLL  PURPOSE : CALCULATES THE EVAPORATION OF PRECIPITATION
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 23/7/92
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE EVP(NPNTS,PRECIP,TEVP,CCA,RHO,DELQ,DELPKM1,EVAP,
     &               BEVAP,IPHASE,AREA_FAC)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS USED IN THIS SUBROUTINE
C-----------------------------------------------------------------------
C
CHBG  include 'C_G.f'
CHBG  include 'DDEVAP.f'
CHBG  include 'DDEVPLQ.f'
CHBG  include 'DDEVPICE.f'
      include 'UKPARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER I                 ! LOOP COUNTER
C
      INTEGER NPNTS             ! IN VECTOR LENGTH
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
C
      REAL DELQ(ln2)            ! IN CHANGE IN HUMIDITY MIXING
                                !    RATIO ACROSS LAYER K (KG/KG)
C
      REAL TEVP(ln2)            ! IN TEMPERATURE OF LAYER K (K)
C
      LOGICAL BEVAP(ln2)        ! IN MASK FOR POINTS WHERE EVAPORATION
                                !    TAKES PLACE
C
      REAL PRECIP(ln2)          ! IN AMOUNT OF PRECIPITATION(KG/M**2/S)
C
      REAL DELPKM1(ln2)         ! IN CHANGE IN PRESSURE ACROSS
                                !    LAYER K-1
C
      REAL CCA(ln2)             ! IN CONVECTIVE CLOUD AMOUNT
C
      REAL RHO(ln2)             ! IN DENSITY OF AIR
C
      INTEGER IPHASE            ! IN INDICATION FOR RAIN (1), OR
                                !    SNOW (2)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL EVAP(ln2)     ! OUT EVAPORATION
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES
C-----------------------------------------------------------------------
C

C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE LOCALLY DEFINED
C-----------------------------------------------------------------------
C
      REAL ECON          ! QUADRATIC TERM
C
      REAL C1            ! CONSTANT
C
      REAL C2            ! CONSTANT
C
      REAL SR_RHO        ! SQUARE ROOT OF DENSITY
C
      REAL LRATE         ! LOCAL RATE OF PRECIPITATION
C
      REAL CA            ! LOCAL CLOUD AREA
C
      REAL AREA_FAC      ! FRACTION OF CONVECTIVE CLOUD AMOUNT TO GIVE
                         ! LOCAL CLOUD AREA
C
C-----------------------------------------------------------------------
C START OF ROUTINE
C-----------------------------------------------------------------------
C
      IF (IPHASE.EQ.1) THEN        ! RAIN
C
      DO I=1,NPNTS
       IF (BEVAP(I)) THEN
         EVAP(I) = 0.0
         IF (PRECIP(I) .GT. 0.0) THEN
           ECON = ((LQ_A*TEVP(I)+LQ_B)*TEVP(I)+LQ_C)
           CA = AREA_FAC*CCA(I)
           LRATE = PRECIP(I)/CA
           SR_RHO = SQRT(RHO(I))
           C1 = RHO_LQA*CA*(LRATE*SR_RHO)**P_LQ1
           C2 = RHO_LQB*CA*(LRATE**P_LQ2)*(RHO(I)**RHO_LQP2)
           EVAP(I) = MIN(ECON*(C1+C2)*DELQ(I)*DELPKM1(I)/G,LRATE)
         END IF
       END IF
      END DO
C
      ELSE IF (IPHASE.EQ.2) THEN        ! SNOW
C
      DO I=1,NPNTS
       IF (BEVAP(I)) THEN
         EVAP(I) = 0.0
         IF (PRECIP(I) .GT. 0.0) THEN
           ECON = ((ICE_A*TEVP(I)+ICE_B)*TEVP(I)+ICE_C)
           CA = AREA_FAC*CCA(I)
           LRATE = PRECIP(I)/CA
           SR_RHO = SQRT(RHO(I))
           C1 = RHO_ICEA*CA*(LRATE*SR_RHO)**P_ICE1
           C2 = RHO_ICEB*CA*(LRATE**P_ICE2)*(RHO(I)**RHO_ICP2)
           EVAP(I) = MIN(ECON*(C1+C2)*DELQ(I)*DELPKM1(I)/G,LRATE)
         END IF
       END IF
      END DO
C
      ENDIF
C
      RETURN
      END
C
CLL  SUBROUTINE FLAG_WET-----------------------------------------------
CLL
CLL  PURPOSE : CALCULATES A MASK FOR WHEN CONDENSATION IS LIQUID
CLL
CLL            IF 0.5 * (TK + TK+1) > TICE THEN ANY CONDENSATION
CLL                                        IN LAYER K+1 IS LIQUID
CLL
CLL            IF 0.5 * (TK + TK+1) < TICE THEN ANY CONDENSATION
CLL                                        IN LAYER K+1 IS ICE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL                  SECTION (2B)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE FLAG_WET (BWATER,TH,EXNER,PSTAR,AKH,BKH,
     &                     NPNTS,NLEV)
C
C-----------------------------------------------------------------------
C   RETURNS 'BWATER' - A BIT VECTOR OF POINTS WHERE CONDENSATE IS WATER
C   RATHER THAN ICE.
C----------------------------------------------- AUTHOR: M FISHER 1987 -
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include 'TICE.f'
CHBG  include 'C_R_CP.f'
      include 'UKPARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS              ! IN VECTOR LENGTH
C
      INTEGER NLEV               ! IN NUMBER OF MODEL LAYERS
C
      INTEGER I,K                ! LOOP COUNTERS
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL TH(ln2,NLEV)             ! IN POTENTIAL TEMPERATURE (K)
C
      REAL EXNER(ln2,NLEV+1)        ! IN EXNER RATIO AT LAYER
                                    ! BOUNDARIES (STARTING WITH THE
                                    ! SURFACE)
C
      REAL PSTAR(ln2)               ! IN Surface pressure
C
      REAL AKH(NLEV+1)              ! IN Hybrid coordinate A at
                                    !    layer boundary
      REAL BKH(NLEV+1)              ! IN Hybrid coordinate B at
                                    !    layer boundary
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      LOGICAL BWATER(ln2,2:NLEV)    ! OUT MASK FOR THOSE POINTS AT
                                    !     WHICH CONDENSATE IS LIQUID
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL EXK                      ! EXNER RATIO FOR LEVEL K
      REAL EXKP1                    ! EXNER RATIO FOR LEVEL K+1
C

      REAL
     &    PU,PL,PU2
      include 'P_EXNERC.f'

C*---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL  NO SIGNIFICANT STRUCTURE
CL---------------------------------------------------------------------
CL
      DO 10 K=1,NLEV-1
       DO 10 I=1,NPNTS
C
        PU2=PSTAR(I)*BKH(K+2) + AKH(K+2)
        PU=PSTAR(I)*BKH(K+1) + AKH(K+1)
        PL=PSTAR(I)*BKH(K) + AKH(K)
        EXK = P_EXNER_C(EXNER(I,K+1),EXNER(I,K),PU,PL,KAPPA)
        EXKP1 = P_EXNER_C(EXNER(I,K+2),EXNER(I,K+1),PU2,PU,KAPPA)
C
        BWATER(I,K+1) = 0.5*(TH(I,K)*EXK + TH(I,K+1)*EXKP1) .GT. TICE
   10 CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE FLX_INIT-----------------------------------------------
CLL
CLL  PURPOSE : CALCULATE INITIAL DOWNDRAUGHT MASSFLUX
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY SUMMER 1992
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL   3.3   23/12/93 : DG020893 : DUE TO CHANGE IN WAY CLOUD TOP IS
CLL                               ESTIMATED BECAUSE OF CHANGES TO THE
CLL                               CALCULATION OF FORCED DETRAINMENT
CLL
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE FLX_INIT (NPNTS,KCT,ICCB,ICCT,FLX,FLX_DD_K,BDDI,
     &                     FLX_STRT)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER I                 ! LOOP COUNTER
C
      INTEGER NPNTS             ! IN NUMBER OF POINTS
C
      INTEGER KCT               ! IN CONVECTIVE CLOUD TOP
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      INTEGER ICCB(ln2)         ! IN CONVECTIVE CLOUD BASE
C
      INTEGER ICCT(ln2)         ! IN CONVECTIVE CLOUD TOP
C
      REAL FLX(ln2,KCT+1)       ! IN CONVECTIVE MASSFLUX (PA/S)
C
      LOGICAL BDDI(ln2)         ! IN MASK FOR THOSE POINTS WHERE
                                !    DOWNDRAUGHT MAY INITIATE
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL FLX_DD_K(ln2)        ! OUT DOWNDRAUGHT MASSFLUX OF LAYER K
                                !     (PA/S)
C
      REAL FLX_STRT(ln2)        ! OUT UPDRAUGHT MASSFLUX AT LEVEL
                                !     DOWNDRAUGHT STARTS (PA/S)
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      INTEGER KDDREF            ! REFERENCE LEVEL FOR DOWNDRAUGHT
                                ! MASSFLUX
C
C----------------------------------------------------------------------
C CALCULATE DOWNDRAUGHT MASSFLUX BASED ON A REFERENCE LEVEL WHICH IS
C 3/4 CLOUD DEPTH
C----------------------------------------------------------------------
C
      DO I=1,NPNTS
       IF (BDDI(I)) THEN
          KDDREF = INT(ICCB(I) + 0.75*(ICCT(I) - ICCB(I)))
          IF (KDDREF .GE. ICCT(I)-1) KDDREF=ICCT(I)-1
          FLX_STRT(I) = FLX(I,KDDREF)
          FLX_DD_K(I) = FLX_STRT(I) * 0.05
       END IF
      END DO
C
      RETURN
      END
C
CLL  SUBROUTINE LATENT_H-----------------------------------------------
CLL
CLL  PURPOSE : CALCULATES A NEW PARCEL TEMPERATURE AFTER
CLL            CONDENSATION HAS OCCURRED
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE LATENT_H (NPNTS,THPKP1,QPKP1,THEKP1,
     &                     QSEKP1,DQSKP1,BGMKP1,BWKP1,EXKP1)
C
C-----------------------------------------------------------------------
C   ADJUSTS PARCEL POTENTIAL TEMPERATURES TO ACCOUNT FOR THE LATENT
C   HEAT RELEASE FROM SUPERSATURATED PARCELS CONDENSING/DEPOSITING OUT
C----------------------------------------------- AUTHOR: M FISHER 1987 -
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include 'C_R_CP.f'
CHBG  include 'C_LHEAT.f'
      include 'UKPARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS            ! VECTOR LENGTHS
C
      INTEGER I                ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL QPKP1(ln2)          ! IN PARCEL MIXING RATIO IN LAYER K+1
                               !    (KG/KG/S)
C
      REAL THEKP1(ln2)         ! IN ENVIRONMENT POTENTIAL FOR LAYER K+1
                               !    (K/S)
C
      REAL QSEKP1(ln2)         ! IN SATURATION MIXING RATIO OF THE
                               !    ENVIRONMENT IN LAYER K+1 (KG/KG/S)
C
      REAL DQSKP1(ln2)         ! IN GRADIENT OF SATURATION MIXING RATIO
                               !    WITH POTENTIAL TEMPERATURE FOR THE
                               !    ENVIRONMENT IN LAYER K+1 (KG/KG/K)
C
      LOGICAL BGMKP1(ln2)      ! IN MASK FOR PARCELS WHICH ARE
                               !    SATURATED LAYER K+1
C
      LOGICAL BWKP1(ln2)       ! IN MASK FOR POINTS IN WHICH
                               ! CONDENSATE IS LIQUID IN LAYER K+1
C
      REAL EXKP1(ln2)          ! IN EXNER RATIO AT MID-POINT OF
                               !    LAYER K+1
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      REAL THPKP1(ln2)         ! INOUT
                               ! IN  INITIAL ESTIMATE OF PARCEL
                               !     POTENTIAL TEMPERATURE IN
                               !     LAYER K+1 (K/S)
                               ! OUT PARCEL POTENTIAL TEMPERATURE
                               !     IN LAYER K+1 AFTER LATENT
                               !     HEATING (K/S)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL EL                  ! LATENT HEAT OF CONDENSATION OR
                               ! (CONDENSATION + FUSION) (J/KG)
C
C*----------------------------------------------------------------------
CL
CL----------------------------------------------------------------------
CL  ADJUST PARCEL POTENTIAL TEMPERATURES TO ACCOUNT FOR LATENT HEATING
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (4), EQUATION(21)
CL----------------------------------------------------------------------
CL
      DO 10 I=1,NPNTS
       IF (BWKP1(I)) THEN
          EL = LC
       ELSE
          EL = LCPLF
       ENDIF
C
       IF (BGMKP1(I)) THEN
         THPKP1(I) = ( THPKP1(I) +
     &   (EL/(EXKP1(I)*CP))*(QPKP1(I)-QSEKP1(I)+THEKP1(I)*DQSKP1(I))
     &              ) / (1.+(EL/(EXKP1(I)*CP))*DQSKP1(I))
       ENDIF
   10  CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE LAYER_CN-----------------------------------------------
CLL
CLL  PURPOSE : CALCULATES LAYER DEPENDENT CONSTANTS FOR LAYER K
CLL            -PRESSURE
CLL            -LAYER THICKNESS
CLL            -ENTRAINMENT COEFFICIENTS
CLL            -DETRAINMENT COEFFICIENTS
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL C.W. , D.G. <- PROGRAMMER OF SOME OR ALL OF PREVIOUS CODE OR CHANGES
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED:P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE LAYER_CN(K,NPNTS,NLEV,EXNER,AK,BK,AKM12,BKM12,
     &                    DELAK,DELBK,PSTAR,PK,PKP1,DELPK,DELPKP1,
     &                    DELPKP12,EKP14,EKP34,AMDETK,EXK,EXKP1,
     &                    DELEXKP1)
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include 'ENTCNST.f'
CHBG  include 'C_R_CP.f'
      include 'UKPARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTER
C----------------------------------------------------------------------
C
      INTEGER NPNTS             ! IN VECTOR LENGTH
C
      INTEGER NLEV              ! IN NUMBER OF MODEL LEVELS
C
      INTEGER K                 ! IN PRESENT MODEL LAYER
C
      INTEGER I                 ! COUNTER FOR DO LOOPS
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL AK(NLEV)             ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BK(NLEV)             ! IN ) MID-LAYER OF LAYER K
C
      REAL AKM12(NLEV+1)        ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BKM12(NLEV+1)        ! IN ) LOWER LAYER BOUNDARY OF LAYER K
C
      REAL DELAK(NLEV)          ! IN ) HYBRID CO-ORDINATE VALUES FOR
      REAL DELBK(NLEV)          ! IN ) FOR THICKNESS OF LAYER K
C
      REAL PSTAR(ln2)           ! IN SURFACE PRESSURE (PA)
C
      REAL EXNER(ln2,NLEV+1)    ! IN EXNER FUNCTION AT LAYER
                                  ! BOUNDARIES STARTING AT LEVEL K-1/2
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL PK(ln2)              ! OUT PRESSURE AT LAYER K (PA)
C
      REAL PKP1(ln2)            ! OUT PRESSURE AT LAYER K+1 (PA)
C
      REAL DELPK(ln2)           ! OUT THICKNESS OF LAYER K (PA)
C
      REAL DELPKP1(ln2)         ! OUT THICHNESS OF LAYER K+1 (PA)
C
      REAL DELPKP12(ln2)        ! OUT THICKNESS BETWEEN LAYER K AND K+1
                                !     (PA)
C
      REAL EKP14(ln2)           ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K+1/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL EKP34(ln2)           ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K+3/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL AMDETK(ln2)          ! OUT MIXING DETRAINMENT COEFFICIENT
                                !     AT LEVEL K MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL EXK(ln2)             ! EXNER FUNCTION AT LEVEL K
C
      REAL EXKP1(ln2)           ! EXNER FUNCTION AT LEVEL K+1
C
      REAL DELEXKP1(ln2)        ! DIFFERENCE IN EXNER FUNCTION
                                ! BETWEEN K+3/2 AND K+1/2
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL AEKP14,AEKP34        ! USED IN CALCULATION OF ENTRAINMENT
                                ! RATE
C

      REAL
     &    PU,PL
      include 'P_EXNERC.f'

C*---------------------------------------------------------------------
C
C----------------------------------------------------------------------
C SET CONSTANT AE USED IN CALCULATION OF ENTARINMENT AND
C DETRAINMENT RATES DEPENDING UPON LEVEL
C----------------------------------------------------------------------
C
      IF(K.EQ.1)THEN
        AEKP14 = AE1
        AEKP34 = AE2
      ELSE
        AEKP14 = AE2
        AEKP34 = AE2
      END IF
C
      DO 10 I=1,NPNTS
CL
CL---------------------------------------------------------------------
CL CALCULATE PK AND DELPK - IF K = 1 (LOWEST MODEL LAYER) THEN
CL VALUES FOR PREVIOUS PASS THROUGH ROUTINE AT (K-1)+1 ARE TAKEN
CL---------------------------------------------------------------------
CL
        IF(K.EQ.1)THEN
          PK(I) = AK(K) + BK(K)*PSTAR(I)
          DELPK(I) = -DELAK(K) - DELBK(K)*PSTAR(I)
        ELSE
          PK(I) = PKP1(I)
          DELPK(I) = DELPKP1(I)
        END IF
CL
CL---------------------------------------------------------------------
CL CALCULATE PKP1, DELPKP1 AND DELPK+1/2
CL---------------------------------------------------------------------
CL
        PKP1(I) = AK(K+1) + BK(K+1)*PSTAR(I)
        DELPKP1(I) = -DELAK(K+1) - DELBK(K+1)*PSTAR(I)
        DELPKP12(I) = PK(I) - PKP1(I)
CL
CL---------------------------------------------------------------------
CL CALCULATE EXNER FUNCTIONS AT MID-LAYES K AND K+1, AND
CL DIFFERENCE OF EXNER FUNCTION ACROSS LAYER K
CL---------------------------------------------------------------------
CL
        IF(K.EQ.1)THEN
          PU=PSTAR(I)*BKM12(K+1) + AKM12(K+1)
          PL=PSTAR(I)*BKM12(K) + AKM12(K)
          EXK(I) = P_EXNER_C(EXNER(I,K+1),EXNER(I,K),PU,PL,KAPPA)
        ELSE
          EXK(I) = EXKP1(I)
        END IF
        PU=PSTAR(I)*BKM12(K+2) + AKM12(K+2)
        PL=PSTAR(I)*BKM12(K+1) + AKM12(K+1)
        EXKP1(I) = P_EXNER_C(EXNER(I,K+2),EXNER(I,K+1),PU,PL,KAPPA)
        DELEXKP1(I) = EXNER(I,K+1)-EXNER(I,K+2)
CL
CL---------------------------------------------------------------------
CL CALCULATE ENTRAINMENT AND MIXING DETRAINMENT COEFFICIENTS
CL---------------------------------------------------------------------
CL
CL
CL---------------------------------------------------------------------
CL CALCULATE ENTRAINMENT COEFFICIENTS MULTIPLIED BY
CL APPROPRIATE LAYER THICKNESS
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (2C), EQUATION(14)
CL---------------------------------------------------------------------
CL
        EKP14(I) = ENTCOEF * AEKP14 * PK(I) *
     &             (PK(I) - AKM12(K+1) - BKM12(K+1)*PSTAR(I)) /
     &             (PSTAR(I) * PSTAR(I))
        EKP34(I) = ENTCOEF * AEKP34 * (AKM12(K+1)+BKM12(K+1)*PSTAR(I)) *
     &             (AKM12(K+1) + BKM12(K+1)*PSTAR(I) - PKP1(I)) /
     &             (PSTAR(I) * PSTAR(I))
CL
CL---------------------------------------------------------------------
CL CALCULATE MIXING DETRAINMENT COEFFICIENT MULTIPLIED BY
CL APPROPRIATE LAYER THICKNESS
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (2C), EQUATION(15)
CL---------------------------------------------------------------------
CL
        IF(K.EQ.1)THEN
          AMDETK(I) = 0.0
        ELSE
          AMDETK(I) = (EKP14(I) + EKP34(I)) * (1.0-1.0/AEKP34)
        END IF
 10   CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE LAYER_DD--------------------------------------------
CLL
CLL  PURPOSE : CALCULATES LAYER DEPENDENT CONSTANTS FOR LAYER K
CLL            -PRESSURE
CLL            -LAYER THICKNESS
CLL            -ENTRAINMENT COEFFICIENTS
CLL            -DETRAINMENT COEFFICIENTS
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT & D.GREGORY SUMMER 1992
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  3.3   23/12/93 : DG060893 : CORRECTION TO PREVENT OVER PREDICTION
CLL                              OF SNOW SHOWERS; REARRANGEMENT OF
CLL                              ENTRAINMENT RATES
CLL
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE LAYER_DD(NPNTS,K,KCT,THE_K,THE_KM1,FLX_STRT,AK,
     &                    BK,AKM12,BKM12,DELAK,DELBK,EXNER_KM12,
     &                    EXNER_KP12,EXNER_KM32,PSTAR,PK,PKM1,DELPK,
     &                    DELPKM1,EXK,EXKM1,AMDETK,EKM14,EKM34,KMIN,
     &                    BDDI)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include 'C_0_DG_C.f'
CHBG  include 'C_R_CP.f'
CHBG  include 'ENTCNST.f'
CHBG  include 'ENTDD.f'
CHBG  include 'DDKMDET.f'
      include 'UKPARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTER
C----------------------------------------------------------------------
C
      INTEGER NPNTS             ! IN VECTOR LENGTH
C
      INTEGER K                 ! IN PRESENT MODEL LAYER
C
      INTEGER I                 ! COUNTER FOR DO LOOPS
C
      INTEGER KCT               ! IN CONVECTIVE CLOUD TOP LAYER
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL AK(K)                ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BK(K)                ! IN ) MID-LAYER OF LAYER K
C
      REAL AKM12(K+1)           ! IN ) HYBRID CO-ORDINATE VALUES AT
      REAL BKM12(K+1)           ! IN ) LOWER LAYER BOUNDARY OF LAYER K
C
      REAL DELAK(K)             ! IN ) HYBRID CO-ORDINATE VALUES FOR
      REAL DELBK(K)             ! IN ) FOR THICKNESS OF LAYER K
C
      REAL PSTAR(ln2)           ! IN SURFACE PRESSURE (PA)
C
      REAL EXNER_KM12(ln2)      ! IN EXNER FUNCTION AT LAYER K-1/2
C
      REAL EXNER_KP12(ln2)      ! IN EXNER FUNCTION AT LAYER K+1/2
C
      REAL EXNER_KM32(ln2)      ! IN EXNER FUNCTION AT LAYER K-3/2
C
      REAL FLX_STRT(ln2)        ! IN UPDRAUGHT MASSFLUX AT LEVEL WHERE
                                !    DOWNDRAUGHT STARTS (PA/S)
C
      REAL THE_K(ln2)           ! IN POTENTIAL TEMPERATURE OF
                                !    ENVIRONMENT IN LAYER K (K)
C
      REAL THE_KM1(ln2)         ! IN POTENTIAL TEMPERATURE OF
                                !    ENVIRONMENT IN LAYER K-1 (K)
C
      LOGICAL BDDI(ln2)         ! IN MASK FOR POINTS WHERE DOWNDRAUGHT
                                !    MAY INITIATE
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C----------------------------------------------------------------------
C
      INTEGER KMIN(ln2)         ! INOUT
                                ! FREEZING LEVEL
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL PK(ln2)              ! OUT PRESSURE AT LAYER K (PA)
C
      REAL PKM1(ln2)            ! OUT PRESSURE AT LAYER K-1 (PA)
C
      REAL DELPK(ln2)           ! OUT THICKNESS OF LAYER K (PA)
C
      REAL DELPKM1(ln2)         ! OUT THICHNESS OF LAYER K-1 (PA)
C
      REAL EKM14(ln2)           ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K-1/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL EKM34(ln2)           ! OUT ENTRAINMENT COEFFICIENT AT
                                !     LEVEL K-3/4 MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL AMDETK(ln2)          ! OUT MIXING DETRAINMENT COEFFICIENT
                                !     AT LEVEL K MULTIPLIED BY
                                !     APPROPRIATE LAYER THICKNESS
C
      REAL EXK(ln2)             ! OUT EXNER FUNCTION AT LEVEL K
C
      REAL EXKM1(ln2)           ! OUT EXNER FUNCTION AT LEVEL K-1
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL TTK                  ! TEMPERATURE STORE AT LAYER K
C
      REAL TTKM1                ! TEMPERATURE STORE AT LAYER K-1
C
      REAL THKM12               ! POTENTIAL TEMPERATURE STORE AT
                                ! LAYER K-1/2
C
      REAL TTKM12               ! TEMPERATURE STORE AT LAYER K-1/2
C
      REAL INCR_FAC             ! INCREMENT FACTOR FOR ENTRAINMENT
                                ! RATES AT FREEZING LEVEL
C
      REAL
     &    PU,PL
      include 'P_EXNERC.f'
 
C----------------------------------------------------------------------
C SET KMIN TO INITIAL VALUE
CL CALCULATE PK, DELPK AND EXNER FUNCTION - IF K = KCT THEN
CL VALUES FOR PREVIOUS PASS THROUGH ROUTINE AT (K-1)+1 ARE TAKEN
C----------------------------------------------------------------------
C
      IF (K.EQ.KCT+1) THEN
       DO I=1,NPNTS
        KMIN(I) = KCT+2
        PK(I) = AK(K) + BK(K)*PSTAR(I)
        DELPK(I) = - DELAK(K) - DELBK(K)*PSTAR(I)
        PU=PSTAR(I)*BKM12(K+1) + AKM12(K+1)
        PL=PSTAR(I)*BKM12(K) + AKM12(K)
        EXK(I) = P_EXNER_C(EXNER_KP12(I),EXNER_KM12(I),PU,PL,KAPPA)
       END DO
      ELSE
       DO I=1,NPNTS
        PK(I) = PKM1(I)
        DELPK(I) = DELPKM1(I)
        EXK(I) = EXKM1(I)
       END DO
      END IF
CL
CL---------------------------------------------------------------------
CL CALCULATE PKM1, DELPKM1
CL CALCULATE EXNER FUNCTIONS AT MID-LAYES K AND K-1, AND
CL DIFFERENCE OF EXNER FUNCTION ACROSS LAYER K
CL---------------------------------------------------------------------
CL
      DO I=1,NPNTS
        PKM1(I) = AK(K-1) + BK(K-1)*PSTAR(I)
        DELPKM1(I) = - DELAK(K-1) - DELBK(K-1)*PSTAR(I)
        PU=PSTAR(I)*BKM12(K) + AKM12(K)
        PL=PSTAR(I)*BKM12(K-1) + AKM12(K-1)
        EXKM1(I) = P_EXNER_C(EXNER_KM12(I),EXNER_KM32(I),PU,PL,KAPPA)
C
CL
CL---------------------------------------------------------------------
CL CALCULATE FREEZING LEVEL : CHECK IF FREEZING LEVEL IN THIS LAYER
CL---------------------------------------------------------------------
CL
       IF (KMIN(I).EQ.KCT+2) THEN
        TTK = THE_K(I)*EXK(I)
        TTKM1 = THE_KM1(I)*EXKM1(I)
        THKM12 = (THE_KM1(I)+THE_K(I))*0.5
        TTKM12 = THKM12*EXNER_KM12(I)
        IF (TTKM12 .GE. TM .AND. TTK .LT. TM) KMIN(I) = K
        IF (TTKM12 .LT. TM .AND. TTKM1 .GE. TM) KMIN(I) = K-1
       END IF
C
CL
CL---------------------------------------------------------------------
CL CALCULATE ENTRAINMENT COEFFICIENTS MULTIPLIED BY
CL APPROPRIATE LAYER THICKNESS
CL
CL CALCULATE MIXING DETRAINMENT COEFFICIENT MULTIPLIED BY
CL APPROPRIATE LAYER THICKNESS
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (2C), EQUATION(14)
CL---------------------------------------------------------------------
CL
      IF (PK(I).LT.PSTAR(I)-DET_LYR) THEN
       EKM14(I) = AE2 * (AKM12(K)+BKM12(K)*PSTAR(I)-PK(I)) / PSTAR(I)
       EKM34(I) = AE2 * (PKM1(I)-AKM12(K)-BKM12(K)*PSTAR(I)) / PSTAR(I)
       AMDETK(I) = (EKM14(I)+EKM34(I)) * (1.0-1.0/AE2)
      ELSE
       EKM14(I) = 0.0
       EKM34(I) = 0.0
       AMDETK(I) = DELPK(I) / (PSTAR(I)*(1.0-BKM12(K+1))-AKM12(K+1))
      END IF
C
      IF (BDDI(I)) THEN
C
      IF (K.EQ.KMIN(I) .AND. PK(I).LT.PSTAR(I)-DET_LYR) THEN
        INCR_FAC = FLX_STRT(I)*DDCOEF1/PSTAR(I)
        IF (INCR_FAC.GT.6.0) INCR_FAC=6.0
        EKM14(I) = EKM14(I)*INCR_FAC
        EKM34(I) = EKM34(I)*INCR_FAC
      ELSE
        EKM14(I) = EKM14(I)*DDCOEF2
        EKM34(I) = EKM34(I)*DDCOEF2
        IF (KMIN(I).NE.KCT+2 .AND. K.LT.KMIN(I) .AND. PK(I).LT.
     & PSTAR(I)-DET_LYR)  AMDETK(I) = AMDETK(I)*DDCOEF2
      END IF
C
      END IF
      END DO
C
      RETURN
      END
C
CLL  SUBROUTINE LIFT_PAR-----------------------------------------------
CLL
CLL  PURPOSE : LIFTS THE PARCEL FROM LAYER K TO K+1
CLL            TAKING ENTRAINEMNT AND MOIST PROCESSES INTO ACOUNT
CLL
CLL            SUBROUTINE LATENT_H CALCULATES THE MOIST PROCESSES
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER NO. ##
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE LIFT_PAR (NPNTS,THPKP1,QPKP1,XSQKP1,BGMKP1,BWKP1,
     &                     THPK,QPK,THEKP1,QEKP1,THEK,QEK,QSEKP1,
     &                     DQSKP1,PKP1,EXKP1,EKP14,EKP34)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER I              ! LOOP COUNTER
C
C----------------------------------------------------------------------
C VARAIBLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL THEK(ln2)         ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K (K)
C
      REAL THEKP1(ln2)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEK(ln2)          ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL QEKP1(ln2)        ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL QSEKP1(ln2)       ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL DQSKP1(ln2)       ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT IN LAYER K+1
                             !    (KG/KG/K)
C
      REAL THPK(ln2)         ! IN PARCEL POTENTIAL TEMPERATURE IN
                             !    LAYER K (K)
C
      REAL QPK(ln2)          ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
C
      LOGICAL BWKP1(ln2)     ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      REAL PKP1(ln2)         ! IN PRESSURE AT LEVEL K+1 (PA)
C
      REAL EXKP1(ln2)        ! IN EXNER RATIO AT MID-POINT OF LAYER K+1
C
      REAL EKP14(ln2)        ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP34(ln2)        ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+3/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL THPKP1(ln2)       ! OUT PARCEL POTENTIAL TEMPERATURE IN
                             !     LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (K)
C
      REAL QPKP1(ln2)        ! OUT PARCEL MIXING RATIO IN LAYER K+1
                             !     AFTER ENTRAINMENT AND LATENT HEATING
                             !     (KG/KG)
C
      REAL XSQKP1(ln2)       ! OUT EXCESS PARCEL WATER AFTER
                             !     LIFTING FROM LAYER K TO K+1
                             !     (KG/KG)
C
      LOGICAL BGMKP1(ln2)    ! OUT MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1 AFTER
                             !     ENTRAINMENT AND LATENT HEATING
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C
C ON THE IBM ARRAYS ARE ALLOCATED USING A PARAMETER STATEMENT
C
C ON THE CRAY ARRAYS ARE DYNAMICALLY ALLOCATED
C----------------------------------------------------------------------
C
      REAL THPKP1T(ln2)      ! INITIAL ESTIMATE OF PARCEL TEMPERATURE
                             ! IN LAYER K+1 AFTER ENTRAINMENT (K)
C
      REAL TT(ln2)           ! TEMPORARY TEMPERATURE USED IN CALCULATION
                             ! OF SATURATION MIXING RATIO (K)
C
      REAL QSPKP1(ln2)       ! SATURATION MIXING RATIO OF PARCEL
                             ! AFTER DRY ASCENT (KG/KG)
C
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL QSATU,LATENT_H
C
C*---------------------------------------------------------------------
C
      DO 10 I=1,NPNTS
CL
CL----------------------------------------------------------------------
CL  LIFT PARCEL MIXING RATIO AND POTENTIAL TEMPERATURE TO THE NEXT LEVEL
CL----------------------------------------------------------------------
CL
CL----------------------------------------------------------------------
CL  INITIAL 'DRY' ASCENT
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (3), EQUATIONS (11B), (12B)
CL----------------------------------------------------------------------
CL
       THPKP1(I) = (  THPK(I)
     &             + EKP14(I)*THEK(I) + EKP34(I)*(1.+EKP14(I))*THEKP1(I)
     &             ) / ((1.+EKP14(I))*(1.+EKP34(I)))
C
       QPKP1(I) = (  QPK(I)
     &             + EKP14(I)*QEK(I) + EKP34(I)*(1.+EKP14(I))*QEKP1(I)
     &             ) / ((1.+EKP14(I))*(1.+EKP34(I)))
C
C-----------------------------------------------------------------------
C   CALCULATE WHERE THE PARCEL IS SUPERSATURATED (IE WHERE GAMMA(K+1)=1
C   SEE DCTN 29 PAGE 123)
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
C PRESSURE OF LAYER K FOR CALCULATION OF SATURATED
C MIXING RATIO
C-----------------------------------------------------------------------
C
       TT(I) = THPKP1(I)*EXKP1(I)
   10 CONTINUE
      CALL QSATU (QSPKP1,TT,PKP1,NPNTS)
C
      DO 20 I=1,NPNTS
       BGMKP1(I) = QPKP1(I) .GT. QSPKP1(I)
CL
CL----------------------------------------------------------------------
CL  CONDENSATION CALCULATION
CL
CL  SUBROUTINE LATENT_H
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (4)
CL----------------------------------------------------------------------
CL
       THPKP1T(I) = THPKP1(I)
   20 CONTINUE
C
      CALL LATENT_H (NPNTS,THPKP1T,QPKP1,THEKP1,QSEKP1,DQSKP1,
     &               BGMKP1,BWKP1,EXKP1)
C
C-----------------------------------------------------------------------
C   CALCULATE A MORE ACCURATE PARCEL SATURATED MIXING RATIO AND CONDENSE
C   OUT ANY EXCESS WATER VAPOUR. STORE THE EXCESS AMOUNTS IN 'XSQKP1'
C   FOR LATER.  SET PARCEL POTENTIAL TEMPERATURES TO THE PROVISIONAL
C   VALUES EXCEPT WHERE THE PARCEL IS NOT SUPERSATURATED WRT THE NEW
C   SATURATED MIXING RATIO. RECALCULATE BIT VECTOR 'BGMKP1'.
C-----------------------------------------------------------------------
C
C
C-----------------------------------------------------------------------
C CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
C PRESSURE OF LAYER K FOR CALCULATION OF SATURATED
C MIXING RATIO
C-----------------------------------------------------------------------
C
      DO 35 I = 1,NPNTS
       TT(I) = THPKP1T(I)*EXKP1(I)
   35 CONTINUE
      CALL QSATU (QSPKP1,TT,PKP1,NPNTS)
C
      DO 40 I=1,NPNTS
       XSQKP1(I) = QPKP1(I) - QSPKP1(I)
C
       IF(XSQKP1(I) .LE. 0.0) THEN
         BGMKP1(I) = .FALSE.
         XSQKP1(I) = 0.0
       ELSE
         BGMKP1(I) = .TRUE.
         THPKP1(I) = THPKP1T(I)
       END IF
C
       QPKP1(I) = QPKP1(I) - XSQKP1(I)
   40 CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE PARCEL-------------------------------------------------
CLL
CLL  PURPOSE : COMPLETES LIFTING OF THE PARCEL FROM LAYER K TO K+1
CLL
CLL            CALL SUBROUTINE DETRAIN, TERM_CON, CLOUD_W
CLL
CLL            AN INITIAL MASS FLUX IS CALCULATED
CLL
CLL            SUBROUTINE DETRAIN CARRIES OUT THE FORCED DETRAINMENT
CLL            CALCULATION
CLL
CLL            SUBROUTINE TERM_CON TESTS FOR ANY CONVECTION WHICH IS
CLL            TERMINATING IN LAYER K+1
CLL
CLL            SUBROUTINE CLOUD_W CARRIES OUT THE CLOUD MICROPHYSICS
CLL            CALCULATION
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL                                   THE RELEVANT POINTS IN CONVECT
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL  3.2  8/07/93 : added convective cloud condensed water diagnostic
CLL               : P Inness
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
      SUBROUTINE PARCEL (K,NPNTS,NLEV,PSTAR,THEKP1,THEK,QEKP1,QEK,
     &                   QSEK,QSEKP1,DQSK,DQSKP1,BLAND,BWKP1,
     &                   DELTAK,FLXK,THPK,QPK,THRK,QRK,
     &                   BTERM,THPKP1,QPKP1,PREKP1,XPK,XPKP1,FLXKP1,
     &                   XSQKP1,THPI,QPI,BGMK,BGMKP1,BLOWST,RBUOY,
     &                   CCA,ICCB,ICCT,TCW,DEPTH,
     &                   EKP14,EKP34,AMDETK,DELPKP1,PK,PKP1,
     &                   EXK,EXKP1,DELEXKP1,CCLWP,CCW,CCWI,CCWL,
     &                   qcloud,ipass,xs,fs,DELPK)
C
      IMPLICIT NONE
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
      include 'PARAMS.f'

CHBG  include 'XSBMIN.f'
CHBG  include 'MASSFC.f'
CHBG  include 'C_G.f'
      include 'UKPARAMS.f'

      integer ntrac
      parameter (ntrac=3)
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER NLEV           ! IN NUMBER OF MODEL LEVELS
C
      INTEGER NDET           ! COMPRESSED VECTOR LENGTH FOR
                             ! FORCED DETRAINMENT CALCULATION
C
      INTEGER K              ! IN PRESENT MODEL LAYER
C
      INTEGER I              ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C VARAIBLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL THEK(ln2)         ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K (K)
C
      REAL THEKP1(ln2)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEK(ln2)          ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL QEKP1(ln2)        ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL QSEKP1(ln2)       ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL DQSKP1(ln2)       ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT IN LAYER K+1
                             !    (KG/KG/K)
C
      REAL PSTAR(ln2)        ! IN SURFACE PRESSURE (PA)
C
      REAL THPK(ln2)         ! IN PARCEL POTENTIAL TEMPERATURE
                             !    IN LAYER K (KG/KG)
C
      REAL QPK(ln2)          ! IN PARCEL MIXING RATIO IN LAYER K (KG/KG)
C
      REAL XSQKP1(ln2)       ! IN EXCESS PARCEL WATER AFER LIFTING FROM
                             !    LAYER K TO K+1 (KG/KG)
C
      REAL RBUOY(ln2)        ! IN PARCEL BUOYANCY IN LAYER K+1 (K)
C
      REAL QSEK(ln2)         ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL DQSK(ln2)         ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    CLOUD ENVIRONMENT OF LAYER K
                             !    (KG/KG/K)
C
      REAL THPI(ln2)         ! IN INITIAL PARCEL POTENTIAL TEMPERATURE
                             !    (K)
C
      REAL QPI(ln2)          ! IN INITIAL PARCEL MIXING RATIO (KG/KG)
C
      REAL XPK(ln2)          ! IN PARCEL CLOUD WATER IN LAYER K (KG/KG)
C
      LOGICAL BWKP1(ln2)     ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      LOGICAL BGMK(ln2)      ! IN MASK FOR PARCELS WHICH ARE
                             !    SATURATED IN LAYER K
C
      LOGICAL BLAND(ln2)     ! IN LAND/SEA MASK
C
      LOGICAL BLOWST(ln2)    ! IN MASK FOR THOSE POINTS AT WHICH
                             !    STABILITY IS LOW ENOUGH FOR
                             !    CONVECTION TO OCCUR
C
      REAL EKP14(ln2)        ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !    K+1/4 MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP34(ln2)        ! IN ENTRAINMENT COEFFICIENT AT LEVEL
                             !   K+3/4 MULTIPLIED BY APPROPRIATE
                             !   LAYER THICKNESS
C
      REAL AMDETK(ln2)       ! IN MIXING DETRAINMENT COEFFICIENT
                             !    AT LEVEL K MULTIPLIED BY
                             !    APPROPORIATE LAYER THICKNESS
C
      REAL DELPKP1(ln2)      ! IN PRESSURE DIFFERENCE ACROSS
                             !    LAYER K+1 (PA)
C
      REAL DELPK(ln2)        ! IN PRESSURE DIFFERENCE ACROSS LAYER K
C
      REAL PK(ln2)           ! IN PRESSURE AT LEVEL K (PA)
C
      REAL PKP1(ln2)         ! IN PRESSURE AT LEVEL K+1 (PA)
C
      REAL EXK(ln2)          ! IN EXNER FUNCTION AT LEVEL K
C
      REAL EXKP1(ln2)        ! IN EXNER FUNCTION AT LEVEL K+1
C
      REAL DELEXKP1(ln2)     ! IN DIFFERENCE IN EXNER FUNCTION ACROSS
                             !    LAYER K+1
C
      logical qcloud         ! IN qcloud = T or F
      integer ipass          ! IN ipass  = 1 or 2 (2=leads)

      real xs(ln2)           ! IN sulfate mixing ratio
C
C
C---------------------------------------------------------------------
C VARAIBLES WHICH ARE BOTH INPUT AND OUTPUT
C---------------------------------------------------------------------
C
      REAL THPKP1(ln2)       ! INOUT
                             ! IN  ESTIMATE OF PARCEL POTENTIAL
                             !     TEMPERATURE IN LAYER K+1 AFTER
                             !     ENTRAINMENT AND LATENT HEATING (K)
                             ! OUT FINAL PARCEL POTENTIAL TEMPERATURE
                             !     IN LAYER K+1 (AFTER FORCED
                             !     DETRAINEMNT) (K)
C
      REAL QPKP1(ln2)        ! INOUT
                             ! IN  ESTIMATE OF PARCEL MIXING RATIO
                             !     IN LAYER K+1 AFTER ENTRAINMENT AND
                             !     LATENT HEATING (KG/KG)
                             ! OUT FINAL PARCEL MIXING RATIO
                             !     IN LAYER K+1 (AFTER FORCED
                             !     DETRAINEMNT) (KG/KG)
C
      REAL FLXK(ln2)         ! INOUT
                             ! IN  PARCEL MASSFLUX IN LAYER K
                             !     (NON-ZERO IF CONVECTION IS NOT
                             !     INITIATED FROM LAYER K) (PA/S)
                             ! OUT PARCEL MASSFLUX IN LAYER K
                             !     (SET IF CONVECTION IS INITIATED
                             !     IN LAYER K) (PA/S)
C
      LOGICAL BGMKP1(ln2)    ! INOUT
                             ! IN  MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1
                             !     CALCULATED ON THE BASIS OF
                             !     INPUT PARCEL POTENTIAL TEMPERATURE
                             !     AND MIXING RATIO
                             ! OUT MASK FOR PARCELS WHICH ARE
                             !     SATURATED IN LAYER K+1 CALCULATED
                             !     FORM PARCEL TEMPERATURE AND
                             !     MIXING RATIO AFTER FORCED
                             !     DETARINMENT CALCULATION
C
      REAL TCW(ln2)          ! INOUT
                             ! IN  TOTAL CONDENSED WATER CONTENT
                             !     SUMMED UPTO LAYER K (KG/M**2/S)
                             ! OUT UPDATED TOTAL CONDENSED WATER
                             !     CONTENT SUMMED UPTO LAYER K+1
                             !     (KG/M**2/S)
C
      REAL DEPTH(ln2)        ! INOUT
                             ! IN  DEPTH OF CONVECTIVE CLOUD TO
                             !     LAYER K (M)
                             ! OUT UPDATED DEPTH OF CONVECTIVE
                             !     CLOUD TO LAYER K+1 (M)
C
      REAL CCLWP(ln2)        ! INOUT
                             ! IN  CONDENSED WATER PATH
                             !     SUMMED UPTO LAYER K (KG/M**2)
                             ! OUT UPDATED CONDENSED WATER PATH
                             !     SUMMED UPTO LAYER K+1 (KG/M**2)

      real fs(ln2,ntrac)     !INOUT Scavenging fraction
C
C
C---------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C---------------------------------------------------------------------
C
      LOGICAL BTERM(ln2)     ! OUT MASK FOR PARCELS WHICH TERMINATE IN
                             !     LAYER K+1
C
      REAL PREKP1(ln2)       ! OUT PRECIPITATION FROM PARCEL AS IT
                             !     RISES FROM LAYER K TO K+1 (KG/M**2/S)
C
      REAL THRK(ln2)         ! OUT PARCEL DETRAINMENT POTENTIAL
                             !     TEMPERATURE IN LAYER K (K)
C
      REAL QRK(ln2)          ! OUT PARCEL DETRAINMENT MIXING RATIO
                             !     IN LAYER K (KG/KG)
C
      REAL XPKP1(ln2)        ! OUT PARCEL CLOUD WATER IN LAYER K+1
                             !     (KG/KG)
C
      REAL FLXKP1(ln2)       ! OUT PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      REAL DELTAK(ln2)       ! OUT PARCEL FORCED DETRAINMENT
                             !     COEFFICIENT IN LAYER K
                             !     MULTIPLIED BY APPROPRIATE
                             !     LAYER THICKNESS
C
      REAL CCA(ln2)          ! OUT CONVECTIVE CLOUD AMOUNT (%)
C
      INTEGER ICCB(ln2)      ! OUT CONVECTIVE CLOUD BASE LEVEL
C
      INTEGER ICCT(ln2)      ! OUT CONVECTIVE CLOUD TOP LEVEL
C
      REAL CCW(ln2)          ! OUT CONVECTIVE CLOUD LIQUID WATER
                             ! (KG/KG) ON MODEL LEVELS
C
      REAL CCWI(ln2)         ! OUT CONVECTIVE CLOUD ICE
                             ! (KG/KG) ON MODEL LEVELS
C
      REAL CCWL(ln2)         ! OUT CONVECTIVE CLOUD LIQUID
                             ! (KG/KG) ON MODEL LEVELS
C
C
C---------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C---------------------------------------------------------------------
C
      REAL THEK_C(ln2)       ! COMPRESSED POTENTIAL TEMPERATURE OF
                             ! CLOUD ENVIRONMENT IN LAYER K (K)
C
      REAL THEKP1_C(ln2)     ! COMPRESSED POTENTIAL TEMPERATURE OF
                             ! CLOUD ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEK_C(ln2)        ! COMPRESSED MIXING RATIO OF CLOUD
                             ! ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL QEKP1_C(ln2)      ! COMPRESSED MIXING RATIO OF CLOUD
                             ! ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL QSEK_C(ln2)       ! COMPRESSED SATURATION MIXING RATIO OF
                             ! CLOUD ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL DQSK_C(ln2)       ! COMPRESSED GRADIENT OF SATURATION MIXING
                             ! RATIO WITH POTENTIAL TEMPERATURE FOR THE
                             ! CLOUD ENVIRONMENT OF LAYER K (KG/KG/K)
C
      REAL QSEKP1_C(ln2)     ! COMPRESSED SATURATION MIXING RATIO OF
                             ! CLOUD ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL DQSKP1_C(ln2)     ! COMPRESSED GRADIENT OF SATURATION MIXING
                             ! RATIO WITH POTENTIAL TEMPERATURE FOR
                             ! THE CLOUD ENVIRONMENT IN LAYER K+1
                             ! (KG/KG/K)
C
      REAL THPK_C(ln2)       ! COMPRESSED PARCEL POTENTIAL
                             ! TEMPERATURE IN LAYER K (K)
C
      REAL QPK_C(ln2)        ! COMPRESSED PARCEL MIXING RATIO IN
                             ! LAYER K (KG/KG)
C
      REAL THPKP1_C(ln2)     ! COMPRESSED PARCEL POTENTIAL
                             ! TEMPERATURE IN LAYER K+1 (K)
C
      REAL QPKP1_C(ln2)      ! COMPRESSED PARCEL MIXING RATIO
                             ! IN LAYER K+1 (KG/KG)
C
      REAL XSQKP1_C(ln2)     ! EXCESS PARCEL WATER AFER LIFTING
                             ! FROM LAYER K TO K+1 (KG/KG)
C
      REAL THRK_C(ln2)       ! COMPRESSED PARCEL DETRAINMENT
                             ! POTENTIAL TEMPERATURE IN LAYER K (K)
C
      REAL QRK_C(ln2)        ! COMPRESSED PARCEL DETRAINMENT MIXING
                             ! RATIO IN LAYER K (KG/KG)
C
      REAL DELTAK_C(ln2)     ! COMPRESSED PARCEL FORCED DETRAINMENT
                             ! COEFFICIENT IN LAYER K
                             ! MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
C
      REAL EKP14_C(ln2)      ! COMPRESSED IN ENTRAINMENT COEFFICIENT AT
                             ! LEVEL K+1/4 MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
C
      REAL EKP34_C(ln2)      ! COMPRESSED ENTRAINMENT COEFFICIENT AT
                             ! LEVEL K+3/4 MULTIPLIED BY APPROPRIATE
                             ! LAYER THICKNESS
C
      REAL PK_C(ln2)         ! COMPRESSED PRESSURE AT LEVEL K (PA)
C
      REAL PKP1_C(ln2)       ! COMPRESSED PRESSURE AT LEVEL K+1 (PA)
C
      REAL EXK_C(ln2)        ! COMPRESSED EXNER FUNCTION AT LEVEL K
C
      REAL EXKP1_C(ln2)      ! COMPRESSED EXNER FUNCTION AT LEVEL K+1
C
      LOGICAL BWKP1_C(ln2)   ! COMPRESSED MASK FOR WHETHER CONDENSATE
                             ! IS LIQUID IN LAYER K+1
C
      LOGICAL BGMK_C(ln2)    ! COMPRESSED MASK FOR PARCELS WHICH ARE
                             ! SATURATED IN LAYER K
C
      LOGICAL BGMKP1_C(ln2)   ! COMPRESSED MASK FOR PARCELS
                              ! WHICH ARESATURATED IN LAYER K+1
C
      INTEGER INDEX1(ln2)    ! INDEX FOR COMPRESS AND EXPAND
C
      LOGICAL BDETK(ln2)     ! MASK FOR POINTS UNDERGOING
                             ! FORCED DETRAINMENT
CHBG
      REAL FRACENV  ! FRACTION OF DETRAINED LIQUID WATER INTO ENVIRON
C
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL DETRAIN,TERM_CON,CLOUD_W
C
C*--------------------------------------------------------------------
C
C
      DO 5 I=1,NPNTS
CL
CL---------------------------------------------------------------------
CL CALCULATE MASK FOR THOSE POINTS UNDERGOING FORCED DETRAINMENT
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (6), EQUATION (23)
CL---------------------------------------------------------------------
CL
       BDETK(I) = RBUOY(I) .LT. XSBMIN
C
   5  CONTINUE
CL
CL----------------------------------------------------------------------
CL  COMPRESS ALL INPUT ARRAYS FOR THE FORCED DETRAINMENT CALCULATIONS
CL----------------------------------------------------------------------
CL
      NDET = 0
      DO 10 I=1,NPNTS
       IF (BDETK(I))THEN
         NDET = NDET + 1
         INDEX1(NDET) = I
       END IF
  10  CONTINUE
C
      IF (NDET .NE. 0) THEN
*vdir nodep
        DO 35 I=1,NDET
          THEK_C(I)  = THEK(INDEX1(I))
          QEK_C(I)   = QEK(INDEX1(I))
          THPK_C(I)  = THPK(INDEX1(I))
          QPK_C(I)   = QPK(INDEX1(I))
          QSEK_C(I)  = QSEK(INDEX1(I))
          DQSK_C(I)  = DQSK(INDEX1(I))
          THEKP1_C(I)= THEKP1(INDEX1(I))
          QEKP1_C(I) = QEKP1(INDEX1(I))
          THPKP1_C(I)= THPKP1(INDEX1(I))
          QPKP1_C(I) = QPKP1(INDEX1(I))
          QSEKP1_C(I)= QSEKP1(INDEX1(I))
          DQSKP1_C(I)= DQSKP1(INDEX1(I))
          XSQKP1_C(I)= XSQKP1(INDEX1(I))
          EKP14_C(I) = EKP14(INDEX1(I))
          EKP34_C(I) = EKP34(INDEX1(I))
          PK_C(I)    = PK(INDEX1(I))
          PKP1_C(I)  = PKP1(INDEX1(I))
          EXK_C(I)   = EXK(INDEX1(I))
          EXKP1_C(I) = EXKP1(INDEX1(I))
C
          BGMK_C(I)  = BGMK(INDEX1(I))
          BGMKP1_C(I)= BGMKP1(INDEX1(I))
          BWKP1_C(I) = BWKP1(INDEX1(I))
  35    CONTINUE
CL
CL-------------------------------------------------------------------
CL DETRAINMENT CALCULATION
CL
CL SUBROUTINE DETRAIN
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (6)
CL-------------------------------------------------------------------
CL
         CALL DETRAIN (NDET,THEK_C,QEK_C,THPK_C,QPK_C,
     &                 QSEK_C,DQSK_C,BGMK_C,THEKP1_C,
     &                 QEKP1_C,THPKP1_C,QPKP1_C,QSEKP1_C,
     &                 DQSKP1_C,BGMKP1_C,BWKP1_C,
     &                 XSQKP1_C,DELTAK_C,
     &                 THRK_C,QRK_C,EKP14_C,EKP34_C,
     &                 PK_C,PKP1_C,EXK_C,EXKP1_C)
C
C-----------------------------------------------------------------------
C   DECOMPRESS/EXPAND OUTPUT ARRAYS FROM THE DETRAINMENT CALCULATIONS
C-----------------------------------------------------------------------
C
C
CDIR$ IVDEP
*vdir nodep
        DO 40 I=1,NDET
          THPKP1(INDEX1(I)) = THPKP1_C(I)
          QPKP1(INDEX1(I))  = QPKP1_C(I)
          XSQKP1(INDEX1(I)) = XSQKP1_C(I)
C
          BGMKP1(INDEX1(I)) = BGMKP1_C(I)
  40    CONTINUE
      ENDIF
C
      DO 45 I=1,NPNTS
        DELTAK(I) = 0.0
        THRK(I) = 0.0
        QRK(I) = 0.0
  45  CONTINUE
C
CDIR$ IVDEP
*vdir nodep
      DO 50 I=1,NDET
        DELTAK(INDEX1(I)) = DELTAK_C(I)
        THRK(INDEX1(I))   = THRK_C(I)
        QRK(INDEX1(I))    = QRK_C(I)
  50  CONTINUE
CL
CL----------------------------------------------------------------------
CL  CALCULATE MASS FLUX AT LEVEL K+1.
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (2B), EQUATION (10A)
CL----------------------------------------------------------------------
CL
      DO 60 I=1,NPNTS
       FLXKP1(I) = FLXK(I)*(1.+EKP14(I))*(1.+EKP34(I))*(1.-DELTAK(I))*
     &                                           (1.-AMDETK(I))
  60  CONTINUE
CL
CL---------------------------------------------------------------------
CL TEST FOR POINTS AT WHICH CONVECTION TERMINATES IN LAYER K+1
CL
CL SUBROUTINE TERM_CON
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (7)
CL---------------------------------------------------------------------
CL
      CALL TERM_CON (NPNTS,NLEV,K,BTERM,BWKP1,FLXKP1,THEKP1,QEKP1,THPI,
     &               QPI,QSEKP1,DELTAK,EXKP1,EKP14,EKP34,PSTAR)
CL
CL----------------------------------------------------------------------
CL CLOUD MICROPHYSICS CALCULATION
CL
CL SUBROUTINE CLOUD_W
CL
CL UM DOCUMENTATION PAPER P27
CL SECTION (8), (9)
CL----------------------------------------------------------------------
CL
CHBG
CHBG PUT A FRACTION OF DETRAINED LIQUID WATER INTO THE ENVIRONMENT
CHBG (REST INTO PRECIPITATION) : UKMO USED FRACENV=1.0
      FRACENV=1.00
      IF(FRACENV.EQ.1.0)THEN
        DO 61 I=1,NPNTS
         PREKP1(I) = 0.0
  61    CONTINUE
      ELSE
        DO 611 I=1,NPNTS
         PREKP1(I) = (1.0-FRACENV)*(FLXK(I)/G)
     &              *XPK(I)*(AMDETK(I)+DELTAK(I)*(1.-AMDETK(I)))
  611   CONTINUE
      ENDIF
CHBG
      CALL CLOUD_W (K,NPNTS,XPKP1,PREKP1,XSQKP1,BLOWST,FLXKP1,
     &              XPK,THEKP1,QEKP1,BWKP1,BLAND,QSEKP1,BGMKP1,
     &              BTERM,CCA,ICCB,ICCT,TCW,DEPTH,EKP14,EKP34,DELEXKP1,
     &              CCLWP,DELPKP1,CCW,CCWI,CCWL,
     &              qcloud,ipass,thpkp1,exkp1,pkp1,xs,fs,DELPK)
C
      RETURN
      END
CLL  SUBROUTINE PEVP_BCB-----------------------------------------------
CLL
CLL  PURPOSE : EVAPORATE RAIN BELOW CLOUD BASE IF NO DOWNDRAUGHT
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4 DATED 5/2/92
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE PEVP_BCB (NPNTS,K,ICCB,TH,PK,Q,DELP,RAIN,SNOW,
     &                     DTHBYDT,DQBYDT,EXK,TIMESTEP,CCA)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C  CONSTANTS
C-----------------------------------------------------------------------
C
CHBG  include 'C_LHEAT.f'
CHBG  include 'C_R_CP.f'
CHBG  include 'C_G.f'
CHBG  include 'CLDAREA.f'
      include 'UKPARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                  ! IN LOOP COUNTER
C
      INTEGER NPNTS              ! VECTOR LENGTH
C
      INTEGER K                  ! IN PRESENT MODEL LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      INTEGER ICCB(ln2)          ! IN CONVECTIVE CLOUD BASE LAYER
C
      REAL PK(ln2)               ! IN PRESSURE (PA)
C
      REAL Q(ln2)                ! IN MIXING RATIO (KG/KG)
C
      REAL TH(ln2)               ! IN POTENTIAL TEMPERATURE (K)
C
      REAL DELP(ln2)             ! IN CHANGE IN PRESSURE ACROSS
                                 !    LAYER K-1 (PA)
C
      REAL EXK(ln2)              ! IN EXNER RATIO OF LAYER K
C
      REAL TIMESTEP              ! IN MODEL TIMESTEP (S)
C
      REAL CCA(ln2)              ! IN CONVECTIVE CLOUD AMOUNT
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT AND OUTPUT
C-----------------------------------------------------------------------
C
      REAL DTHBYDT(ln2)          ! INOUT
                                 ! IN  INCREMENT TO MODEL POTENTIAL
                                 !     TEMPERATURE (K/S)
                                 ! OUT UPDATED INCREMENT TO MODEL
                                 !     POTENTIAL TEMPERATURE (K/S)
C
      REAL DQBYDT(ln2)           ! INOUT
                                 ! IN  INCREMENT TO MODEL MIXING RATIO
                                 !     (KG/KG/S)
                                 ! OUT UPDATED INCREMENT TO MIXING RATIO
                                 !     AFTER EVAPORATION BELOW CLOUD
                                 !     BASE (KG/KG/S)
C
      REAL RAIN(ln2)             ! INOUT
                                 ! IN  AMOUNT OF FALLING RAIN
                                 !     (KG/M**2/S)
                                 ! OUT UPDATED AMOUNT OF FALLING RAIN
                                 !     (KG/M**2/S)
C
      REAL SNOW(ln2)             ! INOUT
                                 ! IN  AMOUNT OF FALLING SNOW
                                 !     (KG/M**2/S)
                                 ! OUT UPDATED AMOUNT OF FALLING SNOW
                                 !     (KG/M**2/S)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
C
      REAL T(ln2)                ! MODEL TEMPERATURE (K)
C
      REAL EVAP_RAIN(ln2)        ! AMOUNT OF EVAPORATION OF RAIN
C
      REAL SUB_SNOW(ln2)         ! AMOUNT OF SNOW SUBLIMATION
C
      REAL QSATE(ln2)            ! SATURATED MIXING RATIO IN
                                 ! ENVIRONMENT (KG/KG)
C
      REAL DELQ(ln2)             ! CHANGE IN MIXING RATIO ACROSS LAYER K
                                 ! (KG/KG)
C
      REAL THS(ln2)              ! SATURATED PARCEL POTENTIAL
                                 ! TEMPERATURE (K)
C
      REAL QS(ln2)               ! SATURATED PARCEL MIXING RATIO
C
      LOGICAL BEVAP(ln2)         ! MASK FOR THOSE POINTS WHERE
                                 ! EVAPORATION OCCURS
C
      REAL DTHBYDT_EVP           ! INCREMENT TO POTENTIAL TEMPERATURE
                                 ! DUE TO EVAPORATION (K)
C
      REAL DQBYDT_EVP            ! INCREMENT TO MIXING RATIO DUE TO
                                 ! EVAPORATION (KG/KG)
C
      REAL DTHBYDT_SAT           ! INCREMENT TO POTENTIAL TEMPERATURE
                                 ! DUE TO SATURATION (K)
C
      REAL FACTOR                ! DTHBYDT_SAT / DTHBYDT_EVP
C
      REAL RHO(ln2)              ! DENSITY OF AIR IN PARCEL
C
      LOGICAL FLAGL
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL QSATU, EVP, SATCAL
C
C-----------------------------------------------------------------------
C EVAPORATE RAIN IN LAYER K IF LAYER K IS BELOW CLOUD BASE
C CALCULATE MOISTURE SUB-SATURATION
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
        T(I) = TH(I)*EXK(I)
        BEVAP(I) = .FALSE.
      END DO
C
      CALL QSATU(QSATE,T,PK,NPNTS)
C
      DO I=1,NPNTS
         DELQ(I) = QSATE(I)-Q(I)
            RHO(I) = PK(I) / (R*T(I))
       IF (K .LT. ICCB(I)) THEN
C
C-----------------------------------------------------------------------
C CHECK IF EVAPORATION POSSIBLE
C-----------------------------------------------------------------------
C
         IF ((RAIN(I).GT.0.0 .OR. SNOW(I).GT.0.0) .AND.
     &        DELQ(I) .GT. 0.0) THEN
C
            BEVAP(I) = .TRUE.
         END IF
       END IF
      END DO
C
      call checkl(BEVAP,NPNTS,FLAGL)
      if(.not.FLAGL)return
C
C-----------------------------------------------------------------------
C CALCULATE EVAPORATION
C-----------------------------------------------------------------------
C
        CALL EVP (NPNTS,RAIN,T,CCA,RHO,DELQ,DELP,EVAP_RAIN,
     &            BEVAP,1,CLDAREA)
C
        CALL EVP (NPNTS,SNOW,T,CCA,RHO,DELQ,DELP,SUB_SNOW,
     &            BEVAP,2,CLDAREA)
C
C-----------------------------------------------------------------------
C CALCULATE TEMPERATURE AND MIXING RATIO IF LAYER BROUGHT TO
C SATURATION BY EVAPORATION AND SUBLIMATION
C-----------------------------------------------------------------------
C
      CALL SATCAL(NPNTS,TH,PK,QS,THS,EXK,Q,TH)
C
C
      DO I=1,NPNTS
        IF (BEVAP(I)) THEN
          DTHBYDT_EVP = -((LC*EVAP_RAIN(I))+((LC+LF)*SUB_SNOW(I)))*G/
     &                   (CP*EXK(I)*DELP(I))
          DQBYDT_EVP = (EVAP_RAIN(I)+SUB_SNOW(I))*G/DELP(I)
C
          DTHBYDT_SAT = (THS(I)-TH(I))/TIMESTEP
C
          IF (DTHBYDT_EVP.LT.DTHBYDT_SAT) THEN
C
C---------------------------------------------------------------------
C  ADJUST EVAPORATION AND SUBLIMATION RATES TO GIVE SATURATION
C---------------------------------------------------------------------
C
            FACTOR = DTHBYDT_SAT/DTHBYDT_EVP
            DTHBYDT_EVP = DTHBYDT_SAT
            DQBYDT_EVP = DQBYDT_EVP*FACTOR
            EVAP_RAIN(I) = EVAP_RAIN(I)*FACTOR
            SUB_SNOW(I) = SUB_SNOW(I)*FACTOR
          END IF
C
C---------------------------------------------------------------------
C  UPDATE INCREMENTS AND RAINFALL AND ADJUST BACK TO GRIDBOX MEANS
C---------------------------------------------------------------------
C
          DTHBYDT(I) = DTHBYDT(I)+DTHBYDT_EVP*CCA(I)*CLDAREA
          DQBYDT(I) = DQBYDT(I)+DQBYDT_EVP*CCA(I)*CLDAREA
          RAIN(I) = RAIN(I)-EVAP_RAIN(I)*CCA(I)*CLDAREA
          SNOW(I) = SNOW(I)-SUB_SNOW(I)*CCA(I)*CLDAREA
        END IF
      END DO
C
      RETURN
      END
C
CLL  SUBROUTINE QSATU--------------------------------------------------
CLL
CLL  PURPOSE : RETURNS A SATURATION MIXING RATIO GIVEN
CLL            A TEMPERATURE AND PRESSURE USING SATURATION
CLL            VAPOUR PRESSURES CALCULATED USING THE
CLL            GOFF-GRATCH FORMULAE, ADOPTED BY THE WMO AS
CLL            TAKEN FROM LANDOLT-BORNSTEIN, 1987 NUMERICAL
CLL            DATA AND FUNCTIONAL RELATIONSHIPS IN SCIENCE
CLL            AND TECHNOLOGY. GROUP V/VOL 4B METEOROLOGY.
CLL            PHYSICAL AND CHEMICAL PROPERTIES OF AIR, P35
CLL
CLL            VALUES IN THE LOOKUP TABLE ARE OVER WATER ABOVE
CLL            0 DEG C AND OVER ICE BELOW THIS TEMPERATURE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED:
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  DOCUMENTATION :
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE QSATU (QS,T,P,NPNTS)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include 'C_EPSLON.f'
      include 'UKPARAMS.f'
C
C----------------------------------------------------------------------
C LOCAL CONSTANTS
C----------------------------------------------------------------------
C
      REAL T_LOW           ! LOWEST TEMPERATURE FOR WHICH LOOK-UP
                           ! TABLE OF SATURATION WATER VAPOUR
                           ! PRESSURE IS VALID (K)
C
      REAL T_HIGH          ! HIGHEST TEMPERATURE FOR WHICH LOOK-UP
                           ! TABLE OF SATURATION WATER VAPOUR
                           ! PRESSURES IS VALID (K)
C
      REAL DELTA_T         ! TEMPERATURE INCREMENT OF THE LOOK-UP
                           ! TABLE OF SATURATION VAPOUR PRESSURES
C
      INTEGER N            ! SIZE OF LOOK-UP TABLE OF SATURATION
                           ! WATER VAPOUR PRESSURES
C
      PARAMETER ( T_LOW = 183.15,
     &            T_HIGH = 338.15,
     &            DELTA_T = 0.1,
     &            N = INT(((T_HIGH - T_LOW + (DELTA_T*0.5))/DELTA_T)
     &                    + 1.0)
     &          )
C
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS        ! VECTOR LENGTH
C
      INTEGER I            ! LOOP COUNTER
C
      INTEGER IES          ! LOOP COUNTER FOR DATA STATEMENT
                           ! LOOK-UP TABLE
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C----------------------------------------------------------------------
C
      REAL T(ln2)          ! IN TEMPERATURE (K)
C
      REAL P(ln2)          ! IN PRESSURE (PA)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL QS(ln2)         ! OUT SATURATION MIXING RATIO AT TEMPERATURE
                           !     T AND PRESSURE P (KG/KG)
C
C
C----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C----------------------------------------------------------------------
C
      REAL ES(N)           ! TABLE OF SATURATION WATER VAPOUR
                           ! PRESSURE (PA) - SET BY DATA STATEMENT
                           ! CALCULATED FROM THE GOFF-GRATCH FORMULAE
                           ! AS TAKEN FROM LANDOLT-BORNSTEIN, 1987
                           ! NUMERICAL DATA AND FUNCTIONAL RELATIONSHIPS
                           ! IN SCIENCE AND TECHNOLOGY. GROUP V/ VOL 4B
                           ! METEOROLOGY. PHYSICAL AND CHEMICAL
                           ! PROPERTIES OF AIR, P35
C
      REAL ATABLE          ! WORK VARIABLES
C
      INTEGER ITABLE       ! WORK VARIABLES
C
C*---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL NO SIGNIFICANT STRUCTURE
CL---------------------------------------------------------------------
CL
C
C----------------------------------------------------------------------
C SATURATION WATER VAPOUR PRESSURE
C
C ABOVE 0 DEG C VALUES ARE OVER WATER
C
C BELOW 0 DEC C VALUES ARE OVER ICE
C----------------------------------------------------------------------
C
      DATA (ES(IES),IES=    1, 95) /
     &0.966483E-02,0.984279E-02,0.100240E-01,0.102082E-01,0.103957E-01,
     &0.105865E-01,0.107803E-01,0.109777E-01,0.111784E-01,0.113825E-01,
     &0.115902E-01,0.118016E-01,0.120164E-01,0.122348E-01,0.124572E-01,
     &0.126831E-01,0.129132E-01,0.131470E-01,0.133846E-01,0.136264E-01,
     &0.138724E-01,0.141225E-01,0.143771E-01,0.146356E-01,0.148985E-01,
     &0.151661E-01,0.154379E-01,0.157145E-01,0.159958E-01,0.162817E-01,
     &0.165725E-01,0.168680E-01,0.171684E-01,0.174742E-01,0.177847E-01,
     &0.181008E-01,0.184216E-01,0.187481E-01,0.190801E-01,0.194175E-01,
     &0.197608E-01,0.201094E-01,0.204637E-01,0.208242E-01,0.211906E-01,
     &0.215631E-01,0.219416E-01,0.223263E-01,0.227172E-01,0.231146E-01,
     &0.235188E-01,0.239296E-01,0.243465E-01,0.247708E-01,0.252019E-01,
     &0.256405E-01,0.260857E-01,0.265385E-01,0.269979E-01,0.274656E-01,
     &0.279405E-01,0.284232E-01,0.289142E-01,0.294124E-01,0.299192E-01,
     &0.304341E-01,0.309571E-01,0.314886E-01,0.320285E-01,0.325769E-01,
     &0.331348E-01,0.337014E-01,0.342771E-01,0.348618E-01,0.354557E-01,
     &0.360598E-01,0.366727E-01,0.372958E-01,0.379289E-01,0.385717E-01,
     &0.392248E-01,0.398889E-01,0.405633E-01,0.412474E-01,0.419430E-01,
     &0.426505E-01,0.433678E-01,0.440974E-01,0.448374E-01,0.455896E-01,
     &0.463545E-01,0.471303E-01,0.479191E-01,0.487190E-01,0.495322E-01/
      DATA (ES(IES),IES= 96,190) /
     &0.503591E-01,0.511977E-01,0.520490E-01,0.529145E-01,0.537931E-01,
     &0.546854E-01,0.555924E-01,0.565119E-01,0.574467E-01,0.583959E-01,
     &0.593592E-01,0.603387E-01,0.613316E-01,0.623409E-01,0.633655E-01,
     &0.644053E-01,0.654624E-01,0.665358E-01,0.676233E-01,0.687302E-01,
     &0.698524E-01,0.709929E-01,0.721490E-01,0.733238E-01,0.745180E-01,
     &0.757281E-01,0.769578E-01,0.782061E-01,0.794728E-01,0.807583E-01,
     &0.820647E-01,0.833905E-01,0.847358E-01,0.861028E-01,0.874882E-01,
     &0.888957E-01,0.903243E-01,0.917736E-01,0.932464E-01,0.947407E-01,
     &0.962571E-01,0.977955E-01,0.993584E-01,0.100942E+00,0.102551E+00,
     &0.104186E+00,0.105842E+00,0.107524E+00,0.109231E+00,0.110963E+00,
     &0.112722E+00,0.114506E+00,0.116317E+00,0.118153E+00,0.120019E+00,
     &0.121911E+00,0.123831E+00,0.125778E+00,0.127755E+00,0.129761E+00,
     &0.131796E+00,0.133863E+00,0.135956E+00,0.138082E+00,0.140241E+00,
     &0.142428E+00,0.144649E+00,0.146902E+00,0.149190E+00,0.151506E+00,
     &0.153859E+00,0.156245E+00,0.158669E+00,0.161126E+00,0.163618E+00,
     &0.166145E+00,0.168711E+00,0.171313E+00,0.173951E+00,0.176626E+00,
     &0.179342E+00,0.182096E+00,0.184893E+00,0.187724E+00,0.190600E+00,
     &0.193518E+00,0.196473E+00,0.199474E+00,0.202516E+00,0.205604E+00,
     &0.208730E+00,0.211905E+00,0.215127E+00,0.218389E+00,0.221701E+00/
      DATA (ES(IES),IES=191,285) /
     &0.225063E+00,0.228466E+00,0.231920E+00,0.235421E+00,0.238976E+00,
     &0.242580E+00,0.246232E+00,0.249933E+00,0.253691E+00,0.257499E+00,
     &0.261359E+00,0.265278E+00,0.269249E+00,0.273274E+00,0.277358E+00,
     &0.281498E+00,0.285694E+00,0.289952E+00,0.294268E+00,0.298641E+00,
     &0.303078E+00,0.307577E+00,0.312135E+00,0.316753E+00,0.321440E+00,
     &0.326196E+00,0.331009E+00,0.335893E+00,0.340842E+00,0.345863E+00,
     &0.350951E+00,0.356106E+00,0.361337E+00,0.366636E+00,0.372006E+00,
     &0.377447E+00,0.382966E+00,0.388567E+00,0.394233E+00,0.399981E+00,
     &0.405806E+00,0.411714E+00,0.417699E+00,0.423772E+00,0.429914E+00,
     &0.436145E+00,0.442468E+00,0.448862E+00,0.455359E+00,0.461930E+00,
     &0.468596E+00,0.475348E+00,0.482186E+00,0.489124E+00,0.496160E+00,
     &0.503278E+00,0.510497E+00,0.517808E+00,0.525224E+00,0.532737E+00,
     &0.540355E+00,0.548059E+00,0.555886E+00,0.563797E+00,0.571825E+00,
     &0.579952E+00,0.588198E+00,0.596545E+00,0.605000E+00,0.613572E+00,
     &0.622255E+00,0.631059E+00,0.639962E+00,0.649003E+00,0.658144E+00,
     &0.667414E+00,0.676815E+00,0.686317E+00,0.695956E+00,0.705728E+00,
     &0.715622E+00,0.725641E+00,0.735799E+00,0.746082E+00,0.756495E+00,
     &0.767052E+00,0.777741E+00,0.788576E+00,0.799549E+00,0.810656E+00,
     &0.821914E+00,0.833314E+00,0.844854E+00,0.856555E+00,0.868415E+00/
      DATA (ES(IES),IES=286,380) /
     &0.880404E+00,0.892575E+00,0.904877E+00,0.917350E+00,0.929974E+00,
     &0.942771E+00,0.955724E+00,0.968837E+00,0.982127E+00,0.995600E+00,
     &0.100921E+01,0.102304E+01,0.103700E+01,0.105116E+01,0.106549E+01,
     &0.108002E+01,0.109471E+01,0.110962E+01,0.112469E+01,0.113995E+01,
     &0.115542E+01,0.117107E+01,0.118693E+01,0.120298E+01,0.121923E+01,
     &0.123569E+01,0.125234E+01,0.126923E+01,0.128631E+01,0.130362E+01,
     &0.132114E+01,0.133887E+01,0.135683E+01,0.137500E+01,0.139342E+01,
     &0.141205E+01,0.143091E+01,0.145000E+01,0.146933E+01,0.148892E+01,
     &0.150874E+01,0.152881E+01,0.154912E+01,0.156970E+01,0.159049E+01,
     &0.161159E+01,0.163293E+01,0.165452E+01,0.167640E+01,0.169852E+01,
     &0.172091E+01,0.174359E+01,0.176653E+01,0.178977E+01,0.181332E+01,
     &0.183709E+01,0.186119E+01,0.188559E+01,0.191028E+01,0.193524E+01,
     &0.196054E+01,0.198616E+01,0.201208E+01,0.203829E+01,0.206485E+01,
     &0.209170E+01,0.211885E+01,0.214637E+01,0.217424E+01,0.220242E+01,
     &0.223092E+01,0.225979E+01,0.228899E+01,0.231855E+01,0.234845E+01,
     &0.237874E+01,0.240937E+01,0.244040E+01,0.247176E+01,0.250349E+01,
     &0.253560E+01,0.256814E+01,0.260099E+01,0.263431E+01,0.266800E+01,
     &0.270207E+01,0.273656E+01,0.277145E+01,0.280671E+01,0.284248E+01,
     &0.287859E+01,0.291516E+01,0.295219E+01,0.298962E+01,0.302746E+01/
      DATA (ES(IES),IES=381,475) /
     &0.306579E+01,0.310454E+01,0.314377E+01,0.318351E+01,0.322360E+01,
     &0.326427E+01,0.330538E+01,0.334694E+01,0.338894E+01,0.343155E+01,
     &0.347456E+01,0.351809E+01,0.356216E+01,0.360673E+01,0.365184E+01,
     &0.369744E+01,0.374352E+01,0.379018E+01,0.383743E+01,0.388518E+01,
     &0.393344E+01,0.398230E+01,0.403177E+01,0.408175E+01,0.413229E+01,
     &0.418343E+01,0.423514E+01,0.428746E+01,0.434034E+01,0.439389E+01,
     &0.444808E+01,0.450276E+01,0.455820E+01,0.461423E+01,0.467084E+01,
     &0.472816E+01,0.478607E+01,0.484468E+01,0.490393E+01,0.496389E+01,
     &0.502446E+01,0.508580E+01,0.514776E+01,0.521047E+01,0.527385E+01,
     &0.533798E+01,0.540279E+01,0.546838E+01,0.553466E+01,0.560173E+01,
     &0.566949E+01,0.573807E+01,0.580750E+01,0.587749E+01,0.594846E+01,
     &0.602017E+01,0.609260E+01,0.616591E+01,0.623995E+01,0.631490E+01,
     &0.639061E+01,0.646723E+01,0.654477E+01,0.662293E+01,0.670220E+01,
     &0.678227E+01,0.686313E+01,0.694495E+01,0.702777E+01,0.711142E+01,
     &0.719592E+01,0.728140E+01,0.736790E+01,0.745527E+01,0.754352E+01,
     &0.763298E+01,0.772316E+01,0.781442E+01,0.790676E+01,0.800001E+01,
     &0.809435E+01,0.818967E+01,0.828606E+01,0.838343E+01,0.848194E+01,
     &0.858144E+01,0.868207E+01,0.878392E+01,0.888673E+01,0.899060E+01,
     &0.909567E+01,0.920172E+01,0.930909E+01,0.941765E+01,0.952730E+01/
      DATA (ES(IES),IES=476,570) /
     &0.963821E+01,0.975022E+01,0.986352E+01,0.997793E+01,0.100937E+02,
     &0.102105E+02,0.103287E+02,0.104481E+02,0.105688E+02,0.106909E+02,
     &0.108143E+02,0.109387E+02,0.110647E+02,0.111921E+02,0.113207E+02,
     &0.114508E+02,0.115821E+02,0.117149E+02,0.118490E+02,0.119847E+02,
     &0.121216E+02,0.122601E+02,0.124002E+02,0.125416E+02,0.126846E+02,
     &0.128290E+02,0.129747E+02,0.131224E+02,0.132712E+02,0.134220E+02,
     &0.135742E+02,0.137278E+02,0.138831E+02,0.140403E+02,0.141989E+02,
     &0.143589E+02,0.145211E+02,0.146845E+02,0.148501E+02,0.150172E+02,
     &0.151858E+02,0.153564E+02,0.155288E+02,0.157029E+02,0.158786E+02,
     &0.160562E+02,0.162358E+02,0.164174E+02,0.166004E+02,0.167858E+02,
     &0.169728E+02,0.171620E+02,0.173528E+02,0.175455E+02,0.177406E+02,
     &0.179372E+02,0.181363E+02,0.183372E+02,0.185400E+02,0.187453E+02,
     &0.189523E+02,0.191613E+02,0.193728E+02,0.195866E+02,0.198024E+02,
     &0.200200E+02,0.202401E+02,0.204626E+02,0.206871E+02,0.209140E+02,
     &0.211430E+02,0.213744E+02,0.216085E+02,0.218446E+02,0.220828E+02,
     &0.223241E+02,0.225671E+02,0.228132E+02,0.230615E+02,0.233120E+02,
     &0.235651E+02,0.238211E+02,0.240794E+02,0.243404E+02,0.246042E+02,
     &0.248704E+02,0.251390E+02,0.254109E+02,0.256847E+02,0.259620E+02,
     &0.262418E+02,0.265240E+02,0.268092E+02,0.270975E+02,0.273883E+02/
      DATA (ES(IES),IES=571,665) /
     &0.276822E+02,0.279792E+02,0.282789E+02,0.285812E+02,0.288867E+02,
     &0.291954E+02,0.295075E+02,0.298222E+02,0.301398E+02,0.304606E+02,
     &0.307848E+02,0.311119E+02,0.314424E+02,0.317763E+02,0.321133E+02,
     &0.324536E+02,0.327971E+02,0.331440E+02,0.334940E+02,0.338475E+02,
     &0.342050E+02,0.345654E+02,0.349295E+02,0.352975E+02,0.356687E+02,
     &0.360430E+02,0.364221E+02,0.368042E+02,0.371896E+02,0.375790E+02,
     &0.379725E+02,0.383692E+02,0.387702E+02,0.391744E+02,0.395839E+02,
     &0.399958E+02,0.404118E+02,0.408325E+02,0.412574E+02,0.416858E+02,
     &0.421188E+02,0.425551E+02,0.429962E+02,0.434407E+02,0.438910E+02,
     &0.443439E+02,0.448024E+02,0.452648E+02,0.457308E+02,0.462018E+02,
     &0.466775E+02,0.471582E+02,0.476428E+02,0.481313E+02,0.486249E+02,
     &0.491235E+02,0.496272E+02,0.501349E+02,0.506479E+02,0.511652E+02,
     &0.516876E+02,0.522142E+02,0.527474E+02,0.532836E+02,0.538266E+02,
     &0.543737E+02,0.549254E+02,0.554839E+02,0.560456E+02,0.566142E+02,
     &0.571872E+02,0.577662E+02,0.583498E+02,0.589392E+02,0.595347E+02,
     &0.601346E+02,0.607410E+02,0.613519E+02,0.619689E+02,0.625922E+02,
     &0.632204E+02,0.638550E+02,0.644959E+02,0.651418E+02,0.657942E+02,
     &0.664516E+02,0.671158E+02,0.677864E+02,0.684624E+02,0.691451E+02,
     &0.698345E+02,0.705293E+02,0.712312E+02,0.719398E+02,0.726542E+02/
      DATA (ES(IES),IES=666,760) /
     &0.733754E+02,0.741022E+02,0.748363E+02,0.755777E+02,0.763247E+02,
     &0.770791E+02,0.778394E+02,0.786088E+02,0.793824E+02,0.801653E+02,
     &0.809542E+02,0.817509E+02,0.825536E+02,0.833643E+02,0.841828E+02,
     &0.850076E+02,0.858405E+02,0.866797E+02,0.875289E+02,0.883827E+02,
     &0.892467E+02,0.901172E+02,0.909962E+02,0.918818E+02,0.927760E+02,
     &0.936790E+02,0.945887E+02,0.955071E+02,0.964346E+02,0.973689E+02,
     &0.983123E+02,0.992648E+02,0.100224E+03,0.101193E+03,0.102169E+03,
     &0.103155E+03,0.104150E+03,0.105152E+03,0.106164E+03,0.107186E+03,
     &0.108217E+03,0.109256E+03,0.110303E+03,0.111362E+03,0.112429E+03,
     &0.113503E+03,0.114588E+03,0.115684E+03,0.116789E+03,0.117903E+03,
     &0.119028E+03,0.120160E+03,0.121306E+03,0.122460E+03,0.123623E+03,
     &0.124796E+03,0.125981E+03,0.127174E+03,0.128381E+03,0.129594E+03,
     &0.130822E+03,0.132058E+03,0.133306E+03,0.134563E+03,0.135828E+03,
     &0.137109E+03,0.138402E+03,0.139700E+03,0.141017E+03,0.142338E+03,
     &0.143676E+03,0.145025E+03,0.146382E+03,0.147753E+03,0.149133E+03,
     &0.150529E+03,0.151935E+03,0.153351E+03,0.154783E+03,0.156222E+03,
     &0.157678E+03,0.159148E+03,0.160624E+03,0.162117E+03,0.163621E+03,
     &0.165142E+03,0.166674E+03,0.168212E+03,0.169772E+03,0.171340E+03,
     &0.172921E+03,0.174522E+03,0.176129E+03,0.177755E+03,0.179388E+03/
      DATA (ES(IES),IES=761,855) /
     &0.181040E+03,0.182707E+03,0.184382E+03,0.186076E+03,0.187782E+03,
     &0.189503E+03,0.191240E+03,0.192989E+03,0.194758E+03,0.196535E+03,
     &0.198332E+03,0.200141E+03,0.201963E+03,0.203805E+03,0.205656E+03,
     &0.207532E+03,0.209416E+03,0.211317E+03,0.213236E+03,0.215167E+03,
     &0.217121E+03,0.219087E+03,0.221067E+03,0.223064E+03,0.225080E+03,
     &0.227113E+03,0.229160E+03,0.231221E+03,0.233305E+03,0.235403E+03,
     &0.237520E+03,0.239655E+03,0.241805E+03,0.243979E+03,0.246163E+03,
     &0.248365E+03,0.250593E+03,0.252830E+03,0.255093E+03,0.257364E+03,
     &0.259667E+03,0.261979E+03,0.264312E+03,0.266666E+03,0.269034E+03,
     &0.271430E+03,0.273841E+03,0.276268E+03,0.278722E+03,0.281185E+03,
     &0.283677E+03,0.286190E+03,0.288714E+03,0.291266E+03,0.293834E+03,
     &0.296431E+03,0.299045E+03,0.301676E+03,0.304329E+03,0.307006E+03,
     &0.309706E+03,0.312423E+03,0.315165E+03,0.317930E+03,0.320705E+03,
     &0.323519E+03,0.326350E+03,0.329199E+03,0.332073E+03,0.334973E+03,
     &0.337897E+03,0.340839E+03,0.343800E+03,0.346794E+03,0.349806E+03,
     &0.352845E+03,0.355918E+03,0.358994E+03,0.362112E+03,0.365242E+03,
     &0.368407E+03,0.371599E+03,0.374802E+03,0.378042E+03,0.381293E+03,
     &0.384588E+03,0.387904E+03,0.391239E+03,0.394604E+03,0.397988E+03,
     &0.401411E+03,0.404862E+03,0.408326E+03,0.411829E+03,0.415352E+03/
      DATA (ES(IES),IES=856,950) /
     &0.418906E+03,0.422490E+03,0.426095E+03,0.429740E+03,0.433398E+03,
     &0.437097E+03,0.440827E+03,0.444570E+03,0.448354E+03,0.452160E+03,
     &0.455999E+03,0.459870E+03,0.463765E+03,0.467702E+03,0.471652E+03,
     &0.475646E+03,0.479674E+03,0.483715E+03,0.487811E+03,0.491911E+03,
     &0.496065E+03,0.500244E+03,0.504448E+03,0.508698E+03,0.512961E+03,
     &0.517282E+03,0.521617E+03,0.525989E+03,0.530397E+03,0.534831E+03,
     &0.539313E+03,0.543821E+03,0.548355E+03,0.552938E+03,0.557549E+03,
     &0.562197E+03,0.566884E+03,0.571598E+03,0.576351E+03,0.581131E+03,
     &0.585963E+03,0.590835E+03,0.595722E+03,0.600663E+03,0.605631E+03,
     &0.610641E+03,0.615151E+03,0.619625E+03,0.624140E+03,0.628671E+03,
     &0.633243E+03,0.637845E+03,0.642465E+03,0.647126E+03,0.651806E+03,
     &0.656527E+03,0.661279E+03,0.666049E+03,0.670861E+03,0.675692E+03,
     &0.680566E+03,0.685471E+03,0.690396E+03,0.695363E+03,0.700350E+03,
     &0.705381E+03,0.710444E+03,0.715527E+03,0.720654E+03,0.725801E+03,
     &0.730994E+03,0.736219E+03,0.741465E+03,0.746756E+03,0.752068E+03,
     &0.757426E+03,0.762819E+03,0.768231E+03,0.773692E+03,0.779172E+03,
     &0.784701E+03,0.790265E+03,0.795849E+03,0.801483E+03,0.807137E+03,
     &0.812842E+03,0.818582E+03,0.824343E+03,0.830153E+03,0.835987E+03,
     &0.841871E+03,0.847791E+03,0.853733E+03,0.859727E+03,0.865743E+03/
      DATA (ES(IES),IES=951,1045) /
     &0.871812E+03,0.877918E+03,0.884046E+03,0.890228E+03,0.896433E+03,
     &0.902690E+03,0.908987E+03,0.915307E+03,0.921681E+03,0.928078E+03,
     &0.934531E+03,0.941023E+03,0.947539E+03,0.954112E+03,0.960708E+03,
     &0.967361E+03,0.974053E+03,0.980771E+03,0.987545E+03,0.994345E+03,
     &0.100120E+04,0.100810E+04,0.101502E+04,0.102201E+04,0.102902E+04,
     &0.103608E+04,0.104320E+04,0.105033E+04,0.105753E+04,0.106475E+04,
     &0.107204E+04,0.107936E+04,0.108672E+04,0.109414E+04,0.110158E+04,
     &0.110908E+04,0.111663E+04,0.112421E+04,0.113185E+04,0.113952E+04,
     &0.114725E+04,0.115503E+04,0.116284E+04,0.117071E+04,0.117861E+04,
     &0.118658E+04,0.119459E+04,0.120264E+04,0.121074E+04,0.121888E+04,
     &0.122709E+04,0.123534E+04,0.124362E+04,0.125198E+04,0.126036E+04,
     &0.126881E+04,0.127731E+04,0.128584E+04,0.129444E+04,0.130307E+04,
     &0.131177E+04,0.132053E+04,0.132931E+04,0.133817E+04,0.134705E+04,
     &0.135602E+04,0.136503E+04,0.137407E+04,0.138319E+04,0.139234E+04,
     &0.140156E+04,0.141084E+04,0.142015E+04,0.142954E+04,0.143896E+04,
     &0.144845E+04,0.145800E+04,0.146759E+04,0.147725E+04,0.148694E+04,
     &0.149672E+04,0.150655E+04,0.151641E+04,0.152635E+04,0.153633E+04,
     &0.154639E+04,0.155650E+04,0.156665E+04,0.157688E+04,0.158715E+04,
     &0.159750E+04,0.160791E+04,0.161836E+04,0.162888E+04,0.163945E+04/
      DATA (ES(IES),IES=1046,1140) /
     &0.165010E+04,0.166081E+04,0.167155E+04,0.168238E+04,0.169325E+04,
     &0.170420E+04,0.171522E+04,0.172627E+04,0.173741E+04,0.174859E+04,
     &0.175986E+04,0.177119E+04,0.178256E+04,0.179402E+04,0.180552E+04,
     &0.181711E+04,0.182877E+04,0.184046E+04,0.185224E+04,0.186407E+04,
     &0.187599E+04,0.188797E+04,0.190000E+04,0.191212E+04,0.192428E+04,
     &0.193653E+04,0.194886E+04,0.196122E+04,0.197368E+04,0.198618E+04,
     &0.199878E+04,0.201145E+04,0.202416E+04,0.203698E+04,0.204983E+04,
     &0.206278E+04,0.207580E+04,0.208887E+04,0.210204E+04,0.211525E+04,
     &0.212856E+04,0.214195E+04,0.215538E+04,0.216892E+04,0.218249E+04,
     &0.219618E+04,0.220994E+04,0.222375E+04,0.223766E+04,0.225161E+04,
     &0.226567E+04,0.227981E+04,0.229399E+04,0.230829E+04,0.232263E+04,
     &0.233708E+04,0.235161E+04,0.236618E+04,0.238087E+04,0.239560E+04,
     &0.241044E+04,0.242538E+04,0.244035E+04,0.245544E+04,0.247057E+04,
     &0.248583E+04,0.250116E+04,0.251654E+04,0.253204E+04,0.254759E+04,
     &0.256325E+04,0.257901E+04,0.259480E+04,0.261073E+04,0.262670E+04,
     &0.264279E+04,0.265896E+04,0.267519E+04,0.269154E+04,0.270794E+04,
     &0.272447E+04,0.274108E+04,0.275774E+04,0.277453E+04,0.279137E+04,
     &0.280834E+04,0.282540E+04,0.284251E+04,0.285975E+04,0.287704E+04,
     &0.289446E+04,0.291198E+04,0.292954E+04,0.294725E+04,0.296499E+04/
      DATA (ES(IES),IES=1141,1235) /
     &0.298288E+04,0.300087E+04,0.301890E+04,0.303707E+04,0.305529E+04,
     &0.307365E+04,0.309211E+04,0.311062E+04,0.312927E+04,0.314798E+04,
     &0.316682E+04,0.318577E+04,0.320477E+04,0.322391E+04,0.324310E+04,
     &0.326245E+04,0.328189E+04,0.330138E+04,0.332103E+04,0.334073E+04,
     &0.336058E+04,0.338053E+04,0.340054E+04,0.342069E+04,0.344090E+04,
     &0.346127E+04,0.348174E+04,0.350227E+04,0.352295E+04,0.354369E+04,
     &0.356458E+04,0.358559E+04,0.360664E+04,0.362787E+04,0.364914E+04,
     &0.367058E+04,0.369212E+04,0.371373E+04,0.373548E+04,0.375731E+04,
     &0.377929E+04,0.380139E+04,0.382355E+04,0.384588E+04,0.386826E+04,
     &0.389081E+04,0.391348E+04,0.393620E+04,0.395910E+04,0.398205E+04,
     &0.400518E+04,0.402843E+04,0.405173E+04,0.407520E+04,0.409875E+04,
     &0.412246E+04,0.414630E+04,0.417019E+04,0.419427E+04,0.421840E+04,
     &0.424272E+04,0.426715E+04,0.429165E+04,0.431634E+04,0.434108E+04,
     &0.436602E+04,0.439107E+04,0.441618E+04,0.444149E+04,0.446685E+04,
     &0.449241E+04,0.451810E+04,0.454385E+04,0.456977E+04,0.459578E+04,
     &0.462197E+04,0.464830E+04,0.467468E+04,0.470127E+04,0.472792E+04,
     &0.475477E+04,0.478175E+04,0.480880E+04,0.483605E+04,0.486336E+04,
     &0.489087E+04,0.491853E+04,0.494623E+04,0.497415E+04,0.500215E+04,
     &0.503034E+04,0.505867E+04,0.508707E+04,0.511568E+04,0.514436E+04/
      DATA (ES(IES),IES=1236,1330) /
     &0.517325E+04,0.520227E+04,0.523137E+04,0.526068E+04,0.529005E+04,
     &0.531965E+04,0.534939E+04,0.537921E+04,0.540923E+04,0.543932E+04,
     &0.546965E+04,0.550011E+04,0.553064E+04,0.556139E+04,0.559223E+04,
     &0.562329E+04,0.565449E+04,0.568577E+04,0.571727E+04,0.574884E+04,
     &0.578064E+04,0.581261E+04,0.584464E+04,0.587692E+04,0.590924E+04,
     &0.594182E+04,0.597455E+04,0.600736E+04,0.604039E+04,0.607350E+04,
     &0.610685E+04,0.614036E+04,0.617394E+04,0.620777E+04,0.624169E+04,
     &0.627584E+04,0.631014E+04,0.634454E+04,0.637918E+04,0.641390E+04,
     &0.644887E+04,0.648400E+04,0.651919E+04,0.655467E+04,0.659021E+04,
     &0.662599E+04,0.666197E+04,0.669800E+04,0.673429E+04,0.677069E+04,
     &0.680735E+04,0.684415E+04,0.688104E+04,0.691819E+04,0.695543E+04,
     &0.699292E+04,0.703061E+04,0.706837E+04,0.710639E+04,0.714451E+04,
     &0.718289E+04,0.722143E+04,0.726009E+04,0.729903E+04,0.733802E+04,
     &0.737729E+04,0.741676E+04,0.745631E+04,0.749612E+04,0.753602E+04,
     &0.757622E+04,0.761659E+04,0.765705E+04,0.769780E+04,0.773863E+04,
     &0.777975E+04,0.782106E+04,0.786246E+04,0.790412E+04,0.794593E+04,
     &0.798802E+04,0.803028E+04,0.807259E+04,0.811525E+04,0.815798E+04,
     &0.820102E+04,0.824427E+04,0.828757E+04,0.833120E+04,0.837493E+04,
     &0.841895E+04,0.846313E+04,0.850744E+04,0.855208E+04,0.859678E+04/
      DATA (ES(IES),IES=1331,1425) /
     &0.864179E+04,0.868705E+04,0.873237E+04,0.877800E+04,0.882374E+04,
     &0.886979E+04,0.891603E+04,0.896237E+04,0.900904E+04,0.905579E+04,
     &0.910288E+04,0.915018E+04,0.919758E+04,0.924529E+04,0.929310E+04,
     &0.934122E+04,0.938959E+04,0.943804E+04,0.948687E+04,0.953575E+04,
     &0.958494E+04,0.963442E+04,0.968395E+04,0.973384E+04,0.978383E+04,
     &0.983412E+04,0.988468E+04,0.993534E+04,0.998630E+04,0.100374E+05,
     &0.100888E+05,0.101406E+05,0.101923E+05,0.102444E+05,0.102966E+05,
     &0.103492E+05,0.104020E+05,0.104550E+05,0.105082E+05,0.105616E+05,
     &0.106153E+05,0.106693E+05,0.107234E+05,0.107779E+05,0.108325E+05,
     &0.108874E+05,0.109425E+05,0.109978E+05,0.110535E+05,0.111092E+05,
     &0.111653E+05,0.112217E+05,0.112782E+05,0.113350E+05,0.113920E+05,
     &0.114493E+05,0.115070E+05,0.115646E+05,0.116228E+05,0.116809E+05,
     &0.117396E+05,0.117984E+05,0.118574E+05,0.119167E+05,0.119762E+05,
     &0.120360E+05,0.120962E+05,0.121564E+05,0.122170E+05,0.122778E+05,
     &0.123389E+05,0.124004E+05,0.124619E+05,0.125238E+05,0.125859E+05,
     &0.126484E+05,0.127111E+05,0.127739E+05,0.128372E+05,0.129006E+05,
     &0.129644E+05,0.130285E+05,0.130927E+05,0.131573E+05,0.132220E+05,
     &0.132872E+05,0.133526E+05,0.134182E+05,0.134842E+05,0.135503E+05,
     &0.136168E+05,0.136836E+05,0.137505E+05,0.138180E+05,0.138854E+05/
      DATA (ES(IES),IES=1426,1520) /
     &0.139534E+05,0.140216E+05,0.140900E+05,0.141588E+05,0.142277E+05,
     &0.142971E+05,0.143668E+05,0.144366E+05,0.145069E+05,0.145773E+05,
     &0.146481E+05,0.147192E+05,0.147905E+05,0.148622E+05,0.149341E+05,
     &0.150064E+05,0.150790E+05,0.151517E+05,0.152250E+05,0.152983E+05,
     &0.153721E+05,0.154462E+05,0.155205E+05,0.155952E+05,0.156701E+05,
     &0.157454E+05,0.158211E+05,0.158969E+05,0.159732E+05,0.160496E+05,
     &0.161265E+05,0.162037E+05,0.162811E+05,0.163589E+05,0.164369E+05,
     &0.165154E+05,0.165942E+05,0.166732E+05,0.167526E+05,0.168322E+05,
     &0.169123E+05,0.169927E+05,0.170733E+05,0.171543E+05,0.172356E+05,
     &0.173173E+05,0.173993E+05,0.174815E+05,0.175643E+05,0.176471E+05,
     &0.177305E+05,0.178143E+05,0.178981E+05,0.179826E+05,0.180671E+05,
     &0.181522E+05,0.182377E+05,0.183232E+05,0.184093E+05,0.184955E+05,
     &0.185823E+05,0.186695E+05,0.187568E+05,0.188447E+05,0.189326E+05,
     &0.190212E+05,0.191101E+05,0.191991E+05,0.192887E+05,0.193785E+05,
     &0.194688E+05,0.195595E+05,0.196503E+05,0.197417E+05,0.198332E+05,
     &0.199253E+05,0.200178E+05,0.201105E+05,0.202036E+05,0.202971E+05,
     &0.203910E+05,0.204853E+05,0.205798E+05,0.206749E+05,0.207701E+05,
     &0.208659E+05,0.209621E+05,0.210584E+05,0.211554E+05,0.212524E+05,
     &0.213501E+05,0.214482E+05,0.215465E+05,0.216452E+05,0.217442E+05/
      DATA (ES(IES),IES=1521,1551) /
     &0.218439E+05,0.219439E+05,0.220440E+05,0.221449E+05,0.222457E+05,
     &0.223473E+05,0.224494E+05,0.225514E+05,0.226542E+05,0.227571E+05,
     &0.228606E+05,0.229646E+05,0.230687E+05,0.231734E+05,0.232783E+05,
     &0.233839E+05,0.234898E+05,0.235960E+05,0.237027E+05,0.238097E+05,
     &0.239173E+05,0.240254E+05,0.241335E+05,0.242424E+05,0.243514E+05,
     &0.244611E+05,0.245712E+05,0.246814E+05,0.247923E+05,0.249034E+05,
     &0.250152E+05/
C
      DO 10 I=1,NPNTS
C
       IF (T(I) .LE. T_LOW) THEN
          QS(I) = EPSILON * ES(1) / P(I)
C
       ELSE IF (T(I) .GE. T_HIGH) THEN
          QS(I) = EPSILON * ES(N) / P(I)
C
       ELSE
C
          ATABLE = (T(I) - T_LOW + DELTA_T) / DELTA_T
          ITABLE = INT(ATABLE)
          ATABLE = ATABLE - ITABLE
C
          QS(I)  = EPSILON * ((1.-ATABLE) * ES(ITABLE)
     &            +     ATABLE  * ES(ITABLE+1)) / P(I)
C
       END IF
C
 10   CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE SATCAL-------------------------------------------------
CLL
CLL  PURPOSE : CALCULATES SATURATED TEMPERATURE
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY AUTUMN 1991
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE SATCAL (NPNTS,TH,PK,QS,THDDS,EXK,Q_K,THE_K)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
CHBG  include 'C_LHEAT.f'
CHBG  include 'C_R_CP.f'
CHBG  include 'C_0_DG_C.f'
      include 'UKPARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
C
      INTEGER I                 ! LOOP COUNTER
C
      INTEGER IC                ! LOOP COUNTER
C
      INTEGER NPNTS             ! VECTOR LENGTH
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL TH(ln2)              ! IN POTENTIAL TEMPERATURE (K)
C
      REAL PK(ln2)              ! IN PRESSURE OF LAYER K (PA)
C
      REAL Q_K(ln2)             ! IN MIXING RATIO OF LAYER K (KG/KG)
C
      REAL EXK(ln2)             ! IN EXNER RATIO OF LAYER K
C
      REAL THE_K(ln2)           ! IN ENVIRONMENTAL POTENTIAL TEMPERATURE
                                !    IN LAYER K
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL QS(ln2)              ! OUT SATURATED SPECIFIC HUMIDITY
                                !     (KG/KG)
C
      REAL THDDS(ln2)           ! OUT SATURATED ENVIRONMENTAL
                                !     POTENTIAL TEMPERATURE (K)
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE LOCALLY DEFINED
C-----------------------------------------------------------------------
C
      REAL T_FG(ln2)            ! TEMPERATURE FIRST GUESS (K)
C
      REAL TH_FG(ln2)           ! POTENTIAL TEMPERATURE FIRST GUESS (K)
C
      REAL DQBYDT(ln2)          ! FIRST GUESS AT MIXING RATIO INCREMENT
                                ! (KG/KG/K)
C
      REAL EL                   ! LATENT HEAT
C
C
C-----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C-----------------------------------------------------------------------
C
      EXTERNAL QSATU, DQS_DTH
C
C-----------------------------------------------------------------------
C SET INITIAL FIRST GUESS TEMPERATURE AND THETA - BASED UPON
C ENVIRONMENTAL TEMPERATURE IN LAYER K
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
       TH_FG(I) = THE_K(I)
       T_FG(I) = TH_FG(I)*EXK(I)
      END DO
C
C----------------------------------------------------------------------
C CALCULATE QSAT FOR INITIAL FIRST GUESS TEMPERATURE
C----------------------------------------------------------------------
C
      CALL QSATU(QS,T_FG,PK,NPNTS)
C
C----------------------------------------------------------------------
C DO TWO ITERATIONS TO FIND SATURATION POINT DUE TO EVAPORATION
C----------------------------------------------------------------------
C
      DO IC=1,2
C
C----------------------------------------------------------------------
C CALCULATE DQSAT/DT FOR FIRST GUESS TEMPERATURE
C----------------------------------------------------------------------
C
       CALL DQS_DTH(DQBYDT,TH_FG,QS,EXK,NPNTS)
C
C----------------------------------------------------------------------
C CALCULATE UPDATED TEMPERATURE AT SATURATION
C----------------------------------------------------------------------
C
       DO I=1,NPNTS
C
        IF (T_FG(I).GT.TM) THEN
         EL = LC
        ELSE
         EL = LCPLF
        END IF
C
        THDDS(I) = (TH(I) - (EL/(CP*EXK(I)))*(QS(I)-Q_K(I)-
     &                  TH_FG(I)*DQBYDT(I))) /
     &                  (1+(EL/(CP*EXK(I)))*DQBYDT(I))
C
C----------------------------------------------------------------------
C CALCULATE TEMPERATURE AT SATURATION AND UPDATE FIRST GUESS
C----------------------------------------------------------------------
C
        TH_FG(I) = THDDS(I)
        T_FG(I) = TH_FG(I)*EXK(I)
C
       END DO
C
C----------------------------------------------------------------------
C CALCULATE REVISED SATURATION MIXING RATIO AT SATURATION
C---------------------------------------------------------------------
C
       CALL QSATU(QS,T_FG,PK,NPNTS)
C
      END DO
C
      RETURN
      END
C
CLL  SUBROUTINE TERM_CON-----------------------------------------------
CLL
CLL  PURPOSE : RETURENS A MASK FOR POINTS AT WHICH CONVECTION
CLL            IS TERMINATING
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED: P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE TERM_CON(NPNTS,NLEV,K,BTERM,BWKP1,FLXKP1,THEKP1,QEKP1,
     &                    THPI,QPI,QSEKP1,DELTAK,EXKP1,EKP14,EKP34,
     &                    PSTAR)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
CHBG  include 'C_R_CP.f'
CHBG  include 'MPARFL.f'
CHBG  include 'XSBMIN.f'
CHBG  include 'C_LHEAT.f'
CHBG  include 'C_EPSLON.f'
CHBG  include 'QSTICE.f'
      include 'UKPARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS          ! IN VECTOR LENGTH
C
      INTEGER NLEV           ! IN NUMBER OF MODEL LAYER
C
      INTEGER K              ! IN PRESENT MODEL LAYER
C
      INTEGER I              ! LOOP COUNTER
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THEKP1(ln2)       ! IN POTENTIAL TEMPERATURE OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (K)
C
      REAL QEKP1(ln2)        ! IN MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL QSEKP1(ln2)       ! IN SATURATION MIXING RATIO OF CLOUD
                             !    ENVIRONMENT IN LAYER K+1 (KG/KG)
C
      REAL THPI(ln2)         ! IN INITIAL PARCEL POTENTIAL TEMPERATURE
                             !    (K)
C
      REAL QPI(ln2)          ! IN INITIAL PARCEL MIXING RATIO (KG/KG)
C
      REAL FLXKP1(ln2)       ! IN PARCEL MASSFLUX IN LAYER K+1 (PA/S)
C
      LOGICAL BWKP1(ln2)     ! IN MASK FOR WHETHER CONDENSATE IS
                             !    LIQUID IN LAYER K+1
C
      REAL DELTAK(ln2)       ! IN FORCED DETRAINMENT IN LAYER K
                             !    MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EXKP1(ln2)        ! IN EXNER RATIO FOR LEVEL K+1
C
      REAL EKP14(ln2)        ! IN ENTRAINMENT RATE FOR LEVEL K+1/4
                             !    MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL EKP34(ln2)        ! IN ENTRAINMENT RATE FOR LEVEL K+3/4
                             !    MULTIPLIED BY APPROPRIATE
                             !    LAYER THICKNESS
C
      REAL PSTAR(ln2)        ! IN SURFACE PRESSURE (PA)
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL BTERM(ln2)     ! OUT MASK OF THOSE POINTS AT WHICH
                             !     CONVECTION IS ENDING
C
C
C-----------------------------------------------------------------------
C VARIABLES THAT ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL EL                ! LATENT HEAT OF CONDENSATION OR
                             ! (CONDENSATION + FUSION) (J/KG)
C
      REAL FLXMIN            ! MINIMUM CONVECTIVE MASSFLUX BELOW
                             ! WHICH TERMINAL DETRAINMENT OCCURS
                             ! (PA/S)
C
      REAL THVUNDI           ! POTENTIAL TEMPERATURE OF AN
                             ! UNDILUTE PARCEL IN LAYER K+1
                             ! FROM THE STARTING LAYER OF
                             ! CONVECTION (K)
C
      REAL THVEKP1           ! VIRTUAL POTENTIAL TEMPERATURE
                             ! OF ENVIRONMENT IN LAYER K+1 (K)
C
C*---------------------------------------------------------------------
C
C----------------------------------------------------------------------
C  CALCULATE MINIMUM MASS FLUX BELOW WHICH CONVECTION IS TERMINATED
C----------------------------------------------------------------------
C
      DO 10 I=1,NPNTS
        FLXMIN = MPARFL*(1.+EKP14(I))*(1.+EKP34(I))*PSTAR(I)
C
C-----------------------------------------------------------------------
C   CREATE A VECTOR OF LATENT HEATS
C-----------------------------------------------------------------------
C
       IF (BWKP1(I)) THEN
          EL = LC
       ELSE
          EL = LCPLF
       ENDIF
CL
CL----------------------------------------------------------------------
CL  PARCELS ARE ONLY CONSIDERED FOR TERMINATION IF THEY ARE DETRAINING
CL  EXCEPT AT THE TOP MODEL LAYER, WHERE ALL CONVECTION TERMINATES
CL
CL  IF THE PARCEL HAS A POTENTIAL TEMPETURE GREATER THAN THE
CL  POTENTIAL TEMPERATURE OF AN UNDILUTE PARCEL FORM THE STARTING
CL  LAYER OF CONVECION IN LAYER K+1 THEN CONVECTION IS TERMINATED
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (7), EQUATION (32)
CL
CL  CONVECTION IS ALSO TERMINATED IF MASS FLUX IN LAYER K+1 IS LESS
CL  IS LESS THAN A MINIMUM VALUE
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (7), EQUATION (33)
CL----------------------------------------------------------------------
CL
       THVUNDI=( THPI(I) + (EL/(EXKP1(I)*CP)) *(QPI(I) - QSEKP1(I))
     &         +((LC-EL)/(EXKP1(I)*CP))*MAX(0.0,(QPI(I)-QSTICE))
     &         )*(1.+C_VIRTUAL*QSEKP1(I))
C
       THVEKP1 = (THEKP1(I)*(1.+C_VIRTUAL*QEKP1(I)) + XSBMIN)
C
       BTERM(I) = (((FLXKP1(I) .LT. FLXMIN) .OR. (THVUNDI .LT. THVEKP1))
     &               .AND. (DELTAK(I).GT.0.0)) .OR. (K+1) .EQ. NLEV
C
  10  CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE TERMDD-------------------------------------------------
CLL
CLL  PURPOSE : CALCULATE WHETHER DOWNDRAUGHT IS ABLE TO CONTINUE
CLL
CLL            CALCULATE BUOYANCY
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE WRITTEN FOR CRAY Y-MP BY S.BETT AND D.GREGORY
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 3
CLL  VERSION NO. 4  DATED 5/2/92
CLL
CLL  SYSTEM TASK : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE TERMDD (NPNTS,BDD_START,THDD_K,QDD_K,THE_K,QE_K,K,
     &                   B_DD_END,BDD_ON)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS USED IN THIS ROUTINE
C-----------------------------------------------------------------------
C
CHBG  include 'C_EPSLON.f'
CHBG  include 'DDKMDET.f'
      include 'UKPARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS                ! IN VECTOR LENGTH
C
      INTEGER I                    ! LOOP COUNTER
C
      INTEGER K                    ! IN PRESENT MODEL LAYER
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THDD_K(ln2)             ! IN MODEL POTENTIAL TEMPERATURE
                                   !    OF DOWNDRAUGHT AT LAYER K (K)
C
      REAL QDD_K(ln2)              ! IN MODEL MIXING RATIO OF
                                   !    DOWNDRAUGHT AT LAYER K
C
      REAL THE_K(ln2)              ! IN POTENTIAL TEMPERATURE OF
                                   !    ENVIRONMENTAL AIR IN LAYER K
C
      REAL QE_K(ln2)               ! IN MODEL MIXING RATIO AT LAYER K
C
      LOGICAL BDD_START(ln2)       ! IN MASK FOR THOSE POINTS WHERE
                                   !    DOWNDRAUGHT MAY OCCUR IN
                                   !    LAYER K-1
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      LOGICAL B_DD_END(ln2)        ! OUT MASK FOR THOSE POINTS WHERE
                                   !     DOWNDRAUGHT IS TERMINATING
C
      LOGICAL BDD_ON(ln2)          ! OUT MASK FOR THOSE POINTS WHERE
                                   !     DOWNDRAUGHT CONTINUES TO LAYER
                                   !     K-1 (AS BDD_START HERE)
C
C-----------------------------------------------------------------------
C VARIABLES WHICH ARE DEFINED LOCALLY
C-----------------------------------------------------------------------
C
      REAL BUOY1                   ! BUOYANCY OF PARCEL
C
      REAL THDD_V                  ! USED IN CALCULATION OF BUOYANCY
C
      REAL THE_V                   ! USED IN CALCULATION OF BUOYANCY
C
C-----------------------------------------------------------------------
C CHECK IF PARCEL STILL NEGATIVELY BUOYANT SUCH THAT DOWNDRAUGHT
C CAN CONTINUE TO NEXT LAYER
C-----------------------------------------------------------------------
C
      DO I=1,NPNTS
         THDD_V = THDD_K(I)*(1.0+C_VIRTUAL*QDD_K(I))
         THE_V = THE_K(I)*(1.0+C_VIRTUAL*QE_K(I))
         BUOY1 = THDD_V - THE_V
C
C-----------------------------------------------------------------------
C CALCULATE STATE OF DOWNDRAUGHT
C-----------------------------------------------------------------------
C
         IF (BDD_START(I) .AND. BUOY1.GT.0.5) THEN
            BDD_ON(I) = .FALSE.
         ELSE IF (BUOY1.GT.0.5 .OR. K.EQ.2) THEN
            B_DD_END(I) = .TRUE.
         END IF
      END DO
C
      RETURN
      END
C
CLL  SUBROUTINE THETAR-------------------------------------------------
CLL
CLL  PURPOSE : CALCULATES THE POTENTIAL TEMPERATURE OF THE DETRAINING
CLL            AIR IN LAYER K AND ALSO THE DIFFERENCE IN THE
CLL            WATER VAPOUR CONTENT OF THE DETRAINING AIR FROM THAT
CLL            OF THE MEAN PARCEL IN LAYER K
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED : P27
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE THETAR (NPNTS,THRK,QRK,XSQR,BGMK,THEK,QEK,QPK,QSEK,
     &                   DQSK,EXK,PK)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C----------------------------------------------------------------------
C MODEL CONSTANTS
C----------------------------------------------------------------------
C
CHBG  include 'C_EPSLON.f'
      include 'UKPARAMS.f'
C
C----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C----------------------------------------------------------------------
C
      INTEGER NPNTS          ! VECTOR LENGTH
C
      INTEGER I              ! LOOP COUNTER
C
C
C----------------------------------------------------------------------
C VARIABLES THAT ARE INPUT
C----------------------------------------------------------------------
C
      REAL THEK(ln2)         ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                             !    IN LAYER K (K)
C
      REAL QEK(ln2)          ! IN ENVIRONMENT MIXING RATIO
                             !    IN LAYER K (KG/KG)
C
      REAL QPK(ln2)          ! IN PARCEL MIXING RATIO IN LAYER K
                             !    (KG/KG)
C
      REAL QSEK(ln2)         ! IN SATURATION MIXING RATIO OF THE
                             !    ENVIRONMENT IN LAYER K (KG/KG)
C
      REAL DQSK(ln2)         ! IN GRADIENT OF SATURATION MIXING RATIO
                             !    WITH POTENTIAL TEMPERATURE FOR THE
                             !    ENVIRONMENT OF LAYER K (KG/KG/K)
C
      LOGICAL BGMK(ln2)      ! IN MASK FOR PARCELS SATURATED IN LAYER K
C
      REAL EXK(ln2)          ! IN EXNER RATIO FOR LEVEL K
C
      REAL PK(ln2)           ! IN PRESSURE AT LEVEL K (PA)
C
C
C----------------------------------------------------------------------
C VARIABLES THAT ARE OUTPUT
C----------------------------------------------------------------------
C
      REAL THRK(ln2)         ! OUT PARCEL DETRAINMENT POTENTIAL
                             !     TEMPERATURE IN LAYER K (K)
C
      REAL QRK(ln2)          ! OUT PARCEL DETRAINMENT MIXING RATIO
                             !     IN LAYER K (KG/KG)
C
      REAL XSQR(ln2)         ! OUT EXCESS WATER VAPOUR OF
                             !     DETRAINING AIR (KG/KG)
C
C
C----------------------------------------------------------------------
C VARIABLES THAT ARE DEFINED LOCALLY
C
C ON THE IBM ARRAYS ARE ALLOCATED USING A PARAMETER STATEMENT
C
C ON THE CRAY ARRAYS ARE DYNAMICALLY ALLOCATED
C----------------------------------------------------------------------
C
      REAL TT(ln2)           ! TEMPORARY TEMPERATURE FOR CALCULATION
                             ! OF SATURATION MIXING RATIO (K)
C
C
C----------------------------------------------------------------------
C EXTERNAL ROUTINES CALLED
C----------------------------------------------------------------------
C
      EXTERNAL QSATU
C
C*----------------------------------------------------------------------
C
      DO 20 I=1,NPNTS
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE POTENTIAL TEMPERATURE OF DETRAINING AIR
CL
CL  UM DOCUMENTATION PAPER P27
CL  SECTION (6), EQUATION (26)
CL----------------------------------------------------------------------
CL
       IF (.NOT.BGMK(I)) THEN
          THRK(I)=THEK(I) * (1. + C_VIRTUAL*QEK(I)) /
     &                 (1. + C_VIRTUAL*QPK(I))
       ELSE
          THRK(I) = THEK(I)*(1.0 + C_VIRTUAL*(QEK(I)-QSEK(I))/
     &                 (1.0 + C_VIRTUAL*THEK(I)*DQSK(I)))
       ENDIF
CL
CL----------------------------------------------------------------------
CL  CALCULATE THE MIXING RATIO OF THE DETRAINING AIR AIR THE
CL  DIFFERENCE BETWEEN THIS AND THE MIXING RATIO OF THE MEAN
CL  PARCEL IN LAYER K
CL
CL  THE MOISTURE DIFFERENCE IS USED TO CALCULATE THE
CL  COND_DET_K TERM OF EQUATION (30), SECTION (6),
CL  UM DOCUMENTATIONM PAPER P27
CL----------------------------------------------------------------------
CL
C
C-----------------------------------------------------------------------
C CONVERT POTENTIAL TEMPERATURE TO TEMPERATURE AND CALCULATE
C PRESSURE OF LAYER K FOR CALCULATION OF SATURATED
C MIXING RATIO
C-----------------------------------------------------------------------
C
       TT(I) = THRK(I)*EXK(I)
   20  CONTINUE
      CALL QSATU (XSQR,TT,PK,NPNTS)
C
      DO 30 I=1,NPNTS
       IF (BGMK(I)) THEN
          QRK(I)  = XSQR(I)
          XSQR(I) = QPK(I) - XSQR(I)
       ELSE
          QRK(I)  = QPK(I)
          XSQR(I) = 0.
       ENDIF
   30  CONTINUE
C
      RETURN
      END
CLL  SUBROUTINE THP_DET------------------------------------------------
CLL
CLL  SUITABLE FOR SINGLE COLUMN MODEL USE
CLL
CLL  CODE REWORKED FOR CRAY Y-MP BY D.GREGORY AUTUMN/WINTER 1989/90
CLL
CLL  MODEL            MODIFICATION HISTORY FROM MODEL VERSION 3.0:
CLL VERSION  DATE
CLL
CLL  PROGRAMMING STANDARDS : UNIFIED MODEL DOCUMENTATION PAPER NO. 4
CLL  VERSION NO. 1
CLL
CLL  LOGICAL COMPONENTS COVERED : P27
CLL
CLL  PURPOSE : CALCULATES POTENTIAL TEMPERATURE OF THE
CLL            PARCEL IN LAYER K+1 AFTER FORCED DETRAINMENT
CLL            IN LAYER K
CLL
CLL  DOCUMENTATION : UNIFIED MODEL DOCUMENTATION PAPER P27
CLL                  SECTION (6), EQUATION (28)
CLL
CLLEND-----------------------------------------------------------------
C
C*L  ARGUMENTS---------------------------------------------------------
C
      SUBROUTINE THP_DET (NPNTS,THPKP1,THEKP1,QPKP1,QEKP1,QSEKP1,
     &                    DQSKP1,BGMKP1,BCALC)
C
      IMPLICIT NONE
      include 'PARAMS.f'
C
C-----------------------------------------------------------------------
C MODEL CONSTANTS
C-----------------------------------------------------------------------
C
CHBG  include 'C_EPSLON.f'
CHBG  include 'XSBMIN.f'
      include 'UKPARAMS.f'
C
C-----------------------------------------------------------------------
C VECTOR LENGTHS AND LOOP COUNTERS
C-----------------------------------------------------------------------
C
      INTEGER NPNTS           ! IN VECTOR LENGTH
C
      INTEGER I               ! LOOP COUNTER
C
C
C-----------------------------------------------------------------------
C VARAIBLES WHICH ARE INPUT
C-----------------------------------------------------------------------
C
      REAL THEKP1(ln2)        ! IN ENVIRONMENT POTENTIAL TEMPERATURE
                              !    IN LAYER K+1 (K)
C
      REAL QPKP1(ln2)         ! IN PARCEL MIXING RATIO IN LAYER K+1
                              !    (KG/KG)
C
      REAL QSEKP1(ln2)        ! IN ENVIRONMENT SATURATED MIXING RATIO
                              !    IN LAYER K+1 (KG/KG)
C
      REAL DQSKP1(ln2)        ! IN GRADIENT OF SATURATION MIXING RATIO
                              !    POTENTIAL TEMPERATURE FOR THE
                              !    ENVIRONMENT IN LAYER K+1 (KG/KG/K)
C
      REAL QEKP1(ln2)         ! IN ENVIRONMENT MIXING RATIO IN
                              !    LAYER K+1 (KG/KG)
C
      LOGICAL BGMKP1(ln2)     ! IN MASK FOR PARCELS WHICH ARE SATURATED
                              !    IN LAYER K+1
C
      LOGICAL BCALC(ln2)      ! IN MASK FOR PARCELS AT WHICH
                              !    CALCULATIONS OF THIS SUBROUTINE ARE
                              !    TO BE CARRIED OUT
C
C
C-----------------------------------------------------------------------
C VARAIBLES WHICH ARE OUTPUT
C-----------------------------------------------------------------------
C
      REAL THPKP1(ln2)        ! OUT PARCEL POTENTIAL TEMPERATURE
                              !     IN LAYER K+1 AFTER FORCED
                              !     DETRAINMENT (K)
C
C*---------------------------------------------------------------------
CL
CL---------------------------------------------------------------------
CL  NO SIGNIFICANT STRUCTURE
CL---------------------------------------------------------------------
CL
C
      DO 10 I=1,NPNTS
        IF (BCALC(I))THEN
         IF (BGMKP1(I)) THEN
           THPKP1(I) = THEKP1(I) +
     &                (C_VIRTUAL*THEKP1(I)*
     &                           (QEKP1(I)-QSEKP1(I)) + XSBMIN)
     &               /( 1. + C_VIRTUAL*THEKP1(I)*DQSKP1(I) )
C
         ELSE
           THPKP1(I) = (THEKP1(I)*(1. + C_VIRTUAL*QEKP1(I))
     &                                                        + XSBMIN)
     &                    /(1. + C_VIRTUAL*QPKP1(I))
         END IF
        END IF
  10  CONTINUE
C
      RETURN
      END
