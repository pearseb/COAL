CUKMO ----------------------------------------------------
CUKMO FILE CLDAREA.f
CUKMO ----------------------------------------------------
      REAL CLDAREA                ! FRACTIONAL CLOUD AREA WHEN NO
                                  ! DD, USED IN EVAPORATION CALC
      PARAMETER (CLDAREA=1.0)
C
CUKMO ----------------------------------------------------
CUKMO FILE CRITDEP.f
CUKMO ----------------------------------------------------
      REAL CRITDSEA  !  CRITICAL DEPTH OF CLOUD FOR THE FORMATION OF
                     !  CONVECTIVE PRECIPITATION OVER SEA (M)
      PARAMETER (CRITDSEA = 1.5E3)
C
      REAL CRITDLND  !  CRITICAL DEPTH OF CLOUD FOR THE FORMATION OF
                     !  CONVECTIVE PRECIPITATION OVER LAND (M)
      PARAMETER (CRITDLND = 4.0E3)
C
      REAL CRITDICE  !  CRITICAL DEPTH OF A GLACIATED CLOUD FOR THE
                     !  FORMATION OF CONVECTIVE PRECIPITATION (M)
      PARAMETER (CRITDICE = 1.0E3)
C
CUKMO ----------------------------------------------------
CUKMO FILE C_0_DG_C.f
CUKMO ----------------------------------------------------
C*L------------------COMDECK C_O_DG_C-----------------------------------
C ZERODEGC IS CONVERSION BETWEEN DEGREES CELSIUS AND KELVIN
C TFS IS TEMPERATURE AT WHICH SEA WATER FREEZES
C TM IS TEMPERATURE AT WHICH FRESH WATER FREEZES AND ICE MELTS
      REAL ZERODEGC,TFS,TM

      PARAMETER(ZERODEGC=273.15,
     &          TFS=271.35,
     &          TM=273.15)
C*----------------------------------------------------------------------

CUKMO ----------------------------------------------------
CUKMO FILE C_EPSLON.f
CUKMO ----------------------------------------------------
C*L------------------COMDECK C_EPSLON-----------------------------------
C EPSILON IS RATIO OF MOLECULAR WEIGHTS OF WATER AND DRY AIR
      REAL EPSILON,C_VIRTUAL

CHBG  PARAMETER(EPSILON=0.62198,
CHBG &          C_VIRTUAL=1./EPSILON-1.)
CHBG Changed to suit CSIRO model
      PARAMETER(EPSILON=0.622,
     &          C_VIRTUAL=1./EPSILON-1.)
C*----------------------------------------------------------------------

CUKMO ----------------------------------------------------
CUKMO FILE C_G.f
CUKMO ----------------------------------------------------
C*L------------------COMDECK C_G----------------------------------------
C G IS MEAN ACCEL DUE TO GRAVITY AT EARTH'S SURFACE
      REAL G
CHBG  PARAMETER(G=9.80665)
      PARAMETER(G=9.80616)
C*----------------------------------------------------------------------
CUKMO ----------------------------------------------------
CUKMO FILE C_LHEAT.f
CUKMO ----------------------------------------------------
C*L------------------COMDECK C_LHEAT------------------------------------
C LC IS LATENT HEAT OF CONDENSATION OF WATER AT 0DEGC
C LF IS LATENT HEAT OF FUSION AT 0DEGC
      REAL LC,LF,LCPLF

CHBG  PARAMETER(LC=2.501E6,
CHBG &          LF=0.334E6)
CHBG Changed to suit CSIRO model : LC=hl=el; LF=hlfusion=hlf
      PARAMETER(LC=2.5E6,LF=3.35E5)
      PARAMETER(LCPLF=LC+LF)
C*----------------------------------------------------------------------
CUKMO ----------------------------------------------------
CUKMO FILE C_R_CP.f
CUKMO ----------------------------------------------------
C*L------------------COMDECK C_R_CP-------------------------------------
C R IS GAS CONSTANT FOR DRY AIR
C CP IS SPECIFIC HEAT OF DRY AIR AT CONSTANT PRESSURE
C PREF IS REFERENCE SURFACE PRESSURE
      REAL R,CP,KAPPA,PREF  ! PREFTOKI

CHBG  PARAMETER(R=287.05,
CHBG &          CP=1005.,
CHBG &          KAPPA=R/CP,
CHBG &          PREF=100000.)
CHBG Changed to suit CSIRO model : KAPPA=cappa; PREF=P00*100.0
      PARAMETER(R=287.04,
     &          CP=1004.64,
     &          KAPPA=R/CP,
     &          PREF=100000.)
C*----------------------------------------------------------------------

CUKMO ----------------------------------------------------
CUKMO FILE DDAREA.f
CUKMO ----------------------------------------------------
      REAL DDCLDFRA               ! FRACTIONAL CLOUD AREA OF DD
      PARAMETER (DDCLDFRA=0.5)
C
CUKMO ----------------------------------------------------
CUKMO FILE DDEVAP.f
CUKMO ----------------------------------------------------
      REAL P_LQ1,P_LQ2   ! EXPONENTS USED IN CALCULATION OF
                         ! EVAPORATION OF LIQUID
      PARAMETER (P_LQ1 = 0.42, P_LQ2 = 0.6)
C
      REAL P_ICE1,P_ICE2 ! EXPONENTS USED IN CALCULATION OF
                         ! EVAPORATION OF ICE
      PARAMETER (P_ICE1 = 0.55, P_ICE2 = 0.76)
C
      REAL RHO_LQP1,RHO_LQP2, ! EXPONENTS AND CONSTANTS ASSOCIATED
     &     RHO_LQA,RHO_LQB    ! WITH DENSITY TERM IN EVAPORATION
                              ! OF LIQUID
      PARAMETER (RHO_LQP1=0.21, RHO_LQP2=0.55, RHO_LQA=67.08,
     &           RHO_LQB=541.16)
C
      REAL RHO_ICP1,RHO_ICP2, ! EXPONENTS AND CONSTANTS ASSOCIATED
     &     RHO_ICEA,RHO_ICEB  ! WITH DENSITY TERM IN EVAPORATION
                              ! OF ICE
      PARAMETER (RHO_ICP1=0.28, RHO_ICP2=0.63, RHO_ICEA=1569.52,
     &           RHO_ICEB=32069.02)
C
CUKMO ----------------------------------------------------
CUKMO FILE DDEVPICE.f
CUKMO ----------------------------------------------------
      REAL ICE_A,ICE_B,ICE_C ! CONSTANTS USED IN QUADRATIC FORMULA
                             ! FOR EVAPORATION OF ICE
      PARAMETER ( ICE_A = -5.2E-9, ICE_B = 2.5332E-6, ICE_C = -2.911E-4)
C
CUKMO ----------------------------------------------------
CUKMO FILE DDEVPLQ.f
CUKMO ----------------------------------------------------
      REAL LQ_A,LQ_B,LQ_C ! CONSTANTS USED IN QUADRATIC FORMULA
                          ! FOR EVAPORATION OF LIQUID
      PARAMETER ( LQ_A = 2.008E-9, LQ_B = -1.385E-6, LQ_C = 2.424E-4)
C
CUKMO ----------------------------------------------------
CUKMO FILE DDKMDET.f
CUKMO ----------------------------------------------------
      REAL DET_LYR  ! THICKNESS LEVEL USED IN CALCULATION OF MIXING
                    ! DETRAINMENT FOR DOWNDRAUGHT  (PA)
      PARAMETER (DET_LYR = 10000.0)
C
CUKMO ----------------------------------------------------
CUKMO FILE DELTHST.f
CUKMO ----------------------------------------------------
      REAL DELTHST  !  DIFFERENCE IN POTENTIAL TEMPERATURE BETWEEN
                    !  LEVELS ABOVE WHICH THE ATMOSPHERE IF ASSUMED
                    !  TO BE TOO STABLE TO CONVECT (K)
      PARAMETER (DELTHST = 1.5)
C
CUKMO ----------------------------------------------------
CUKMO FILE ENTCNST.f
CUKMO ----------------------------------------------------
      REAL AE1,AE2,       ! COEFFICIENTS USED IN CALCULATION
     &     ENTCOEF        ! OF ENTRAINMENT RATE
      PARAMETER(AE1=1.0,AE2=1.5)
      PARAMETER (ENTCOEF = 3.0)
C
CUKMO ----------------------------------------------------
CUKMO FILE ENTDD.f
CUKMO ----------------------------------------------------
      REAL DDCOEF1, ! COEFFICIENTS USED IN CALCULATION OF DOWNDRAUGHT
     &     DDCOEF2  ! ENTRAINMENT RATES
      PARAMETER (DDCOEF1 = 1.8E6, DDCOEF2 = 3.0)
C
CUKMO ----------------------------------------------------
CUKMO FILE MASSFC.f
CUKMO ----------------------------------------------------
      REAL C,D  ! CONSTANTS USED TO DETERMINE THE INITIAL CONVECTIVE
                ! MASS FLUX FROM PARCEL BUOYANCY
                ! MI = 1.E-3 * (D + C*BUOYANCY/DELP)
       PARAMETER ( C = 5.17E-4, D = 0.0 ) ! UKMO value
       REAL C_R21,C_T63 ! C value for 2 CSIRO model resolutions
c       PARAMETER ( C_R21 = 1.00E-4 )
       PARAMETER ( C_R21 = 3.00E-4 )
       PARAMETER ( C_T63 = 3.00E-4 )
C
CUKMO ----------------------------------------------------
CUKMO FILE MPARB.f
CUKMO ----------------------------------------------------
      REAL MPARB  !  MINIMUM (PARCEL BUOYANCY/LAYER THICKNESS) (K/PA)
      PARAMETER (MPARB = 1.0)
C
CUKMO ----------------------------------------------------
CUKMO FILE MPARFL.f
CUKMO ----------------------------------------------------
      REAL MPARFL  !  MINIMUM PARCEL MASS FLUX
                   !  = 1E-3 * MINIMUM PARCEL BUOYANCY *
                   !              MASS FLUX PARAMETER C
      PARAMETER (MPARFL = 1.0E-3 * 1.0 * 3.33E-4)
C
CUKMO ----------------------------------------------------
CUKMO FILE MPARWTR.f
CUKMO ----------------------------------------------------
      REAL MPARWTR  !  MINIMUM PARCEL WATER IN GRAMS PER KILOGRAM
                    !  BEFORE PRECIPITATION IS ALLOWED (KG/KG)
CHBG  PARAMETER (MPARWTR = 1.0E-3)
      PARAMETER (MPARWTR = 0.5E-3)
C
CUKMO ----------------------------------------------------
CUKMO FILE PARXS.f
CUKMO ----------------------------------------------------
      REAL THPIXS, QPIXS ! INITIAL EXCESS POTENTIAL TEMPERATURE (K)
                         ! AND MIXING RATIO (KG/KG)
       PARAMETER ( THPIXS = 0.2, QPIXS = 0.0 )
C
CUKMO ----------------------------------------------------
CUKMO FILE QSTICE.f
CUKMO ----------------------------------------------------
      REAL QSTICE   !  APPROXIMATION TO SATURATION MIXING RATIO
                    !  AT TEMPERATURE AT WHICH LIQUID WATER TURNS TO
                    !  ICE (SEE COMDECK TICE) (KG/KG)
c      PARAMETER (QSTICE = 3.5E-3)
      PARAMETER (QSTICE = 1.1E-3) !This is for Tice = 258.15 K
C
CUKMO ----------------------------------------------------
CUKMO FILE RV.f
CUKMO ----------------------------------------------------
      REAL RV  !  GAS CONSTANT FOR WATER VAPOUR (J/KG/K)
CHBG  PARAMETER ( RV = 461.1 )
CHBG Changed to suit CSIRO model : RV=ars
      PARAMETER ( RV = 461.0 )
C
CUKMO ----------------------------------------------------
CUKMO FILE TICE.f
CUKMO ----------------------------------------------------
      REAL TICE  ! TEMPERATURE IN KELVIN AT WHICH LIQUID WATER TURNS
                 ! TO ICE IN THE CONVECTIVE DOWNDRAUGHT SCHEME (K)
c      PARAMETER (TICE = 273.15)
      PARAMETER (TICE = 258.15) !If changing this, also change TICEuk in progcld.f (LDR, 3/01)
C
CUKMO ----------------------------------------------------
CUKMO FILE XSBMIN.f
CUKMO ----------------------------------------------------
      REAL XSBMIN !  MINIMUM EXCESS BUOYANCY TO CONTINUE PARCEL ASCENT
                  !  (K)
      PARAMETER (XSBMIN = 0.2)
C
