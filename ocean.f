# 1 "ocean.F"
c Added bulk_force to contrl namelist group
c Pearse Buchanan 16/07/2018
c
c (1) Removing the ITM=1 option from the ocean model source code.
c (2) Updating the default configuration of the high-resolution version of the
c     ocean model.
c SJP 2009/06/23
c
c End-of-month work for the ocean model transferred from OCEND to a new
c subroutine, OCEMON.
c SJP 2009/05/11
c
c Further modifications to the ocean model output routines to improve memory
c usage. The features that were causing stack overflows have now been resolved.
c SJP 2009/05/11
c
c Modifying the ocean model output routines for improved memory usage.
c SJP 2009/05/06
c
c (1) Adding density as an ocean model statistic.
c (2) Enhancing user control over the coupling between the atmosphere and
c     ocean.
c SJP 2009/04/21
c
c Enabling user control over the ocean model statistics that are saved to file.
c SJP 2009/04/19
c
c Incorporating the code which calculates the meridional overturning
c streamfunctions into the core model source code.
c SJP 2009/04/17
c
c Add calls to NINT to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/15
c
c Implementing the equation of state of McDougall et al (2003).
c SJP 2009/04/08
c
c (1) Implementing a Robert time filter.
c (2) Adding code to check for conservation of tracers within the ocean.
c SJP 2008/12/17
c
c Increase the default value of MXSCAN for the new, high-resolution version of
c the model from 80 -> 1000.
c SJP 2008/03/09
c
c (1) Modified for the conversion of the ocean model restart file to netCDF.
c (2) Removed one array bounds violation.
c SJP 2008/03/08
c
c (1) Add the NAMELIST variable COS_AMF. When this is TRUE, the horizontal
c     viscosity varies as the cosine of the latitude. Otherwise, a constant
c     value is used (except when using the old ocean model grid, in which case
c     the values are reduced over the high-latitude Arctic Ocean).
c (2) Add the NAMELIST variable NO_GM_ARCTIC. When this is TRUE,
c     Gent-McWilliams eddy diffusion is replaced with horizontal diffusion
c     within the Fourier-filtered region of the Arctic Ocean. Otherwise, eddy
c     diffusion is used at all latitudes.
c (3) Add default values for all NAMELIST variables, for both the old and new
c     ocean model grids.
c (4) Add code which displays the values of the isopycnal, eddy and horizontal
c     diffusivities on each model level, as well as the value of the horizontal
c     viscosity at each latitude.
c (5) Fix a bug in the section of code which checks the values of the
c     permuting disk indicators.
c SJP 2008/02/27
c
c Modify such that, when using the new ocean model grid, the horizontal
c viscosity now varies as the cosine of the latitude across the entire globe.
c SJP 2007/11/25
c
c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Rationalise the diagnostic information written to standard output.
c SJP 2007/11/21
c
c (1) Tidied up the sections of code which set the values of AHIFAC and AMFAC.
c     When using the new ocean model grid, the horizontal viscosity now varies
c     as the cosine of the latitude over the Arctic Ocean.
c (2) Added code which prints the parameter values specified at compile time.
c SJP 2007/07/06
c
c Subroutine OCCOUP1 split off into a separate source file.
c SJP 2007/06/20
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/18
c
c Modified to use OCEAN_LOW preprocessor macro, rather than OCEAN_DBL.
c SJP 2007/06/04
c
c (1) Modified to enable the new ocean model grid.
c (2) Removed redundant references to the array HSUM and the header files
c     ICEOC.f and DEPL.f.
c As a result of these changes, the model now only reads from the topography
c file sttop.bot_ind when it is initialised from scratch. Otherwise, it obtains
c the topography information from the restart file.
c SJP 2007/06/01
c
c Added the line "include 'PARAMS.f'" to subroutines OCEAN and OCCOUP1, as this
c line is no longer included via the header file OPARAMS.f.
c SJP 2007/05/31
c
c Added IMPLICIT NONE statements, plus variable declarations as necessary.
c SJP 2007/05/30
c
c Transferred COMMON block declarations to separate header files.
c SJP 2007/05/29
c
c Modified to enable relaxation of the coupled model SSTs and SSSs towards
c prescribed values.
c SJP 2006/01/05
c
c Added TRELAX to /PARMS/, which contains the relaxation timescale for the
c stand-alone OGCM in days. Initialise TRELAX to 20.0, so that the default
c behaviour of the model is identical to previously. Also create a new COMMON
c block /TRELAX/, to hold TRELAX.
c SJP 2005/03/19
c
c Re-inserted the zeroing of wind stress adjustment arrays TXCOR and TYCOR.
c SJP 2004/11/18
c
c Modify AHEFAC and AHHFAC again, so that eddy diffusion extends all the way to
c the North Pole. There is now no change to horizontal diffusion whatsoever.
c SJP 2004/09/13
c
c Modify AHEFAC and AHHFAC in the same manner as Dave Bi, causing the latitude
c at which eddy diffusion is replaced with horizontal diffusion to be shifted
c northward from J=53 [74.8 degN] to J=55 [81.2 degN].
c 2004/03/30
c
c (1) Removed redundant code.
c (2) Changed default value of NA from 0 to 1.
c (3) Added default values for IGM, ITM, IYES, NDIAG, ITDB, URAT, DIFRAT,
c     ITSET and DTXF.
c SJP 2004/02/14
c
c (1) Redundant code removed from OCCOUP1.
c (2) Remove unnecessary variable declarations.
c SJP 2004/01/06
c
c Modified so that the Gent-McWilliams eddy parameterisation scheme is
c always used, enabling some optimisation. The parameter IGM is thus rendered
c irrelevant.
c SJP 2003/12/18
c
c (1) Modified for changes to /orestart/.
c (2) Commented out zeroing of wind stress adjustment arrays TXCOR and TYCOR.
c (3) Writes to Fortran unit 32 by coupled model commented out.
c SJP 2003/09/04
c
c Replaced calls to OGET/OPUT with data access via COMMON block /orestart/.
c SJP 2003/09/02
c
c Fix up calls to OGET and OPUT.
c SJP 2003/05/01
c
c Calls to OSTART removed, as they did nothing. Parameter definitions moved
c to include file OPARAMS.f.
c SJP 2003/04/29
c
c Removed write to fort.1, as this file is now redundant.
c SJP 2002/02/15
c
c Correction to calculations of KNITD and NLAST, to ensure that rounding
c does not result in them having the wrong value.
c SJP 2002/01/11
c
c $Log: ocean.f,v $
c Revision 1.21  2000/06/20 02:08:33  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.20  1998/12/10 00:55:42  ldr
c HBG changes to V5-1-21
c
c Revision 1.19  1998/05/27  02:10:17  ldr
c Merge TIE and ACH changes.
c
c Revision 1.18  1998/05/26  05:10:49  ldr
c Final 10k run changes (mainly ACH salinity stuff) merged into V5-1-9.
c
c Revision 1.17.1.1  1998/05/27  02:07:37  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.17  1997/12/23  01:03:49  ldr
c REmove spurious <<<< signs from merging.
c
c Revision 1.16  1997/12/23  00:23:36  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.15.1.1  1998/05/26  04:48:56  ldr
c Final changes (inc. ACH salinity stuff) for 10k run.
c
c Revision 1.15  1997/12/19  01:25:38  ldr
c Changes from ACH for 21 OGCM levels and optinal GM mixing scheme...
c
c Include array for storing depth of convective mixed layer
c Insert code to initialize arrays for uedd, vedd, wedd accumulation
c Include new array edplt(i,j,k,3)
c Change to 21 levels in ocean model, and insertion of
c eddy-induced transport (major changes delineated)
c
c Revision 1.14.1.1  1997/12/19  02:03:15  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.14  1997/03/06  23:33:59  ldr
c Corrections from HBG to the recent tidy-ups to ocean routines.
c
c Revision 1.13  1996/03/21  03:18:56  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.12  1995/05/04  04:10:14  ldr
c Initialize iacnt to 0 in executable code rather than in data statement
c which upsets the Fujitsu.
c
c Revision 1.11  1994/08/08  17:21:45  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.10  94/03/30  12:34:37  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c
c Revision 1.9  93/12/17  15:33:11  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.8  93/11/29  11:38:32  ldr
c Changes to V4-4-32l from HBG for coupled model
c
c Revision 1.7  93/11/03  11:44:20  ldr
c Changes to V4-4-24l from HBG for coupled model
c
c Revision 1.6  93/10/05  13:06:41  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c
c Revision 1.5  93/08/19  15:08:44  ldr
c Minor cosmetic changes.
c
c Revision 1.4  93/08/13  14:43:28  ldr
c Minor changes to get coupled model to run on SGI.
c
c Revision 1.3  93/08/10  15:27:35  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c
c Revision 1.2  93/07/12  14:14:22  ldr
c Minor changes from SPO for coupled model V4-4.
c
      subroutine ocean(lcouple)
C
C=======================================================================
C                                                                    ===
C  OCEAN IS THE PRIMARY CALLING ROUTINE.  IT PERFORMS ALL            ===
C        OPERATIONS WHICH NEED BE DONE ONLY ONCE AT THE              ===
C        BEGINNING OF EACH RUN OF THE MODEL, CALLS STEP              ===
C        ONCE PER TIMESTEP, AND ATTENDS TO OPERATIONS                ===
C        WHICH MUST BE DONE ONLY AT THE END OF EACH RUN.             ===
C                                                                    ===
C  Isopycnal routine includes code to calculate the                  ===
C  eddy ("bolus") transport velocity of GM (J.P.O. 1990)             ===
C  (uedd, vedd, wedd).                                               ===
C                                                                    ===
C        THE CALLING SEQUENCE OF THE CODE IS AS FOLLOWS:             ===
C                                                                    ===
C                             ---> CLINIC - STATE/M2003              ===
C                           /                                        ===
C                          /                                         ===
C    -->  OCEAN --> STEP -- - ---> TRACER - STATEC/M2003             ===
C                   /      \                                         ===
C                MATRIX     \                                        ===
C                             ---> RELAX                             ===
C                                                                    ===
C      Note : main.f calls                                           ===
C                  ocstart = OCEAN                                   ===
C                      occoup1 (startup if coupled)                  ===
C                  ocstep  = STEP                                    ===
C                      ocinit (coupled atmos interface)              ===
C                      ocend (housekeeping per step if coupled)      ===
C                  ocemon (end of month jobs)                        ===
C                  ocfinal (end of run job)                          ===
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
      include 'FIELDS.f'
      include 'WORKSP.f'
      include 'COEFFS.f'
      include 'ACHEXT.f'
      include 'ETRANS2.f'
      include 'TTFC.f'  !agcm
      include 'KTSIN.f'
      include 'ATM2OC.f'        !agcm
      include 'TIMEX.f'
      include 'CHANGEM.f'
      include 'CRELAX.f'
      include 'AOGCM.f'
      include 'OHIST.f'
C
C---------------------------------------------------------------------
C  DEFINE AND EQUIVALENCE LOCAL DATA; DEFINE NAMELIST INPUT
C---------------------------------------------------------------------
C
      integer kpr, i, j, k, l, m, ierr, lsegp, ibk, isp, iept,
     &        iepu, jrev, kz, ndk, ndka, ib
      real fkmp, fkmq, fkmz, fx, fxa, fxb
      DIMENSION FKMP(IMT,JMT),FKMQ(IMT,JMT),FKMZ(IMT,JMT)
      EQUIVALENCE (P,FKMP),(PB,FKMZ),(ZTD,FKMQ)
      DIMENSION KPR(IMT)
      INTEGER JTOP(JMT)
      logical lcouple
      REAL SUMDY
      REAL TINITF(KM,NT)
      NAMELIST /EDDY2/AHH1F,AHH2F,AHH3F
      NAMELIST /EDDY/ cos_amf, AMF,FKPMF,AHI1F,AHI2F,AHI3F,SLMXRF
      NAMELIST /TSTEPS/ DTTSF,DTUVF,DTSFF
      NAMELIST /PARMS/ ACORF,MXSCAN,SORF,CRITF,TRELAX
      NAMELIST /CONTRL/ NFIRST,NNERGY,NMIX,NTSI,NA, robert_time_filter,
     &                  pnu, check_conservation, m2003_eos, bulk_force,
     &                  ncar_bulk_method, mit_bulk_method
      NAMELIST /TSPROF/ TINITF
      NAMELIST /IBOX/ ISIS,IEIS,JSIS,JEIS
      NAMELIST /ACCEL/ DTXF
      NAMELIST /PLTG/ JPLOT
      namelist /etrans/ no_gm_arctic, ahe1f,ahe2f,ahe3f
      NAMELIST /COEFS/ CDRAG,ITSET
      NAMELIST /ICPLE/ IOCYR,IOCMN
      namelist /osave/ save_smfzon, save_smfmer, save_stfht,
     &  save_stfsal, save_temp, save_sal, save_rho, save_u, save_v,
     &  save_w, save_uedd, save_vedd, save_wedd, save_res,
     &  save_cdepthm, save_over

      CHARACTER*1 DOT,BLK,ABT
      data slopemx/2*500.,430.,350.,290.,220.,150.,100.,70.,12*50./
      DATA DOT/'.'/,BLK/' '/
      DIMENSION ABT(IMT)
      integer ii, jj, kk
      character title*60

C
C---------------------------------------------------------------------
C  BEGIN EXECUTABLE CODE
C---------------------------------------------------------------------
C
      entry ocstart(lcouple)
C
C=======================================================================
C  BEGIN INTRODUCTORY SECTION WHICH IS NEEDED FOR EACH RUN  ============
C       (INCLUDING RESTARTS)                                ============
C=======================================================================
C

c...  Print the parameter values specified at compile time
      write (*, *)
      write (*, *) "The ocean model was compiled with the following ",
     &             "parameter values:"
      write (*, *)
      write (*, *) "IMT   = ", imt
      write (*, *) "JMT   = ", jmt
      write (*, *) "KM    = ", km
      write (*, *) "NT    = ", nt
      write (*, *) "LSEG  = ", lseg
      write (*, *) "NISLE = ", nisle
      write (*, *) "JFRST = ", jfrst
      write (*, *) "JFT0  = ", jft0
      write (*, *) "JFT1  = ", jft1
      write (*, *) "JFT2  = ", jft2
      write (*, *) "JFU0  = ", jfu0
      write (*, *) "JFU1  = ", jfu1
      write (*, *) "JFU2  = ", jfu2
      write (*, *)

c...  Initialise the ocean model history arrays
      do j = 2, jmtm1
        do i = 2, imtm1
          ohist_smfzon(i, j) = 0.0
          ohist_smfmer(i, j) = 0.0
          ohist_stfht(i, j) = 0.0
          ohist_stfsal(i, j) = 0.0
          ohist_res(i, j) = 0.0
        end do
      end do

      do k = 1, km
        do j = 2, jmtm1
           do i = 2, imtm1
            ohist_temp(i, j, k) = 0.0
            ohist_sal(i, j, k) = 0.0
            ohist_rho(i, j, k) = 0.0
            ohist_u(i, j, k) = 0.0
            ohist_v(i, j, k) = 0.0
            ohist_w(i, j, k) = 0.0
            ohist_uedd(i, j, k) = 0.0
            ohist_vedd(i, j, k) = 0.0
            ohist_wedd(i, j, k) = 0.0
          end do
        end do
      end do

c...  Initialise the arrays which contain the relaxation surface fluxes within
c...  the coupled model
      do j = 1, jmt
        do i = 1, imt
          crelax_qnf(i, j) = 0.0
          crelax_salf(i, j) = 0.0
          crelax_stfht(i, j) = 0.0
          crelax_stfsal(i, j) = 0.0
        end do
      end do

C
C---------------------------------------------------------------------
C  SET THE TYPE OF MIXING TIMESTEP.  IF EB=.TRUE., AN EULER
C  BACKWARD STEP IS DONE; IF EB=.FALSE., A FORWARD STEP IS DONE.
C---------------------------------------------------------------------
C
      EB=.TRUE.
C
C---------------------------------------------------------------------
C  INITIALIZE VARIOUS QUANTITIES
C---------------------------------------------------------------------
C

c...  Specify default values for the NAMELIST parameters
c...
c...  &CONTRL
# 447

      nfirst = 0
      nnergy = 87600
      nmix = 19
      ntsi = 720
      na = 1
      robert_time_filter = .true.
      pnu = 0.005
      check_conservation = .false.
      m2003_eos = .true.
      bulk_force = .false.
      ncar_bulk_method = .true.
      mit_bulk_method = .false.


c...  &EDDY
# 471

      cos_amf = .true.
      amf = 3.2e9
      fkpmf = 20.0
      ahi1f = 0.6e7
      ahi2f = 0.6e7
      ahi3f = 5.0e4
      slmxrf = 100.0


c...  &EDDY2
# 486

      ahh1f = 0.6e7
      ahh2f = 0.6e7
      ahh3f = 5.0e4


c...  &ETRANS
# 498

      no_gm_arctic = .true.
      ahe1f = 0.6e7
      ahe2f = 0.6e7
      ahe3f = 5.0e4


c...  &TSTEPS
      dttsf = 3600.0
      dtuvf = 3600.0
      dtsff = 3600.0

c...  &PARMS
# 517

      acorf = 0.5
      mxscan = 1000
      sorf = 1.4
      critf = 1.0e8
      trelax = 20.0


c...  &ICPLE
      iocyr = 50
      iocmn = 12

c...  &PLTG
# 532

      jplot = 25


c...  &COEFS
      cdrag = 2.6e-3
      itset = 1

c...  &ACCEL
      do i = 1, km
        dtxf(i) = 1.0
      end do

c...  &OSAVE
      save_smfzon  = .true.
      save_smfmer  = .true.
      save_stfht   = .true.
      save_stfsal  = .true.
      save_temp    = .true.
      save_sal     = .true.
      save_rho     = .false.
      save_u       = .true.
      save_v       = .true.
      save_w       = .true.
      save_uedd    = .true.
      save_vedd    = .true.
      save_wedd    = .true.
      save_res     = .true.
      save_cdepthm = .true.
      save_over    = .true.

c When using the new ocean model grid, ensure that the ocean and atmosphere
c models use the same values for pi and the radius of the Earth. For the
c purposes of backwards compatibility, use the original values when using the
c old ocean model grid.
c
c In due course, the model should be modified to ensure consistent use of
c physical parameters throughout.
c
c SJP 2007/06/01

# 578

      pi = 3.14159265
      omega = pi/43082.0
      radius = 6.37122e8
      radian = 180.0/pi


      GRAV=980.6
C
C---------------------------------------------------------------------
C  SET THE LATITUDE (IN DEGREES) OF THE SOUTHERN WALL
C---------------------------------------------------------------------
C

# 594

      swldeg = -89.1331111924962


C
C---------------------------------------------------------------------
C  SET Y DIMENSION OF BOXES IN DEGREES AND CONVERT TO CENTIMETERS
C---------------------------------------------------------------------
C

# 661

        dyt(  2) =  1.5717758068763
        dyt(  3) =  1.5717758068763
        dyt(  4) =  1.5873125395736
        dyt(  5) =  1.5873125395736
        dyt(  6) =  1.5904318409297
        dyt(  7) =  1.5904318409297
        dyt(  8) =  1.5915276112852
        dyt(  9) =  1.5915276112852
        dyt( 10) =  1.5920305756336
        dyt( 11) =  1.5920305756336
        dyt( 12) =  1.5923010877064
        dyt( 13) =  1.5923010877064
        dyt( 14) =  1.5924626726494
        dyt( 15) =  1.5924626726494
        dyt( 16) =  1.5925666534402
        dyt( 17) =  1.5925666534402
        dyt( 18) =  1.5926373800807
        dyt( 19) =  1.5926373800807
        dyt( 20) =  1.5926875876849
        dyt( 21) =  1.5926875876849
        dyt( 22) =  1.5927244533195
        dyt( 23) =  1.5927244533195
        dyt( 24) =  1.5927522680130
        dyt( 25) =  1.5927522680130
        dyt( 26) =  1.5927737232154
        dyt( 27) =  1.5927737232154
        dyt( 28) =  1.5927905752292
        dyt( 29) =  1.5927905752292
        dyt( 30) =  1.5928040087625
        dyt( 31) =  1.5928040087625
        dyt( 32) =  1.5928148456915
        dyt( 33) =  1.5928148456915
        dyt( 34) =  1.5928236699541
        dyt( 35) =  1.5928236699541
        dyt( 36) =  1.5928309049419
        dyt( 37) =  1.5928309049419
        dyt( 38) =  1.5928368629406
        dyt( 39) =  1.5928368629406
        dyt( 40) =  1.5928417775583
        dyt( 41) =  1.5928417775583
        dyt( 42) =  1.5928458254894
        dyt( 43) =  1.5928458254894
        dyt( 44) =  1.5928491414079
        dyt( 45) =  1.5928491414079
        dyt( 46) =  1.5928518283288
        dyt( 47) =  1.5928518283288
        dyt( 48) =  1.5928539649014
        dyt( 49) =  1.5928539649014
        dyt( 50) =  1.5928556105851
        dyt( 51) =  1.5928556105851
        dyt( 52) =  1.5928568093168
        dyt( 53) =  1.5928568093168
        dyt( 54) =  1.5928575920821
        dyt( 55) =  1.5928575920821
        dyt( 56) =  1.5928579786506
        dyt( 57) =  1.5928579786506
        dyt( 58) =  1.5928579786506
        dyt( 59) =  1.5928579786506
        dyt( 60) =  1.5928575920821
        dyt( 61) =  1.5928575920821
        dyt( 62) =  1.5928568093168
        dyt( 63) =  1.5928568093168
        dyt( 64) =  1.5928556105851
        dyt( 65) =  1.5928556105851
        dyt( 66) =  1.5928539649014
        dyt( 67) =  1.5928539649014
        dyt( 68) =  1.5928518283288
        dyt( 69) =  1.5928518283288
        dyt( 70) =  1.5928491414079
        dyt( 71) =  1.5928491414079
        dyt( 72) =  1.5928458254894
        dyt( 73) =  1.5928458254894
        dyt( 74) =  1.5928417775583
        dyt( 75) =  1.5928417775583
        dyt( 76) =  1.5928368629406
        dyt( 77) =  1.5928368629406
        dyt( 78) =  1.5928309049419
        dyt( 79) =  1.5928309049419
        dyt( 80) =  1.5928236699541
        dyt( 81) =  1.5928236699541
        dyt( 82) =  1.5928148456915
        dyt( 83) =  1.5928148456915
        dyt( 84) =  1.5928040087625
        dyt( 85) =  1.5928040087625
        dyt( 86) =  1.5927905752292
        dyt( 87) =  1.5927905752292
        dyt( 88) =  1.5927737232154
        dyt( 89) =  1.5927737232154
        dyt( 90) =  1.5927522680130
        dyt( 91) =  1.5927522680130
        dyt( 92) =  1.5927244533195
        dyt( 93) =  1.5927244533195
        dyt( 94) =  1.5926875876849
        dyt( 95) =  1.5926875876849
        dyt( 96) =  1.5926373800807
        dyt( 97) =  1.5926373800807
        dyt( 98) =  1.5925666534402
        dyt( 99) =  1.5925666534402
        dyt(100) =  1.5924626726494
        dyt(101) =  1.5924626726494
        dyt(102) =  1.5923010877064
        dyt(103) =  1.5923010877064
        dyt(104) =  1.5920305756336
        dyt(105) =  1.5920305756336
        dyt(106) =  1.5915276112852
        dyt(107) =  1.5915276112852
        dyt(108) =  1.5904318409297
        dyt(109) =  1.5904318409297
        dyt(110) =  1.5873125395736
        dyt(111) =  1.5873125395736
        dyt(112) =  1.5717758068763
        dyt(113) =  1.5717758068763


      DO 52 J=2,JMTM1
        DYT(J)=DYT(J)*RADIUS/RADIAN
 52   CONTINUE
      DYT(1)=DYT(2)
      DYT(JMT)=DYT(JMTM1)
C
C---------------------------------------------------------------------
C  SET X DIMENSION OF BOXES IN DEGREES AND CONVERT TO CENTIMETERS
C---------------------------------------------------------------------
C
      DO 57 I=2,IMTM1
        DXT(I)=360./FLOAT(IMT-2)
        DXT(I)=DXT(I)*RADIUS/RADIAN
 57   CONTINUE
      DXT(1)=DXT(2)
      DXT(IMT)=DXT(IMTM1)
C
C  SET CYCLIC CONDITIONS ON DXT
C
      DXT(1)=DXT(IMTM1)
      DXT(IMT)=DXT(2)
C
C---------------------------------------------------------------------
C  SET Z DIMENSION OF BOXES (IN CENTIMETERS)
C---------------------------------------------------------------------
C
        DZ(1) =  25.E2       
        DZ(2) =  25.E2      
        DZ(3) =  30.E2     
        DZ(4) =  37.E2    
        DZ(5) =  43.E2   
        DZ(6) =  50.E2  
        DZ(7) =  60.E2 
        DZ(8) =  80.E2
        DZ(9) = 120.E2      
        DZ(10) = 150.E2    
        DZ(11) = 180.E2   
        DZ(12) = 210.E2  
        DZ(13) = 240.E2 
        DZ(14) = 290.E2
        DZ(15) = 360.E2     
        DZ(16) = 450.E2    
        DZ(17) = 450.E2   
        DZ(18) = 450.E2  
        DZ(19) = 450.E2 
        DZ(20) = 450.E2
        DZ(21) = 450.E2

        dzlev1=DZ(1)*0.01 ! Depth of layer 1 in m (for ice model)
C
C---------------------------------------------------------------------
C  READ IN RUN PARAMETERS
C---------------------------------------------------------------------
C
      READ (5,CONTRL)
      WRITE(6,CONTRL)
      READ (5,EDDY)
      WRITE(6,EDDY)
      READ(5,EDDY2)
      WRITE(6,EDDY2)
      read(5,etrans)
      write(6,etrans)
      READ (5,TSTEPS)
      WRITE(6,TSTEPS)
      READ (5,PARMS)
      WRITE(6,PARMS)
      READ(5,ICPLE)
      WRITE(6,ICPLE)
      READ(5,PLTG)
      WRITE(6,PLTG)
      READ(5,COEFS)
      WRITE(6,COEFS)
C  ************** SET ACCEL FACTORS FOR 21-LEVEL MODEL  ****************
C  ***** N.B. MUST SET ALL FACTORS TO 1.0 TO GET TRUE SEASONAL BEHAVIOUR
      READ(5,ACCEL)
      PRINT*,'SPEEDUP FACTORS:',(DTXF(K),K=1,KM)
      read (5, osave)
      write (6, osave)

c...  Check for internal consistency between run parameters
      if (save_rho .and. .not. m2003_eos) then
        write (*, *)
        write (*, *) "***  WARNING: Inconsistency in run parameters"
        write (*, *) "***"
        write (*, *) "***  M2003_EOS = F"
        write (*, *) "***  SAVE_RHO = T"
        write (*, *) "***"
        write (*, *) "***  Density can only be saved if M2003_EOS = T"
        write (*, *) "***"
        write (*, *) "***  Setting SAVE_RHO = F for this simulation ..."
        write (*, *)
        save_rho = .false.
      end if

      knitd = int(24.0 * 3600.0 / DTTSF + 0.1)

      pnu2m = 1.0 - 2.0 * pnu
C
C---------------------------------------------------------------------
C  COMPUTE AUXILIARY ARRAYS BASED UPON THE SPACING SPECIFIED ABOVE
C---------------------------------------------------------------------
C
      DO 100 K=1,KM
        C2DZ(K)=2.0*DZ(K)
        DZ2R(K)=1.0/C2DZ(K)
 100  CONTINUE
      DZZ(1)=0.5*DZ(1)
      ZDZ(1)=DZ(1)
      DO 110 K=2,KM
        DZZ(K)=0.5*(DZ(K-1)+DZ(K))
        ZDZ(K)=ZDZ(K-1)+DZ(K)
 110  CONTINUE
      DZZ(KM+1)=0.5*DZ(KM)
      DZZ2R(KMP1)=0.5/DZZ(KMP1)
      ZDZZ(1)=DZZ(1)
      DO 120 K=1,KM
        DZZ2R(K)=.5/DZZ(K)
        ZDZZ(K+1)=ZDZZ(K)+DZZ(K+1)
        EEM(K)=FKPMF/(DZ(K)*DZZ(K))
        FFM(K)=FKPMF/(DZ(K)*DZZ(K+1))

c...  ZDB contains the depth of each model level in metres, and is used as an
c...  approximation for the pressure in db when calculating the density.

        zdb(k) = 0.01 * zdzz(k)

C  **************   SET SCALING FOR HIGH LATITUDE AHI ******************
        KAR(K)=K
        AHI(K)=AHI2F+(AHI1F-AHI2F)*EXP(-ZDZZ(K)/AHI3F)
        AHH(K)=AHH2F+(AHH1F-AHH2F)*EXP(-ZDZZ(K)/AHH3F)

c                set vertical profile of ahe
c         Phase in ahe from top in such a way as to give roughly uniform
c          zonally averaged velocity profile in upper 270 m, in Southern Ocean
c        ahe(2), ahe(3) chosen experimentally to give roughly uniform profiles
c       in levels 1--3 despite strong differences in vertical diffusivity
c       between these levels. In constant v. diff. zone (levels 4--7),
c       differences between the ahe are proportional to the grid box
c       thicknesses.

       if(k.gt.1) ahe(k)=ahe2f+(ahe1f-ahe2f)*exp(-zdz(k-1)/ahe3f)

       if(k.eq.7)then
         ahe(1) = 0.00
         ahe(2) = 0.07e7
         ahe(3) = 0.18e7
         ahe(4) = 0.29e7
         ahe(5) = 0.42e7
         if(ahe(6).gt.0.58e7)ahe(6) = 0.58e7
         if(ahe(7).gt.0.77e7)ahe(7) = 0.77e7
       endif

c               set change from eddy trans. to horiz. diff at high N. lats
c
c Note that this feature is currently disabled for both the old and new ocean
c model grids, with Gent-McWilliams eddy diffusion extending all the way to the
c North Pole in both cases.
c
c This feature was originally implemented because of numerical problems arising
c from the use of eddy diffusion at high northern latitudes (Tony Hirst, pers.
c comm.). However, it no longer appears to be necessary.
c
c SJP 2007/06/01

      do j = 1,jmt
        ahifac(j, k) = 1.0
        ahhfac(j, k) = 0.0
        ahefac(j, k) = 1.0
        if (no_gm_arctic) then
          if (j .eq. jft2-1) then
            ahhfac(j, k) = 0.5
            ahefac(j, k) = 0.5
          end if
          if (j .ge. jft2) then
            ahhfac(j, k) = 1.0
            ahefac(j, k) = 0.0
          end if
        end if
      end do
            
c...  The following section of code merely has the effect of ensuring a maximum
c...  value of 1.0e7 cm^2/s for the isopycnal tracer diffusivity at the North
c...  Pole, with a smooth transition towards that value over the Arctic Ocean.
c...
c...  Disable this feature when using the new ocean model grid, as a value
c...  higher than 1.0e7 cm^2/s should never be specified for the isopycnal
c...  tracer diffusivity.
c...
c...  SJP 2007/07/04

# 975


 120  CONTINUE

      write (*, *)
      write (*, '(a)') " Ocean isopycnal diffusivity:"
      write (*, *)
      write (*, '(a)') "   k       Depth      Diffusivity"
      write (*, '(a)') "            (m)         (cm^2/s)"
      write (*, *)
      do k = 1, km
        write (*, '(i4, f12.2, 1pe17.6)') k, 0.01*zdzz(k), ahi(k)
      end do
      write (*, *)

      write (*, *)
      write (*, '(a)') " Ocean eddy diffusivity:"
      write (*, *)
      write (*, '(a)') "   k       Depth      Diffusivity"
      write (*, '(a)') "            (m)         (cm^2/s)"
      write (*, *)
      do k = 1, km
        write (*, '(i4, f12.2, 1pe17.6)') k, 0.01*zdzz(k), ahe(k)
      end do
      write (*, *)

      if (no_gm_arctic) then
        write (*, *)
        write (*, '(a)') " Ocean horizontal diffusivity:"
        write (*, *)
        write (*, '(a)') "   k       Depth      Diffusivity"
        write (*, '(a)') "            (m)         (cm^2/s)"
        write (*, *)
        do k = 1, km
          write (*, '(i4, f12.2, 1pe17.6)') k, 0.01*zdzz(k), ahh(k)
        end do
        write (*, *)
      end if

      PHI(1)=SWLDEG/RADIAN
      PHIT(1)=PHI(1)-.5*DYT(1)/RADIUS
      SUMDY=PHI(1)
      DYU(JMT)=DYT(JMT)
      DO 130 J=1,JMT
        IF(J.NE.JMT) DYU(J)=.5*(DYT(J)+DYT(J+1))
        DYTR(J)=1./DYT(J)
        DYT2R(J)=.5/DYT(J)
        DYT4R(J)=.25/DYT(J)
        DYUR(J)=1./DYU(J)
        DYU2R(J)=.5/DYU(J)
        DYU4R(J)=.25/DYU(J)
        IF(J.NE.JMT) SUMDY=SUMDY+DYT(J+1)/RADIUS
        IF(J.NE.JMT) PHI(J+1)=SUMDY
        IF(J.NE.1) PHIT(J)=.5*(PHI(J-1)+PHI(J))
        CST(J)=COS(PHIT(J))
        CS (J)=COS(PHI (J))
        SINE(J)=SIN(PHI(J))
        CSTR(J)=1.0/CST(J)
        CSR(J)=1.0/CS(J)
        TNG(J)=SINE(J)/CS(J)
 130  CONTINUE
      DXU(IMT)=DXT(IMT)
      DXU(IMT)=.5*(DXT(2)+DXT(3))
      DO 140  I=1,IMT
        IF(I.NE.IMT) DXU(I)=.5*(DXT(I)+DXT(I+1))
        DXTR(I)=1./DXT(I)
        DXT2R(I)=.5/DXT(I)
        DXT4R(I)=.25/DXT(I)
        DXUR(I)=1./DXU(I)
        DXU2R(I)=.5/DXU(I)
        DXU4R(I)=.25/DXU(I)
  140 CONTINUE
C
C---------------------------------------------------------------------
C  COMPUTE SIN AND COS VALUES FOR VECTOR CORRECTION BEFORE FILTER
C---------------------------------------------------------------------
C
      FX=1.0E-10
      FXA=DXT(1)/RADIUS
      DO 670 I=2,IMUM1
        FXB=FXA*FLOAT(I-2)
        SPSIN(I)=SIN(FXB)
        SPCOS(I)=COS(FXB)
        IF(ABS(SPSIN(I)).LT.FX)SPSIN(I)=0.0
        IF(ABS(SPCOS(I)).LT.FX)SPCOS(I)=0.0
  670 CONTINUE
      SPSIN(1)=0.0
      SPCOS(1)=0.0
      SPSIN(IMU)=0.0
      SPCOS(IMU)=0.0

      DO 674 J = 1,JMT
        AMFAC(J)   = 1.
  674 CONTINUE

      if (cos_amf) then

        do j = 2, jmtm1
          amfac(j) = cs(j)
          if (j .ge. jfu2) amfac(j) = amfac(jfu0)
        end do

      else

c...  For backwards compatibility, retain this section of code for the old,
c...  low-resolution version of the model. SJP 2007/02/27

# 1095


      end if

      write (*, *)
      write (*, '(a)') " Ocean horizontal viscosity:"
      write (*, *)
      write (*, '(a)') "    j       Latitude        Viscosity"
      write (*, '(a)') "              (deg)          (cm^2/s)"
      write (*, *)
      do j = 2, jmtm1
        write (*, '(i5, f15.6, 1pe17.6)') j, radian*phi(j),
     &                                    amfac(j)*amf
      end do
      write (*, *)

C
C---------------------------------------------------------------------
C  PRINT GRID GEOMETRY ARRAYS
C---------------------------------------------------------------------
C
c     PRINT 9701
c9701 FORMAT(50H0 GRID BOX THICKNESS  'DZ'                        )
c     PRINT 970, DZ
c     PRINT 9702
c9702 FORMAT(50H0 GRID POINT SEPARATION  'DZZ'                    )
c     PRINT 970, DZZ
c     PRINT 9703
c9703 FORMAT(50H0 DEPTH OF BOX BOTTOM  'ZDZ'                      )
c     PRINT 970, ZDZ
c     PRINT 9704
c9704 FORMAT(50H0 DEPTH OF GRID POINT  'ZDZZ'                     )
c     PRINT 970, ZDZZ
c     PRINT 9705
c9705 FORMAT(50H0 LATITUDE OF T,S POINTS (RADIANS)  'PHIT'        )
c     PRINT 970, PHIT
c     PRINT 9706
c9706 FORMAT(50H0 LATITUDE OF U,V POINTS (RADIANS)  'PHI'         )
c     PRINT 970, PHI
c     PRINT 9707
c9707 FORMAT(50H0 COSINE OF T,S LATITUDE  'CST'                   )
c     PRINT 970, CST
c     PRINT 9708
c9708 FORMAT(50H0 COSINE OF U,V LATITUDE  'CS'                    )
c     PRINT 970, CS
c     PRINT 9709
c9709 FORMAT(50H0 SINE OF U,V LATITUDE  'SINE'                    )
c     PRINT 970, SINE
c 970 FORMAT(1X,6E13.5)
c
C---------------------------------------------------------------------
C  OPEN THE DISK DATASETS
C---------------------------------------------------------------------
C

c...  Read the data from the ocean model restart file, unless the model is
c...  being initialised from scratch
      if (nfirst .eq. 0) then
        write (*, *)
        write (*, *) "Reading data from ocean model restart file..."
        write (*, *)
        call orest_read
      end if

c...  Within the coupled model, set the volume of the ocean for the purposes of
c...  freshwater conservation
      if (lcouple) then
        ovolume = 1.0e-6 * volume
        write (*, *)
        write (*, *) "Setting OVOLUME = ", ovolume, " m^3 ..."
        write (*, *)
      end if

c...  Derive a logical mask which indicates whether or not each latitude row
c...  on the tracer grid contains any ocean gridpoints. This mask is used to
c...  avoid performing unnecessary density calculations.
      do k = 1, km
        do j = 1, jmt
          omask(j, k) = .false.
          do i = 1, imt
            if (kmt(i, j) .ge. k) omask(j, k) = .true.
          end do
        end do
      end do

      if (save_over) then

c.....  Read the ocean basin masks. These are used when calculating the
c.....  meridional overturning streamfunctions.
        write (*, *)
        write (*, *) "Reading ocean basin masks :"
        write (*, *)
        write (*, *) "File     =  bsnmask.nc"
        write (*, *) "Variable =  mask"

        call read_int_2d("bsnmask.nc", "mask", imt-2, jmt-2, basin,
     &                   title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

c.....  Derive logical masks which indicate whether or not each latitude row in
c.....  each ocean basin contains any large-scale velocity or any eddy-induced
c.....  velocity gridpoints. These masks are used when calculating the
c.....  meridional overturning streamfunctions. Note that there are NBASIN-1
c.....  basins; basin number NBASIN represents the world ocean.
        do ib = 1, nbasin
          do k = 1, km
            do j = 2, jmtm1
              vmask(j, k, ib) = .false.
              veddmask(j, k, ib) = .false.
            end do
          end do
        end do
        do k = 1, km
          do j = 2, jmtm1
            do i = 2, imtm1
              if (kmu(i, j) .ge. k) then
                vmask(j, k, nbasin) = .true.
                do ib = 1, nbasin-1
                  if (basin(i, j) .eq. ib) vmask(j, k, ib) = .true.
                end do
              end if
              if ((kmt(i, j) .ge. k) .and. (kmt(i, j+1) .ge. k)) then
                veddmask(j, k, nbasin) = .true.
                do ib = 1, nbasin-1
                  if (basin(i, j) .eq. ib) veddmask(j, k, ib) =  .true.
                end do
              end if
            end do
          end do
        end do

      end if

c...  Initialise the array CDEPTHM, which contains the convective layer depth
      do j = 2, jmtm1
        do i = 2, imtm1
          ohist_cdepthm(i, j) = zdz(1)
        end do
      end do

      RHOCW = 4.1E6
C=======================================================================
C  END INTRODUCTORY SECTION  ===========================================
C=======================================================================
C
      print*, "INTRODUCTORY SECTION COMPLETE"
      print*, " "
C=======================================================================
C  BEGIN SECTION WHICH IS EXECUTED ONLY WHEN STARTING  =================
C       A RUN FROM SCRATCH                             =================
C=======================================================================
C
C   This section must be activated (i.e. NFIRST==1) if the ocean grid is
C   changed.

      IF (NFIRST.EQ.0) GO TO 160
C
C---------------------------------------------------------------------
C  SET MAXIMUM LEVEL INDICATORS FOR TOPOGRAPHY
C---------------------------------------------------------------------
C
      DO 690 J=1,JMT
      DO 690 I=1,IMT
        FKMP(I,J)=0
        FKMQ(I,J)=0
        FKMZ(I,J)=0
 690  CONTINUE
C
C  1ST, SET NUMBER OF VERTICAL LEVELS FOR T POINTS
C
C *********************    READ IN TOPOGRAPHY    ********************
        OPEN(25,FILE='sttop.bot_ind',STATUS='old',iostat=ierr)
        call errcheck(ierr,'sttop.bot_ind    ','ocean     ')
        DO 700 I = 2,IMTM1
          READ(25,705)(JTOP(J),J=1,JMTM2)
  705     FORMAT(56I2)
          DO 700 J=2,JMTM1
            FKMP(I,J) = JTOP(J-1)
C           IF(FKMP(I,J).GE.1)FKMP(I,J) = 21
  700   CONTINUE
        CLOSE(25)
C
C  SET CYCLIC BOUNDARY CONDITIONS
C
      DO 728 J=1,JMT
        FKMP(  1,J)=FKMP(IMTM1,J)
        FKMP(IMT,J)=FKMP(    2,J)
 728  CONTINUE
C
C  2ND, COMPUTE NUMBER OF VERTICAL LEVELS AT EACH U,V POINT
C
      DO 730 J=1,JMTM1
      DO 730 I=1,IMTM1
       FKMQ(I,J)=MIN(FKMP(I,J),FKMP(I+1,J),FKMP(I,J+1),FKMP(I+1,J+1))
 730  CONTINUE
C
C  SET CYCLIC CONDITIONS
C
      DO 732 J=1,JMT
        FKMQ(IMT,J)=FKMQ(2,J)
 732  CONTINUE
C
C  3RD, COMPUTE AN ARRAY TO INDICATE "INTERIOR" GRID BOXES
C
      ! this variable (FKMZ) must be altered to accomodate different
      ! ocean grids, otherwise later calculations will be told to occur
      ! in the wrong places (i.e. over land)
      DO 740 J=2,JMTM1
      DO 740 I=2,IMU
        FKMZ(I,J)=MIN(FKMQ(I-1,J-1),FKMQ(I,J-1),FKMQ(I-1,J),FKMQ(I,J))
 740  CONTINUE
C
C  SET CYCLIC CONDITIONS
C
      DO 742 J=1,JMT
        FKMZ(1,J)=FKMZ(IMTM1,J)
 742  CONTINUE
C
C---------------------------------------------------------------------
C  COMPUTE START & END INDICES FOR STREAM FUNCTION CALCULATIONS
C---------------------------------------------------------------------
C

      do j = 1, lseg
        do i = 1, jmt
          isz(i, j) = 0
          iez(i, j) = 0
        end do
      end do
      do i = 1, nisle
        isis(i) = 0
        ieis(i) = 0
        jsis(i) = 0
        jeis(i) = 0
      end do
      do k = 1, km
        do j = 1, lsegf
          do i = 1, njtbft
            istf(i, j, k) = 0
            ietf(i, j, k) = 0
          end do
        end do
      end do
      do k = 1, km
        do j = 1, lsegf
          do i = 1, njtbfu
            isuf(i, j, k) = 0
            ieuf(i, j, k) = 0
          end do
        end do
      end do
      do j = 1, lsegf
        do i = 1, njtbfu
          iszf(i, j) = 0
          iezf(i, j) = 0
        end do
      end do
      LSEGP=LSEG+1
      DO 780 J=3,JSCAN
        L=1
      DO 780 I=2,IMUM1
        IF(FKMZ(I-1,J).EQ.0. .AND. FKMZ(I,J).NE.0.) ISZ(J,L)=I
        IF(I.EQ.2 .AND. FKMZ(IMT,J).NE.0.) ISZ(J,L)=2
        IF(FKMZ(I,J).NE.0. .AND. FKMZ(I+1,J).EQ.0.) IEZ(J,L)=I
        IF(I.EQ.IMUM1 .AND. FKMZ(I+1,J).NE.0.) IEZ(J,L)=I
        IF(FKMZ(I,J).NE.0. .AND. FKMZ(I+1,J).EQ.0.) L=L+1
        IF(L.GT.LSEGP) STOP 780
 780  CONTINUE
C
C---------------------------------------------------------------------
C  FIND AND PRINT START & END INDICES FOR FILTERING
C---------------------------------------------------------------------
C
      PRINT 833
      IF (LSEGF.GT.11) PRINT 834
      PRINT 835
      CALL FINDEX(FKMP,NJTBFT,KM,JFT1,JFT2,IMT,ISTF,IETF)
      PRINT 836
      CALL FINDEX(FKMQ,NJTBFU,KM,JFU1,JFU2,IMU,ISUF,IEUF)
      PRINT 837
      CALL FINDEX(FKMZ,NJTBFU, 1,JFU1,JFU2,IMT,ISZF,IEZF)
  833 FORMAT (1H1,'START AND END INDICES FOR FOURIER FILTERING:'/)
  834 FORMAT (1X,'ONLY 11 SETS OF INDICES FIT ACCROSS THE PAGE.',
     &       '  OTHERS WILL NOT BE PRINTED.'/)
  835 FORMAT (///1X,'FILTERING INDICES FOR T,S:')
  836 FORMAT (///1X,'FILTERING INDICES FOR U,V:')
  837 FORMAT (///1X,'FILTERING INDICES FOR STREAM FUNCTION:')
C
C---------------------------------------------------------------------
C  COMPUTE FIELD OF RECIPROCAL DEPTH
C---------------------------------------------------------------------
C
      DO 790 J=1,JMT
      DO 790 I=1,IMT
        HR(I,J)=0.0
        IF(FKMQ(I,J).NE.0) HR(I,J)=1./ZDZ(INT(FKMQ(I,J)))
 790  CONTINUE
C
C---------------------------------------------------------------------
C  COMPUTE THE SURFACE AREA AND VOLUME OF THE OCEAN
C---------------------------------------------------------------------
C
      AREA=0.0
      VOLUME=0.0
      DO 800 J=2,JMTM1
      DO 800 I=2,IMTM1
        IF(FKMP(I,J).GT.0) THEN
          AREA=AREA+CST(J)*DXT(I)*DYT(J)
          VOLUME=VOLUME+CST(J)*DXT(I)*DYT(J)*ZDZ(INT(FKMP(I,J)))
          print*, "Ocean area =", AREA
          print*, "Ocean volume =", VOLUME
        ENDIF
 800  CONTINUE
C
C---------------------------------------------------------------------
C  PRINT TOPOGRAPHY MAP
C  (..NOTE.. THE NUMBER OF LEVELS ARE PRINTED IN HEX;
C            A DOT SUPERIMPOSED ===> ADD AN ADDITIONAL 16)
C---------------------------------------------------------------------
C
      PRINT 950
 950  FORMAT(50H1 NUMBER OF LEVELS AT T,S POINTS AND U,V POINTS   )
      DO 810 IBK=1,IMT,65
        PRINT 960
 960    FORMAT(/)
        ISP=IBK
        IEPT=IBK+64
        IEPU=IBK+64
        IF(IEPT.GT.IMT) IEPT=IMT
        IF(IEPU.GT.IMU) IEPU=IMU
      DO 810 JREV=1,JMT
        J=JMT-JREV+1
        IF(J.NE.JMT) THEN
          DO 968 I=1,IMT
            KPR(I)=NINT(FKMQ(I,J))
 968      CONTINUE
          PRINT 972, (KPR(I),I=ISP,IEPU)
 972      FORMAT(2X,65(1X,Z1))
          DO 969 I=1,IMT
            ABT(I)=BLK
 969      CONTINUE
          DO 952 I=ISP,IEPU
            IF(KPR(I).GT.15)ABT(I)=DOT
 952      CONTINUE
          PRINT 971,(ABT(I),I=ISP,IEPU)
 971      FORMAT(2H+ ,65(1X,A1))
        ENDIF
        DO 953 I=1,IMT
          KPR(I)=NINT(FKMP(I,J))
 953    CONTINUE
        PRINT 982, (KPR(I),I=ISP,IEPT)
 982    FORMAT(1X,65(1X,Z1))
        DO 979 I=1,IMT
          ABT(I)=BLK
 979    CONTINUE
        DO 954 I=ISP,IEPT
          IF(KPR(I).GT.15)ABT(I)=DOT
 954    CONTINUE
        PRINT 981,(ABT(I),I=ISP,IEPT)
 981    FORMAT(1H+,65(1X,A1))
 810  CONTINUE
C
C---------------------------------------------------------------------
C  PRINT AREA AND VOLUME OF THE OCEAN, AS WELL AS START & END
C  INDICES FOR THE STREAM FUNCTION CALCULATION
C---------------------------------------------------------------------
C
      PRINT 940, AREA,VOLUME
 940  FORMAT(//,15H SURFACE AREA =,1PE13.6,5X,9H VOLUME =,1PE13.6)
      PRINT 9502
 9502 FORMAT(43H1 START AND END INDICES FOR STREAM FUNCTION)
      DO 830 JREV=1,JMT
        J=JMT-JREV+1
        PRINT 930,J,(ISZ(J,L),IEZ(J,L),L=1,LSEG)
 930    FORMAT(' J=',I3,5X,5(2I5,10X))
 830  CONTINUE
C
C---------------------------------------------------------------------
C  READ IN INITIAL TRACER VALUES
C  AND ISLAND BOX CORNER POINT INDICES
C---------------------------------------------------------------------
C
      READ (5,TSPROF)
      WRITE(6,TSPROF)
      DO 832 M=1,NT
      DO 832 K=1,KM
        TINIT(K,M)=TINITF(K,M)
 832  CONTINUE
      READ (5,IBOX)
      WRITE(6,IBOX)
C
C---------------------------------------------------------------------
C  INITIALIZE SLAB DATA ON DISK
C---------------------------------------------------------------------
C
      DO 880 J=1,JMT
        DO 840 I=1,IMT
C
C  SET WIND STRESS TO SPECIFIED DISTRIBUTION
C
        WSX(I,J) = 0.0
        WSY(I,J) = 0.0
C
C  SET MAXIMUM LEVEL INDICATORS TO VALUES COMPUTED ABOVE
C
          kmt(i, j) = int(fkmp(i, j) + 0.1)
          kmu(i, j) = int(fkmq(i, j) + 0.1)
 840    CONTINUE
        DO 842 K=1,KM
        DO 842 I=1,IMT
C
C  SET INTERNAL MODE VELOCITIES TO ZERO
C
          UB(I,K)=0.0
          U (I,K)=0.0
          VB(I,K)=0.0
          V (I,K)=0.0
        DO 842 M=1,NT
          TB(I,K,M)=0.0
          T (I,K,M)=0.0
 842    CONTINUE
C
C  SET TRACERS OVER OCEAN POINTS TO SPECIFIED VALUES
C
        DO 870 I=1,IMT
          KZ=NINT(FKMP(I,J))
          IF(KZ.NE.0) THEN
            DO 860 K=1,KZ
            DO 860 M=1,NT
              TB(I,K,M)=TINIT(K,M)
              IF(M.EQ.2) TB(I,K,M)=TINIT(K,M)-0.035
              T(I,K,M)=TB(I,K,M)
 860        CONTINUE
          ENDIF
 870    CONTINUE
C
C  SEND THE INITIAL SLABS TO DISK
C
        do kk = 1, nt
          do jj = 1, km
            do ii = 1, imt
              odam_t(ii, j, jj, kk, 1) = tb(ii, jj, kk)
            end do
          end do
        end do
        do jj = 1, km
          do ii = 1, imt
            odam_u(ii, j, jj, 1) = ub(ii, jj)
          end do
        end do
        do jj = 1, km
          do ii = 1, imt
            odam_v(ii, j, jj, 1) = vb(ii, jj)
          end do
        end do
        do kk = 1, nt
          do jj = 1, km
            do ii = 1, imt
              odam_t(ii, j, jj, kk, 2) = t(ii, jj, kk)
            end do
          end do
        end do
        do jj = 1, km
          do ii = 1, imt
            odam_u(ii, j, jj, 2) = u(ii, jj)
          end do
        end do
        do jj = 1, km
          do ii = 1, imt
            odam_v(ii, j, jj, 2) = v(ii, jj)
          end do
        end do
 880  CONTINUE
C
C---------------------------------------------------------------------
C   INITIALIZE REMAINDER OF DISK
C---------------------------------------------------------------------
C
C  SET INITIAL STREAM FUNCTION TO ZERO
C
      DO 890 J=1,JMT
      DO 890 I=1,IMT
        PB(I,J)=0.0
        P (I,J)=0.0
 890  CONTINUE
C
C  Set relaxation solutions to zero
C
      do kk = 1, 2
        do j = 1, jmt
          do i = 1, imt
            ptd2(i, j, k) = 0.0
          end do
        end do
      end do
C
C  SET TIMESTEP COUNTER AND TOTAL ELAPSED TIME TO ZERO
C
      ITT=0
      TTSEC=0.0
C
C=======================================================================
C  END SECTION TO START FROM SCRATCH  ==================================
C=======================================================================
C
  160 CONTINUE
C=======================================================================
C  BEGIN SECTION TO TIMESTEP THE MODEL  ================================
C=======================================================================
C

c rjm...
      if (nt.gt.2) then
# 1612

      print*,"rjm call bio_i"
      call bio_i

      endif
c rjm

C************** PRINT START AND END TIMES AND RENUMBER RESTART FILE ****
C************** AFTER RESETTING TIME COUNTERS IF DESIRED  **************


C************** PRINT START AND END TIMES ******************************
C************** AFTER RESETTING TIME COUNTERS IF DESIRED  **************
      If (.not.lcouple) Then
      IF(ITSET.EQ.1)THEN
        ITT = 0
        TTSEC = 0.
      ENDIF
      End If
C
C  COMPUTE PERMUTING DISC INDICATORS
C
      NDISK =MOD(ITT+0,2)+1
      NDISKA=MOD(ITT+1,2)+1
      ittx=0

      if(lcouple)then
c****************************************************************
C***** For coupled model only
c****************************************************************
c----
c---- Here the counter for the ocean model is made to coincide
c---- with the atmospheric model counter in coupled mode.
c----
        ittx=ndays*knitd
        if(ittx.ne.ITT)then
          print *,'****************** Warning *********************'
          print *,'Coupled model mode :'
          print *,'Changing ITT in ocean model to match atmos model'
          print *,'ITT was ',ITT,' , now changed to ',ittx
          print *,'****************** Warning *********************'
          ITT=ittx
          ndk=MOD(ITT+0,2)+1
          ndka=MOD(ITT+1,2)+1
          if((ndk.ne.NDISK).or.(ndka.ne.NDISKA))then
            print *,'Curious - permuting indicators wrong - stop'
            stop
          end if
          TTSEC=ndays*3600.0*24.0
        end if
c****************************************************************
c****************************************************************
c****************************************************************
      end if

c.......set itfs= restart timestep.....
        itfs=ITT
      WRITE(*,1317)ITT,itfs
 1317 FORMAT(1X,'ITT=',I10,' itfs=',i10)
C
C---------------------------------------------------------------------
C  INITIALIZE SEVERAL VARIABLES TO ZERO TO AVOID AN "UNINITIALIZED
C       VARIABLE" TYPE OF ERROR LATER WHERE, FOR PURPOSES OF VEC-
C       TORIZATION, THE COMPUTATION PROCEEDS ACROSS LAND POINTS
C---------------------------------------------------------------------
C
      DO 168 I=1,IMT
        UUNDER(I)=0.0
        VUNDER(I)=0.0
        DO 166 J=1,JMT
          ZTD(I,J)=0.0
 166    CONTINUE
        DO 167 K=1,KMP1
          TEMPA(I,K)=0.0
          TEMPB(I,K)=0.0
 167    CONTINUE
      DO 168 M=1,NTMIN2
      DO 168 K=1,KMP2
        TDIF(I,K,M)=0.0
 168  CONTINUE

c---- Set indicator for startup data read (see subroutine STEP)
      lmon1o=-1
      lmon1n=0

      print*, " "
      print*, "OCEAN SETUP COMPLETE"
      print*, " "

      RETURN
      END
