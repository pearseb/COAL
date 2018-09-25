C     *==========================================================*
C     | bulkf_const.h
C     | o Header file for bulk formula:
C     |   - basic parameter ( I/O frequency, etc ...)
C     |   - physical constants
C     | ***Adapted from MITgcm source code***
C     *==========================================================*

C ---------------- CONSTANTS ----------------
C.. densities
C     rhoA      ::  density of air [kg/m^3]
C     rhoFW     ::  density of fresh water [kg/m^3]
C     p0        ::  reference sea-level atmospheric pressure [mb]
C.. specific heats
C     cpair     ::  specific heat of air [J/kg/K]
C.. latent heat
C     Lvap      ::  latent heat of vaporization at 0.oC [J/kg]
C     Lfresh    ::  latent heat of melting of pure ice [J/kg]
C.. constants
C     Rgas      ::  gas constant for dry air   [J/kg/K]
C     karman      ::  von Karman constant  [-]
C     sboltz    ::  Stefan-Boltzmann constant [W/m^2/K^4]
C.. for transfer coefficient
C     zref      :: reference height [m] for transfer coefficient
C     zrou      :: roughness length scale = 0.0005
C     zwd       :: height [m] of near-surface wind-speed input data
C     zth       :: height [m] of near-surface air-temp. & air-humid. input
C     cdrag_[n] ::  n = 1,2,3 coefficients used to evaluate drag coefficient
C     cStantonS,U :: coefficients used to evaluate Stanton number (for
C                 sensib. Heat Flx), under Stable / Unstable stratification
C     cDalton   :: coefficient used to evaluate Dalton number (for Evap)
C.. for bulk formula
C     umin      :: minimum wind speed used in bulk-formulae [m/s]
C     fhumid    :: dry-air - water-vapor molecular mass ratio (minus one)
C                    (used to calculate virtual temp.) = 0.606
C     saltQs    :: reduction of sat. vapor pressure over salty water
C     gamma_blk :: adiabatic lapse rate
C     ssq[0-2]  :: constants for calculation of saturated specific humidity
C     qs[1-2]w  :: constants for calculation of saturated specific humidity
C                   over open water
C     qs[1-2]i  :: constants for calculation of saturated specific humidity
C                   over ice 
C.. Albedo for shortwave radiation
C     alb_oce   :: ocean surface albedo [0-1]
C     alb_ice   :: ice surface albedo [0-1]
C.. for Long-Wave Radiation
C     emi_atm   ::
C     emi_oce   ::
C     emi_ice   ::
C.. Sea ice related constants, 
C     rdn_i     :: neutral transfer coefficients over ice (Cd = Ce = Ch) 
C     ocefrz    :: temperature of freezing point of seawater

      real :: rhoA = 1.22, 
     +        rhoFW = 1025.0,
     +        p0 = 1013.0,
     +        cpair = 1.005e3,
     +        Lvap = 2.5e6,
     +        Lfresh = 3.337e5,
     +        Rgas = 287.058,
     +        karman = 0.4,
     +        sboltz = 5.67e-8,
     +        gravity = 9.8,
     +        zref = 10.0, 
     +        zrou = 0.0005,
     +        zwd = 10.0,
     +        zth = 10.0,
     +        cdrag_1 = 2.7e-3,
     +        cdrag_2 = 0.142e-3,
     +        cdrag_3 = 0.0764e-3,
     +        cStantonS = 18.0e-3,
     +        cStantonU = 32.7e-3,
     +        cDalton = 34.6e-3,
     +        umin = 1.0,
     +        fhumid = 0.606,
     +        saltQs = 0.98,
     +        gamma_blk = 0.01,
     +        ssq0 = 3.797915,
     +        ssq1 = 7.93252e-6,
     +        ssq2 = 2.166847e-3,
     +        qs1w = 640.38e3,
     +        qs2w = -5107.4,
     +        qs1i = 11637.8,
     +        qs2i = -5897.8,
     +        alb_oce = 0.1,
     +        alb_ice = 0.7,
     +        emi_atm = 0.90,
     +        emi_oce = 0.985,
     +        emi_ice = 0.950,
     +        rdn_i = 1.63e-3,
     +        ocefrz = -1.85
     
C    niter_bulk :: Number of iterations to find turbulent transfer coeff.
      integer :: niter_bulk = 5



