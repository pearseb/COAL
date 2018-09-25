! cable_types.f90
!
! This file declares derived type variables for CABLE and
! contains the interface that allocates and deallocates them.
!
! Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach, 
! Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!
! Fortran-95 recoding by Harvey Davies, Gab Abramowitz, Martin Dix & Bernard Pak
!
! Please send bug reports to Bernard.Pak@csiro.au
!
! This file has the module define_types containing subroutines:
!      alloc_*_type, and
!      dealloc_*_type
! used via the interface alloc_cbm_var and dealloc_cbm_var respectively.
!
MODULE define_types
  ! Contains all variables which are not subroutine-internal
  USE define_dimensions, ONLY: r_1,r_2,i_d,ms,msn,ncp,ncs,nrb,mf
  IMPLICIT NONE
  PUBLIC
  PRIVATE r_1,r_2,i_d,ms,msn,ncp,ncs,nrb,mf
  ! Energy and water balance variables:
  TYPE balances_type 
    REAL(r_1), DIMENSION(:), POINTER :: drybal ! energy balance for dry canopy
    REAL(r_1), DIMENSION(:), POINTER :: ebal   ! energy balance per time step (W/m^2)
    REAL(r_1), DIMENSION(:), POINTER :: ebal_tot ! cumulative energy balance (W/m^2)
    REAL(r_1), DIMENSION(:), POINTER :: evap_tot ! cumulative evapotranspiration (mm/dels)
    REAL(r_1), DIMENSION(:), POINTER :: osnowd0  ! snow depth, first time step
    REAL(r_1), DIMENSION(:), POINTER :: precip_tot ! cumulative precipitation (mm/dels)
    REAL(r_1), DIMENSION(:), POINTER :: rnoff_tot  ! cumulative runoff (mm/dels)
    REAL(r_1), DIMENSION(:), POINTER :: wbal   ! water balance per time step (mm/dels)
    REAL(r_1), DIMENSION(:), POINTER :: wbal_tot ! cumulative water balance (mm/dels)
    REAL(r_1), DIMENSION(:), POINTER :: wbtot0 ! total soil water (mm), first time step
    REAL(r_1), DIMENSION(:), POINTER :: wetbal ! energy balance for wet canopy
  END TYPE balances_type
  ! Soil parameters:
  TYPE soil_parameter_type 
    REAL(r_1), DIMENSION(:,:), POINTER :: albsoil ! soil reflectance (second dimension, BP 21Oct2009)
    REAL(r_1), DIMENSION(:), POINTER :: bch  ! parameter b in Campbell equation
    REAL(r_1), DIMENSION(:), POINTER :: c3   ! c3 drainage coeff (fraction) (EK nov 2007)
    REAL(r_1), DIMENSION(:), POINTER :: clay ! fraction of soil which is clay
    REAL(r_1), DIMENSION(:), POINTER :: cnsd ! thermal conductivity of dry soil [W/m/K]
    REAL(r_1), DIMENSION(:), POINTER :: css  ! soil specific heat capacity [kJ/kg/K]
    REAL(r_1), DIMENSION(:), POINTER :: hsbh  ! difsat * etasat (=hyds*abs(sucs)*bch)
    REAL(r_1), DIMENSION(:), POINTER :: hyds  ! hydraulic conductivity @ saturation [m/s], Ksat
    INTEGER(i_d), DIMENSION(:), POINTER :: i2bp3  ! parameter one in K vis suction (=nint(bch)+2)
    INTEGER(i_d), DIMENSION(:), POINTER :: ibp2   ! parameter two in K vis suction (function of pbch)
    INTEGER(i_d), DIMENSION(:), POINTER :: isoilm ! integer soil type
    REAL(r_1), DIMENSION(:), POINTER :: rhosoil ! soil density [kg/m3]
    REAL(r_1), DIMENSION(:), POINTER :: rs20  ! soil respiration at 20 C [mol m-2 s-1]
    REAL(r_1), DIMENSION(:), POINTER :: sand  ! fraction of soil which is sand
    REAL(r_1), DIMENSION(:), POINTER :: sfc   ! vol H2O @ field capacity
    REAL(r_1), DIMENSION(:), POINTER :: silt  ! fraction of soil which is silt
    REAL(r_1), DIMENSION(:), POINTER :: ssat  ! vol H2O @ saturation
    REAL(r_1), DIMENSION(:), POINTER :: sucs  ! suction at saturation (m)
    REAL(r_1), DIMENSION(:), POINTER :: swilt ! vol H2O @ wilting
    REAL(r_1), DIMENSION(ms) :: zse   ! thickness of each soil layer (1=top) in m
    REAL(r_1), DIMENSION(ms+1) :: zshh ! distance between consecutive layer midpoints (m)
  END TYPE soil_parameter_type
  ! Soil and snow variables:
  TYPE soil_snow_type 
    REAL(r_1), DIMENSION(:,:), POINTER :: albsoilsn ! soil + snow reflectance
    REAL(r_1), DIMENSION(:), POINTER :: cls     ! factor for latent heat
    REAL(r_1), DIMENSION(:), POINTER :: dfn_dtg ! d(canopy%fns)/d(ssoil%tgg)
    REAL(r_1), DIMENSION(:), POINTER :: dfh_dtg ! d(canopy%fhs)/d(ssoil%tgg)
    REAL(r_1), DIMENSION(:), POINTER :: dfe_ddq ! d(canopy%fes)/d(dq)
    REAL(r_1), DIMENSION(:), POINTER :: ddq_dtg ! d(dq)/d(ssoil%tgg)
    REAL(r_1), DIMENSION(:), POINTER :: evapsn  ! snow evaporation
    REAL(r_1), DIMENSION(:), POINTER :: fwtop   ! water flux to the soil (EK nov 2007)
    REAL(r_2), DIMENSION(:,:), POINTER :: gammzz ! heat capacity for each soil layer
    INTEGER(i_d), DIMENSION(:), POINTER :: isflag ! 0 => no snow 1 => snow
    REAL(r_1), DIMENSION(:), POINTER :: osnowd  ! snow depth from previous time step
    REAL(r_1), DIMENSION(:), POINTER :: potev   ! potential evapotranspiration
    REAL(r_2), DIMENSION(:), POINTER :: pwb_min
    REAL(r_1), DIMENSION(:), POINTER :: runoff  ! total runoff (mm/dels)
    REAL(r_1), DIMENSION(:), POINTER :: rnof1   ! surface runoff (mm/dels)
    REAL(r_1), DIMENSION(:), POINTER :: rnof2   ! deep drainage (mm/dels)
    REAL(r_1), DIMENSION(:), POINTER :: rtsoil  ! turbulent resistance for soil
    REAL(r_1), DIMENSION(:,:), POINTER :: sconds ! EK nov 2007
    REAL(r_1), DIMENSION(:,:), POINTER :: sdepth ! snow depth
    REAL(r_1), DIMENSION(:,:), POINTER  :: smass ! snow mass
    REAL(r_1), DIMENSION(:), POINTER :: snage   ! snow age
    REAL(r_1), DIMENSION(:), POINTER :: snowd   ! snow depth (liquid water)
    REAL(r_1), DIMENSION(:), POINTER :: smelt   ! snow melt (EK nov 2007)
    REAL(r_1), DIMENSION(:,:), POINTER  :: ssdn ! snow densities
    REAL(r_1), DIMENSION(:), POINTER :: ssdnn   ! average snow density
    REAL(r_1), DIMENSION(:,:), POINTER :: tgg   ! soil temperature in K
    REAL(r_1), DIMENSION(:,:), POINTER  :: tggsn ! snow temperature in K
    REAL(r_1), DIMENSION(:), POINTER :: tss     ! surface temperature (weighted soil, snow)
    REAL(r_2), DIMENSION(:,:), POINTER :: wb    ! volumetric soil moisture (solid+liq)
    REAL(r_1), DIMENSION(:,:), POINTER :: wbfice ! fraction of ssat that is ice
    REAL(r_2), DIMENSION(:,:), POINTER :: wbice  ! volumentric soil ice
    REAL(r_2), DIMENSION(:,:), POINTER :: wblf  ! fraction of ssat that is liquid
    REAL(r_1), DIMENSION(:), POINTER :: wbtot   ! total soil water (mm)
    REAL(r_1), DIMENSION(:), POINTER :: wetfac ! surface wetness fact. at current time step
    REAL(r_1), DIMENSION(:), POINTER :: owetfac ! surface wetness fact. at previous time step
  END TYPE soil_snow_type
  ! Vegetation parameters:
  TYPE veg_parameter_type
    REAL(r_1), DIMENSION(:), POINTER :: canst1 ! max intercepted water by canopy (mm/LAI)
    REAL(r_1), DIMENSION(:), POINTER :: dleaf  ! chararacteristc legnth of leaf (m)
    REAL(r_1), DIMENSION(:), POINTER :: ejmax  ! max pot. electron transport rate top leaf(mol/m2/s)
    REAL(r_1), DIMENSION(:), POINTER :: frac4  ! fraction of c4 plants
    REAL(r_1), DIMENSION(:,:), POINTER :: froot  ! fraction of root in each soil layer
    REAL(r_1), DIMENSION(:), POINTER :: hc     ! roughness height of canopy
    INTEGER(i_d),DIMENSION(:), POINTER :: iveg ! vegetation type
    INTEGER(i_d),DIMENSION(:), POINTER :: meth ! method for calculation of canopy fluxes and temp.
    REAL(r_1), DIMENSION(:), POINTER :: rp20   ! plant respiration coefficient at 20 C
    REAL(r_1), DIMENSION(:), POINTER :: rpcoef ! temperature coef nonleaf plant respiration (1/C)
    REAL(r_1), DIMENSION(:), POINTER :: shelrb ! sheltering factor (dimensionless) 
    REAL(r_1), DIMENSION(:), POINTER :: wai   ! wood area index (stem+branches+twigs)
    REAL(r_1), DIMENSION(:), POINTER :: vegcf  ! biome-specific soil respiration rate
    REAL(r_1), DIMENSION(:), POINTER :: tminvj ! min temperature of the start of photosynthesis
    REAL(r_1), DIMENSION(:), POINTER :: tmaxvj ! max temperature of the start of photosynthesis
    REAL(r_1), DIMENSION(:), POINTER :: vbeta  ! stomatal sensitivity to soil water
    REAL(r_1), DIMENSION(:), POINTER :: xalbnir ! modifier for albedo in near ir band
!     REAL(r_1), DIMENSION(:), POINTER :: rootbeta  ! parameter for estimating
!                                    ! vertical root mass distribution (froot)
    REAL(r_1), DIMENSION(:), POINTER :: vcmax  ! maximum RuBP carboxylation rate top leaf (mol/m2/s)
    REAL(r_1), DIMENSION(:), POINTER :: vlai   ! leaf area index
    REAL(r_1), DIMENSION(:), POINTER :: xfang  ! leaf angle PARAMETER
    REAL(r_1), DIMENSION(:), POINTER :: extkn  ! extinction coef for vertical
                                                ! nitrogen profile in canopy(-)
!   rml 22/10/07
    LOGICAL,   DIMENSION(:), POINTER :: deciduous ! flag used for phenology fix
  END TYPE veg_parameter_type
  ! Canopy/vegetation variables:
  TYPE canopy_type
    REAL(r_1), DIMENSION(:), POINTER :: cansto ! canopy water storage (mm)
    REAL(r_1), DIMENSION(:), POINTER :: cduv  ! drag coefficient for momentum
    REAL(r_1), DIMENSION(:), POINTER :: delwc ! change in canopy water store (mm/dels)
    REAL(r_1), DIMENSION(:), POINTER :: dewmm ! dewfall (mm)
    REAL(r_2), DIMENSION(:), POINTER :: dgdtg ! derivative of gflux wrt soil temp
    REAL(r_1), DIMENSION(:), POINTER :: fe    ! total latent heat (W/m2)
    REAL(r_1), DIMENSION(:), POINTER :: fh    ! total sensible heat (W/m2)
    REAL(r_1), DIMENSION(:), POINTER :: fpn   ! plant photosynthesis (g C s-1)
    REAL(r_1), DIMENSION(:), POINTER :: frp   ! plant respiration (g C m-2 s-1)
    REAL(r_1), DIMENSION(:), POINTER :: frpw  ! plant respiration (g C m-2 s-1)???
    REAL(r_1), DIMENSION(:), POINTER :: frpr  ! plant respiration (g C m-2 s-1)???
    REAL(r_1), DIMENSION(:), POINTER :: frs   ! soil respiration (g C m-2 s-1)
    REAL(r_1), DIMENSION(:), POINTER :: fnee  ! net carbon flux (g C m-2 s-1)
    REAL(r_1), DIMENSION(:), POINTER :: frday ! daytime leaf resp
    REAL(r_1), DIMENSION(:), POINTER :: fnv   ! net rad. avail. to canopy (W/m2)
    REAL(r_1), DIMENSION(:), POINTER :: fev   ! latent hf from canopy (W/m2)
    REAL(r_2), DIMENSION(:), POINTER :: fevc  ! dry canopy transpiration (W/m2)
    REAL(r_2), DIMENSION(:), POINTER :: fevw  ! lat heat fl wet canopy (W/m2)
    REAL(r_2), DIMENSION(:), POINTER :: potev_c ! canopy potential evapotranspitation (YP & Mao, jun08)
    REAL(r_1), DIMENSION(:), POINTER :: fhv   ! sens heatfl from canopy (W/m2)
    REAL(r_2), DIMENSION(:), POINTER :: fhvw  ! sens heatfl from wet canopy (W/m2)
    REAL(r_1), DIMENSION(:), POINTER :: fns   ! net rad avail to soil (W/m2)
    REAL(r_1), DIMENSION(:), POINTER :: fes   ! latent heatfl from soil (W/m2)
    REAL(r_1), DIMENSION(:), POINTER :: fhs   ! sensible heat flux from soil
    REAL(r_1), DIMENSION(:), POINTER :: fwet   ! fraction of canopy wet
    REAL(r_1), DIMENSION(:), POINTER :: ga    ! ground heat flux (W/m2) ???
    REAL(r_1), DIMENSION(:), POINTER :: ghflux  ! ground heat flux (W/m2) ???
    REAL(r_1), DIMENSION(:), POINTER :: precis! throughfall to soil, after snow (mm)
    REAL(r_1), DIMENSION(:), POINTER :: qscrn ! specific humudity at screen height (g/g)
    REAL(r_1), DIMENSION(:), POINTER :: rnet  ! net radiation absorbed by surface (W/m2)
    REAL(r_1), DIMENSION(:), POINTER :: segg  ! latent heatfl from soil mm (EK nov 2007)
    REAL(r_1), DIMENSION(:), POINTER :: sghflux ! ground heat flux (W/m2) ???
    REAL(r_1), DIMENSION(:), POINTER :: spill ! can.storage excess after dewfall (mm)
    REAL(r_1), DIMENSION(:), POINTER :: through ! canopy throughfall (mm)
    REAL(r_1), DIMENSION(:), POINTER :: tscrn ! air temperature at screen height (oC)
    REAL(r_1), DIMENSION(:), POINTER :: tv    ! vegetation temp (K)
    REAL(r_1), DIMENSION(:), POINTER :: us    ! friction velocity
    REAL(r_1), DIMENSION(:), POINTER :: uscrn ! wind speed at screen height (m/s)
    REAL(r_1), DIMENSION(:), POINTER :: vlaiw  ! lai adjusted for snow depth for calculation of resistances
    REAL(r_1), DIMENSION(:), POINTER :: wcint ! canopy rainfall interception (mm)
!! ypw (26oct2010)
    REAL(r_1), DIMENSION(:,:), POINTER :: rwater
    REAL(r_2), DIMENSION(:,:), POINTER :: evapfbl
!! ypw (26oct2010)
  END TYPE canopy_type
  ! Radiation variables:
  TYPE radiation_type
    REAL(r_1), DIMENSION(:,:), POINTER :: albedo ! canopy+soil albedo
    REAL(r_1), DIMENSION(:), POINTER     :: extkb  ! beam radiation extinction coeff
    REAL(r_1), DIMENSION(:), POINTER     :: extkd2 ! diffuse 2D radiation extinction coeff
    REAL(r_1), DIMENSION(:), POINTER ::  extkd ! diffuse radiation extinction coeff (-)
    REAL(r_1), DIMENSION(:), POINTER :: flws   ! soil long-wave radiation
    REAL(r_1), DIMENSION(:,:), POINTER  :: fvlai  ! leaf area index of big leaf
    REAL(r_2), DIMENSION(:,:), POINTER  :: gradis ! radiative conductance
    REAL(r_1), DIMENSION(:), POINTER :: latitude  ! latitude
    REAL(r_1), DIMENSION(:), POINTER :: lwabv ! long wave absorbed by vegetation
    REAL(r_1), DIMENSION(:,:,:), POINTER :: qcan ! absorbed radiation for canopy (W/m^2)
    REAL(r_1), DIMENSION(:), POINTER     :: qssabs ! absorbed short-wave radiation for soil
    REAL(r_1), DIMENSION(:,:), POINTER ::  rhocdf ! canopy diffuse reflectance (-)
    REAL(r_1), DIMENSION(:,:), POINTER  :: rniso  !  sum(rad%qcan, 3) total abs by canopy (W/m2)
    REAL(r_1), DIMENSION(:,:), POINTER  :: scalex ! scaling PARAMETER for big leaf
    REAL(r_1), DIMENSION(:), POINTER     :: transd ! fraction SW diffuse transmitted through canopy
    REAL(r_1), DIMENSION(:), POINTER ::  trad  !  radiative temperature (soil and veg)
    ! new types, ypw 11/july/2008
    REAL(r_1),DIMENSION(:,:), POINTER  :: reffdf  !effective conopy diffuse reflectance
    REAL(r_1),DIMENSION(:,:), POINTER  :: reffbm  !effective conopy beam reflectance
    REAL(r_1), DIMENSION(:,:), POINTER :: extkbm  !modified k beam(6.20)(for leaf scattering)
    REAL(r_1), DIMENSION(:,:), POINTER :: extkdm  !modified k diffuse(6.20)(for leaf scattering)
    REAL(r_1), DIMENSION(:), POINTER   :: fbeam   !beam fraction
    REAL(r_1), DIMENSION(:,:),POINTER  :: cexpkbm ! canopy beam transmittance
    REAL(r_1), DIMENSION(:,:),POINTER  :: cexpkdm ! canopy diffuse transmittance
    ! ********* ypw 11/july/2008
  END TYPE radiation_type
  ! Roughness variables:
  TYPE roughness_type
    ! "coexp": coefficient in exponential in-canopy wind profile
    ! U(z) = U(h)*exp(coexp*(z/h-1)), found by gradient-matching
    ! canopy and roughness-sublayer U(z) at z=h
    REAL(r_1), DIMENSION(:), POINTER   :: coexp ! Extinction coefficient for wind profile in canopy
    REAL(r_1), DIMENSION(:), POINTER   :: disp  ! zero-plane displacement
    REAL(r_1), DIMENSION(:), POINTER   :: hruff ! canopy height above snow level
    REAL(r_1), DIMENSION(:), POINTER   :: hruff_grmx ! grid maximum of hruff
    REAL(r_1), DIMENSION(:), POINTER   :: rt0us ! eq. 3.54, SCAM manual (CSIRO tech report 132)
    REAL(r_1), DIMENSION(:), POINTER   :: rt1usa ! resistance from disp to hruf
    REAL(r_1), DIMENSION(:), POINTER   :: rt1usb ! resistance from hruf to zruffs (or zref if zref<zruffs)
    REAL(r_1), DIMENSION(:), POINTER   :: rt1 ! 1/aerodynamic conductance
    REAL(r_1), DIMENSION(:), POINTER   :: term2, term3, term5, term6 ! for aerodynamic resistance calc.
    ! "usuh": us/uh (us=friction velocity, uh = mean velocity at z=h)
    REAL(r_1), DIMENSION(:), POINTER   :: usuh ! Friction velocity/windspeed at canopy height
!    REAL(r_1), DIMENSION(:), POINTER   :: za   ! level of lowest atmospheric model layer
    REAL(r_1), DIMENSION(:), POINTER   :: za_uv ! level of lowest uv grid
    REAL(r_1), DIMENSION(:), POINTER   :: za_tq ! level of lowest tq grid
    REAL(r_1), DIMENSION(:), POINTER   :: z0m  ! roughness length
!    REAL(r_1), DIMENSION(:), POINTER   :: zref ! Reference height for met forcing
    REAL(r_1), DIMENSION(:), POINTER   :: zref_uv ! Ref height wrt uv grid
    REAL(r_1), DIMENSION(:), POINTER   :: zref_tq ! Ref height wrt tq grid
    REAL(r_1), DIMENSION(:), POINTER   :: zruffs ! SCALAR Roughness sublayer depth (ground=origin)
    REAL(r_1), DIMENSION(:), POINTER   :: z0soilsn ! roughness length of bare soil surface
    REAL(r_1), DIMENSION(:), POINTER   :: z0soil ! roughness length of bare soil surface
  END TYPE roughness_type
  ! Air variables:
  TYPE air_type
    REAL(r_1), DIMENSION(:), POINTER :: rho  ! dry air density (kg m-3)
    REAL(r_1), DIMENSION(:), POINTER :: volm ! molar volume (m3 mol-1)
    REAL(r_1), DIMENSION(:), POINTER :: rlam ! latent heat for water (j/kg)
    REAL(r_1), DIMENSION(:), POINTER :: qsat ! saturation specific humidity
    REAL(r_1), DIMENSION(:), POINTER :: epsi ! d(qsat)/dT ((kg/kg)/K)
    REAL(r_1), DIMENSION(:), POINTER :: visc ! air kinematic viscosity (m2/s)
    REAL(r_1), DIMENSION(:), POINTER :: psyc ! psychrometric constant
    REAL(r_1), DIMENSION(:), POINTER :: dsatdk ! d(es)/dT (mb/K)
    REAL(r_1), DIMENSION(:), POINTER :: cmolar ! conv. from m/s to mol/m2/s
  END TYPE air_type
  ! Flux variables from a homogenous patch (subset of a grid cell)
  TYPE patch_flux_type
    REAL(r_1), DIMENSION(:), POINTER :: fe ! latent heat (W/m2)
    REAL(r_1), DIMENSION(:), POINTER :: fh ! sensible heat (W/m2)
  END TYPE patch_flux_type
  ! Meterological data:
  TYPE met_type
    REAL(r_1), DIMENSION(:), POINTER :: ca   ! CO2 concentration (mol/mol)
    INTEGER(i_d), DIMENSION(:), POINTER :: year ! local time year AD 
    INTEGER(i_d), DIMENSION(:), POINTER :: moy  ! local time month of year 
    REAL(r_1), DIMENSION(:), POINTER :: doy  ! local time day of year = days since 0 hr 1st Jan 
    REAL(r_1), DIMENSION(:), POINTER :: hod  ! local hour of day
    REAL(r_1), DIMENSION(:), POINTER :: fsd  ! downward short-wave radiation (W/m2)
    REAL(r_1), DIMENSION(:), POINTER :: fld  ! downward long-wave radiation (W/m2)
    REAL(r_1), DIMENSION(:), POINTER :: precip  ! rainfall (liquid+solid)(mm/dels)
    REAL(r_1), DIMENSION(:), POINTER :: precip_s ! solid preipitation only (mm/dels) (EK nov 2007)
!    REAL(r_1), DIMENSION(:), POINTER :: tc    ! surface air temperature (oC)
    REAL(r_1), DIMENSION(:), POINTER :: tk    ! surface air temperature (oK)
    REAL(r_1), DIMENSION(:), POINTER :: tvair ! within canopy air temperature (oK)
    REAL(r_1), DIMENSION(:), POINTER :: tvrad ! radiative vegetation temperature (K)
    REAL(r_1), DIMENSION(:), POINTER :: pmb   ! surface air pressure (mbar)
    REAL(r_1), DIMENSION(:), POINTER :: ua    ! surface wind speed (m/s)
    REAL(r_1), DIMENSION(:), POINTER :: qv    ! surface specific humidity (g/g)
    REAL(r_1), DIMENSION(:), POINTER :: qvair ! within canopy specific humidity (g/g)
    REAL(r_1), DIMENSION(:), POINTER :: da    ! water vap pressure deficit at ref height (Pa)
    REAL(r_1), DIMENSION(:), POINTER :: dva   ! in canopy water vap pressure deficit (Pa)
    REAL(r_1), DIMENSION(:), POINTER :: coszen  ! cos(zenith angle of sun)
  END TYPE met_type
  ! Cumulative flux variables:
  TYPE sum_flux_type
    REAL(r_1), DIMENSION(:), POINTER :: sumpn  ! sum of canopy photosynthesis (g C m-2)
    REAL(r_1), DIMENSION(:), POINTER :: sumrp  ! sum of plant respiration (g C m-2)
    REAL(r_1), DIMENSION(:), POINTER :: sumrpw ! sum of plant respiration (g C m-2)
    REAL(r_1), DIMENSION(:), POINTER :: sumrpr ! sum of plant respiration (g C m-2)
    REAL(r_1), DIMENSION(:), POINTER :: sumrs  ! sum of soil respiration (g C m-2)
    REAL(r_1), DIMENSION(:), POINTER :: sumrd  ! sum of daytime respiration (g C m-2)
    REAL(r_1), DIMENSION(:), POINTER :: dsumpn ! daily sumpn
    REAL(r_1), DIMENSION(:), POINTER :: dsumrp ! daily sumrp
    REAL(r_1), DIMENSION(:), POINTER :: dsumrs ! daily sumrs
    REAL(r_1), DIMENSION(:), POINTER :: dsumrd ! daily sumrd
    REAL(r_1), DIMENSION(:), POINTER :: sumxrp ! sum plant resp. modifier
    REAL(r_1), DIMENSION(:), POINTER :: sumxrs ! sum soil resp. modifier
  END TYPE sum_flux_type
  TYPE bgc_pool_type
    REAL(r_1), DIMENSION(:,:), POINTER :: cplant ! plant carbon (g C/m2))
    REAL(r_1), DIMENSION(:,:), POINTER :: csoil  ! soil carbon (g C/m2)
    REAL(r_1), DIMENSION(ncp) :: ratecp ! plant carbon rate constant (1/year)
    REAL(r_1), DIMENSION(ncs) :: ratecs ! soil carbon rate constant (1/year)
  END TYPE bgc_pool_type

  ! Functions for allocating these types
  ! All overloaded so code only needs to call alloc_cbm_var
  ! Alloc routines could all initialise to NaN or zero for debugging?
  ! Don't need the mland/mp argument here as it's a module variable.
  PUBLIC :: alloc_cbm_var
  PRIVATE :: alloc_bgc_pool_type, dealloc_bgc_pool_type
  INTERFACE alloc_cbm_var
     MODULE PROCEDURE alloc_balances_type,            &
          alloc_soil_parameter_type,      &
          alloc_soil_snow_type,           &
          alloc_veg_parameter_type,       &
          alloc_canopy_type,              &
          alloc_radiation_type,           &
          alloc_roughness_type,           &
          alloc_air_type,                 &
          alloc_met_type,                 &
          alloc_sum_flux_type,            &
          alloc_bgc_pool_type
  END INTERFACE
  INTERFACE dealloc_cbm_var
     MODULE PROCEDURE dealloc_balances_type,            &
          dealloc_soil_parameter_type,      &
          dealloc_soil_snow_type,           &
          dealloc_veg_parameter_type,       &
          dealloc_canopy_type,              &
          dealloc_radiation_type,           &
          dealloc_roughness_type,           &
          dealloc_air_type,                 &
          dealloc_met_type,                 &
          dealloc_sum_flux_type,            &
          dealloc_bgc_pool_type
  END INTERFACE
CONTAINS
  
  SUBROUTINE alloc_balances_type(var, landunits)
    TYPE(balances_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var % drybal(landunits) )
    ALLOCATE ( var % ebal(landunits) )
    ALLOCATE ( var % ebal_tot(landunits) )
    ALLOCATE ( var % evap_tot(landunits) )
    ALLOCATE ( var % osnowd0(landunits) )
    ALLOCATE ( var % precip_tot(landunits) )
    ALLOCATE ( var % rnoff_tot(landunits) )
    ALLOCATE ( var % wbal(landunits) )
    ALLOCATE ( var % wbal_tot(landunits) )
    ALLOCATE ( var % wbtot0(landunits) )
    ALLOCATE ( var % wetbal(landunits) )
  END SUBROUTINE alloc_balances_type

  SUBROUTINE alloc_soil_parameter_type(var, landunits)
    TYPE(soil_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var % albsoil(landunits,nrb) )  ! (BP 21Oct2009)
    ALLOCATE ( var % bch(landunits) )
    ALLOCATE ( var % c3(landunits) )
    ALLOCATE ( var % clay(landunits) )
    ALLOCATE ( var % cnsd(landunits) )
    ALLOCATE ( var % css(landunits) )
    ALLOCATE ( var % hsbh(landunits) )
    ALLOCATE ( var % hyds(landunits) )
    ALLOCATE ( var % i2bp3(landunits) )
    ALLOCATE ( var % ibp2(landunits) )
    ALLOCATE ( var % isoilm(landunits) )
    ALLOCATE ( var % rhosoil(landunits) )
    ALLOCATE ( var % rs20(landunits) )
    ALLOCATE ( var % sand(landunits) )
    ALLOCATE ( var % sfc(landunits) )
    ALLOCATE ( var % silt(landunits) )
    ALLOCATE ( var % ssat(landunits) )
    ALLOCATE ( var % sucs(landunits) )
    ALLOCATE ( var % swilt(landunits) )
  END SUBROUTINE alloc_soil_parameter_type
 
  SUBROUTINE alloc_soil_snow_type(var, landunits)
    TYPE(soil_snow_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var % albsoilsn(landunits,nrb) )
    ALLOCATE ( var % cls(landunits) )
    ALLOCATE ( var % dfn_dtg(landunits) )
    ALLOCATE ( var % dfh_dtg(landunits) )
    ALLOCATE ( var % dfe_ddq(landunits) )
    ALLOCATE ( var % ddq_dtg(landunits) )
    ALLOCATE ( var % evapsn(landunits) )
    ALLOCATE ( var % fwtop(landunits) )
    ALLOCATE ( var % gammzz(landunits,ms) )
    ALLOCATE ( var % isflag(landunits) )
    ALLOCATE ( var % osnowd(landunits) )
    ALLOCATE ( var % potev(landunits) )
    ALLOCATE ( var % pwb_min(landunits) )
    ALLOCATE ( var % runoff(landunits) )
    ALLOCATE ( var % rnof1(landunits) )
    ALLOCATE ( var % rnof2(landunits) )
    ALLOCATE ( var % rtsoil(landunits) )
    ALLOCATE ( var % sconds(landunits,msn) )
    ALLOCATE ( var % sdepth(landunits,msn) )
    ALLOCATE ( var % smass(landunits,msn) )
    ALLOCATE ( var % snage(landunits) )
    ALLOCATE ( var % snowd(landunits) )
    ALLOCATE ( var % smelt(landunits) )
    ALLOCATE ( var % ssdn(landunits,msn) )
    ALLOCATE ( var % ssdnn(landunits) )
    ALLOCATE ( var % tgg(landunits,ms) )
    ALLOCATE ( var % tggsn(landunits,msn) )
    ALLOCATE ( var % tss(landunits) )
    ALLOCATE ( var % wb(landunits,ms) )
    ALLOCATE ( var % wbfice(landunits,ms) )
    ALLOCATE ( var % wbice(landunits,ms) )
    ALLOCATE ( var % wblf(landunits,ms) )
    ALLOCATE ( var % wbtot(landunits) )
    ALLOCATE ( var % wetfac(landunits) )
    ALLOCATE ( var % owetfac(landunits) )
  END SUBROUTINE alloc_soil_snow_type
   
  SUBROUTINE alloc_veg_parameter_type(var, landunits)
    TYPE(veg_parameter_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var % iveg(landunits) )
    ALLOCATE ( var % meth(landunits) )
    ALLOCATE ( var % vlai(landunits) )
    ALLOCATE ( var % froot(landunits,ms) )
    ALLOCATE ( var % canst1(landunits) )
    ALLOCATE ( var % ejmax(landunits) )
    ALLOCATE ( var % frac4(landunits) )
    ALLOCATE ( var % wai(landunits) )  ! new addition in Oct 2007 (YP)
    ALLOCATE ( var % vegcf(landunits) )  ! new addition in Oct 2007 (YP)
    ALLOCATE ( var % tminvj(landunits) )
    ALLOCATE ( var % tmaxvj(landunits) )
    ALLOCATE ( var % vbeta(landunits) )
    ALLOCATE ( var % xalbnir(landunits) )
!    ALLOCATE ( var % rootbeta(landunits) )  ! new addition in Oct 2007 (YP)
    ALLOCATE ( var % hc(landunits) )
    ALLOCATE ( var % shelrb(landunits) )
    ALLOCATE ( var % vcmax(landunits) )
    ALLOCATE ( var % xfang(landunits) )
    ALLOCATE ( var % dleaf(landunits) )
    ALLOCATE ( var % rp20(landunits) )
    ALLOCATE ( var % rpcoef(landunits) )
    ALLOCATE ( var % extkn(landunits) )  ! new addition in Oct 2007 (YP)
    ALLOCATE ( var % deciduous(landunits) )   ! rml addition 22/10/07
  END SUBROUTINE alloc_veg_parameter_type
   
  SUBROUTINE alloc_canopy_type(var, landunits)
    TYPE(canopy_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var % cansto(landunits) )
    ALLOCATE ( var % delwc(landunits) )
    ALLOCATE ( var % dewmm(landunits) )
    ALLOCATE ( var % fe(landunits) )
    ALLOCATE ( var % fh(landunits) )
    ALLOCATE ( var % fpn(landunits) )
    ALLOCATE ( var % frp(landunits) )
    ALLOCATE ( var % frpw(landunits) )
    ALLOCATE ( var % frpr(landunits) )
    ALLOCATE ( var % frs(landunits) )
    ALLOCATE ( var % fnee(landunits) )
    ALLOCATE ( var % frday(landunits) )
    ALLOCATE ( var % fnv(landunits) )
    ALLOCATE ( var % fev(landunits) )
    ALLOCATE ( var % fevc(landunits) )
    ALLOCATE ( var % fevw(landunits) )
    ALLOCATE ( var % fhv(landunits) )
    ALLOCATE ( var % fhvw(landunits) )
    ALLOCATE ( var % fns(landunits) )
    ALLOCATE ( var % fes(landunits) )
    ALLOCATE ( var % fhs(landunits) )
    ALLOCATE ( var % vlaiw(landunits) )
    ALLOCATE ( var % fwet(landunits) )
    ALLOCATE ( var % tv(landunits) )
    ALLOCATE ( var % ga(landunits) )
    ALLOCATE ( var % ghflux(landunits) )
    ALLOCATE ( var % segg(landunits) )
    ALLOCATE ( var % sghflux(landunits) )
    ALLOCATE ( var % dgdtg(landunits) )
    ALLOCATE ( var % through(landunits) )
    ALLOCATE ( var % precis(landunits) )
    ALLOCATE ( var % rnet(landunits) )
    ALLOCATE ( var % spill(landunits) )
    ALLOCATE ( var % wcint(landunits) )
    ALLOCATE ( var % us(landunits) )
    ALLOCATE ( var % tscrn(landunits) )
    ALLOCATE ( var % qscrn(landunits) )
    ALLOCATE ( var % uscrn(landunits) )
    ALLOCATE ( var % cduv(landunits) )
    ALLOCATE ( var % potev_c(landunits) )
!! ypw (26oct2010)
    ALLOCATE ( var % rwater(landunits,ms) )
    ALLOCATE ( var % evapfbl(landunits,ms) )
!! ypw (26oct2010)
  END SUBROUTINE alloc_canopy_type
   
  SUBROUTINE alloc_radiation_type(var, landunits)
    TYPE(radiation_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var % albedo(landunits,nrb) )
    ALLOCATE ( var % extkb(landunits) )
    ALLOCATE ( var % extkd2(landunits) )
    ALLOCATE ( var % extkd(landunits) )
    ALLOCATE ( var % flws(landunits) )
    ALLOCATE ( var % fvlai(landunits,mf) )
    ALLOCATE ( var % gradis(landunits,mf) )
    ALLOCATE ( var % latitude(landunits) )
    ALLOCATE ( var % lwabv(landunits) )
    ALLOCATE ( var % qcan(landunits,mf,nrb) )
    ALLOCATE ( var % qssabs(landunits) )
    ALLOCATE ( var % rhocdf(landunits,nrb) )
    ALLOCATE ( var % rniso(landunits,mf) )
    ALLOCATE ( var % scalex(landunits,mf) )
    ALLOCATE ( var % transd(landunits) )
    ALLOCATE ( var % trad(landunits) )
    ! new types, ypw 11/july/2008
    ALLOCATE ( var % reffdf(landunits,nrb) )
    ALLOCATE ( var % reffbm(landunits,nrb) )
    ALLOCATE ( var % extkbm(landunits,nrb) )
    ALLOCATE ( var % extkdm(landunits,nrb) )
    ALLOCATE ( var % fbeam(landunits) )
    ALLOCATE ( var % cexpkbm(landunits,nrb) )
    ALLOCATE ( var % cexpkdm(landunits,nrb) )
    ! ********** ypw 11/july/2008
  END SUBROUTINE alloc_radiation_type
   
  SUBROUTINE alloc_roughness_type(var, landunits)
    TYPE(roughness_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var % coexp(landunits) )
    ALLOCATE ( var % disp(landunits) )
    ALLOCATE ( var % hruff(landunits) )
    ALLOCATE ( var % hruff_grmx(landunits) )
    ALLOCATE ( var % rt0us(landunits) )
    ALLOCATE ( var % rt1usa(landunits) )
    ALLOCATE ( var % rt1usb(landunits) )
    ALLOCATE ( var % rt1(landunits) )
    ALLOCATE ( var % term2(landunits) )
    ALLOCATE ( var % term3(landunits) )
    ALLOCATE ( var % term5(landunits) )
    ALLOCATE ( var % term6(landunits) )
    ALLOCATE ( var % usuh(landunits) )
!    ALLOCATE ( var % za(landunits) )
    ALLOCATE ( var % za_uv(landunits) )
    ALLOCATE ( var % za_tq(landunits) )
    ALLOCATE ( var % z0m(landunits) )
!    ALLOCATE ( var % zref(landunits) )
    ALLOCATE ( var % zref_uv(landunits) )
    ALLOCATE ( var % zref_tq(landunits) )
    ALLOCATE ( var % zruffs(landunits) )
    ALLOCATE ( var % z0soilsn(landunits) )
    ALLOCATE ( var % z0soil(landunits) )
  END SUBROUTINE alloc_roughness_type
   
  SUBROUTINE alloc_air_type(var, landunits)
    TYPE(air_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var % rho(landunits) )
    ALLOCATE ( var % volm(landunits) )
    ALLOCATE ( var % rlam(landunits) )
    ALLOCATE ( var % qsat(landunits) )
    ALLOCATE ( var % epsi(landunits) )
    ALLOCATE ( var % visc(landunits) )
    ALLOCATE ( var % psyc(landunits) )
    ALLOCATE ( var % dsatdk(landunits) )
    ALLOCATE ( var % cmolar(landunits) )
  END SUBROUTINE alloc_air_type
   
  SUBROUTINE alloc_met_type(var, landunits)
    TYPE(met_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var % ca(landunits) )
    ALLOCATE ( var % year(landunits) )
    ALLOCATE ( var % moy(landunits) )
    ALLOCATE ( var % doy(landunits) )
    ALLOCATE ( var % hod(landunits) )
    ALLOCATE ( var % fsd(landunits) )
    ALLOCATE ( var % fld(landunits) )
    ALLOCATE ( var % precip(landunits) )
    ALLOCATE ( var % precip_s(landunits) )
!    ALLOCATE ( var % tc(landunits) )
    ALLOCATE ( var % tk(landunits) )
    ALLOCATE ( var % tvair(landunits) )
    ALLOCATE ( var % tvrad(landunits) )
    ALLOCATE ( var % pmb(landunits) )
    ALLOCATE ( var % ua(landunits) )
    ALLOCATE ( var % qv(landunits) )
    ALLOCATE ( var % qvair(landunits) )
    ALLOCATE ( var % da(landunits) )
    ALLOCATE ( var % dva(landunits) )
    ALLOCATE ( var % coszen(landunits) )
  END SUBROUTINE alloc_met_type
   
  SUBROUTINE alloc_sum_flux_type(var, landunits)
    TYPE(sum_flux_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var % sumpn(landunits) )
    ALLOCATE ( var % sumrp(landunits) )
    ALLOCATE ( var % sumrpw(landunits) )
    ALLOCATE ( var % sumrpr(landunits) )
    ALLOCATE ( var % sumrs(landunits) )
    ALLOCATE ( var % sumrd(landunits) )
    ALLOCATE ( var % dsumpn(landunits) )
    ALLOCATE ( var % dsumrp(landunits) )
    ALLOCATE ( var % dsumrs(landunits) )
    ALLOCATE ( var % dsumrd(landunits) )
    ALLOCATE ( var % sumxrp(landunits) )
    ALLOCATE ( var % sumxrs(landunits) )
  END SUBROUTINE alloc_sum_flux_type

  SUBROUTINE alloc_bgc_pool_type(var, landunits)
    TYPE(bgc_pool_type), INTENT(inout) :: var
    INTEGER, INTENT(in) :: landunits
    ALLOCATE ( var % cplant(landunits,ncp) )
    ALLOCATE ( var % csoil(landunits,ncs) )
  END SUBROUTINE alloc_bgc_pool_type

  ! Begin deallocation routines:
   SUBROUTINE dealloc_balances_type(var)
    TYPE(balances_type), INTENT(inout) :: var
    DEALLOCATE ( var % drybal )
    DEALLOCATE ( var % ebal )
    DEALLOCATE ( var % ebal_tot )
    DEALLOCATE ( var % evap_tot )
    DEALLOCATE ( var % osnowd0 )
    DEALLOCATE ( var % precip_tot )
    DEALLOCATE ( var % rnoff_tot )
    DEALLOCATE ( var % wbal )
    DEALLOCATE ( var % wbal_tot )
    DEALLOCATE ( var % wbtot0 )
    DEALLOCATE ( var % wetbal )
  END SUBROUTINE dealloc_balances_type

  SUBROUTINE dealloc_soil_parameter_type(var)
    TYPE(soil_parameter_type), INTENT(inout) :: var
    DEALLOCATE ( var % albsoil )
    DEALLOCATE ( var % bch )
    DEALLOCATE ( var % c3 )
    DEALLOCATE ( var % clay )
    DEALLOCATE ( var % cnsd )
    DEALLOCATE ( var % css )
    DEALLOCATE ( var % hsbh )
    DEALLOCATE ( var % hyds )
    DEALLOCATE ( var % i2bp3 )
    DEALLOCATE ( var % ibp2 )
    DEALLOCATE ( var % isoilm )
    DEALLOCATE ( var % rhosoil )
    DEALLOCATE ( var % rs20 )
    DEALLOCATE ( var % sand )
    DEALLOCATE ( var % sfc )
    DEALLOCATE ( var % silt )
    DEALLOCATE ( var % ssat )
    DEALLOCATE ( var % sucs )
    DEALLOCATE ( var % swilt )
  END SUBROUTINE dealloc_soil_parameter_type
 
  SUBROUTINE dealloc_soil_snow_type(var)
    TYPE(soil_snow_type), INTENT(inout) :: var
    DEALLOCATE ( var % albsoilsn )
    DEALLOCATE ( var % cls )
    DEALLOCATE ( var % dfn_dtg )
    DEALLOCATE ( var % dfh_dtg )
    DEALLOCATE ( var % dfe_ddq )
    DEALLOCATE ( var % ddq_dtg )
    DEALLOCATE ( var % evapsn )
    DEALLOCATE ( var % fwtop )
    DEALLOCATE ( var % gammzz )
    DEALLOCATE ( var % isflag )
    DEALLOCATE ( var % osnowd )
    DEALLOCATE ( var % potev )
    DEALLOCATE ( var % pwb_min )
    DEALLOCATE ( var % runoff )
    DEALLOCATE ( var % rnof1 )
    DEALLOCATE ( var % rnof2 )
    DEALLOCATE ( var % rtsoil )
    DEALLOCATE ( var % sconds )
    DEALLOCATE ( var % sdepth )
    DEALLOCATE ( var % smass )
    DEALLOCATE ( var % snage )
    DEALLOCATE ( var % snowd )
    DEALLOCATE ( var % smelt )
    DEALLOCATE ( var % ssdn )
    DEALLOCATE ( var % ssdnn )
    DEALLOCATE ( var % tgg )
    DEALLOCATE ( var % tggsn )
    DEALLOCATE ( var % tss )
    DEALLOCATE ( var % wb )
    DEALLOCATE ( var % wbfice )
    DEALLOCATE ( var % wbice )
    DEALLOCATE ( var % wblf )
    DEALLOCATE ( var % wbtot )
    DEALLOCATE ( var % wetfac )
    DEALLOCATE ( var % owetfac )
  END SUBROUTINE dealloc_soil_snow_type
   
  SUBROUTINE dealloc_veg_parameter_type(var)
    TYPE(veg_parameter_type), INTENT(inout) :: var
    DEALLOCATE ( var % iveg )
    DEALLOCATE ( var % meth )
    DEALLOCATE ( var % vlai )
    DEALLOCATE ( var % froot )
    DEALLOCATE ( var % canst1 )
    DEALLOCATE ( var % ejmax )
    DEALLOCATE ( var % frac4 )
    DEALLOCATE ( var % wai )  ! new addition in Oct 2007 (YP)
    DEALLOCATE ( var % vegcf )  ! new addition in Oct 2007 (YP)
    DEALLOCATE ( var % tminvj )
    DEALLOCATE ( var % tmaxvj )
    DEALLOCATE ( var % vbeta )
    DEALLOCATE ( var % xalbnir )
!    DEALLOCATE ( var % rootbeta )  ! new addition in Oct 2007 (YP)
    DEALLOCATE ( var % hc )
    DEALLOCATE ( var % shelrb )
    DEALLOCATE ( var % vcmax )
    DEALLOCATE ( var % xfang )
    DEALLOCATE ( var % dleaf )
    DEALLOCATE ( var % rp20 )
    DEALLOCATE ( var % rpcoef )
    DEALLOCATE ( var % extkn )  ! new addition in Oct 2007 (YP)
    DEALLOCATE ( var % deciduous )  ! rml addition 22/10/07
  END SUBROUTINE dealloc_veg_parameter_type
   
  SUBROUTINE dealloc_canopy_type(var)
    TYPE(canopy_type), INTENT(inout) :: var
    DEALLOCATE ( var % cansto )
    DEALLOCATE ( var % delwc )
    DEALLOCATE ( var % dewmm )
    DEALLOCATE ( var % fe )
    DEALLOCATE ( var % fh )
    DEALLOCATE ( var % fpn )
    DEALLOCATE ( var % frp )
    DEALLOCATE ( var % frpw )
    DEALLOCATE ( var % frpr )
    DEALLOCATE ( var % frs )
    DEALLOCATE ( var % fnee )
    DEALLOCATE ( var % frday )
    DEALLOCATE ( var % fnv )
    DEALLOCATE ( var % fev )
    DEALLOCATE ( var % fevc )
    DEALLOCATE ( var % fevw )
    DEALLOCATE ( var % fhv )
    DEALLOCATE ( var % fhvw )
    DEALLOCATE ( var % fns )
    DEALLOCATE ( var % fes )
    DEALLOCATE ( var % fhs )
    DEALLOCATE ( var % vlaiw )
    DEALLOCATE ( var % fwet )
    DEALLOCATE ( var % tv )
    DEALLOCATE ( var % ga )
    DEALLOCATE ( var % ghflux )
    DEALLOCATE ( var % segg )
    DEALLOCATE ( var % sghflux )
    DEALLOCATE ( var % dgdtg )
    DEALLOCATE ( var % through )
    DEALLOCATE ( var % precis )
    DEALLOCATE ( var % rnet )
    DEALLOCATE ( var % spill )
    DEALLOCATE ( var % wcint )
    DEALLOCATE ( var % us )
    DEALLOCATE ( var % tscrn )
    DEALLOCATE ( var % qscrn )
    DEALLOCATE ( var % uscrn )
    DEALLOCATE ( var % cduv )
    DEALLOCATE ( var % potev_c )
!! ypw (26oct2010)
    DEALLOCATE ( var % rwater )
    DEALLOCATE ( var % evapfbl )
!! ypw (26oct2010)
  END SUBROUTINE dealloc_canopy_type
   
  SUBROUTINE dealloc_radiation_type(var)
    TYPE(radiation_type), INTENT(inout) :: var
    DEALLOCATE ( var % albedo )
    DEALLOCATE ( var % extkb )
    DEALLOCATE ( var % extkd2 )
    DEALLOCATE ( var % extkd )
    DEALLOCATE ( var % flws )
    DEALLOCATE ( var % fvlai )
    DEALLOCATE ( var % gradis )
    DEALLOCATE ( var % latitude )
    DEALLOCATE ( var % lwabv )
    DEALLOCATE ( var % qcan )
    DEALLOCATE ( var % qssabs )
    DEALLOCATE ( var % rhocdf )
    DEALLOCATE ( var % rniso )
    DEALLOCATE ( var % scalex )
    DEALLOCATE ( var % transd )
    DEALLOCATE ( var % trad )
    ! new types, ypw 11/july/2008
    DEALLOCATE ( var % reffdf )
    DEALLOCATE ( var % reffbm )
    DEALLOCATE ( var % extkbm )
    DEALLOCATE ( var % extkdm )
    DEALLOCATE ( var % fbeam )
    DEALLOCATE ( var % cexpkbm )
    DEALLOCATE ( var % cexpkdm )
    ! ********** ypw 11/july/2008
  END SUBROUTINE dealloc_radiation_type
   
  SUBROUTINE dealloc_roughness_type(var)
    TYPE(roughness_type), INTENT(inout) :: var
    DEALLOCATE ( var % coexp )
    DEALLOCATE ( var % disp )
    DEALLOCATE ( var % hruff )
    DEALLOCATE ( var % hruff_grmx )
    DEALLOCATE ( var % rt0us )
    DEALLOCATE ( var % rt1usa )
    DEALLOCATE ( var % rt1usb )
    DEALLOCATE ( var % rt1 )
    DEALLOCATE ( var % term2 )
    DEALLOCATE ( var % term3 )
    DEALLOCATE ( var % term5 )
    DEALLOCATE ( var % term6 )
    DEALLOCATE ( var % usuh )
!    DEALLOCATE ( var % za )
    DEALLOCATE ( var % za_uv )
    DEALLOCATE ( var % za_tq )
    DEALLOCATE ( var % z0m )
!    DEALLOCATE ( var % zref )
    DEALLOCATE ( var % zref_uv )
    DEALLOCATE ( var % zref_tq )
    DEALLOCATE ( var % zruffs )
    DEALLOCATE ( var % z0soilsn )
    DEALLOCATE ( var % z0soil )
  END SUBROUTINE dealloc_roughness_type
   
  SUBROUTINE dealloc_air_type(var)
    TYPE(air_type), INTENT(inout) :: var
    DEALLOCATE ( var % rho )
    DEALLOCATE ( var % volm )
    DEALLOCATE ( var % rlam )
    DEALLOCATE ( var % qsat )
    DEALLOCATE ( var % epsi )
    DEALLOCATE ( var % visc )
    DEALLOCATE ( var % psyc )
    DEALLOCATE ( var % dsatdk )
    DEALLOCATE ( var % cmolar )
  END SUBROUTINE dealloc_air_type
   
  SUBROUTINE dealloc_met_type(var)
    TYPE(met_type), INTENT(inout) :: var
    DEALLOCATE ( var % ca )
    DEALLOCATE ( var % year )
    DEALLOCATE ( var % moy )
    DEALLOCATE ( var % doy )
    DEALLOCATE ( var % hod )
    DEALLOCATE ( var % fsd )
    DEALLOCATE ( var % fld )
    DEALLOCATE ( var % precip )
    DEALLOCATE ( var % precip_s )
!    DEALLOCATE ( var % tc )
    DEALLOCATE ( var % tk )
    DEALLOCATE ( var % tvair )
    DEALLOCATE ( var % tvrad )
    DEALLOCATE ( var % pmb )
    DEALLOCATE ( var % ua )
    DEALLOCATE ( var % qv )
    DEALLOCATE ( var % qvair )
    DEALLOCATE ( var % da )
    DEALLOCATE ( var % dva )
    DEALLOCATE ( var % coszen )
  END SUBROUTINE dealloc_met_type
   
  SUBROUTINE dealloc_sum_flux_type(var)
    TYPE(sum_flux_type), INTENT(inout) :: var
    DEALLOCATE ( var % sumpn )
    DEALLOCATE ( var % sumrp )
    DEALLOCATE ( var % sumrpw )
    DEALLOCATE ( var % sumrpr )
    DEALLOCATE ( var % sumrs )
    DEALLOCATE ( var % sumrd )
    DEALLOCATE ( var % dsumpn )
    DEALLOCATE ( var % dsumrp )
    DEALLOCATE ( var % dsumrs )
    DEALLOCATE ( var % dsumrd )
    DEALLOCATE ( var % sumxrp )
    DEALLOCATE ( var % sumxrs )
  END SUBROUTINE dealloc_sum_flux_type

  SUBROUTINE dealloc_bgc_pool_type(var)
    TYPE(bgc_pool_type), INTENT(inout) :: var
    DEALLOCATE ( var % cplant )
    DEALLOCATE ( var % csoil )
  END SUBROUTINE dealloc_bgc_pool_type

END MODULE define_types
