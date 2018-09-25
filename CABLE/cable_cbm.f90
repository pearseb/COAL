! cable_cbm.f90
!
! Source file containing main routine execution order for CABLE
!
! Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach, 
! Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!
! Fortran-95 recoding by Harvey Davies, Gab Abramowitz, Martin Dix & Bernard Pak
!
! Please send bug reports to Bernard.Pak@csiro.au
!
! This file contains modules:
!   cbm_module 
! The subroutines included are:
!   cbm
!
MODULE cbm_module
  USE canopy_module
  USE carbon_module
  USE soil_snow_module
  USE define_types
  USE define_dimensions, ONLY: i_d,r_1, mp
  USE physical_constants, ONLY: emsoil, sboltz
  USE roughness_module
  USE radiation_module
  USE albedo_module
  IMPLICIT NONE
  PRIVATE
  PUBLIC cbm 
CONTAINS

  SUBROUTINE cbm(ktau, kstart, kend, dels, air, bgc, canopy, met, &
       bal, rad, rough, soil, ssoil, sum_flux, veg)
! added mvtype and mstype to define_dimensions (BP sep2010)
!  SUBROUTINE cbm(ktau, kstart, kend, dels, air, bgc, canopy, met, &
!       bal, rad, rough, soil, ssoil, sum_flux, veg, mvtype, mstype)
    ! BP added mvtype and mstype to the list (dec 2007)
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number
    INTEGER(i_d), INTENT(IN)            :: kstart ! starting value of ktau
    INTEGER(i_d), INTENT(IN)            :: kend ! total # timesteps in run
    REAL(r_1), INTENT(IN)               :: dels ! time setp size (s)
    TYPE (air_type), INTENT(INOUT)      :: air
    TYPE (bgc_pool_type), INTENT(INOUT) :: bgc
    TYPE (canopy_type), INTENT(INOUT)   :: canopy
    TYPE (met_type), INTENT(INOUT)      :: met
    TYPE (balances_type), INTENT(INOUT) :: bal
    TYPE (radiation_type), INTENT(INOUT)      :: rad
    TYPE (roughness_type), INTENT(INOUT)      :: rough
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE (sum_flux_type), INTENT(INOUT)       :: sum_flux
    TYPE (veg_parameter_type), INTENT(INOUT)  :: veg
!    INTEGER(i_d), INTENT(IN)            :: mvtype  ! Number of vegetation types
!    INTEGER(i_d), INTENT(IN)            :: mstype ! Number of soil types
!    INTEGER(i_d) :: j, tmp   ! counter for output

    ! Fix in-canopy turbulence scheme globally:
    veg%meth = 1

    CALL ruff_resist(veg, rough, ssoil, canopy)
    CALL init_radiation(met, rad, veg, canopy)
    CALL surface_albedo(ssoil, veg, met, rad, canopy)

    ! Calculate canopy variables:
    CALL define_canopy(ktau,bal,rad,rough,air,met,dels,ssoil,soil,veg, &
                     & bgc,canopy)

    ! owetfac is updated here instead of in canopy_module
    ! because the UM driver will call define_canopy twice
    ssoil%owetfac = ssoil%wetfac

    ! Calculate soil and snow variables:
    CALL soil_snow(dels, ktau, soil, ssoil, veg, canopy, met)
    !	need to adjust fe after soilsnow
    canopy%fev = REAL(canopy%fevc + canopy%fevw, r_1)
    ! Calculate total latent heat flux:
    canopy%fe = canopy%fev + canopy%fes
    ! Calculate net radiation absorbed by soil + veg
    canopy%rnet = canopy%fns + canopy%fnv
    
    ! Calculate radiative/skin temperature:
    rad%trad = ( (1.-rad%transd)*canopy%tv**4 &
         & + rad%transd * ssoil%tss**4 )**0.25
    
!    rad%flws = sboltz*emsoil* ssoil%tss **4   ! YP Nov 2009

! sumcflux now called in CABLE_run or cable_driver (BP apr2010)
!    ! Call carbon routines:
!    CALL sumcflux(ktau, kstart, kend, dels, bgc, canopy,  &
!         & soil, ssoil, sum_flux, veg, met, mvtype, mstype)

  END SUBROUTINE cbm

END MODULE cbm_module
