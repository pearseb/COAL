! cable_albedo.f90
!
! Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach, 
! Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!
! Fortran-95 recoding by Harvey Davies, Gab Abramowitz, Martin Dix & Bernard Pak
!
! Please send bug reports to Bernard.Pak@csiro.au
!
! This file contains modules:
!   albedo_module
! The subroutines included are:
!   surface_albedo
!
MODULE albedo_module
  USE define_dimensions, ONLY: r_1,i_d,mp,nrb
  USE other_constants, ONLY: refl,taul
  USE define_types
  IMPLICIT NONE
  ! This module contains the following subroutines:
  PRIVATE
  PUBLIC surface_albedo
CONTAINS
  !-----------------------------------------------------------------------------
  SUBROUTINE surface_albedo(ssoil, veg, met, rad, canopy)
    TYPE (soil_snow_type),INTENT(INOUT)   :: ssoil
    TYPE (veg_parameter_type),INTENT(IN)  :: veg
    TYPE (met_type),INTENT(INOUT)         :: met
    TYPE (radiation_type),INTENT(INOUT)   :: rad
    TYPE (canopy_type),INTENT(INOUT)      :: canopy
    REAL(r_1), DIMENSION(nrb) :: c1 ! sqrt(1. - taul - refl)
    LOGICAL, DIMENSION(mp)   :: mask ! select points for calculation
    REAL(r_1), DIMENSION(mp,nrb) :: rhocbm ! modified canopy beam reflectance (6.21)
    INTEGER(i_d) :: b ! rad. band 1=visible, 2=near-infrared, 3=long-wave
    REAL(r_1), DIMENSION(nrb) :: rhoch ! canopy reflection black horizontal leaves(6.19)

    rad%cexpkbm = 0.0
    rad%extkbm  = 0.0
    rad%cexpkdm = 0.0
    rad%extkdm  = 0.0
    rhocbm      = 0.0

    ! Initialise effective conopy beam reflectance:
    rad%reffbm = ssoil%albsoilsn
    rad%reffdf = ssoil%albsoilsn
    rad%albedo = ssoil%albsoilsn

    ! Define vegetation mask:
    mask = canopy%vlaiw > 1.0e-2 .AND. met%fsd > 1.0e-2
    
    c1 = SQRT(1. - taul - refl)
    ! Define canopy reflection black horizontal leaves(6.19)
    rhoch = (1.0 - c1) / (1.0 + c1)
    ! Update extinction coefficients and fractional transmittance for 
    ! leaf transmittance and reflection (ie. NOT black leaves):
    DO b = 1, 2 ! 1 = visible, 2 = nir radiaition
       
       rad%extkdm(:,b) = rad%extkd * c1(b)
       ! Define canopy diffuse transmittance (fraction):
       rad%cexpkdm(:,b) = EXP(-rad%extkdm(:,b) * canopy%vlaiw)
       ! Calculate effective diffuse reflectance (fraction):
       WHERE (canopy%vlaiw > 1.0e-2) &
          rad%reffdf(:,b) = rad%rhocdf(:,b) + (ssoil%albsoilsn(:,b) &
                           - rad%rhocdf(:,b)) * rad%cexpkdm(:,b)**2     
       WHERE (mask) ! i.e. vegetation and sunlight are present
          rad%extkbm(:,b) = rad%extkb * c1(b)
          ! Canopy reflection (6.21) beam:
          rhocbm(:,b) = 2.*rad%extkb/(rad%extkb+rad%extkd)*rhoch(b)
          ! Canopy beam transmittance (fraction):
          rad%cexpkbm(:,b) = EXP(-rad%extkbm(:,b)*canopy%vlaiw)
          ! Calculate effective beam reflectance (fraction):
          rad%reffbm(:,b) = rhocbm(:,b) + (ssoil%albsoilsn(:,b) &
               - rhocbm(:,b))*rad%cexpkbm(:,b)*rad%cexpkbm(:,b)
       END WHERE
       ! Define albedo:
       WHERE (canopy%vlaiw > 1.0e-2) &
          rad%albedo(:,b) = (1.0-rad%fbeam)*rad%reffdf(:,b) &
                          + rad%fbeam*rad%reffbm(:,b)
    END DO
!    ! Modify near IR albedo (YP Apr08)
!    WHERE (mask) rad%reffbm(:,2) = rad%reffbm(:,2) * veg%xalbnir(:)
!    WHERE (canopy%vlaiw > 1.0e-2)
!       rad%reffdf(:,2) = rad%reffdf(:,2) * veg%xalbnir(:)
!       rad%albedo(:,2) = rad%albedo(:,2) * veg%xalbnir(:)
!    END WHERE
!    ! Define IR albedo - CURRENTLY NOT USED elsewhere
!    rad%albedo(:,3) = 0.05
!    WRITE(48,'(a7,200f5.2)') 'vlaiw  ', canopy%vlaiw
!    WRITE(48,'(a7,200f5.2)') 'extkdm1', rad%extkdm(:,1)
!    WRITE(48,'(a7,200f5.2)') 'cexpkdm', rad%cexpkdm(:,1)
!    WRITE(48,'(a7,200f5.2)') 'rhocdf1', rad%rhocdf(:,1)
!    WRITE(48,'(a7,200f5.1)') 'xalbnir', veg%xalbnir
!    WRITE(48,'(a7,200f5.2)') 'reffdf1', rad%reffdf(:,1)
!    WRITE(48,'(a7,200f5.2)') 'reffdf2', rad%reffdf(:,2)
!    WRITE(48,'(a7,200f5.2)') 'reffbm1', rad%reffbm(:,1)
!    WRITE(48,'(a7,200f5.2)') 'reffbm2', rad%reffbm(:,2)
!    WRITE(48,'(a7,200f5.2)') 'fbeam  ', rad%fbeam

  END SUBROUTINE surface_albedo

END MODULE albedo_module
