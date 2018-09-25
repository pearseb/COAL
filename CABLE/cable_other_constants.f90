! cable_other_constants.f90
!
! Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach, 
! Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!
! Fortran-95 recoding by Harvey Davies, Gab Abramowitz, Martin Dix & Bernard Pak
!
! Please send bug reports to Bernard.Pak@csiro.au
!
! This file contains modules:
!  other_constants

MODULE other_constants
  USE define_dimensions, ONLY : r_1, nrb
  IMPLICIT NONE
  PRIVATE r_1, nrb
  ! Gaussian integ. weights
  REAL(r_1), PARAMETER, DIMENSION(nrb) :: gauss_w=(/0.308,0.514,0.178/)
  ! leaf reflectance
  REAL(r_1), PARAMETER, DIMENSION(nrb) :: refl = (/ 0.07, 0.425, 0.00 /) ! YP nov2009
!  REAL(r_1), PARAMETER, DIMENSION(nrb) :: refl = (/ 0.07, 0.425, 0.02 /)
  ! leaf transmittance
  REAL(r_1), PARAMETER, DIMENSION(nrb) :: taul = (/ 0.07, 0.425, 0.00 /) ! YP nov2009
!  REAL(r_1), PARAMETER, DIMENSION(nrb) :: taul = (/ 0.07, 0.425, 0.02 /)
END MODULE other_constants
