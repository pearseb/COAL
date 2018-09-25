! cable_math_constants.f90
!
! Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach, 
! Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!
! Fortran-95 recoding by Harvey Davies, Gab Abramowitz, Martin Dix & Bernard Pak
!
! Please send bug reports to Bernard.Pak@csiro.au
!
! This file contains modules:
!  math_constants

MODULE math_constants
  USE define_dimensions, ONLY : i_d, r_1
  IMPLICIT NONE
  PRIVATE i_d, r_1
  REAL(r_1), PARAMETER :: pi_c = 3.1415927
  REAL(r_1), PARAMETER :: pi180 = pi_c / 180.0 ! radians / degree
  REAL(r_1), PARAMETER :: two_pi = 2.0 * pi_c
END MODULE math_constants
