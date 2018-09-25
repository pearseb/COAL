! cable_dimensions.f90
!
! Source file containing dimensions variables for CABLE
!
! Development by Ying-Ping Wang, Eva Kowalczyk, Ray Leuning
! Gab Abramowitz, Martin Dix, Harvey Davies, Mike Raupach
!
! bugs to bernard.pak@csiro.au

MODULE define_dimensions
  INTEGER :: mland    ! # land grid cells
  INTEGER :: mp       ! # patches(>=mland)
! moved max_vegpatches to MODULE io_variables (BP sep2010)
!  INTEGER :: max_vegpatches ! The maximum # of patches in any grid cell
! added mvtype and mstype (BP sep2010)
  INTEGER            :: mvtype  ! total # vegetation types,   from input
  INTEGER            :: mstype  ! total # soil types,         from input
  INTEGER, PARAMETER :: mf = 2  ! # leaves (sunlit, shaded)
  INTEGER, PARAMETER :: nrb = 3 ! # radiation bands
  INTEGER, PARAMETER :: ms = 6  ! # soil layers
  INTEGER, PARAMETER :: msn = 3 ! max # snow layers
  INTEGER, PARAMETER :: ncp = 3 ! # vegetation carbon stores
  INTEGER, PARAMETER :: ncs = 2 ! # soil carbon stores
  ! i_d is default kind for representing integer values.
  INTEGER, PARAMETER :: i_d = KIND(9)
  ! r_1 is default kind for representing REAL values (typically 32 or 64 bits).
  INTEGER, PARAMETER :: r_1  = KIND(1.0)
  ! r_2 is kind for representing REAL values with at least 10-digit precision
  ! (typically 64 bits).
  INTEGER, PARAMETER :: r_2  = SELECTED_REAL_KIND(12, 50)

  !$OMP THREADPRIVATE(mland,mp)

END MODULE define_dimensions
