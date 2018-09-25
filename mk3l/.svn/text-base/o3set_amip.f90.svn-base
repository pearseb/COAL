      subroutine o3set_amip ( alat, mins, prh, ihem, qo3 )

!     Interpolate the AMIP2 ozone data in time, latitude and pressure to 
!     produce ozone mixing ratio for model radiation code. 

!     This version takes an ihem argument that specifies which half of 
!     the qo3 array to fill (double row physics version).

!     Fixed format f90

!     Modified to enable 100,000-year runs.
!     SJP 2008/02/28
!
!     Fixing syntax errors identified by the g95 Fortran compiler.
!     SJP 2009/04/14

      implicit none
      include 'PHYSPARAMS.f'
      include 'PARAMS.f'
      include 'O3AMIP.f'
      real, intent(in) :: alat     ! Latitude
      integer(kind=8), intent(in) :: mins  ! Time
      integer, intent(in) :: ihem
      real, intent(in),  dimension(ln2,nlp) :: prh  ! Half level pressure
      real, intent(out), dimension(ln2,nl)  :: qo3  ! Ozone mixing ratio
!     Day no of middle of month
      real, parameter, dimension(12) :: monmid =  
     &  (/ 15.5, 45.0, 74.5, 105.0, 135.5, 166.0, 196.5, 227.5, 258.0, 
     &     288.5, 319.0, 349.5 /)
      real, dimension(kg) :: oz    ! Column for this date and lat.
      real, dimension(lg) :: ozcol ! Integrated column
      real, dimension(lon,nlp) :: qo3p
      integer, dimension(lon) :: k1
      integer :: i, m1, m2, j1, j2, k, kk, k1min, ioff
      real :: fac1, fac2, tfac1, tfac2, date, theta

      ioff = (ihem-1) * lon  ! Offset for double row arrays.
      theta = alat*180.0/pi

!     Factors for interpolation in latitude, 
!     note that grid latitudes run S to N.
      if ( theta <= glat(1) ) then
         j1 = 1
         j2 = 1
         fac1 = 1.0
         fac2 = 0.0
      else if ( theta >= glat(jg) ) then
         j1 = jg
         j2 = jg
         fac1 = 1.0
         fac2 = 0.0
      else
!     Input data isn't exactly equally spaced, so search to find the spanning
!     latitudes. Want to find first j such that glat(j) >= theta
         do j2=2,jg
            if ( glat(j2) >= theta ) exit
         end do
         j1 = j2 - 1
         fac1 = ( glat(j2) - theta ) / ( glat(j2) - glat(j1) )
         fac2 = 1.0 - fac1
      end if

!     Time interpolation factors (assume year of exactly 365 days)
      date = real(modulo(mins,525600_8)) / 1440.0
      if ( date <= monmid(1) ) then  ! First half of Jan
         m1 = 12
         m2 = 1
         tfac1 = (monmid(1) - date)/31.0
         tfac2 = 1.0 - tfac1
      else if ( date >= monmid(12) ) then
         m1 = 12
         m2 = 1
         tfac2 = (date - monmid(12))/31.0
         tfac1 = 1 - tfac2
      else
!        Search for bracketing dates, i such that monmid(i) >= date
         do m2=2,12
            if ( monmid(m2) >= date ) exit
         end do
         m1 = m2 - 1
         tfac1 = ( monmid(m2) - date ) / ( monmid(m2) - monmid(m1) )
         tfac2 = 1.0 - tfac1
      end if

!     Now interpolate in latitude and time.
      oz = tfac1 * ( fac1*gdat(j1,:,m1) + fac2*gdat(j2,:,m1 ) ) + 
     &     tfac2 * ( fac1*gdat(j1,:,m2) + fac2*gdat(j2,:,m2 ) ) 

!     Each longitude has same vertical profile as a function of pressure,
!     but model levels depend on surface pressure so each must be 
!     interpolated separately. To ensure the column amount is maintained
!     correctly, first calculate ozone column between every half level and
!     TOA.
      ozcol(1) = 0.0
      do k=1,kg
         ozcol(k+1) = ozcol(k) + oz(k)*dp(k)
      end do

!     To calculate model amounts, find the data half level immediately above
!     the given model half level. Note that model levels start at surface while
!     data levels start at the top.
!
      k1(:) = 1
      qo3p(:,nlp) = 0.0 ! TOA
      do k=nl,1,-1  ! Start at the top model level

!        Find largest k1 such that gpri(k1) <= prh(k)
!        Minimise work by starting from min value of k1 on the previous level
         k1min = minval(k1)
         do kk = k1min, kg
            do i=1,lon
               if ( gpri(kk) <= prh(i+ioff,k) ) then
                  k1(i) = kk
               end if
            end do
         end do

         do i=1,lon
            qo3p(i,k) = ozcol(k1(i)) +
     $                  (prh(i+ioff,k)-gpri(k1(i)))*oz(k1(i))
         end do

      end do

!     Finally get the model ozone mixing ratios by differencing the half level
!     path amounts. The vertical order of the qo3 array is reversed to 
!     match the radiation code.
      do k=1,nl
         do i=1, lon
            qo3(i+ioff,nl+1-k) = ( qo3p(i,k+1) - qo3p(i,k) ) / 
     &                      ( prh(i+ioff,k+1) - prh(i+ioff,k) )
         end do
      end do


      end subroutine o3set_amip
