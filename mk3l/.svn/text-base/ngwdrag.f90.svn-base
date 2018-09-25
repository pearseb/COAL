! (1) Fixing syntax errors identified by the g95 Fortran compiler.
! (2) Minor fixes to resolve warnings issued by the g95 Fortran compiler.
! SJP 2009/04/14
!
! $Log: ngwdrag.f90,v $
! Revision 1.6  2001/10/12 02:06:56  rot032
! HBG changes, chiefly for new river-routing scheme and conservation
!
! Revision 1.1  2001/02/18 23:35:53  rot032
! Initial revision
!
! Revision 1.5  1999/05/20 06:23:48  rot032
! HBG changes to V5-2
!
! Revision 1.4  1998/12/10  00:55:28  ldr
! HBG changes to V5-1-21
!
! Revision 1.3  1997/12/17  23:22:44  ldr
! Changes from MRD for parallel execution on NEC.
!
! Revision 1.2  1997/10/03  05:29:55  mrd
! Add a check on full level winds to see if there's a lower critical level.
!
! Revision 1.1  1996/12/23  03:54:05  mrd
! Initial revision
!

      subroutine ngwdrag(lg,tdt,pg,tg,ttg,u,v,sd,gamma,theta,slope,
     &                   prf,prh,dudt,dvdt,dkegw)

!  This routine implements a new gravity wave drag scheme, following the
!  ECMWF scheme described in 
!  Baines, P.G. and T. N. Palmer, 1990: Rationale for a new physically-based
!    parameterization of subgridscale orographic effects. ECMWF Tech. Memo 169
!  Lott, F. and M. Miller: A new sub-grid scale orographic drag 
!    parametrization: Its formulation and testing. QJRMS (submitted)

!  The code is fixed format Fortran 90. Fixed format is required to use the
!  standard f77 include files.

      implicit none
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

      include 'CNSTA.f'
      include 'SI.f'

!  Subroutine arguments

      integer, intent(in) ::  lg             ! Latitude index
      real, intent(in)    ::  tdt            ! Time step
      real, intent(in)    ::  pg(ln2)        ! Surface pressure
      real, intent(in)    ::  ttg(ln2,nl)    ! Atmospheric temperature
      real, intent(in)    ::  tg(ln2)        ! Surface temperature
      real, intent(in)    ::  u(ln2,nl)      ! Zonal wind
      real, intent(in)    ::  v(ln2,nl)      ! Meridional wind
      real, intent(in)    ::  sd(ln2)        ! Orographic standard deviation
      real, intent(in)    ::  gamma(ln2)     ! Orographic anisotropy factor
      real, intent(in)    ::  theta(ln2)     ! Direction of orographic 
                                             ! principal axis
      real, intent(in)    ::  slope(ln2)     ! RMS slope of the orography
      real, intent(in)    ::  prf(ln2,nl)
      real, intent(in)    ::  prh(ln2,nlp)
      real, intent(out)    :: dudt(ln2,nl)   ! Zonal wind tendency
      real, intent(out)    :: dvdt(ln2,nl)   ! Meridional wind tendency
      real, intent(out)    :: dkegw(ln2,nl)  ! Change in KE

!  Working arrays
      real     press(ln2,nl)  ! Full level pressure
      real     press2(ln2,nlp)! Half level pressure
      real     dp(ln2,nl)     ! Thickness of half-layers
      real     geop(ln2,nl)   ! Level geopotential heights
      real     ulow(ln2)      ! Average low level wind
      real     vlow(ln2)
      real     norm(ln2)      ! Low level wind speed
      real     taux(ln2)      ! x-component of GWD unit vector
      real     tauy(ln2)      ! y-component of GWD unit vector
      real     vph(ln2,nlp)   ! Half level wind projected in plane of 
                              ! surface stress
      real     vpf(ln2,nl)    ! Full level wind projected in plane of
                              ! surface stress
      real     deltainv
      real     phi(ln2)       ! Surface wind direction
      real     psi(ln2)       ! Orographic axis wrt wind direction.
      real     alpha(ln2)     ! Direction of GW stress
      real     tau(ln2,nlp)   ! Gravity wave stress (on half levels)
      real     dwind(ln2) 
      real     stab(ln2,nlp)  ! Stability
      real     rho(ln2,nlp)   ! Density
      logical  lo(ln2)
      logical  lo1(ln2,nl)
      real     hcrit(ln2,nl)
      real     ncrit(ln2,nl)
      real     ri(ln2,nlp)     ! Richardson no (half levels)
      integer  kcrith(ln2)
      integer  khlim(ln2)
      real     nu(ln2)         ! = 2 mu N/U
      real     sqr(ln2)
      real     dz2n(ln2)
      integer  icrit(ln2)
      real, parameter :: hmax = 1.e6  ! Height limit (geopotential)
      real     dz2(ln2)
      real     kdrag(ln2)
      real     alfa(ln2)
      real     riw(ln2)
      real     zr(ln2)
      real     tl(ln2)
      real     tfr(ln2)
      real     dels
      integer  kh

      !  GWD scheme parameters
!hbg  real, parameter :: gfac = 0.1     ! G/pi^2 approx 0.1
      real, parameter :: gfac = 1.0     ! G/pi^2 approx 1.0 (G is not gravity)
      real, parameter :: ldrag = 4.     ! Drag multiplier for supercritical stress
      real, parameter :: rahilo = 0.30  ! Ratio of high to low level stress
      !   real, parameter :: sigt = 0.93   ! Top of layer for stress computation
      !  Changed for 9L model
      real, parameter :: sigt = 0.90    ! Top of layer for stress computation
      real, parameter :: sigcr = 0.80   ! Top of layer for low level drag
      real, parameter :: rcrit = 0.25   ! Critical Richardson number for onset
                                        ! of wave turbulence
      real, parameter :: vcrit = 0.0    ! Lower limit for low level wind used 
                                        ! for GWD

      !  Security parameters
      real, parameter :: vsec = 0.10
      real, parameter :: tsec = 1.0e-7
      real, parameter :: ssec = 1.0e-12
      real, parameter :: sdmin = 1.e-6

      integer  kcrit          ! Index of highest level below sigcr
      integer  ktop           ! Index of highest level below sigt

      !  Loop index variables
      integer j, k, mg

      real cons1, cons2, cons5

      cons1 = 1. / rdry
      cons2 = grav**2 / cp
      cons5 = 1.5 * pi

!  Compute basic state variables

!  Define low level wind, project winds in plane of low level wind 
!  and set indicator for critical levels.

!  Define top of low level flow
      do k=1,nl
	 if ( sigh(k) >= sigcr ) kcrit = k
      end do

!  Define top of subcritical low level drag
      do k=1,nl
	 if ( sig(k) >= sigt ) ktop = k
      end do
      if ( ktop == 1 ) then
	 print*, ' Error - only one level in surface region '
	 stop
      end if
!
!  Full level pressure (Pa)
      press  = 100.*prf
      press2 = 100.*prh
      do k=2,nl
	 dp(:,k) = press(:,k-1) - press(:,k)
      end do

!  Calculate average low level wind speed
      ulow = u(:,1)*dsk(1)
      vlow = v(:,1)*dsk(1)
      do k=2,ktop
	 ulow = ulow + u(:,k)*dsk(k)
	 vlow = vlow + v(:,k)*dsk(k)
      end do

!  Thickness of the averaged levels
      deltainv = 1./(1. - sigh(ktop+1))
      ulow = ulow * deltainv
      vlow = vlow * deltainv
      norm = max(sqrt(ulow**2+vlow**2), vsec)
      phi = atan2(vlow, ulow)

!  Calculate GWD unit vectors 
      psi = theta - phi  ! Angle of orography wrt wind
      alpha = atan ( tan(psi) * (b(gamma)-c(gamma)) /
     &             ( b(gamma) + c(gamma)*tan(psi)**2) ) 
      taux = cos(alpha+phi)
      tauy = sin(alpha+phi)
      where ( sd > sdmin )
         kdrag = gfac * slope/sd *
     &         sqrt( (b(gamma)*cos(psi))**2 + (c(gamma)*sin(psi))**2  )
      elsewhere 
!        Use a non-zero value to avoid numerical problems later
	 kdrag = 1.e-10 
      endwhere

!  Project wind into plane of stress
      do k=1,nl
	 vpf(:,k) = taux*u(:,k) + tauy*v(:,k)
      end do
      tau = 0.0
      hcrit = 0.0
      ncrit = 0.0
      lo1 = .false.

      icrit = nlp
      vph(:,1) = taux*ulow + tauy*vlow
!  Half level winds.
!  The loop counts down from the top to look for the lowest critical
!  level
      do k=nl,2,-1
	 vph(:,k) = ( vpf(:,k-1)*dsk(k) + vpf(:,k)*dsk(k-1) ) /
     &                  ( dsk(k-1) + dsk(k) )
	 where ( vph(:,k) .lt. vsec ) 
	    icrit = k
	    vph(:,k) = vsec
	 endwhere
      end do

!  Now check full level winds to see if there's a lower critical level. This 
!  can happen on odd cases. Without this check the drag can accelerate upper
!  level easterlies.
!  This sets the critical level to the half level below.
      do k=nl,1,-1
	 where ( vpf(:,k) .lt. vsec )
	    icrit = k
	 endwhere
      end do
            
!  Brunt-Vaisala frequency and density at half levels.
      do k=2,nl
	 rho(:,k) = 2.*cons1*press2(:,k)/(ttg(:,k-1)+ttg(:,k))
	 stab(:,k) = 2.*cons2 / (ttg(:,k-1)+ttg(:,k)) *
     &        ( 1. - cp*rho(:,k)*( ttg(:,k-1) - ttg(:,k)) / dp(:,k) )
      end do
      stab(:,nlp) = 0.

!  Surface values are an average of lowest layers
      stab(:,1) = 0.
      rho (:,1) = 0.
      do k=2,ktop
	 stab(:,1) = stab(:,1) + stab(:,k)*(sig(k)-sig(k-1))
	 rho (:,1) = rho (:,1) + rho (:,k)*(sig(k)-sig(k-1))
      end do
      stab(:,1) = stab(:,1) / (sig(ktop)-sig(1))
      rho (:,1) = rho (:,1) / (sig(ktop)-sig(1))
      stab = max(stab, ssec)
!      print*, ' N ', (sqrt(stab(1,k)),k=2,nl)


!  Mean flow Richardson number
      do k=2,nl
	 dwind = max(abs(vpf(:,k)-vpf(:,k-1)), vsec)
	 ri(:,k) = stab(:,k)*(dp(:,k)/(grav*rho(:,k)*dwind))**2
	 ri(:,k) = max(ri(:,k), rcrit)
      end do
            

!  The GWD routine requires the heights of the model levels. These
!  are the height above the surface, not the usual height above msl.
      geop = 0.
      do k=1,nl
         do j=1,nl
            geop(:,k) = geop(:,k)+am(k,j)*ttg(:,j)
         end do
      end do

!  Calculate critical height for Froude layer
!  Use hcrit to hold the wind integral temporarily.
      hcrit(:,1) = 0.
      ncrit(:,1) = 0.
      do k=2,nl
	 hcrit(:,k) = hcrit(:,k-1) + vph(:,k)*dp(:,k)
	 ncrit(:,k) = ncrit(:,k-1) + sqrt(stab(:,k))*dp(:,k)
      end do

      do k=2,nl
!        Convert hcrit to geopotential to match geop
	 hcrit(:,k) = grav*cons5*hcrit(:,k)/ncrit(:,k) ! cons5 = 3pi/2
	 lo1(:,k) = geop(:,k) > hcrit(:,k)
!         print*, ' ZCRIT ', k, hcrit(1,k), geop(1,k)
!        Note the use of xor here is not standard. It works on Cray and SGI.
!	 where ( lo1(:,k) .xor. lo1(:,k-1) ) 
!           The index of the half-level between levels k and k-1 is k
        where ( ( lo1(:,k) .or. lo1(:,k-1) ) .and.
     &              .not. ( lo1(:,k) .and. lo1(:,k-1) ) )  ! Standard version
	    kcrith = k
	 endwhere
	 where ( geop(:,k) <  hmax )
	    khlim =  k
	 endwhere
      end do
!      print*, ' kcrith ', kcrith

!  Lowest level stress
      nu = 2.*sd*sqrt(stab(:,1))/vph(:,1)
      tau(:,1) = rho(:,1)*kdrag*sqrt(stab(:,1))*sd**2*vph(:,1)
!  Hydraulic case
      where ( nu > 1 ) tau(:,1) = tau(:,1) / nu**2
!  If the drag or wind are below limits, or if there is a
!  critical level in the surface layer set the drag to zero.
      lo = (tau(:,1) < tsec) .or. (icrit <= kcrit) .or.
     &     (vph(:,1) < vcrit)
      where ( lo ) 
	 icrit = 1
	 tau(:,1) = 0.
      end where
      tau(:,kcrit) = rahilo*tau(:,1)
!  Limit the supercritical drag to the region below hmax
      kcrith = min(kcrith, khlim)
!  Limit the supercritical drag to the low level drag region defined by 
!  kcrit, except if there is a critical level above this.
!  Is this really the intention of Peltier and Clark?
      kcrith = min(kcrith, max(icrit, kcrit))

      tl = ldrag*gamma(:)*(nu-1)*tau(:,1)
      where (nu > 1.) 
!        Supercritical low level stress
	 tfr = tl
      elsewhere
	 tfr = 0.0
      endwhere

!  Compute stress profile in the upper region (above kcrit)
      do k=kcrit+1,nl

!     Wave displacement at next level
         norm = rho(:,k)*kdrag*sqrt(stab(:,k))*vph(:,k)
         dz2(:) = 2.0*tau(:,k-1)/norm
!     Wave Richardson number, new wave displacement and stress
         sqr = sqrt(ri(:,k))
	 alfa = sqrt(stab(:,k)*dz2)/vph(:,k)
	 riw = ri(:,k)*(1.-alfa)/(1.+alfa*sqr)**2
	 zr = 2.+1./sqr
	 dz2n = (vph(:,k)*(2.*sqrt(zr)-zr))**2/stab(:,k)
	 where ( riw < rcrit ) 
	    tau(:,k) = norm*dz2n*0.5
	 elsewhere
	    tau(:,k) = tau(:,k-1)
	 endwhere
!        Wave is completely absorbed at a critical level.
	 where ( (tau(:,k-1) < tsec) .or. (k >= icrit) ) 
	    tau(:,k) = 0.
	 endwhere
	 tau(:,k) = min(tau(:,k), tau(:,k-1))

      end do ! Loop over k

!  Low level stress profile (subcritical component)
      do k=1,kcrit
	 tau(:,k) = ( tau(:,kcrit)*(sigh(1)-sigh(k)) +
     &                tau(:,1)*(sigh(k)-sigh(kcrit)) ) /
     &                 (sigh(1)-sigh(kcrit))
      end do
!      print*,  ' Std stress ', (tau(1,k),k=1,nl)

!  Addition of low level stress profile (supercritical)
!  kcrith defines the top of the supercritical layer.
      do k=1,nl
	 do mg=1,ln2
	    kh = kcrith(mg)
	    if ( k.lt.kh ) then
	       dels = sigh(k)-sigh(kh)
	    else
	       dels = 0.
	    end if
!           Indirect addressing here can't be put in array form?
	    tau(mg,k) = tau(mg,k) + tfr(mg)*dels/(sigh(1)-sigh(kh))
	 end do
      end do
!      print*, ' Std + supercritical ', (tau(1,k),k=1,nl)

!  Top layer stress set to zero to ensure all the wave drag is absorbed 
!  somewhere.
      tau (:,nlp) = 0.
      stab(:,nlp) = 0.
      ri  (:,nlp) = 0.
      vph (:,nlp) = 0.

!  Tendencies
      do k=1,nl
	 dudt(:,k) = -grav*taux*(tau(:,k+1)-tau(:,k)) /
     &                          (press2(:,k+1)-press2(:,k))
	 dvdt(:,k) = -grav*tauy*(tau(:,k+1)-tau(:,k)) /
     &      	                (press2(:,k+1)-press2(:,k))
         dkegw(:,k)=u(:,k)*dudt(:,k)+v(:,k)*dvdt(:,k)
      end do

      return

      contains

!     Internal functions for B(gamma) and C(gamma)
!     Use simple polynomial form for these and let the optimiser do its best.
      function b(gamma)
	 real, intent(in), dimension(ln2) :: gamma(ln2)
	 real, dimension(ln2) :: b
	 b = 1. - 0.18 * gamma - 0.04 * gamma*gamma
      end function b
      function c(gamma)
	 real, intent(in), dimension(ln2) :: gamma(ln2)
	 real, dimension(ln2) :: c
	 c = 0.48 * gamma + 0.3 * gamma*gamma
      end function c

      end subroutine ngwdrag
