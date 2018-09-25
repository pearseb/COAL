! soil_snow.f90
!
! Science development by Eva Kowalczyk, CSIRO Marine and Atmospheric Research
! 
! Fortran-95 recoding by Harvey Davies, Gab Abramowitz, Martin Dix & Bernard Pak
!
! Please send bug reports to Bernard.Pak@csiro.au
!
! This file contains the soil_snow_module only, with the following subroutines:
!   trimb,
!   smoisturev,
!   snowdensity,
!   snow_melting,
!   snow_accum,
!   surfbv,
!   snow_albedo,
!   stempv,
!   stempvsn,
!   snowcheck,
!   snowl_adjust,
!   soilfreeze,
!   remove_trans, and
!   soil_snow
!
MODULE soil_snow_module
  USE define_dimensions, ONLY: r_1,r_2,i_d,mp,ms
  USE physical_constants, ONLY: tfrz 
  USE define_types       
  USE io_variables, ONLY: landpt
  IMPLICIT NONE
  PRIVATE
  REAL(r_1), PARAMETER :: cgsnow = 2090.0 ! specific heat capacity for snow
  REAL(r_1), PARAMETER :: csice = 2.100e3 ! specific heat capacity for ice
  REAL(r_1), PARAMETER :: cswat = 4.218e3 ! specific heat capacity for water
  REAL(r_1), PARAMETER :: hl = 2.5104e6   ! latent heat of evaporation
  REAL(r_1), PARAMETER :: hlf = 0.335e6   ! latent heat of fusion
  REAL(r_1), PARAMETER :: cp = 1004.64    ! specific heat capacity for air
  REAL(r_1), PARAMETER :: rhowat = 1000.0 ! density of water
  REAL(r_1), PARAMETER :: snmin = 0.11    ! for 3-layer;

  ! This module contains the following subroutines:
  PUBLIC soil_snow ! must be available outside this module
  PRIVATE trimb, smoisturev, snow_accum, stempv
  PRIVATE snowdensity, snow_melting, snowcheck, snowl_adjust 
  PRIVATE soilfreeze, remove_trans, snow_albedo
  INTEGER(i_d), PARAMETER, PRIVATE  :: idjd = 1

CONTAINS

  !----------------------------------------------------------------------
  ! SUBROUTINE trimb
  !
  !      this routine solves the system
  !	   a(k)*u(k-1)+b(k)*u(k)+c(k)*u(k+1)=rhs(k)    for k=2,kmax-1
  !	   with  b(k)*u(k)+c(k)*u(k+1)=rhs(k)	       for k=1
  !	   and	 a(k)*u(k-1)+b(k)*u(k)=rhs(k)	       for k=kmax
  !
  !	 the Thomas algorithm is used for solving sets of linear equation
  !	 rhs initially contains rhs; leaves with answer (jlm)
  !	 n.b. this one does not assume b = 1-a-c
  !
!  SUBROUTINE trimb (a, b, c, rhs, kmax)
!    REAL(r_2), DIMENSION(:,:), INTENT(IN)     :: a ! coef "A" in finite diff eq
!    REAL(r_2), DIMENSION(:,:), INTENT(IN)     :: b ! coef "B" in finite diff eq
!    REAL(r_2), DIMENSION(:,:), INTENT(IN)     :: c ! coef "C" in finite diff eq
!    REAL(r_2), DIMENSION(:,:), INTENT(INOUT)  :: rhs ! right hand side of eq
!    INTEGER(i_d), INTENT(IN)                  :: kmax ! no. of discrete layers
!    INTEGER(i_d)                              :: k ! do lloop counter
!    REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) :: e 
!    REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) :: g 
!    REAL(r_2), DIMENSION(SIZE(a,1),SIZE(a,2)) :: temp 
!
!    WRITE(46,'(i2,1x,a5,200e12.3)') 1,' b = ',b(:,1)
!    WRITE(46,'(i2,1x,a5,200e12.3)') 1,' c = ',c(:,1)
!    e(:,1) = c(:,1) / b(:,1)
!    WRITE(46,'(i2,1x,a5,200e12.3)') 1,' e = ',e(:,1)
!    DO k = 2, kmax - 1
!        temp(:,k) = 1.0 / (b(:,k) - a(:,k) * e(:,k-1) )
!        e(:,k) = c(:,k) * temp(:,k)
!    WRITE(46,'(i2,1x,a5,200e12.3)') k,' e = ',e(:,k)
!    END DO
!   
!    g(:,1) = rhs(:,1) / b(:,1)
!   PRINT *, 'In the middle of trimb'    
!    DO k = 2, kmax - 1
!      g(:,k) = (rhs(:,k) - a(:,k) * g(:,k-1) ) * temp(:,k)
!    END DO
!    WRITE(46,'(i2,1x,a5,200e12.3)') kmax,' b = ',b(:,kmax)
!    WRITE(46,'(i2,1x,a5,200e12.3)') kmax,' a = ',a(:,kmax)
!    WRITE(46,'(i2,1x,a5,200e12.3)') kmax,' a*e ',(a(:,kmax) * e(:,kmax-1))
!   PRINT *, '(b(:,kmax) - a(:,kmax) * e(:,kmax-1)) = '
!   PRINT *, (b(:,kmax) - a(:,kmax) * e(:,kmax-1))
!
!    ! do back substitution to give answer now
!    rhs(:,kmax) = (rhs(:,kmax) - &
!       & a(:,kmax) * g(:,kmax-1)) / (b(:,kmax) - a(:,kmax) * e(:,kmax-1))
!    DO k = kmax - 1, 1, - 1
!      rhs(:,k) = g(:,k) - e(:,k) * rhs(:,k + 1)
!    END DO
!    
!  END SUBROUTINE trimb
  SUBROUTINE trimb (a, b, c, rhs, kmax)
    INTEGER(i_d), INTENT(IN)  :: kmax ! no. of discrete layers
    REAL(r_2), INTENT(IN)     :: a(mp,kmax) ! coef "A" in finite diff eq
    REAL(r_2), INTENT(IN)     :: b(mp,kmax) ! coef "B" in finite diff eq
    REAL(r_2), INTENT(IN)     :: c(mp,kmax) ! coef "C" in finite diff eq
    REAL(r_2), INTENT(INOUT)  :: rhs(mp,kmax) ! right hand side of eq
    INTEGER(i_d)              :: k ! do lloop counter
    REAL(r_2)                 :: e(mp,kmax) 
    REAL(r_2)                 :: g(mp,kmax)
    REAL(r_2)                 :: temp(mp,kmax)

!    WRITE(46,'(i2,1x,a5,e12.3)') 1,' b = ',b(mp-1,1)
!    WRITE(46,'(i2,1x,a5,e12.3)') 1,' c = ',c(mp-1,1)
    e(:,1) = c(:,1) / b(:,1)
!    WRITE(46,'(i2,1x,a5,e12.3)') 1,' e = ',e(mp-1,1)
    DO k = 2, kmax - 1
        temp(:,k) = 1.0 / (b(:,k) - a(:,k) * e(:,k-1) )
        e(:,k) = c(:,k) * temp(:,k) 
!    WRITE(46,'(i2,1x,a5,e12.3)') k,' e = ',e(mp-1,k)
    END DO
    
    g(:,1) = rhs(:,1) / b(:,1)
!   PRINT *, 'In the middle of trimb'
    DO k = 2, kmax - 1
      g(:,k) = (rhs(:,k) - a(:,k) * g(:,k-1) ) * temp(:,k)
    END DO
!    WRITE(46,'(i2,1x,a5,e12.3)') kmax,' b = ',b(mp-1,kmax)
!    WRITE(46,'(i2,1x,a5,e12.3)') kmax,' a = ',a(mp-1,kmax)
!    WRITE(46,'(i2,1x,a5,e12.3)') kmax,' a*e ',(a(mp-1,kmax) * e(mp-1,kmax-1))
!   PRINT *, '(b(:,kmax) - a(:,kmax) * e(:,kmax-1)) = '
!   PRINT *, (b(mp-1,kmax) - a(mp-1,kmax) * e(mp-1,kmax-1))
   
    ! do back substitution to give answer now
    rhs(:,kmax) = (rhs(:,kmax) - &
       & a(:,kmax) * g(:,kmax-1)) / (b(:,kmax) - a(:,kmax) * e(:,kmax-1))
    DO k = kmax - 1, 1, - 1
      rhs(:,k) = g(:,k) - e(:,k) * rhs(:,k + 1)
    END DO
    
  END SUBROUTINE trimb

  !-------------------------------------------------------------------------
  ! SUBROUTINE smoisturev (fwtop,dels,ktau,ssoil,soil)
  !      Solves implicit soil moisture equation
  !      Science development by Eva Kowalczyk and John McGregor, CMAR
  !
  SUBROUTINE smoisturev (dels,ktau,ssoil,soil)
    REAL(r_1), INTENT(IN)                     :: dels    ! time step size (s)
    INTEGER(i_d), INTENT(IN)                  :: ktau  ! integration step number
    TYPE (soil_snow_type), INTENT(INOUT)      :: ssoil ! soil and snow variables
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil  ! soil parameters
    INTEGER(i_d), PARAMETER                   :: ntest = 0 ! 2 for funny pre-set
                                                           ! for idjd
    ! nmeth selects the solution method
    INTEGER(i_d), PARAMETER                   :: nmeth = -1 ! preferred method
    !                                  Values as follows:
    !                                   -1 for simple implicit D
    !                                    1 for fully implicit solution
    !                                    2 for simpler implicit
    !                                    3 for simple implicit D, explicit K 
    !                                    4 for simple implicit D, implicit K
    !                                    0 for simple implicit D, new jlm TVD K
! change dimension of at,bt,ct from 3*ms to ms (BP Jun2010)
    REAL(r_2), DIMENSION(mp,ms)   :: at     ! coef "A" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms)   :: bt     ! coef "B" in finite diff eq
    REAL(r_2), DIMENSION(mp,ms)   :: ct     ! coef "C" in finite diff eq
    REAL(r_2), DIMENSION(mp)      :: fact
    REAL(r_2), DIMENSION(mp)      :: fact2
    REAL(r_2), DIMENSION(mp)      :: fluxhi
    REAL(r_2), DIMENSION(mp)      :: fluxlo
    REAL(r_2), DIMENSION(mp)      :: hydss  ! hydraulic
                                                ! conductivity adjusted for ice
    INTEGER(i_d)                                 :: k
    REAL(r_2), DIMENSION(mp)      :: phi
    REAL(r_2), DIMENSION(mp)      :: pwb
    REAL(r_2), DIMENSION(mp)      :: rat
    REAL(r_2), DIMENSION(mp)      :: speed_k
    REAL(r_2), DIMENSION(mp)      :: ssatcurr_k
    REAL(r_1), DIMENSION(mp)      :: wblfmn
    REAL(r_1), DIMENSION(mp)      :: wblfmx
    REAL(r_2), DIMENSION(mp,ms+1) :: wbh
    REAL(r_2), DIMENSION(mp,ms+1) :: z1
    REAL(r_2), DIMENSION(mp,ms+1) :: z2
    REAL(r_2), DIMENSION(mp,ms+1) :: z3
    REAL(r_1), DIMENSION(mp,ms+1) :: z1mult
    REAL(r_2), DIMENSION(mp,0:ms) :: fluxh
    REAL(r_2), DIMENSION(mp,0:ms) :: delt
    REAL(r_2), DIMENSION(mp,0:ms) :: dtt
    REAL(r_2), DIMENSION(mp)      :: pwb_wbh
    REAL(r_2), DIMENSION(mp,ms)   :: ssatcurr
    REAL(r_1), DIMENSION(mp)      :: totwba ! diagnostic
    REAL(r_1), DIMENSION(mp)      :: totwbb
    REAL(r_1), DIMENSION(mp)      :: totwbc
    REAL(r_1), DIMENSION(mp)      :: totwblb
    REAL(r_1), DIMENSION(mp)      :: totwblc
    REAL(r_1), DIMENSION(mp)      :: wbficemx
    REAL(r_2), DIMENSION(mp)      :: wbh_k
    REAL(r_2), DIMENSION(mp)      :: wbl_k
    REAL(r_2), DIMENSION(mp)      :: wbl_kp
    REAL(r_2), DIMENSION(mp)      :: wh
    REAL(r_2), DIMENSION(mp)      :: z3_k
    REAL(r_1), DIMENSION(mp)      :: zsetot
    LOGICAL :: is_open     ! Is file open?
    INTEGER(i_d) :: u           ! I/O unit
    !
!    IF (ktau == 1) THEN
       ! For some sorts of nested models, this might have been allocated on
       ! another nest
       !IF ( ALLOCATED(ssoil%pwb_min) ) DEALLOCATE(ssoil%pwb_min)
       !ALLOCATE(ssoil%pwb_min(mp))
       ! Set working variable:
       ssoil%pwb_min = (soil%swilt / soil%ssat ) **soil%ibp2
!    END IF
    ! Diagnostic block below:
!!$    IF (ntest > 0) THEN
!!$     ! Print variables for all active veg/soil patches in idjd grid point:
!!$     PRINT * , 'entering smoisturev fwtop,i2bp3,swilt,sfc,ssat: ', &
!!$     ssoil%fwtop(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$     soil%i2bp3(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$     soil%swilt(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$     soil%sfc(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$     soil%ssat(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1))
!!$     u = 97
!!$     inquire (u, opened=is_open)
!!$     IF (.NOT. is_open) THEN
!!$       open (u, file='f97.txt', status='replace')
!!$       write (u, *) 'ktau', ' fwtop', ' i2bp3', ' swilt', ' sfc', ' ssat'
!!$     END IF
!!$     ! Write variables for all active veg/soil patches in idjd grid point:
!!$     write (u, *) ktau,&
!!$     ssoil%fwtop(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$     soil%i2bp3(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$     soil%swilt(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$     soil%sfc(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$     soil%ssat(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1))
!!$     IF (ntest == 2) THEN ! just to test conservation
!!$       IF (ktau == 1) ssoil%wb(:,ms) = soil%swilt
!!$       ssoil%fwtop = 0.0
!!$     END IF
!!$     WRITE (6, " ('wb   ', 6f8.3) ") (ssoil%wb(landpt(idjd)%cstart: &
!!$              & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k) , k = 1, ms)
!!$     WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(landpt(idjd)%cstart: &
!!$              & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k) , k = 1, ms)
!!$     totwba = 0.0
!!$     DO k = 1, ms
!!$       totwba = totwba + soil%zse(k)*REAL(ssoil%wb(:,k),r_1)
!!$     END DO
!!$    END IF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! preset to allow for non-land & snow points in trimb
!    PRINT *, 'mp in smoisturev = ', mp
    at = 0.0
    bt = 1.0
    ct = 0.0
    z1mult(:,1) = 0.0       ! corresponds to 2b+3
    z1mult(:,ms+1) = 0.0    ! corresponds to 2b+3
    z1(:,1) = 0.0           ! i.e. K(.5),    value at surface
    z1(:,ms+1) = 0.0        ! i.e. K(ms+.5), value at bottom
    ! nmeth: equation solution technique
    IF (nmeth <= 0) THEN
       ! jlm split TVD version
       ! all land points
       delt(:,0) = 0.0
       fluxh(:,0) = 0.0
       fluxh(:,ms) = 0.0
       DO k = 1, ms-1
          ! Calculate amount of liquid soil water:
          wbl_k = MAX(  0.01_r_2, ssoil%wb(:,k)   - ssoil%wbice(:,k)   )
          wbl_kp = MAX( 0.01_r_2, ssoil%wb(:,k+1) - ssoil%wbice(:,k+1) )
          ! Calculate difference in liq soil water b/w consecutive layers:
          delt(:,k) = wbl_kp - wbl_k
          ! especially to allow for isolated frozen layers, use min speed
          wh = MIN(wbl_k, wbl_kp)
          ! with 50% wbice, reduce hyds by 1.e-5
          ! Calculate hyd conductivity adjusted for ice:
          hydss = soil%hyds * (1.0 - MIN(2.0_r_2 * ssoil%wbice(:,k) &
              & / MAX(0.01_r_2,  ssoil%wb(:,k) ), 0.99999_r_2) )
          speed_k = hydss * (wh / soil%ssat ) ** (soil%i2bp3 - 1)
          ! update wb by TVD method
          rat = delt(:,k - 1) / (delt(:,k)+SIGN(REAL(1.0e-20,r_2), delt(:,k)))
          phi = MAX(0.0_r_2, MIN(1.0_r_2, 2.0_r_2 * rat), &
               MIN(2.0_r_2, rat) ) ! 0 for -ve rat
          fluxhi = wh
          fluxlo = wbl_k
!!$          ! Block below for testing purposes only:
!!$          IF (ntest > 0) THEN
!!$             ! Print variables for all active veg/soil patches in idjd grid point:
!!$             PRINT * , 'in TVD for k= ', k
!!$             PRINT * , 'wbl,wh,hydss ', &
!!$                  wbl_k(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$                  wh(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$                  hydss(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1))
!!$             PRINT * , 'speeda,speedb,fluxhi,fluxlo,delt,rat,phi ', &
!!$                  speed_k(landpt(idjd)%cstart: &
!!$                  (landpt(idjd)%cstart+landpt(idjd)%nap-1)),0.5*soil%zse(k)/dels, &
!!$                  fluxhi(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$                  fluxlo(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$                  delt(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1),k),&
!!$                  rat(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
!!$                  phi(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1))
!!$          END IF
          ! scale speed to grid lengths per dels & limit speed for stability
          ! 1. OK too for stability
          speed_k = MIN(speed_k, REAL(0.5 * soil%zse(k) / dels , r_2))
          fluxh(:,k) = speed_k * (fluxlo + phi * (fluxhi - fluxlo) )
      END DO
!    PRINT *, 'Before calculating drainage'
      ! calculate drainage (this code replaces the code in the surfb)
      k = ms 
      WHERE( ssoil%wb(:,ms) > soil%sfc(:))
         wbl_k = MAX(0.01_r_2, ssoil%wb(:,ms) - ssoil%wbice(:,ms) )
         wbl_kp = MAX(0.01_r_2, soil%ssat(:) - ssoil%wbice(:,ms) )
         wh = 0.9*wbl_k + 0.1*wbl_kp
         ! Calculate hyd conductivity adjusted for ice:
         hydss = soil%hyds * ( 1.0 - MIN( 2.0_r_2 * ssoil%wbice(:,ms) &
              / MAX(0.01_r_2, ssoil%wb(:,ms)), 0.99999_r_2 ) )
         speed_k = hydss * (wh / soil%ssat ) ** (soil%i2bp3 - 1)
         fluxlo = wbl_k
         ! scale speed to grid lengths per dels & limit speed for stability
         speed_k = MIN(speed_k, REAL(0.5 * soil%zse(ms) / dels , r_2))
         fluxh(:,ms) = MAX(0.0_r_2,speed_k * fluxlo )
      END WHERE
!     PRINT *, 'Before update wb'
      ! update wb by TVD method
      DO k = ms, 1, -1
        IF (nmeth == -1) THEN ! each new wb constrained by ssat
          fluxh(:,k-1) = MIN(fluxh (:,k-1), (soil%ssat &
               - ssoil%wb(:,k) ) * soil%zse(k) / dels + fluxh(:,k) )
        END IF
        ! fluxh (:,ms) is drainage
        ssoil%wb(:,k) = ssoil%wb(:,k) + dels * (fluxh(:,k-1) - fluxh(:,k)) &
             & / soil%zse(k)
        ! re-calculate wblf
        ssatcurr_k = soil%ssat - ssoil%wbice(:,k)
        dtt(:,k) = dels / (soil%zse(k) * ssatcurr_k)
        ! this defn of wblf has different meaning from previous one in surfbv
        ! N.B. are imposing wbice<wb, so wblf <1
        ssoil%wblf(:,k) = (ssoil%wb(:,k) - ssoil%wbice(:,k) ) / ssatcurr_k
     END DO
!    PRINT *, 'Before rnof2'
     ssoil%rnof2 = dels * REAL(fluxh(:,ms),r_1) * 1000.0
     ! wbh_k represents wblf(k-.5)
     DO k = 2, ms
        ssatcurr_k = REAL(soil%ssat,r_2) - ssoil%wbice(:,k)
        wbh_k = ( soil%zse(k) * ssoil%wblf(:,k-1) + soil%zse(k-1) &
             & * ssoil%wblf(:,k) ) / ( soil%zse(k) + soil%zse(k-1) )
        ! i.e. wbh**(bch+1)
        fact = wbh_k** (soil%ibp2 - 1)
        ! with 50% wbice, reduce hbsh by 1.e-5
        pwb_wbh = (soil%hsbh * (1.0 - MIN(2.0 * MIN(0.99_r_2, MAX( &
             & ssoil%wbice(:,k-1) / MAX(0.01_r_2, ssoil%wb(:,k-1)), &
             & ssoil%wbice(:,k)   / MAX(0.01_r_2, ssoil%wb(:,k)  ) )) &
             & , 0.99999_r_2) )) &
             & * MAX(ssoil%pwb_min, wbh_k * fact)
        ! moisture diffusivity (D) is  wbh*pwb; hsbh includes b
        ! i.e. D(k-.5)/soil%zshh(k)
        z3_k = pwb_wbh / soil%zshh (k)
        ! PRINT * , 'z3_k ', z3_k
        ! where dtt=dels/(soil%zse(k)*ssatcurr_k)
        at (:,k) = - dtt(:,k) * z3_k
        ct (:,k-1) = - dtt(:,k-1) * z3_k
      END DO
      bt = 1.0 - at - ct
      ! Block below for testing purposes only:
!!$      IF (ntest > 0) THEN
!!$        ! Print variables for all active veg/soil patches in idjd grid point:
!!$        !        PRINT * , 'midway through nmeth<=0'
!!$        !        PRINT * , 'fluxh ', (fluxh(idjd,k) , k = 1, ms)
!!$        WRITE (6, " ('wb   ', 6f8.3) ") (ssoil%wb(landpt(idjd)%cstart: &
!!$                  & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1,ms)
!!$        WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(landpt(idjd)%cstart: &
!!$                  & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1,ms)
!!$        totwbb = 0.
!!$        totwblb = 0.
!!$        DO k = 1, ms
!!$           totwbb = totwbb + soil%zse(k) * REAL(ssoil%wb(:,k),r_1) ! diagnostic
!!$           totwblb=totwblb + soil%zse(k) * REAL(ssoil%wblf(:,k),r_1) ! diag.
!!$        END DO
!!$        ! PRINT *, 'nmeth, b+2, 2b+3: ',nmeth, soil%ibp2(idjd), soil%i2bp3(idjd)
!!$        WRITE (6, " ('wb   ', 6f8.3) ") (ssoil%wb(landpt(idjd)%cstart: &
!!$             & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k=1,ms)
!!$        WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(landpt(idjd)%cstart: &
!!$             & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k=1,ms)
!!$        WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(landpt(idjd)%cstart: &
!!$             & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k=1,ms)
!!$        !        PRINT * , 'zse ', soil%zse
!!$        !        PRINT * , 'zshh ', soil%zshh
!!$        !        PRINT * , 'dtt ', (dtt(idjd,k) , k = 1, ms)
!!$        !        PRINT * , 'at ', (at(idjd,k) , k = 1, ms)
!!$        !        PRINT * , 'bt ', (bt(idjd,k) , k = 1, ms)
!!$        !        PRINT * , 'ct ', (ct(idjd,k) , k = 1, ms)
!!$     END IF
     ssoil%wblf(:,1) = ssoil%wblf(:,1) + dtt(:,1) &
          & * ssoil%fwtop / rhowat
  END IF

  IF (nmeth > 0) THEN
     wbficemx = 0.0
     DO k = 1, ms
        ssatcurr(:,k) = REAL(soil%ssat,r_2) - ssoil%wbice(:,k)
        ! this defn of wblf has different meaning from previous one in surfbv
        ! N.B. are imposing wbice<wb, so wblf <1
        ssoil%wblf(:,k) = (ssoil%wb(:,k) - ssoil%wbice(:,k) ) / ssatcurr(:,k)
        ssoil%wbfice(:,k) = REAL(ssoil%wbice(:,k),r_1) / soil%ssat
        wbficemx = MAX(wbficemx, ssoil%wbfice(:,k) )
        dtt(:,k) = dels / (soil%zse(k) * ssatcurr(:,k) )
     END DO

     IF (nmeth == 1) THEN ! full implicit method
        DO k = 2, ms
           ! wbh(k)=MIN(1.,ww(k)*wblf(:,k-1)+(1.-ww(k))*wblf(:,k))
           ! jlm: this is same as:
           wbh(:,k) = (soil%zse(k) * ssoil%wblf(:,k-1) + soil%zse(k-1) &
                * ssoil%wblf(:,k) ) / (soil%zse(k) + soil%zse(k-1) )
           fact = wbh(:,k) ** (soil%ibp2 - 1) ! i.e. wbh**(bch+1)
           fact2 = fact * fact
           pwb = soil%hsbh * fact
           ! moisture diffusivity (D) is  wbh*pwb
           ! other term (K) is wbh*soil%hyds*fact2
           z1(:,k) = wbh(:,k) * ( (soil%i2bp3 - 1) * soil%hyds * fact2 &
                - soil%ibp2 * pwb * (ssoil%wblf(:,k) - ssoil%wblf(:,k-1) ) &
                / soil%zshh (k) )
           z2(:,k) = - soil%i2bp3 * soil%hyds * fact2 + soil%ibp2 * pwb &
                * (ssoil%wblf(:,k) - ssoil%wblf(:,k-1) ) / soil%zshh (k)
           z3(:,k) = pwb * wbh(:,k) / soil%zshh (k)
           at(:,k) = dtt(:,k) * (z2(:,k) * 0.5 * soil%zse(k) / soil%zshh (k) &
                - z3(:,k) )
        END DO
        DO k = 1, ms - 1
           ct(:,k) = dtt(:,k) * ( - z2(:,k+1) * 0.5 * soil%zse(k) &
                / soil%zshh (k+1) - z3(:,k+1) )
           bt(:,k) = 1.0 + dtt(:,k) * ( - z2(:,k+1) * 0.5 * soil%zse(k+1) &
                / soil%zshh (k+1) + z2(:,k) * 0.5 * soil%zse(MAX(k-1,1)) &
                / soil%zshh (k) + z3(:,k+1) + z3(:,k) )
        END DO
        bt(:,ms) = 1.0 + dtt(:,ms) * (z2(:,ms) * 0.5 &
             * soil%zse(ms) / soil%zshh (ms) + z3(:,ms) )
        DO k = 1, ms
           ssoil%wblf(:,k) = ssoil%wblf(:,k) &
                + dtt(:,k) * (z1(:,k+1) - z1(:,k) )
        END DO
     END IF ! (nmeth == 1)
     
     IF (nmeth >= 2) THEN ! part implicit method
        DO k = 2, ms
           z1mult(:,k) = soil%i2bp3 ! corresponds to 2b+3
        END DO
        DO k = 2, ms ! wbh(k) represents wblf(k-.5)
           wbh(:,k) = (soil%zse(k) * ssoil%wblf(:,k-1) + soil%zse(k-1) &
                * ssoil%wblf(:,k) ) / (soil%zse(k) + soil%zse(k-1) )
           fact = wbh(:,k) ** (soil%ibp2 - 1) ! i.e. wbh**(bch+1)
           IF (nmeth == 2) pwb_wbh = soil%hsbh * wbh(:,k) * fact
           IF (nmeth >= 3) & 
                pwb_wbh = soil%hsbh * MAX(ssoil%pwb_min, wbh(:,k) * fact)
           fact2 = fact * fact
           ! moisture diffusivity (D) is  wbh*pwb
           ! other term (K) is wbh*soil%hyds*fact2
           z1(:,k) = soil%hyds * fact2 !  i.e. K(k-.5)/wbh(:,k)
           z3(:,k) = pwb_wbh / soil%zshh(k) !  i.e. D(k-.5)/soil%zshh(k)
           at(:,k) = - dtt(:,k) * z3(:,k)
           ct(:,k-1) = - dtt(:,k-1) * z3(:,k)
        END DO
        bt = 1.0 - at - ct
        IF (nmeth == 4) THEN ! for simple implicit D, implicit K
           bt(:,1) = bt(:,1) + dtt(:,1) * z1mult(:,1+1) &
                * z1(:,1+1) * soil%zse(1+1) / (soil%zse(1) + soil%zse(1+1) )
           DO k = 2, ms
              at(:,k)   = at(:,k) - dtt(:,k) * z1mult(:,k) * z1(:,k) &
                   * soil%zse(k) / (soil%zse(k) + soil%zse(k-1) )
              ct(:,k-1) = ct(:,k-1) + dtt(:,k-1) * z1mult(:,k) * z1(:,k) &
                   * soil%zse(k-1) / (soil%zse(k) + soil%zse(k-1) )
              bt(:,k)   = bt(:,k) - dtt(:,k) * z1mult(:,k) * z1(:,k) &
                   * soil%zse(k-1) / (soil%zse(k) + soil%zse(k-1) ) &
                   + dtt(:,k) * z1mult(:,k+1) * z1(:,k+1) &
                   * soil%zse(k+1) / (soil%zse(k) + soil%zse(k+1) )
           END DO
        END IF ! (nmeth == 4)
        DO k = 2, ms
          ! i.e. now K(k-.5)
          z1(:,k) = wbh(:,k) * z1(:,k)
        END DO
        ! the following top & bottom b.c.'s will preserve a uniform column
        !     z1(1) =z1(2)   ! simple dk/dz=0
        !     z1(ms+1)=z1(ms) ! simple dk/dz=0
        ! N.B. z1 are here +ve
        z1(:,1) = MIN(z1(:,2), z1(:,ms) )
        z1(:,ms + 1) = z1(:,1)
        ! no gravit. term if too much ice 11/12/00
        DO k = 1, ms
           IF (nmeth == 4) THEN
              WHERE (wbficemx < 0.75)
                 ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k) &
                      * ( (z1mult(:,k+1) - 1.0) * z1(:,k+1) &
                      - (z1mult(:,k) - 1.0) * z1(:,k) )
              END WHERE
           ELSE
              WHERE (wbficemx < 0.75)
                 ssoil%wblf(:,k) = ssoil%wblf(:,k) + dtt(:,k) &
                      * (z1(:,k) - z1(:,k+1) )
              END WHERE
           END IF
        END DO
     END IF ! (nmeth >= 2)

      ! Block below for testing purposes only:
!!$      IF (ntest > 0) THEN
!!$         wblfmx = MAXVAL(REAL(ssoil%wblf,r_1), 2)
!!$         wblfmn = MINVAL(REAL(ssoil%wblf,r_1), 2)
!!$         totwbb = 0.0
!!$         totwblb = 0.0
!!$         DO k = 1, ms
!!$            totwbb = totwbb + soil%zse(k) * REAL(ssoil%wb(:,k),r_1) ! diag.
!!$            totwblb = totwblb + soil%zse(k) * REAL(ssoil%wblf(:,k),r_1) ! diag.
!!$         END DO
!!$         ! Print variables for all active veg/soil patches in idjd grid point:
!!$         ! PRINT *, 'nmeth, b+2, 2b+3: ',nmeth, soil%ibp2(idjd), soil%i2bp3(idjd)
!!$         WRITE (6, " ('wb   ', 6f8.3) ") (ssoil%wb(landpt(idjd)%cstart: &
!!$              & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
!!$         WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(landpt(idjd)%cstart: &
!!$              & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
!!$         WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(landpt(idjd)%cstart: &
!!$              & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
!!$         WRITE (6, " ('wbh  ', 7f8.3) ") wbh(landpt(idjd)%cstart: &
!!$              & (landpt(idjd)%cstart+landpt(idjd)%nap-1),:)
!!$         WRITE (6, " ('ssatcurr', 6f8.3) ") ssatcurr(landpt(idjd)%cstart: &
!!$              & (landpt(idjd)%cstart+landpt(idjd)%nap-1),:)
!!$         ! PRINT * , 'pwb_wbh,ssoil%pwb_min* for ms ', pwb_wbh(idjd), &
!!$         !                          & soil%hsbh(idjd) * ssoil%pwb_min(idjd)
!!$         ! PRINT * , 'wblfmx,wblfmn', wblfmx(idjd), wblfmn(idjd)
!!$         ! PRINT * , 'zse ', soil%zse
!!$         ! PRINT * , 'zshh ', soil%zshh
!!$         ! PRINT * , 'at ', (at(idjd,k), k = 1, ms)
!!$         ! PRINT * , 'bt ', (bt(idjd,k), k = 1, ms)
!!$         ! PRINT * , 'ct ', (ct(idjd,k), k = 1, ms)
!!$      END IF
      IF (nmeth == 3) THEN
         ! artificial fix applied here for safety (explicit nmeth only)
         DO k = 1, ms
            ssoil%wblf(:,k) = &
                 & MAX(0.0_r_2, MIN(ssoil%wblf(:,k), 1.0_r_2) )
         END DO
      END IF ! (nmeth == 3)
      ssoil%wblf(:,1) = ssoil%wblf(:,1) &
           & + dtt(:,1) * ssoil%fwtop / rhowat
   END IF ! (nmeth > 0)

!  PRINT *, 'Before trimb'
! DO k = 1, ms
!  WRITE(46,'(i2,1x,a5,200e12.3)') k,'at = ', at(:,k)
!  WRITE(46,'(i2,1x,a5,200e12.3)') k,'bt = ', bt(:,k)
!  WRITE(46,'(i2,1x,a5,200e12.3)') k,'ct = ', ct(:,k)
! END DO
!  PRINT *, 'ssoil%wblf = ', ssoil%wblf
   CALL trimb(at, bt, ct, ssoil%wblf, ms)
!  PRINT *, 'After trimb'
 
   DO k = 1, ms
      ssatcurr(:,k) = soil%ssat - ssoil%wbice(:,k)
      ssoil%wb(:,k) = ssoil%wblf(:,k) * ssatcurr(:,k) + ssoil%wbice(:,k)
      ssoil%wbice(:,k) = MIN(ssoil%wbice(:,k), 0.99 * ssoil%wb(:,k) )
   END DO

   ! Block below for testing purposes only:
!!$   IF (ntest > 0) THEN
!!$      ! Print variables for all active veg/soil patches in idjd grid point:
!!$      !  PRINT * , 'at end of smoisturev,fwtop ', ssoil%fwtop(idjd)
!!$      !  PRINT * , 'tgg ', (ssoil%tgg(idjd,k), k = 1, ms)
!!$      WRITE (6, " ('wb   ', 6f8.3) ") (ssoil%wb(landpt(idjd)%cstart: &
!!$           (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
!!$      WRITE (6, " ('wbice', 6f8.3) ") (ssoil%wbice(landpt(idjd)%cstart: &
!!$           (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
!!$      WRITE (6, " ('wblf ', 6f8.3) ") (ssoil%wblf(landpt(idjd)%cstart: &
!!$           (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
!!$      totwbc = 0.
!!$      totwblc = 0.
!!$      zsetot = 0.
!!$      DO k = 1, ms
!!$         totwbc = totwbc + soil%zse(k) * REAL(ssoil%wb(:,k),r_1)
!!$         totwblc = totwblc + soil%zse(k) * REAL(ssoil%wblf(:,k),r_1)
!!$         zsetot = zsetot + soil%zse(k)
!!$      END DO
!!$      !      PRINT *, 'totwba,totwbb,totwbc ', totwba(idjd), totwbb(idjd),totwbc(idjd)
!!$      !      PRINT *, 'totwblb,totwblc ', totwblb(idjd), totwblc(idjd)
!!$      !      PRINT *, 'with totwbc/zsetot: ', totwbc(idjd) / zsetot(idjd)
!!$   END IF
 END SUBROUTINE smoisturev

  !-------------------------------------------------------------------------
  SUBROUTINE snowdensity (dels, ssoil)
    REAL(r_1), INTENT(IN)   :: dels   ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT) :: ssoil  ! soil+snow variables
    
    WHERE (ssoil%snowd > 0.1 .AND. ssoil%isflag == 0)
       ssoil%ssdn(:,1) = MIN(400.0,MAX(120.0, ssoil%ssdn(:,1) + dels &
            * ssoil%ssdn(:,1) * 3.1e-6 * EXP( -0.03 &
            * (273.15 - MIN(tfrz, ssoil%tgg(:,1) )) &
            - MERGE(0.046, 0.0, ssoil%ssdn(:,1) >= 150.0) &
            * (ssoil%ssdn(:,1) - 150.0) ) ))
       ssoil%ssdn(:,1) = MIN(400.0,ssoil%ssdn(:,1) + dels * 9.806 &
            * ssoil%ssdn(:,1) * 0.75 * ssoil%snowd &
            / (3.0e7 * EXP(0.021 * ssoil%ssdn(:,1) + 0.081 &
            * (273.15 - MIN(tfrz, ssoil%tgg(:,1))))))
       ssoil%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,1) ** 2 &
            + 0.074, 1.0) )
       ssoil%ssdnn = ssoil%ssdn(:,1)
    END WHERE
    WHERE (ssoil%isflag == 1)
       ssoil%ssdn(:,1) = ssoil%ssdn(:,1) + dels * ssoil%ssdn(:,1) * 3.1e-6 &
            * EXP( -0.03 * (273.15 - MIN(tfrz, ssoil%tggsn(:,1))) &
            - MERGE(0.046, 0.0, ssoil%ssdn(:,1) >= 150.0) &
            * (ssoil%ssdn(:,1) - 150.0) )
       ssoil%ssdn(:,2) = ssoil%ssdn(:,2) + dels * ssoil%ssdn(:,2) * 3.1e-6 &
            * EXP( -0.03 * (273.15 - MIN(tfrz, ssoil%tggsn(:,2))) &
            - MERGE(0.046, 0.0, ssoil%ssdn(:,2) >= 150.0) &
            * (ssoil%ssdn(:,2) - 150.0) )
       ssoil%ssdn(:,3) = ssoil%ssdn(:,3) + dels * ssoil%ssdn(:,3) * 3.1e-6 &
            * EXP( -0.03 * (273.15 - MIN(tfrz, ssoil%tggsn(:,3))) &
            - MERGE(0.046, 0.0, ssoil%ssdn(:,3) >= 150.0) &
            * (ssoil%ssdn(:,3) - 150.0) )
       ssoil%ssdn(:,1) = ssoil%ssdn(:,1) + dels * 9.806 * ssoil%ssdn(:,1) &
            * 0.05*ssoil%ssdn(:,1) &
            / (3.0e7 * EXP(.021 * ssoil%ssdn(:,1) + 0.081 &
            * (273.15 - MIN(tfrz, ssoil%tggsn(:,1)))))
       ssoil%ssdn(:,2) = ssoil%ssdn(:,2) + dels * 9.806 * ssoil%ssdn(:,2) &
            * (0.05 * ssoil%ssdn(:,1) + 0.5 * ssoil%smass(:,2) ) &
            / (3.0e7 * EXP(.021 * ssoil%ssdn(:,2) + 0.081 &
            * (273.15 - MIN(tfrz, ssoil%tggsn(:,2)))))
       ssoil%ssdn(:,3) = ssoil%ssdn(:,3) + dels * 9.806 * ssoil%ssdn(:,3) &
            * (0.05*ssoil%ssdn(:,1) + ssoil%smass(:,2) &
            + 0.5*ssoil%smass(:,3)) &
            / (3.0e7 * EXP(.021 * ssoil%ssdn(:,3) + 0.081 &
            * (273.15 - MIN(tfrz, ssoil%tggsn(:,3)))))
       ssoil%sdepth(:,1) =  ssoil%smass(:,1) / ssoil%ssdn(:,1) 
       ssoil%sdepth(:,2) =  ssoil%smass(:,2) / ssoil%ssdn(:,2) 
       ssoil%sdepth(:,3) =  ssoil%smass(:,3) / ssoil%ssdn(:,3) 
       ssoil%ssdnn = (ssoil%ssdn(:,1) * ssoil%smass(:,1) + ssoil%ssdn(:,2) &
            * ssoil%smass(:,2) + ssoil%ssdn(:,3) * ssoil%smass(:,3) ) &
            / ssoil%snowd
       ssoil%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,1) ** 2 &
            + 0.074, 1.0) )
       ssoil%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,2) ** 2 &
            + 0.074, 1.0) )
       ssoil%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,3) ** 2 &
            + 0.074, 1.0) )
    END WHERE
    
  END SUBROUTINE snowdensity

  !-------------------------------------------------------------------------
  SUBROUTINE snow_melting (dels, snowmlt, ktau, ssoil )
    REAL(r_1), INTENT(IN)                 :: dels   ! integration time step (s)
    REAL(r_1), DIMENSION(mp), INTENT(OUT) :: snowmlt ! snow melt   
    INTEGER(i_d), INTENT(IN)              :: ktau ! integration step number
    TYPE(soil_snow_type), INTENT(INOUT)   :: ssoil  ! soil+snow variables
    INTEGER(i_d), PARAMETER      :: ntest = 0 ! for snow diag prints
    INTEGER(i_d)                 :: k
    REAL(r_1), DIMENSION(mp)     :: osm
    REAL(r_1), DIMENSION(mp)     :: sgamm
    REAL(r_1), DIMENSION(mp,0:3) :: smelt1
    REAL(r_1), DIMENSION(mp)     :: snowflx
   
    snowmlt= 0.0
    smelt1 = 0.0
    WHERE (ssoil%snowd > 0.0 .AND. ssoil%isflag == 0 &
         .AND. ssoil%tgg(:,1) >= tfrz )
       ! snow covered land
       ! following done in sflux  via  ga= ... +cls*egg + ...
       ! ** land,snow,melting
       snowflx = REAL((ssoil%tgg(:,1) - tfrz) * ssoil%gammzz(:,1),r_1)
       ! prevent snow depth going negative
       snowmlt = MIN(snowflx / hlf, ssoil%snowd )
       ssoil%snowd = ssoil%snowd - snowmlt
       ssoil%tgg(:,1) = &
            & REAL(ssoil%tgg(:,1) - snowmlt * hlf / ssoil%gammzz(:,1),r_1)
    END WHERE
 
    smelt1(:,0) = 0.0
    DO k = 1, 3
       WHERE (ssoil%snowd > 0.0 .AND. ssoil%isflag > 0)
          sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
          ! snow melt refreezing
          snowflx = smelt1(:,k-1) * hlf / dels
          ssoil%tggsn(:,k) = ssoil%tggsn(:,k) + snowflx * dels / sgamm
          ! increase density due to snowmelt
          osm = ssoil%smass(:,k)
          ssoil%smass(:,k) = ssoil%smass(:,k) + smelt1(:,k-1)
          ssoil%ssdn(:,k) = MAX(120.0,MIN(ssoil%ssdn(:,k)*osm/ssoil%smass(:,k) &
                  + rhowat*(1.0-osm/ssoil%smass(:,k)), 400.0))
          ssoil%sdepth(:,k) = ssoil%smass(:,k) / ssoil%ssdn(:,k)
          sgamm = ssoil%smass(:,k) * cgsnow
          smelt1(:,k-1) = 0.0
          smelt1(:,k) = 0.0
          ! snow melting
          WHERE (ssoil%tggsn(:,k) > tfrz)
             snowflx = (ssoil%tggsn(:,k) - tfrz) * sgamm
             smelt1(:,k) = MIN(snowflx / hlf, 0.9 * ssoil%smass(:,k) )
             osm = ssoil%smass(:,k)
             ssoil%smass(:,k) = ssoil%smass(:,k) - smelt1(:,k)
             ssoil%tggsn(:,k) = MIN(ssoil%tggsn(:,k) - smelt1(:,k)*hlf/sgamm, &
                  & tfrz)
             ssoil%sdepth(:,k) = ssoil%smass(:,k) / ssoil%ssdn(:,k)
          END WHERE
       END WHERE
    END DO
    WHERE (ssoil%snowd > 0.0 .AND. ssoil%isflag > 0)
       snowmlt = smelt1(:,1) + smelt1(:,2) + smelt1(:,3)
       ssoil%snowd = ssoil%snowd - snowmlt
    END WHERE
    
  END SUBROUTINE snow_melting

  !-------------------------------------------------------------------------
  SUBROUTINE snow_accum (dels,  ktau, canopy, met, ssoil, soil)
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    INTEGER(i_d), INTENT(IN)                 :: ktau ! integration step number
    TYPE(canopy_type), INTENT(INOUT)         :: canopy ! vegetation variables
    TYPE(met_type), INTENT(INOUT)            :: met   ! all met forcing
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    INTEGER(i_d), PARAMETER  :: ntest = 0 ! for snow diag prints
    INTEGER(i_d), PARAMETER  :: nglacier = 0 ! 0 original, 1 off, 2 new Eva
!    REAL(r_1), DIMENSION(mp) :: evapsn ! EAK aug08
    INTEGER(i_d)             :: k
    REAL(r_1), DIMENSION(mp) :: osm
    REAL(r_1), DIMENSION(mp) :: sgamm
    REAL(r_1), DIMENSION(mp) :: snowmlt
    REAL(r_1), DIMENSION(mp) :: xxx
    
    ! Diagnostic block below:
    IF (ntest > 0) THEN
      ! Print variables for all active veg/soil patches in idjd grid point:
      PRINT *, 'surfb1'
      PRINT *, 'entering snow_accum', &
            & canopy%precis(landpt(idjd)%cstart: &
            &              (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'osnowd,snowd,isflag', &
            & ssoil%osnowd(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & ssoil%snowd(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
            & ssoil%isflag(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'tggsn ', (ssoil%tggsn(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, 3)
      PRINT *, 'tgg ', (ssoil%tgg(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
      PRINT *, 'wb ', (ssoil%wb(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
      PRINT *, 'wbice ', (ssoil%wbice(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
      PRINT *, 'gammzz ', (ssoil%gammzz(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
      PRINT *, 'albnew ', (ssoil%albsoilsn(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, 2)
    END IF
    
    snowmlt =0.0
   
    ! canopy%precis is both liquid and snow
    WHERE (canopy%precis > 0.0 .AND. ssoil%isflag == 0)
       ssoil%snowd = MAX(ssoil%snowd + met%precip_s, 0.0) ! accumulate solid part
       canopy%precis = canopy%precis - met%precip_s
       ssoil%ssdn(:,1) = MAX(120.0, ssoil%ssdn(:,1) &
            * ssoil%osnowd / MAX(0.01, ssoil%snowd) &
            + 120.0 * met%precip_s / MAX(0.01, ssoil%snowd))
       ssoil%ssdnn = ssoil%ssdn(:,1)
       WHERE (canopy%precis > 0.0 .AND. ssoil%tgg(:,1) < tfrz)
          ssoil%snowd = MAX(ssoil%snowd + canopy%precis, 0.0)
          ssoil%tgg(:,1) = ssoil%tgg(:,1) + canopy%precis * hlf &
               / REAL(ssoil%gammzz(:,1),r_1)
          ! change density due to water being added 
          ssoil%ssdn(:,1) = MIN(400.0, MAX(120.0, ssoil%ssdn(:,1) &
               * ssoil%osnowd / MAX(0.01, ssoil%snowd) &
               + rhowat * canopy%precis / MAX(0.01, ssoil%snowd)))
          ! ssoil%tgg(:,1) = ssoil%tgg(:,1) + canopy%precis * hlf  &
          !                / (ssoil%gammzz(:,1) + cswat * canopy%precis)
          canopy%precis = 0.0
          ssoil%ssdnn = ssoil%ssdn(:,1)
       END WHERE
    END WHERE ! (canopy%precis > 0. .and. ssoil%isflag == 0) 
   
    WHERE (canopy%precis > 0.0 .AND. ssoil%isflag > 0)
       ! add solid precip
       ssoil%snowd = MAX(ssoil%snowd + met%precip_s, 0.0)
       canopy%precis = canopy%precis - met%precip_s  ! remaining liquid precip
       ! update top snow layer with fresh snow
       osm = ssoil%smass(:,1)
       ssoil%smass(:,1) = ssoil%smass(:,1) + met%precip_s
       ssoil%ssdn(:,1) = MAX(120.0,ssoil%ssdn(:,1)*osm/ssoil%smass(:,1) &
            + 120.0 * met%precip_s/ssoil%smass(:,1))
       ssoil%sdepth(:,1) = MAX(0.02,ssoil%smass(:,1) / ssoil%ssdn(:,1))
       ! add liquid precip
       WHERE (canopy%precis > 0.0)
          ssoil%snowd = MAX(ssoil%snowd + canopy%precis, 0.0)
          sgamm = ssoil%ssdn(:,1) * cgsnow * ssoil%sdepth(:,1)
          osm = ssoil%smass(:,1)
          ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + canopy%precis * hlf &
               * osm / (sgamm * ssoil%osnowd )
          !                        * ssoil%smass(:,1) / (sgamm * ssoil%osnowd )
          ssoil%smass(:,1) = ssoil%smass(:,1) + canopy%precis &
               * osm/ssoil%osnowd
          !                        * ssoil%smass(:,1)/ssoil%osnowd
          ssoil%ssdn(:,1) = MAX(120.0,MIN(ssoil%ssdn(:,1)*osm/ssoil%smass(:,1) &
               +  rhowat*(1.0-osm/ssoil%smass(:,1)), 400.0))
          ssoil%sdepth(:,1) = ssoil%smass(:,1)/ssoil%ssdn(:,1)
          
          sgamm = ssoil%ssdn(:,2) * cgsnow * ssoil%sdepth(:,2)
          osm = ssoil%smass(:,2)
          ssoil%tggsn(:,2) = ssoil%tggsn(:,2) + canopy%precis * hlf &
               * osm / (sgamm * ssoil%osnowd )
          !                        * ssoil%smass(:,2) / (sgamm * ssoil%osnowd )
          ssoil%smass(:,2) = ssoil%smass(:,2) + canopy%precis &
               * osm/ssoil%osnowd
          !                        ssoil%smass(:,2)/ssoil%osnowd
          ssoil%ssdn(:,2) = MAX(120.0,MIN(ssoil%ssdn(:,2)*osm/ssoil%smass(:,2) &
               + rhowat*(1.0-osm/ssoil%smass(:,2)), 400.0))
          ssoil%sdepth(:,2) = ssoil%smass(:,2) / ssoil%ssdn(:,2)
          
          sgamm = ssoil%ssdn(:,3) * cgsnow * ssoil%sdepth(:,3)
          osm = ssoil%smass(:,3)
          ssoil%tggsn(:,3) = ssoil%tggsn(:,3) + canopy%precis * hlf &
               * osm / (sgamm * ssoil%osnowd )
          !                        * ssoil%smass(:,3) / (sgamm * ssoil%osnowd )
          ssoil%smass(:,3) = ssoil%smass(:,3) + canopy%precis &
               * osm/ssoil%osnowd
          !                        * ssoil%smass(:,3)/ssoil%osnowd
          ssoil%ssdn(:,3) = MAX(120.0,MIN(ssoil%ssdn(:,3)*osm/ssoil%smass(:,3) &
               + rhowat*(1.0-osm/ssoil%smass(:,3)), 400.0))
          ssoil%sdepth(:,3) = ssoil%smass(:,3)/ssoil%ssdn(:,3)
          
          canopy%precis = 0.0
          
       END WHERE
    END WHERE
       
    WHERE (ssoil%snowd < 0.1 .AND. canopy%fes > 0.0)
       canopy%fes = MIN(canopy%fes, &
            MAX(0.0,(REAL(ssoil%wb(:,1),r_1)-soil%swilt))* soil%zse(1) &
            * 1000.0 * hl / dels)
       canopy%fes = MIN(canopy%fes, &
            REAL((ssoil%wb(:,1)-ssoil%wbice(:,1)),r_1) * soil%zse(1) &
            * 1000.0 * hl / dels)
    END WHERE
    ! Calculate soil evaporation total in mm (from W/m2):
    canopy%segg = canopy%fes / hl
    WHERE (ssoil%snowd > 0.1) canopy%segg = canopy%fes / (hl+hlf) ! EAK aug08
    ! Initialise snow evaporation:
    ssoil%evapsn = 0
    ! Snow evaporation and melting
    WHERE (ssoil%snowd > 0.1)
       ssoil%evapsn = dels * canopy%fes / ( hl + hlf ) ! EAK aug08
       xxx = ssoil%evapsn
       WHERE ( ssoil%isflag == 0 .AND. canopy%fes.gt.0.0) &
            & ssoil%evapsn = MIN(ssoil%snowd, xxx ) 
       WHERE ( ssoil%isflag  > 0 .AND. canopy%fes.gt.0.0) &
            & ssoil%evapsn = MIN(0.9*ssoil%smass(:,1), xxx )
       ssoil%snowd = ssoil%snowd - ssoil%evapsn
       WHERE ( ssoil%isflag > 0 )
          ssoil%smass(:,1) = ssoil%smass(:,1)  - ssoil%evapsn
          ssoil%sdepth(:,1) = MAX(0.02,ssoil%smass(:,1) / ssoil%ssdn(:,1))
       END WHERE
       canopy%segg = 0.0
       canopy%fes = ssoil%evapsn * (hl + hlf) / dels ! return for hyd. balance
    END WHERE

  END SUBROUTINE snow_accum 

  !-------------------------------------------------------------------------
  SUBROUTINE surfbv (dels, ktau, met, ssoil, soil )
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    INTEGER(i_d), INTENT(IN)                 :: ktau ! integration step number
    TYPE(met_type), INTENT(INOUT)            :: met    ! all met forcing
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    INTEGER(i_d), PARAMETER      :: ntest = 0 ! for snow diag prints
    INTEGER(i_d), PARAMETER      :: nglacier = 2 ! 0 original, 1 off, 2 new Eva
    INTEGER(i_d)                 :: k
    REAL(r_1), DIMENSION(mp)     :: rnof5
    REAL(r_1), DIMENSION(mp)     :: sgamm
    REAL(r_1), DIMENSION(mp)     :: smasstot
    REAL(r_1), DIMENSION(mp,0:3) :: smelt1
    REAL(r_2), DIMENSION(mp)     :: xxx

    CALL smoisturev ( dels, ktau, ssoil, soil)
!   PRINT *, 'After smoisturev'

    ! Diagnostic block below:
    IF (ntest > 0) THEN
    ! Print variables for all active veg/soil patches in idjd grid point:
    PRINT *, 'in surfbv after smoisturev '
    PRINT *, 'osnowd,snowd,isflag,ssat,runoff', &
    ssoil%osnowd(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
    ssoil%snowd(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
    ssoil%isflag(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
    soil%ssat(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
    ssoil%runoff(landpt(idjd)%cstart:(landpt(idjd)%cstart+landpt(idjd)%nap-1))
    PRINT *, 'tggsn_d ', (ssoil%tggsn(landpt(idjd)%cstart: &
          & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, 3)
    END IF
    DO k = 1, ms
       xxx = REAL(soil%ssat,r_2)
       ssoil%rnof1 = ssoil%rnof1 + REAL((MAX(ssoil%wb(:,k) - xxx, 0.0_r_2) &
            & * 1000.0),r_1) * soil%zse(k)
       ssoil%wb(:,k) = MIN( ssoil%wb(:,k), xxx )
       !      ssoil%wb(:,k) = MIN( ssoil%wb(:,k), soil%ssat )
    END DO
    ! for deep runoff use wb-sfc, but this value not to exceed .99*wb-wbice
    ! account for soil/ice cracking
    !   fracm = MIN(0.2, 1. - MIN(1., ssoil%wb(:,ms) / soil%sfc ) )
    !   ssoil%wb(:,ms) = ssoil%wb(:,ms) &
    !                  + fracm*ssoil%rnof1/(1000.0*soil%zse(ms))
    !   ssoil%rnof1 = (1. - fracm) * ssoil%rnof1 
    ! the code below is replaced, see subroutine smoistv 
    !   tmp = MAX(MIN(ssoil%wb(:,ms) - soil%sfc, .99 * ssoil%wb(:,ms) &
    !       - ssoil%wbice(:,ms) ) * soil%c3 / 86400., 0.)
    !   ssoil%rnof2 = soil%zse(ms) * 1000. * tmp * dels
    !   ssoil%wb(:,ms) = ssoil%wb(:,ms) - tmp * dels
    ssoil%runoff = ssoil%rnof1 + ssoil%rnof2
    ssoil%wbtot = 0.0
    DO k = 1, ms
       ssoil%wbtot = ssoil%wbtot &
            & + REAL(ssoil%wb(:,k),r_1) * 1000.0 * soil%zse(k)
    END DO
    !
    !---  glacier formation
    IF (nglacier == 2) THEN
       WHERE (ssoil%snowd > 1000.0)
          rnof5 = ssoil%snowd - 1000.0
          ssoil%runoff = ssoil%runoff + rnof5
          !---- change local tg to account for energy - clearly not best method
          WHERE (ssoil%isflag == 0)
             smasstot = 0.0
             ssoil%tgg(:,1) = ssoil%tgg(:,1) - rnof5 * hlf &
                  & / REAL(ssoil%gammzz(:,1),r_1)
             ssoil%snowd = 1000.0
          ELSEWHERE
             smasstot = ssoil%smass(:,1) + ssoil%smass(:,2) + ssoil%smass(:,3)
          END WHERE
       END WHERE
       
       DO k = 1, 3
          WHERE (ssoil%snowd > 1000.0 .AND. ssoil%isflag > 0)
             sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
             smelt1(:,k) = MIN(rnof5 * ssoil%smass(:,k) / smasstot, &
                  & 0.9 * ssoil%smass(:,k) )
             ssoil%smass(:,k) = ssoil%smass(:,k) - smelt1(:,k)
             ssoil%snowd = ssoil%snowd - smelt1(:,k)
             ssoil%tggsn(:,k) = ssoil%tggsn(:,k) - smelt1(:,k) * hlf / sgamm
          END WHERE
       END DO
    END IF

    ! Diagnostic block below:
    IF (ntest > 0) THEN
       ! Print variables for all active veg/soil patches in idjd grid point:
       PRINT *, 'end surfbv  rnof1,runoff ',&
            & ssoil%rnof1(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
            & ssoil%runoff(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1))

      sgamm = ssoil%ssdn(:,1) * cgsnow * ssoil%sdepth(:,1)

      PRINT *, 'snowd,isflag,sgamm ', &
            & ssoil%snowd(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
            & ssoil%isflag(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
            & sgamm(landpt(idjd)%cstart: &
            &      (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'tggsn_d ', (ssoil%tggsn(landpt(idjd)%cstart: &
            (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, 3)
      PRINT *, 'tgg ', (ssoil%tgg(landpt(idjd)%cstart: &
            (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
      PRINT *, 'wb ', (ssoil%wb(landpt(idjd)%cstart: &
            (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
    END IF
  END SUBROUTINE surfbv

  !-------------------------------------------------------------------------
  SUBROUTINE snow_albedo (dels, ktau, met, ssoil, soil )
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    INTEGER(i_d), INTENT(IN)                 :: ktau ! integration step number
    TYPE(met_type), INTENT(INOUT)            :: met    ! all met forcing
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil  ! soil+snow variables
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil ! soil parameters
    INTEGER(i_d), PARAMETER      :: ntest = 0 ! for snow diag prints
    REAL(r_1), DIMENSION(mp)     :: alv ! Snow albedo for visible
    REAL(r_1), DIMENSION(mp)     :: alir ! Snow albedo for near infra-red (NIR)
    REAL(r_1), PARAMETER         :: alvo  = 0.95 ! albedo for vis. on new snow
    REAL(r_1), PARAMETER         :: aliro = 0.65 ! albedo for NIR on new snow
    REAL(r_1), DIMENSION(mp)     :: ar1 ! crystal growth  (-ve)
    REAL(r_1), DIMENSION(mp)     :: ar2 ! freezing of melt water
    REAL(r_1), DIMENSION(mp)     :: ar3
    REAL(r_1), DIMENSION(mp)     :: dnsnow ! new snow albedo
    REAL(r_1), DIMENSION(mp)     :: dtau
    REAL(r_1), DIMENSION(mp)     :: fage !age factor
    REAL(r_1), DIMENSION(mp)     :: fzenm
    REAL(r_1), DIMENSION(mp)     :: sfact
    REAL(r_1), DIMENSION(mp)     :: snr
    REAL(r_1), DIMENSION(mp)     :: snrat
    REAL(r_1), DIMENSION(mp)     :: talb ! snow albedo
    REAL(r_1), DIMENSION(mp)     :: tmp ! temporary value

!$$ Removed because soil%albsoil is now in radiation bands (BP 21Oct2009)
!    !	 calculate soil/snow albedo
!    !	  cuvrf(i,1) = albsav(iq) ! use surface albedo from indata
!    sfact = 0.68
!    WHERE (soil%albsoil <= 0.14)
!       sfact = 0.5
!    ELSEWHERE (soil%albsoil > 0.14 .AND. soil%albsoil <= 0.20)
!       sfact = 0.62
!    END WHERE
!    ssoil%albsoilsn(:,2) = 2.0 * soil%albsoil / (1.0 + sfact)
!    ssoil%albsoilsn(:,1) = sfact * ssoil%albsoilsn(:,2)
!$$ End of removal (BP 21Oct2009)
    ! new snow albedo (needs osnowd from the previous dels)
    dnsnow = MIN( 1.0, 0.1 * MAX(0.0,ssoil%snowd-ssoil%osnowd) ) ! new snow in cm
    ! Snow age depends on snow crystal growth, freezing of melt water,
    ! accumulation of dirt and amount of new snow.
    tmp = ssoil%isflag*ssoil%tggsn(:,1) + (1 - ssoil%isflag )*ssoil%tgg(:,1)
    tmp = MIN(tmp, 273.15)
    ar1 = 5000.0 * (1.0 / 273.15 - 1.0 / tmp) ! crystal growth  (-ve)
    ar2 = 10.0 * ar1 ! freezing of melt water
    snr = ssoil%snowd / MAX(ssoil%ssdnn, 100.0)
    ! fixes for Arctic & Antarctic
    WHERE (soil%isoilm == 9)
       ar3 = 0.0005
       dnsnow = MAX(dnsnow, 0.002) !increase refreshing of snow in Antarctic
       snrat = MIN(1.0, snr / (snr + 0.001) )
    ELSEWHERE
       ! accumulation of dirt
       ar3 = 0.1
       snrat = MIN(1.0, snr / (snr + 0.01) )
    END WHERE
    dtau = 1.0e-6 * (EXP(ar1) + EXP(ar2) + ar3) * dels
    WHERE (ssoil%snowd <= 1.0)
       ssoil%snage = 0.0
    ELSEWHERE
       ssoil%snage = MAX(0.0, (ssoil%snage + dtau) * (1.0 - dnsnow) )
    END WHERE
    fage = 1.0 - 1.0 / (1.0 + ssoil%snage ) !age factor
    !
    ! Snow albedo is dependent on zenith angle and  snow age.
    ! albedo zenith dependence
    ! alvd = alvo * (1.0-cs*fage); alird = aliro * (1.-cn*fage)
    ! where cs = 0.2, cn = 0.5, b = 2.0
    tmp = MAX(0.17365, met%coszen )
    tmp = MAX(0.01, met%coszen )
    fzenm = MAX(merge(0.0, (1.0+0.5)/(1.0+4.0*tmp)-0.5, tmp > 0.5), 0.0)
    tmp = alvo * (1.0 - 0.2 * fage)
    alv = 0.4 * fzenm * (1.0 - tmp) + tmp
    tmp = aliro * (1.0 - 0.5 * fage)
    alir = 0.4 * fzenm * (1.0 - tmp) + tmp
    ! talb = 0.5 * (alv + alir) ! snow albedo
    ! alss = (1. - snrat)*soil%albsoil + snrat*talb ! canopy free surf albedo
!$$ Changed because soil%albsoil is now in radiation bands (BP 21Oct2009)
    ssoil%albsoilsn(:,2) = (1.0 - snrat) * soil%albsoil(:,2) + snrat * alir
    ssoil%albsoilsn(:,1) = (1.0 - snrat) * soil%albsoil(:,1) + snrat * alv
!    ssoil%albsoilsn(:,2) = (1.0 - snrat) * ssoil%albsoilsn(:,2) + snrat * alir
!    ssoil%albsoilsn(:,1) = (1.0 - snrat) * ssoil%albsoilsn(:,1) + snrat * alv
    
  END SUBROUTINE snow_albedo 

  !-------------------------------------------------------------------------
  ! SUBROUTINE stempv
  !	 calculates temperatures of the soil
  !	 tgg - new soil/snow temperature
  !	 ga - heat flux from the atmosphere (ground heat flux)
  !	 ccnsw - soil thermal conductivity, including water/ice
  !
  SUBROUTINE stempv(dels, canopy, ssoil, soil)
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    INTEGER(i_d), PARAMETER          :: ntest = 0
    REAL(r_2), DIMENSION(mp, -2:ms)  :: at
    REAL(r_2), DIMENSION(mp, -2:ms)  :: bt
    REAL(r_2), DIMENSION(mp, -2:ms)  :: ct
    REAL(r_2), DIMENSION(mp,ms)      :: ccnsw  ! soil thermal
                                               ! conductivity (incl water/ice)
    REAL(r_1), DIMENSION(mp)         :: coefa
    REAL(r_1), DIMENSION(mp)         :: coefb
    REAL(r_2), DIMENSION(mp)         :: dtg
    REAL(r_2), DIMENSION(mp)         :: ew
    REAL(r_2), DIMENSION(mp,-2:ms+1) :: coeff
    INTEGER(i_d)                     :: k
    REAL(r_1), DIMENSION(mp,-2:ms)   :: rhs
    REAL(r_1), DIMENSION(mp)         :: sgamm
    REAL(r_2), DIMENSION(mp,ms+3)    :: tmp_mat ! temp. matrix for tggsn & tgg
    REAL(r_2), DIMENSION(mp)         :: xx
    REAL(r_2), DIMENSION(mp)         :: wblfsp 

    !
    at = 0.0
    bt = 1.0
    ct = 0.0
    coeff = 0.0
    DO k = 1, ms
       WHERE (soil%isoilm == 9)
          ccnsw(:,k) = 1.5
       ELSEWHERE
          ew = ssoil%wblf(:,k) * soil%ssat
          ccnsw(:,k) = MIN(soil%cnsd * EXP( ew * LOG(60.0) + ssoil%wbfice(:,k) &
               & * soil%ssat * LOG(250.0) ), 1.5_r_2) &
               & * MAX(1.0_r_2, SQRT( MIN( 2.0_r_2, 0.5 * soil%ssat &
               & / MIN(ew, 0.5_r_2 * soil%ssat )) ) )
       END WHERE
    END DO
   
    WHERE (ssoil%isflag == 0)
       xx = MAX(0., ssoil%snowd / ssoil%ssdnn )
       ccnsw(:,1) = (ccnsw(:,1) - 0.2) * (soil%zse(1)/(soil%zse(1) + xx)) + 0.2
    END WHERE
    
    DO k = 3, ms
       WHERE (ssoil%isflag == 0)
          coeff(:,k) = 2.0 / (soil%zse(k-1)/ccnsw(:,k-1)+soil%zse(k)/ccnsw(:,k))
       END WHERE
    END DO
    k = 1
    WHERE (ssoil%isflag == 0)
       coeff(:,2) = 2.0 / ((soil%zse(1)+xx)/ccnsw(:,1)+soil%zse(2)/ccnsw(:,2))
       coefa = 0.0
       coefb = REAL(coeff(:,2),r_1)
       wblfsp = ssoil%wblf(:,k)
       xx=soil%css * soil%rhosoil
       ssoil%gammzz(:,k) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
            & + soil%ssat * (wblfsp * cswat * rhowat + ssoil%wbfice(:,k) &
            & * csice * rhowat * 0.9), xx ) * soil%zse(k)
       ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
       dtg = dels / ssoil%gammzz(:,k)
       at(:,k) = - dtg * coeff(:,k)
       ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
       bt(:,k) = 1.0 - at(:,k) - ct(:,k)
    END WHERE
   
    DO k = 2, ms
       WHERE (ssoil%isflag == 0)
          wblfsp = ssoil%wblf(:,k)
          xx=soil%css * soil%rhosoil
          ssoil%gammzz(:,k) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
                  & + soil%ssat * (wblfsp * cswat * rhowat + ssoil%wbfice(:,k) &
                  & * csice * rhowat * 0.9), xx ) * soil%zse(k)
          dtg = dels / ssoil%gammzz(:,k)
          at(:,k) = - dtg * coeff(:,k)
          ct(:,k) = - dtg * coeff(:,k+1) ! c3(ms)=0 & not really used
          bt(:,k) = 1.0 - at(:,k) - ct(:,k)
       END WHERE
    END DO
    WHERE (ssoil%isflag == 0)
       bt(:,1) = bt(:,1) - canopy%dgdtg * dels / ssoil%gammzz(:,1)
       ssoil%tgg(:,1) = ssoil%tgg(:,1) + (canopy%ga - ssoil%tgg(:,1) &
            & * REAL(canopy%dgdtg,r_1)) * dels / REAL(ssoil%gammzz(:,1),r_1)
    END WHERE
   
    ! Diagnostic block:
    IF (ntest > 0) THEN
      ! Print variables for all active veg/soil patches in idjd grid point:
      PRINT *, 'tgg1,ga,gammzz ', &
            & ssoil%tgg(landpt(idjd)%cstart: &
            &          (landpt(idjd)%cstart+landpt(idjd)%nap-1),1), &
            & canopy%ga(landpt(idjd)%cstart: &
            &          (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & ssoil%gammzz(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1),1)
      PRINT *, 'dgdtg,degdt,dfgdt ', &
            & canopy%dgdtg(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'ssat,css,rhos,cswat,rhowat,csice ', &
            & soil%ssat(landpt(idjd)%cstart: &
            &          (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & soil%css(landpt(idjd)%cstart: &
            &         (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & soil%rhosoil(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
            & cswat, rhowat, csice
      PRINT *, 'wblf1,wbfice1,zse1,cgsnow ', &
            & ssoil%wblf(landpt(idjd)%cstart: &
            &           (landpt(idjd)%cstart+landpt(idjd)%nap-1),1), &
            & ssoil%wbfice(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1),1),&
            & soil%zse(1), cgsnow
      PRINT *, 'at ', (at(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k) , k = 1, ms)
      PRINT *, 'bt ', (bt(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k) , k = 1, ms)
      PRINT *, 'ct ', (ct(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k) , k = 1, ms)
      PRINT *, 'rhs ', (ssoil%tgg(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k) , k = 1, ms)
    END IF
    coeff(:,1-3) = 0.0
    ! 3-layer snow points done here
    WHERE (ssoil%isflag /= 0)
       ssoil%sconds(:,1) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,1)**2 &
            & + 0.074, 1.0) )
       ssoil%sconds(:,2) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,2)**2 &
            & + 0.074, 1.0) )
       ssoil%sconds(:,3) = MAX(0.2, MIN(2.876e-6 * ssoil%ssdn(:,3)**2 &
            & + 0.074, 1.0) )
       coeff(:,-1) = 2.0 / (ssoil%sdepth(:,1) / ssoil%sconds(:,1) &
            & + ssoil%sdepth(:,2) / ssoil%sconds(:,2) )
       coeff(:,0) = 2.0 / (ssoil%sdepth(:,2) / ssoil%sconds(:,2) &
            & + ssoil%sdepth(:,3) / ssoil%sconds(:,3) )
       coeff(:,1) = 2.0 / (ssoil%sdepth(:,3) / ssoil%sconds(:,3) &
            & + soil%zse(1) / ccnsw (:,1) )
    END WHERE
    
    DO k = 2, ms
       WHERE (ssoil%isflag /= 0)
          coeff(:,k) = 2.0 / (soil%zse(k-1)/ccnsw(:,k-1)+soil%zse(k)/ccnsw(:,k))
       END WHERE
    END DO
    WHERE (ssoil%isflag /= 0)
       coefa = REAL(coeff (:,-1),r_1)
       coefb = REAL(coeff (:,1),r_1)
    END WHERE
    DO k = 1, 3
       WHERE (ssoil%isflag /= 0)
          sgamm = ssoil%ssdn(:,k) * cgsnow * ssoil%sdepth(:,k)
          dtg = dels / sgamm
          ! rhs(k-3) = ssoil%tggsn(:,k)  ! A
          at(:,k-3) = - dtg * coeff(:,k-3)
          ct(:,k-3) = - dtg * coeff(:,k-2)
          bt(:,k-3) = 1.0 - at(:,k-3) - ct(:,k-3)
       END WHERE
    END DO
    DO k = 1, ms
       WHERE (ssoil%isflag /= 0)
          wblfsp = ssoil%wblf(:,k)
          xx=soil%css * soil%rhosoil
          ssoil%gammzz(:,k) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
               & + soil%ssat * (wblfsp * cswat * rhowat + ssoil%wbfice(:,k) &
               & * csice * rhowat * 0.9), xx ) * soil%zse(k)
          dtg = dels / ssoil%gammzz(:,k)
          at(:,k) = - dtg * coeff(:,k)
          ct(:,k) = - dtg * coeff(:,k + 1) ! c3(ms)=0 & not really used
          bt(:,k) = 1.0 - at(:,k) - ct(:,k)
       END WHERE
    END DO 
    WHERE (ssoil%isflag /= 0)
       sgamm = ssoil%ssdn(:,1) * cgsnow * ssoil%sdepth(:,1)
       ! rhs(1-3) = rhs(1-3)+canopy%ga*dels/sgamm
       ! new code
       bt(:,-2) = bt(:,-2) - canopy%dgdtg * dels / sgamm
       ssoil%tggsn(:,1) = ssoil%tggsn(:,1) + (canopy%ga - ssoil%tggsn(:,1) &
            & * REAL(canopy%dgdtg,r_1) ) * dels / sgamm
       rhs(:,1-3) = ssoil%tggsn(:,1)
    END WHERE
 
    ! Diagnostic block:
    IF (ssoil%isflag(idjd) /= 0 .AND.  ntest > 0) THEN
      ! Print variables for all active veg/soil patches in idjd grid point:
      PRINT *, 'in stempv 3-layer snow code '
      PRINT *, 'ccnsw ', (ccnsw(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, ms)
      PRINT *, 'sdepth d ', (ssoil%sdepth(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, 3)
      PRINT *, 'sconds ', (ssoil%sconds(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1, 3)
      PRINT *, 'coeff ', coeff
      PRINT *, 'at ', (at(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = - 2, ms)
      PRINT *, 'bt ', (bt(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = - 2, ms)
      PRINT *, 'ct ', (ct(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = - 2, ms)
      PRINT *, 'rhs(tggsn,tgg) ', &
            & (ssoil%tggsn(landpt(idjd)%cstart: &
            &           (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1,3), &
            & (ssoil%tgg(landpt(idjd)%cstart: &
            &           (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), k = 1,ms)
      PRINT *, 'tggsn,ga,sgamm ', &
            & ssoil%tggsn(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1),k), &
            & canopy%ga(landpt(idjd)%cstart: &
            &          (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & sgamm(landpt(idjd)%cstart: &
            &      (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'dgdtg ', canopy%dgdtg(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1))
    END IF
    !
    !     note in the following that tgg and tggsn are processed together
    tmp_mat(:,1:3) = REAL(ssoil%tggsn,r_2) 
    tmp_mat(:,4:(ms+3)) = REAL(ssoil%tgg,r_2)

    CALL trimb (at, bt, ct, tmp_mat, ms + 3)
   
    ssoil%tggsn = REAL(tmp_mat(:,:3),r_1)
    ssoil%tgg   = REAL(tmp_mat(:,4:(ms+3)),r_1)
    canopy%sghflux = coefa * (ssoil%tggsn(:,1) - ssoil%tggsn(:,2) )
    !      canopy%sghflux = coefb * (ssoil%tggsn(:,3) - ssoil%tgg(:,1) )
    canopy%ghflux = coefb * (ssoil%tgg(:,1) - ssoil%tgg(:,2) ) ! +ve downwards
    
  END SUBROUTINE stempv
  !-------------------------------------------------------------------------
  SUBROUTINE snowcheck(dels, ktau, ssoil )
    REAL(r_1), INTENT(IN)               :: dels ! integration time step (s)
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number
    TYPE(soil_snow_type), INTENT(INOUT) :: ssoil
    INTEGER(i_d), PARAMETER :: ntest = 0 !  for prints
    INTEGER(i_d)            :: k

    ! Diagnostic block:
    IF (ntest > 0) THEN
      ! Print variables for all active veg/soil patches in idjd grid point:
      PRINT *, 'in snowcheck ktau= ', ktau, ssoil%isflag(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1)), idjd
      PRINT *, 'dels,ssdn ', dels, (ssoil%ssdn(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,3)
      PRINT *, 'osnowd,snowd', &
            & ssoil%osnowd(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
            & ssoil%snowd(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'tgg ', (ssoil%tgg(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'tggsn ', (ssoil%tggsn(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,3)
    END IF

    ! The following initialization could be commented out as it is already
    ! done in subroutines write_default_params, get_default_inits and
    ! get_restart_data (BP may2010) 
    IF (ktau <= 1) THEN ! initialisations for first time step
      ssoil%ssdn = 120.0
      ssoil%ssdnn = 120.0
      ssoil%tggsn = tfrz
      ssoil%isflag = 0
    END IF

    WHERE (ssoil%snowd <= 0.0) ! i.e. no snow:
       ssoil%isflag = 0
       ssoil%ssdn(:,1) = 120.0
       ssoil%ssdn(:,2) = 120.0
       ssoil%ssdn(:,3) = 120.0
       ssoil%ssdnn = 120.0
       ssoil%tggsn(:,1) = tfrz
       ssoil%tggsn(:,2) = tfrz
       ssoil%tggsn(:,3) = tfrz
       ssoil%sdepth(:,1) = ssoil%snowd / ssoil%ssdn(:,1)
       ssoil%sdepth(:,2) = ssoil%snowd / ssoil%ssdn(:,2)
       ssoil%sdepth(:,3) = ssoil%snowd / ssoil%ssdn(:,3)
       ssoil%smass(:,1) = ssoil%snowd
       ssoil%smass(:,2) = 0.0     ! EK to fix -ve sdepth 21Dec2007
       ssoil%smass(:,3) = 0.0     ! EK to fix -ve sdepth 21Dec2007
       !      ssoil%smass(:,2) = ssoil%snowd
       !      ssoil%smass(:,3) = ssoil%snowd
    ELSEWHERE (ssoil%snowd < snmin * ssoil%ssdnn)
       WHERE (ssoil%isflag == 1)
          ssoil%ssdn(:,1) = ssoil%ssdnn
          ssoil%tgg(:,1) = ssoil%tggsn(:,1)
       END WHERE
       ssoil%ssdnn = MIN( 400.0, MAX(120.0, ssoil%ssdn(:,1)) )
       ssoil%isflag = 0
       ssoil%tggsn(:,1) = tfrz
       ssoil%tggsn(:,2) = tfrz
       ssoil%tggsn(:,3) = tfrz
       ssoil%sdepth(:,1) = ssoil%snowd / ssoil%ssdn(:,1)
       ssoil%sdepth(:,2) = 0.0     ! EK to fix -ve sdepth 21Dec2007
       ssoil%sdepth(:,3) = 0.0     ! EK to fix -ve sdepth 21Dec2007
       ssoil%smass(:,1) = ssoil%snowd     ! EK to fix -ve sdepth 21Dec2007
       ssoil%smass(:,2) = 0.0     ! EK to fix -ve sdepth 21Dec2007
       ssoil%smass(:,3) = 0.0     ! EK to fix -ve sdepth 21Dec2007
    ELSEWHERE ! sufficient snow now
       WHERE (ssoil%isflag == 0)
          ssoil%tggsn(:,1) = ssoil%tgg(:,1)
          ssoil%tggsn(:,2) = ssoil%tgg(:,1)
          ssoil%tggsn(:,3) = ssoil%tgg(:,1)
          ssoil%ssdn(:,2) = ssoil%ssdn(:,1)
          ssoil%ssdn(:,3) = ssoil%ssdn(:,1)
          ssoil%sdepth(:,1) = 0.05
          ! next 5 lines replaced to fix -ve sdepth (EK 21Dec2007)
          ssoil%smass(:,1)  = 0.05 * ssoil%ssdn(:,1)
          ssoil%smass(:,2)  = (ssoil%snowd - ssoil%smass(:,1)) * 0.4
          ssoil%sdepth(:,2) = ssoil%smass(:,2)/ssoil%ssdn(:,2)
          ssoil%smass(:,3)  = (ssoil%snowd - ssoil%smass(:,1)) * 0.6
          ssoil%sdepth(:,3) = ssoil%smass(:,3)/ssoil%ssdn(:,3)
          !        ssoil%sdepth(:,2) =  (ssoil%snowd / ssoil%ssdn(:,1) - 0.05) &
          !             & * merge(0.4, 0.4, ssoil%snowd > 20.0)
          !        ssoil%sdepth(:,3) = (ssoil%snowd / ssoil%ssdn(:,1) - 0.05) &
          !             & * merge(0.6, 0.6, ssoil%snowd > 20.0)
          !        ssoil%smass(:,1) = 0.05 * ssoil%ssdn(:,1)
          !        ssoil%smass(:,2) = ssoil%sdepth(:,2) * ssoil%ssdn(:,2)
          !        ssoil%smass(:,3) = ssoil%sdepth(:,3) * ssoil%ssdn(:,3)
       END WHERE
       ssoil%isflag = 1
    END WHERE
  
  END SUBROUTINE snowcheck 
  !******************************************************************
  SUBROUTINE snowl_adjust(dels, ktau,  ssoil, canopy )
    REAL(r_1), INTENT(IN)               :: dels ! integration time step (s)
    INTEGER(i_d), INTENT(IN)            :: ktau ! integration step number
    TYPE(soil_snow_type), INTENT(INOUT) :: ssoil
    TYPE(canopy_type), INTENT(INOUT)    :: canopy
    INTEGER(i_d), PARAMETER  :: ntest = 0 !  for prints
    INTEGER(i_d)             :: k
    REAL(r_2), DIMENSION(mp) :: excd
    REAL(r_2), DIMENSION(mp) :: excm
    REAL(r_2), DIMENSION(mp) :: frac 
    REAL(r_2), DIMENSION(mp) :: xfrac 
    REAL(r_1), DIMENSION(mp) :: osm
! AJA ##############################################
    INTEGER(i_d) :: api ! active patch counter
! ##################################################
    !						~.11 to turn on 3-layer snow
    ! Diagnostic block:
    IF (ntest > 0) THEN
      ! Print variables for all active veg/soil patches in idjd grid point:
      PRINT *, 'in soilsnowv before stempv,  ktau= ', &
            & ktau,ssoil%isflag(landpt(idjd)%cstart: &
            &                  (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'soilsnowv ', &
            & canopy%dgdtg(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & canopy%fevc(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & canopy%fes(landpt(idjd)%cstart: &
            &           (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & canopy%precis(landpt(idjd)%cstart: &
            &              (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'ga,dels,ssdn ', &
            & canopy%ga(landpt(idjd)%cstart: &
            &          (landpt(idjd)%cstart+landpt(idjd)%nap-1)), dels, &
            & (ssoil%ssdn(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,3)
      PRINT *, 'osnowd,snowd,isflag', &
            & ssoil%osnowd(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & ssoil%snowd(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & ssoil%isflag(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'tgg ', (ssoil%tgg(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'tggsn ', (ssoil%tggsn(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,3)
      PRINT *, 'wb ', (ssoil%wb(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'wbice ', (ssoil%wbice(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
    END IF

    
    ! adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...
    WHERE (ssoil%isflag > 0)
       WHERE ( ssoil%sdepth(:,1) > 0.05 )
          excd = ssoil%sdepth(:,1) - 0.05
          excm = excd * ssoil%ssdn(:,1)
          ssoil%sdepth(:,1) = ssoil%sdepth(:,1) - REAL(excd,r_1)
          !        ssoil%sdepth(:,1) = 0.05
          osm = ssoil%smass(:,1)
          ssoil%smass(:,1) = ssoil%smass(:,1) - REAL(excm,r_1)
             
          osm = ssoil%smass(:,2)
          ssoil%smass(:,2) = MAX(0.01, ssoil%smass(:,2) + REAL(excm,r_1))
          ssoil%ssdn(:,2) = REAL(MAX(120.0_r_2, MIN(500.0_r_2, ssoil%ssdn(:,2) &
                        & * osm/ssoil%smass(:,2) + ssoil%ssdn(:,1) * excm &
                        & / ssoil%smass(:,2))),r_1)
          ssoil%sdepth(:,2) =  ssoil%smass(:,2) / ssoil%ssdn(:,2) 
          ssoil%tggsn(:,2) = REAL(ssoil%tggsn(:,2) * osm / ssoil%smass(:,2) &
                         & + ssoil%tggsn(:,1) * excm/ ssoil%smass(:,2),r_1)
          ! following line changed to fix -ve sdepth (EK 21Dec2007)
          ssoil%smass(:,3) = MAX(0.01, &
                         & ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2))
          ! ssoil%smass(:,3) = ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2)
       ELSEWHERE ! ssoil%sdepth(:,1) < 0.05
          ! 1st layer
          excd = 0.05 - ssoil%sdepth(:,1)
          excm = excd * ssoil%ssdn(:,2)
          osm = ssoil%smass(:,1)
          ssoil%smass(:,1) = ssoil%smass(:,1) + REAL(excm,r_1)
          ssoil%sdepth(:,1) = 0.05
          ssoil%ssdn(:,1) = REAL(MAX(120.0_r_2, MIN(500.0_r_2, ssoil%ssdn(:,1) &
                        & * osm/ssoil%smass(:,1) + ssoil%ssdn(:,2) * excm &
                        & / ssoil%smass(:,1))),r_1)
          ssoil%tggsn(:,1) = REAL(ssoil%tggsn(:,1) * osm / ssoil%smass(:,1) &
                         & + ssoil%tggsn(:,2) * excm/ ssoil%smass(:,1),r_1)
          ! 2nd layer
          ssoil%smass(:,2) = MAX(0.01, ssoil%smass(:,2) - REAL(excm,r_1))
          ssoil%sdepth(:,2) = ssoil%smass(:,2)/ssoil%ssdn(:,2)
          ! following line changed to fix -ve sdepth (EK 21Dec2007)
          ssoil%smass(:,3) = MAX(0.01, &
                         & ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2))
          ! ssoil%smass(:,3) = ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2)
        END WHERE
   ! AJA REPLACING THE WHERE LOOP #########################  
      END WHERE 
    ! AJA END WHERE

    DO  api=1,mp
        IF(ssoil%isflag(api).gt.0)THEN
          frac(api) = ssoil%smass(api,2) / MAX(0.02, ssoil%smass(api,3))
          ! if frac > 0.6 or frac < 0.74 do nothing HOW TO translate this to xfrac
          xfrac(api) = 2.0/3.0/ frac(api)
          IF(xfrac(api) > 1.0 )THEN
           excm(api) = (xfrac(api) - 1.0) * ssoil%smass(api,2)
           osm(api) = ssoil%smass(api,2)
           ! changed 0.02 to 0.01 to fix -ve sdepth (EK 21Dec2007)
           ssoil%smass(api,2) = MAX(0.01, ssoil%smass(api,2) + REAL(excm(api),r_1))
           ssoil%tggsn(api,2) = ssoil%tggsn(api,2) * osm(api) / ssoil%smass(api,2) +  &
                         & ssoil%tggsn(api,3) * REAL(excm(api),r_1)/ ssoil%smass(api,2)
           ssoil%ssdn(api,2) = MAX(120.0, MIN(500.0, ssoil%ssdn(api,2)* &
              & osm(api)/ssoil%smass(api,2) + ssoil%ssdn(api,3) &
              * REAL(excm(api),r_1) /ssoil%smass(api,2)) )
           ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
           ssoil%smass(api,3) = MAX(0.01, &
                         & ssoil%snowd(api) - ssoil%smass(api,1) - ssoil%smass(api,2))
           ssoil%sdepth(api,3) = MAX(0.02, ssoil%smass(api,3) / ssoil%ssdn(api,3) )
         ELSE! xfrac < 1
           excm(api) = (1 - xfrac(api)) * ssoil%smass(api,2)
           ssoil%smass(api,2) = MAX(0.01, ssoil%smass(api,2) - REAL(excm(api),r_1))
           ssoil%sdepth(api,2) = MAX(0.02, ssoil%smass(api,2) / ssoil%ssdn(api,2) )

           osm(api) = ssoil%smass(api,3)
           ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
           ssoil%smass(api,3) = MAX(0.0, &
                         & ssoil%snowd(api) - ssoil%smass(api,1) - ssoil%smass(api,2))
           ssoil%tggsn(api,3) = ssoil%tggsn(api,3) * osm(api) / ssoil%smass(api,3) +  &
                         & ssoil%tggsn(api,2) * REAL(excm(api),r_1)/ ssoil%smass(api,3)
           ssoil%ssdn(api,3) = MAX(120.0, MIN(500.0, ssoil%ssdn(api,3)* &
                osm(api)/ssoil%smass(api,3) + ssoil%ssdn(api,2) &
                * REAL(excm(api),r_1) / ssoil%smass(api,3)) )
           ssoil%sdepth(api,3) = ssoil%smass(api,3) /  ssoil%ssdn(api,3)
         END IF
           ssoil%isflag(api) = 1
           ssoil%ssdnn(api) = (ssoil%ssdn(api,1) * ssoil%sdepth(api,1) + ssoil%ssdn(api,2) &
             & * ssoil%sdepth(api,2) + ssoil%ssdn(api,3) * ssoil%sdepth(api,3) ) &
             & / (ssoil%sdepth(api,1) + ssoil%sdepth(api,2) + ssoil%sdepth(api,3))
         END IF
   END DO
     


! AJA COMMENTED FROM THE ORIGINAL   #########################################################
!        frac = ssoil%smass(:,2) / MAX(0.02, ssoil%smass(:,3))
!        ! if frac > 0.6 or frac < 0.74 do nothing HOW TO translate this to xfrac
!        xfrac = 2.0/3.0/ frac
!        WHERE ( xfrac > 1.0 )
!          excm = (xfrac - 1.0) * ssoil%smass(:,2)
!          osm = ssoil%smass(:,2)
!          ! changed 0.02 to 0.01 to fix -ve sdepth (EK 21Dec2007)
!          ssoil%smass(:,2) = MAX(0.01, ssoil%smass(:,2) + REAL(excm,r_1))
!          ssoil%tggsn(:,2) = REAL(ssoil%tggsn(:,2) * osm / ssoil%smass(:,2) &
!                         & + ssoil%tggsn(:,3) * excm/ ssoil%smass(:,2),r_1)
!          ssoil%ssdn(:,2) = REAL(MAX(120.0_r_2, MIN(500.0_r_2, ssoil%ssdn(:,2) &
!                        & * osm/ssoil%smass(:,2) + ssoil%ssdn(:,3) * excm &
!                        & / ssoil%smass(:,2))),r_1)
!          ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
!          ssoil%smass(:,3) = MAX(0.01, &
!                         & ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2))
!          ssoil%sdepth(:,3) = MAX(0.02, ssoil%smass(:,3) / ssoil%ssdn(:,3) )
!        ELSEWHERE ! xfrac < 1
!          excm = (1 - xfrac) * ssoil%smass(:,2)
!          ssoil%smass(:,2) = MAX(0.01, ssoil%smass(:,2) - REAL(excm,r_1))
!          ssoil%sdepth(:,2) = MAX(0.02, ssoil%smass(:,2) / ssoil%ssdn(:,2) )
!          osm = ssoil%smass(:,3)
!          ! following line added MAX function to fix -ve sdepth (EK 21Dec2007)
!          ssoil%smass(:,3) = MAX(0.0, &
!                         & ssoil%snowd - ssoil%smass(:,1) - ssoil%smass(:,2))
!          ssoil%tggsn(:,3) = REAL(ssoil%tggsn(:,3) * osm / ssoil%smass(:,3) &
!                         & + ssoil%tggsn(:,2) * excm/ ssoil%smass(:,3),r_1)
!          ssoil%ssdn(:,3) = REAL(MAX(120.0_r_2, MIN(500.0_r_2, ssoil%ssdn(:,3) &
!                        & * osm/ssoil%smass(:,3) + ssoil%ssdn(:,2) * excm &
!                        & / ssoil%smass(:,3)) ),r_1)
!          ssoil%sdepth(:,3) = ssoil%smass(:,3) /  ssoil%ssdn(:,3)
!        END WHERE
!        ssoil%isflag = 1
!        ssoil%ssdnn = (ssoil%ssdn(:,1) * ssoil%sdepth(:,1) + ssoil%ssdn(:,2) &
!               * ssoil%sdepth(:,2) + ssoil%ssdn(:,3) * ssoil%sdepth(:,3) ) &
!               / (ssoil%sdepth(:,1) + ssoil%sdepth(:,2) + ssoil%sdepth(:,3))
!      END WHERE
!    END WHERE

    ! Diagnostic block:
    IF (ntest > 0) THEN
      ! Print variables for all active veg/soil patches in idjd grid point:
      PRINT *, 'in soilsnowv before stempv,  ktau= ',ktau
      PRINT *, 'ga,dels,ssdn ', &
            & canopy%ga(landpt(idjd)%cstart: &
            &          (landpt(idjd)%cstart+landpt(idjd)%nap-1)), dels, &
            & (ssoil%ssdn(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,3)
      PRINT *, 'osnowd,snowd,isflag', &
            & ssoil%osnowd(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & ssoil%snowd(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & ssoil%isflag(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'tgg ', (ssoil%tgg(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'tggsn ', (ssoil%tggsn(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,3)
      PRINT *, 'wb ', (ssoil%wb(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'wbice ', (ssoil%wbice(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'wblf ', (ssoil%wblf(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'wbfice ', (ssoil%wbfice(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
    END IF
    IF (ntest == 1) THEN
      PRINT *, 'in soilsnow printing wbfice_max'
      PRINT *, 'sdepth c2 ', (ssoil%sdepth(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,3)
    END IF

  END SUBROUTINE snowl_adjust


  !******************************************************************
  SUBROUTINE soilfreeze(dels, soil, ssoil)
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    REAL(r_2), DIMENSION(mp)           :: sicefreeze
    REAL(r_2), DIMENSION(mp)           :: sicemelt
    REAL(r_1), DIMENSION(mp)           :: xx
    INTEGER(i_d) k
    
    ! Always allow for some water (< 2%) to remain unfrozen
    DO k = 1, ms
       WHERE (ssoil%tgg(:,k) < tfrz &
            & .AND. 0.98 * ssoil%wb(:,k) - ssoil%wbice(:,k) > .001)
          sicefreeze = MIN( MAX(0.0_r_2,(0.98*ssoil%wb(:,k)-ssoil%wbice(:,k))) &
               & * soil%zse(k) * 1000.0, &
               & (tfrz - ssoil%tgg(:,k) ) * ssoil%gammzz(:,k) / hlf )
          ssoil%wbice(:,k) = MIN( ssoil%wbice(:,k) + sicefreeze / (soil%zse(k) &
               * 1000.0), 0.98 * ssoil%wb(:,k) )
          xx=soil%css * soil%rhosoil
          ssoil%gammzz(:,k) = MAX( &
               REAL((1.0 - soil%ssat) * soil%css * soil%rhosoil ,r_2) &
               + (ssoil%wb(:,k) - ssoil%wbice(:,k)) * REAL(cswat * rhowat,r_2) &
               + ssoil%wbice(:,k) * REAL(csice * rhowat * 0.9,r_2), REAL(xx,r_2)) &
               * REAL(soil%zse(k),r_2)
          WHERE (k == 1 .AND. ssoil%isflag == 0)
             ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
          END WHERE
          ssoil%tgg(:,k) = ssoil%tgg(:,k) + REAL(sicefreeze,r_1) &
               * hlf / REAL(ssoil%gammzz(:,k),r_1)
       ELSEWHERE (ssoil%tgg(:,k) > tfrz .AND. ssoil%wbice(:,k) > 0.)
          sicemelt = MIN(ssoil%wbice(:,k) * soil%zse(k) * 1000.0, &
               & (ssoil%tgg(:,k) - tfrz) * ssoil%gammzz(:,k) / hlf)
          ssoil%wbice(:,k) = MAX(0.0_r_2, ssoil%wbice(:,k) - sicemelt &
               / (soil%zse(k) * 1000.0) )
          xx = soil%css * soil%rhosoil
          ssoil%gammzz(:,k) = MAX( &
               REAL((1.0-soil%ssat) * soil%css * soil%rhosoil,r_2) &
               + (ssoil%wb(:,k) - ssoil%wbice(:,k)) * REAL(cswat*rhowat,r_2) &
               + ssoil%wbice(:,k) * REAL(csice * rhowat * 0.9,r_2), &
               REAL(xx,r_2) ) * REAL(soil%zse(k),r_2)
          WHERE (k == 1 .AND. ssoil%isflag == 0)
             ssoil%gammzz(:,k) = ssoil%gammzz(:,k) + cgsnow * ssoil%snowd
          END WHERE
          ssoil%tgg(:,k) = ssoil%tgg(:,k) - REAL(sicemelt,r_1) &
               * hlf / REAL(ssoil%gammzz(:,k),r_1)
       END WHERE
    END DO
  END SUBROUTINE soilfreeze
  
  !******************************************************************
  SUBROUTINE remove_trans(dels, soil, ssoil, canopy, veg)
    ! Removes transpiration water from soil.
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    REAL(r_2), DIMENSION(mp,ms)   :: evapfbl
    REAL(r_2), DIMENSION(mp,0:ms) :: diff 
    REAL(r_2), DIMENSION(mp)      :: xx
    INTEGER(i_d) k

! ypw (26oct2010) replacing the whole subroutine
!                 as canopy%evapfbl has already been calculated in define_canopy
    DO k = 1,ms
      WHERE (canopy%fevc > 0.0)
        ssoil%wb(:,k) = ssoil%wb(:,k)-canopy%evapfbl(:,k)/(soil%zse(k)*1000.0)
      END WHERE
    END DO

!    diff(:,:) = 0.0  ! (Bug fixed, BP may08)
!!    diff(:,0) = 0.0
!    evapfbl(:,:) = 0.0
!    DO k = 1,ms
!       ! Removing transpiration from soil:
!       WHERE (canopy%fevc > 0.0)     ! convert to mm/dels
!          ! Calculate the amount (perhaps moisture/ice limited)
!          ! which can be removed:
!          xx = canopy%fevc * dels / hl * veg%froot(:,k) + diff(:,k-1)
!          diff(:,k) = (MAX( 0.0_r_2, MIN( ssoil%wb(:,k) - &
!               REAL(soil%swilt,r_2), ssoil%wb(:,k) - ssoil%wbice(:,k) ) ) &
!               * REAL(soil%zse(k),r_2) * 1000.0_r_2 - xx)  &
!               / REAL(soil%zse(k) * 1000.0,r_2)
!          WHERE ( diff(:,k) > 0.0 )
!            ssoil%wb(:,k) = ssoil%wb(:,k) - xx / (soil%zse(k) * 1000.0)
!            diff(:,k) = 0.0
!            evapfbl(:,k) = xx / (soil%zse(k) * 1000.0)
!          ELSEWHERE
!            diff(:,k) = xx
!          ENDWHERE
!          !    evapfbl(:,k) = (MIN(canopy%fevc * dels / hl * veg%froot(:,k), &
!          !         MAX(0._r_2, MIN(ssoil%wb(:,k) &
!          !         - soil%swilt,ssoil%wb(:,k)-ssoil%wbice(:,k))) &
!          !         * soil%zse(k) * 1000.)) / (soil%zse(k) * 1000.)
!          !      ! Remove this amount from  the soil:
!          !      ssoil%wb(:,k) = ssoil%wb(:,k) -   evapfbl(:,k)
!        END WHERE
!     END DO
!
!!    PRINT *, 'Before adjusting fevc'
!     ! Adjust fevc
!     WHERE (canopy%fevc > 0.) ! convert to mm/dels
!        canopy%fevc = (evapfbl(:,1)*soil%zse(1)+evapfbl(:,2)*soil%zse(2) &
!            & +evapfbl(:,3)*soil%zse(3)+evapfbl(:,4)*soil%zse(4)+evapfbl(:,5) &
!            & *soil%zse(5)+evapfbl(:,6)*soil%zse(6))*1000.*hl/dels
!     END WHERE

  END SUBROUTINE remove_trans 

  !----------------------------------------------------------------------------
  ! SUBROUTINE soil_snow
  !
  ! Replaces following
  ! SUBROUTINE soilsnow (dt_in, ktau_in, ga, dgdtg, condxpr, scondxpr, fev, &
  !                      fes, t, coszen)
  !	    for snow diag prints set ntest to 1 throughout
  !	    or, usefully can edit 'ntest > 0' to 'ktau > nnn'
  !----------------------------------------------------------------------
  ! Inputs:
  !	 dt_in - time step in sec
  !	 ktau_in - time step no.
  !	 ga	 - ground heat flux W/m^2
  !	 dgdtg	 -
  !	 condxpr - total precip reaching the ground (liquid and solid)
  !	 scondxpr - precip (solid only)
  !	 fev   - transpiration (W/m2)
  !	 fes   - soil evaporation (W/m2)
  !	 isoil - soil type
  !	 ivegt - vegetation type
  ! Output
  !	 ssoil

  SUBROUTINE soil_snow(dels, ktau, soil, ssoil, veg, canopy, met)
    REAL(r_1), INTENT(IN)                    :: dels ! integration time step (s)
    INTEGER(i_d), INTENT(IN)                 :: ktau ! integration step number
    TYPE(soil_parameter_type), INTENT(INOUT) :: soil
    TYPE(soil_snow_type), INTENT(INOUT)      :: ssoil
    TYPE(canopy_type), INTENT(INOUT)         :: canopy
    TYPE(veg_parameter_type), INTENT(INOUT)  :: veg
    TYPE(met_type), INTENT(INOUT)            :: met ! all met forcing
    INTEGER(i_d), PARAMETER  :: ntest = 0 !  for prints
    INTEGER(i_d)             :: k
!    REAL(r_2), DIMENSION(mp):: sicefreeze
!    REAL(r_2), DIMENSION(mp):: sicemelt
    REAL(r_1), DIMENSION(mp) :: snowmlt
    REAL(r_1), DIMENSION(mp) :: totwet
    REAL(r_1), DIMENSION(mp) :: weting
    REAL(r_2), DIMENSION(mp) :: xxx
    REAL(r_2), DIMENSION(mp) :: xx
    ! Diagnostic block:
    IF(ntest>0) THEN
      ! Print variables for all active veg/soil patches in idjd grid point:
      PRINT *, 'in soilsnowv before stempv,  ktau= ', ktau, &
            & ssoil%isflag(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'soilsnowv ', &
            & canopy%dgdtg(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1)),&
            & canopy%fevc(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & canopy%fes(landpt(idjd)%cstart: &
            &           (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & canopy%precis(landpt(idjd)%cstart: &
            &              (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'ga,dels,ssdn ', &
            & canopy%ga(landpt(idjd)%cstart: &
            &          (landpt(idjd)%cstart+landpt(idjd)%nap-1)), dels, &
            & (ssoil%ssdn(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,3)
      PRINT *, 'osnowd,snowd,isflag', &
            & ssoil%osnowd(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & ssoil%snowd(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & ssoil%isflag(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'tgg ', (ssoil%tgg(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'tggsn ', (ssoil%tggsn(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,3)
      PRINT *, 'wb ', (ssoil%wb(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'wbice ', (ssoil%wbice(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
    END IF
    
    ssoil%fwtop = 0.0
    ssoil%runoff = 0.0 ! initialise total runoff
    ssoil%rnof1 = 0.0 ! initialise surface runoff
    ssoil%rnof2 = 0.0 ! initialise deep drainage
    ssoil%smelt = 0.0 ! initialise snowmelt
    ssoil%osnowd = ssoil%snowd

!    IF (ktau <= 1) THEN ! initialisations if first timestep
!      canopy%ga = 0.0
!      ! N.B. snmin should exceed sum of layer depths, i.e. .11 m
!      !     IF (ntest == 3) snmin = .11 ! to force 3-layer snow for testing
!      ssoil%wbtot = 0.0
!      DO k = 1, ms
!          ssoil%wbtot = ssoil%wbtot + REAL(ssoil%wb(:,k),r_1)*1000.0*soil%zse(k)
!          ! Update soil ice:
!          WHERE (ssoil%tgg(:,k) <= tfrz .AND. ssoil%wbice(:,k) <= 0.01)
!            ssoil%wbice(:,k) = 0.1 * ssoil%wb(:,k)
!          END WHERE
!          WHERE (ssoil%tgg(:,k) < tfrz)
!            ssoil%wbice(:,k) = MIN(0.98_r_2 * ssoil%wb(:,k), &
!                 REAL(soil%ssat,r_2))
!          END WHERE
!      END DO
!      xx=soil%css * soil%rhosoil
!      ssoil%gammzz(:,1) = MAX( (1.0 - soil%ssat) * soil%css * soil%rhosoil &
!           + (ssoil%wb(:,1) - ssoil%wbice(:,1) ) * cswat * rhowat &
!           + ssoil%wbice(:,1) * csice * rhowat * .9, xx ) * soil%zse(1)
!    END IF
    
    DO k = 1, ms ! for stempv
       ! Set liquid soil water fraction (fraction of saturation value):
       ssoil%wblf(:,k) = MAX( 0.01_r_2, (ssoil%wb(:,k) - ssoil%wbice(:,k)) ) &
            & / REAL(soil%ssat,r_2)
       ! Set ice soil water fraction (fraction of saturation value):
       ssoil%wbfice(:,k) = REAL(ssoil%wbice(:,k),r_1) / soil%ssat
    END DO

    CALL snowcheck (dels, ktau, ssoil)
    CALL snowdensity (dels, ssoil)
!   PRINT *, 'Before snow_accum'
    CALL snow_accum (dels, ktau, canopy, met, ssoil, soil)
    CALL snow_melting (dels, snowmlt, ktau, ssoil)

    ! Add snow melt to global snow melt variable:
    ssoil%smelt = snowmlt
!   PRINT *, 'Before snowl_adjust'
    ! Adjust levels in the snowpack due to snow accumulation/melting,
    ! snow aging etc...
    CALL snowl_adjust(dels, ktau, ssoil, canopy )

    ! Diagnostic block:
    IF (ntest > 0) THEN
      ! Print variables for all active veg/soil patches in idjd grid point:
      PRINT *, 'in soilsnowv before stempv,  ktau= ',ktau
      PRINT *, 'ga,dels,ssdn ', &
            & canopy%ga(landpt(idjd)%cstart: &
            &          (landpt(idjd)%cstart+landpt(idjd)%nap-1)), dels, &
            & (ssoil%ssdn(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,3)
      PRINT *, 'osnowd,snowd,isflag', &
            & ssoil%osnowd(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & ssoil%snowd(landpt(idjd)%cstart: &
            &            (landpt(idjd)%cstart+landpt(idjd)%nap-1)), &
            & ssoil%isflag(landpt(idjd)%cstart: &
            &             (landpt(idjd)%cstart+landpt(idjd)%nap-1))
      PRINT *, 'tgg ', (ssoil%tgg(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'tggsn ', (ssoil%tggsn(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,3)
      PRINT *, 'wb ', (ssoil%wb(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'wbice ', (ssoil%wbice(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'wblf ', (ssoil%wblf(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
      PRINT *, 'wbfice ', (ssoil%wbfice(landpt(idjd)%cstart: &
            & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
    END IF

!   PRINT *, 'Before stempv'
    CALL stempv(dels, canopy, ssoil, soil)
    CALL snow_melting (dels, snowmlt, ktau, ssoil )

    ! Add new snow melt to global snow melt variable: 
    ssoil%smelt = ssoil%smelt + snowmlt    
    ! Define total water into soil (throughfall+snowmelt):
    totwet = canopy%precis + ssoil%smelt
    ! Move any excessive precipitation/melting (> 300mm/day)
    ! to surface runoff:
    ssoil%rnof1 = MAX(0., totwet - 300.0*dels / 86400.0)
    ! Define water available to go into soil:
    weting = totwet - ssoil%rnof1
    ! Define available pore space in top soil layer:
    xxx=soil%ssat - ssoil%wb(:,1)
    ! Add to runoff any additional available soil water that doesn't 
    ! fit into top layer:
    ssoil%rnof1 = ssoil%rnof1 + MAX(0.0, REAL(weting - &
         0.99 * MIN( REAL(xxx,r_1) * soil%zse(1) * rhowat, weting),r_1))
    ! Redefine water into soil to reflect this:
    weting = totwet - ssoil%rnof1
    ! Define flux of water into top soil layer (weting - evaporation):
    ssoil%fwtop = weting / dels - canopy%segg

    CALL soilfreeze(dels, soil, ssoil)

!   PRINT *, 'Before remove_trans'
    ! Take transpiration water from soil:
    CALL remove_trans(dels, soil, ssoil, canopy, veg)

!   PRINT *, 'Before surfbv'
    CALL surfbv(dels, ktau,  met, ssoil, soil )

!    WRITE(47,*) 'albsoilsn Before snow_albedo'
!    WRITE(47,'(i5,200f5.2)') ktau, ssoil%albsoilsn(:,1)
!    WRITE(47,'(i5,200f5.2)') ktau, ssoil%albsoilsn(:,2)

!   PRINT *, 'Before snow_albedo'
    CALL snow_albedo(dels, ktau, met, ssoil, soil )

!    WRITE(47,*) 'albsoilsn After snow_albedo'
!    WRITE(47,'(i5,200f5.2)') ktau, ssoil%albsoilsn(:,1)
!    WRITE(47,'(i5,200f5.2)') ktau, ssoil%albsoilsn(:,2)

    ! Diagnostic block:
    IF (ntest > 0) THEN
       ! Print variables for all active veg/soil patches in idjd grid point:
       PRINT *, 'after surfbv,isflag ', ssoil%isflag(landpt(idjd)%cstart: &
             & (landpt(idjd)%cstart+landpt(idjd)%nap-1))
       PRINT *, 'tgg ', (ssoil%tgg(landpt(idjd)%cstart: &
             & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
       PRINT *, 'wb ', (ssoil%wb(landpt(idjd)%cstart: &
             & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
       PRINT *, 'wblf ', (ssoil%wblf(landpt(idjd)%cstart: &
             & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
       PRINT *, 'wbfice ', (ssoil%wbfice(landpt(idjd)%cstart: &
             & (landpt(idjd)%cstart+landpt(idjd)%nap-1),k),k=1,ms)
    END IF

    ! Set weighted soil/snow surface temperature:
    ssoil%tss = (1-ssoil%isflag)*ssoil%tgg(:,1) &
         & + ssoil%isflag*ssoil%tggsn(:,1)
    ssoil%wbtot = 0.0
    DO k = 1, ms
       ssoil%wbtot = ssoil%wbtot + REAL(ssoil%wb(:,k)*1000.0*soil%zse(k),r_1)
    END DO

  END SUBROUTINE soil_snow

END MODULE soil_snow_module
