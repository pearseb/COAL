! cable_canopy.f90
!
! Science development by Ying-Ping Wang, Eva Kowalczyk, Mike Raupach, 
! Ray Leuning et al at CSIRO Marine and Atmospheric Research.
!
! Fortran-95 recoding by Harvey Davies, Gab Abramowitz, Martin Dix & Bernard Pak
!
! Please send bug reports to Bernard.Pak@csiro.au
!
! This file contains the canopy_module, with subroutines define_canopy, dryLeaf,
!                                                    wetLeaf and photosynthesis
! The functions included are:
!   qsatf,
!   ej3x,
!   ej4x,
!   xvcmxt4,
!   xvcmxt3,
!   xejmxt3,
!   psim,
!   psis
!
MODULE canopy_module
  USE define_dimensions, ONLY: r_1,r_2,mp,mf,ms,i_d
  USE define_types
  USE photosynthetic_constants
  USE radiation_module, ONLY: radiation, sinbet
  USE air_module, ONLY: define_air
  USE physical_constants
  IMPLICIT NONE
  REAL(r_1), DIMENSION(:), POINTER :: cansat ! max canopy intercept. (mm)
  REAL(r_2), DIMENSION(:), POINTER :: ghwet  ! cond for heat for a wet canopy
  REAL(r_1), DIMENSION(:), POINTER :: dsx ! leaf surface vpd
  REAL(r_1), DIMENSION(:), POINTER :: fwsoil ! soil water modifier of stom. cond
  REAL(r_1), DIMENSION(:), POINTER :: tlfx ! leaf temp prev. iter (K)
  REAL(r_1), DIMENSION(:), POINTER :: tlfy ! leaf temp (K)
  REAL(r_2), DIMENSION(:), POINTER :: ecy ! lat heat fl dry big leaf
  REAL(r_2), DIMENSION(:), POINTER :: hcy ! veg. sens heat
  REAL(r_2), DIMENSION(:), POINTER :: rny ! net rad
  REAL(r_2), DIMENSION(:,:), POINTER :: gbhu ! forcedConvectionBndryLayerCond
  REAL(r_2), DIMENSION(:,:), POINTER :: gbhf ! freeConvectionBndryLayerCond
                                             ! mol/m2/s
  REAL(r_1), DIMENSION(:,:), POINTER :: gswmin ! min stomatal conductance
  REAL(r_2), DIMENSION(:,:), POINTER :: gswx ! stom cond for water
  REAL(r_2), DIMENSION(:,:), POINTER :: csx ! leaf surface CO2 concentration
  PRIVATE
  PUBLIC define_canopy

  !$OMP THREADPRIVATE(cansat,ghwet,dsx,fwsoil,tlfx,tlfy,ecy,hcy,rny)
  !$OMP THREADPRIVATE(gbhu,gbhf,gswmin,gswx,csx)

CONTAINS
  
  !--------------------------------------------------------------------------
  SUBROUTINE define_canopy(ktau,bal,rad,rough,air,met,dels, &
                           ssoil,soil,veg,bgc,canopy)
    TYPE (balances_type),       INTENT(INOUT) :: bal
    TYPE (radiation_type),      INTENT(INOUT) :: rad
    TYPE (roughness_type),      INTENT(INOUT) :: rough
    TYPE (air_type),            INTENT(INOUT) :: air
    TYPE (met_type),            INTENT(INOUT) :: met
    TYPE (soil_snow_type),      INTENT(INOUT) :: ssoil
    TYPE (bgc_pool_type),       INTENT(IN)    :: bgc
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (veg_parameter_type),  INTENT(INOUT) :: veg
    TYPE (canopy_type),         INTENT(INOUT) :: canopy
    INTEGER(i_d),               INTENT(IN)    :: ktau ! integration step number
    REAL(r_1),                  INTENT(IN)    :: dels ! integration timesetp (s)
    REAL(r_2), DIMENSION(mp) :: gbvtop ! bnd layer cond. top leaf
    REAL(r_1), DIMENSION(mp) :: rt0 ! turbulent resistance
    REAL(r_1), DIMENSION(mp) :: ortsoil ! turbulent resistance, prev time step
    REAL(r_1), DIMENSION(mp) :: rt1usc ! eq. 3.53, SCAM manual, 1997
! ypw testing new formulation (26oct2010), move rwater to canopy_type
    REAL(r_1), DIMENSION(mp) :: rwater ! soil water availability
    REAL(r_1), DIMENSION(mp) :: oldcansto ! prev t step canopy storage
    REAL(r_1), DIMENSION(mp) :: cc ! limitation term for canopy interception per timestep		   
    REAL(r_1), DIMENSION(mp) :: denom ! denominator in calculating screen temperature, humidity etc
    REAL(r_1), DIMENSION(mp) :: tstar ! 
    REAL(r_1), DIMENSION(mp) :: zscrn !
    REAL(r_1), DIMENSION(mp) :: qstar !
    REAL(r_1), DIMENSION(mp) :: rsts  !
    REAL(r_1), DIMENSION(mp) :: qsurf !
    REAL(r_1), DIMENSION(mp) :: qtgnet !
    REAL(r_1), DIMENSION(mp) :: poolcoef1 ! non-leaf carbon turnover rate * non-leaf pool size
    REAL(r_1), DIMENSION(mp) :: poolcoef1w ! wood carbon turnover rate * wood pool size
    REAL(r_1), DIMENSION(mp) :: poolcoef1r ! root carbon turnover rate * root pool size
    REAL(r_1), DIMENSION(mp) :: rbw ! leaf boundary layer resistance for water
    REAL(r_1), DIMENSION(mp) :: rsw ! stomatal resistance for water
    REAL(r_1), DIMENSION(mp) :: dmah ! A_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp) :: dmbh ! B_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp) :: dmch ! C_{H} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp) :: dmae ! A_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp) :: dmbe ! B_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp) :: dmce ! C_{E} in eq. 3.41 in SCAM, CSIRO tech report 132
    REAL(r_1), DIMENSION(mp) :: tss4 ! soil/snow temperature**4
    REAL(r_1), DIMENSION(mp) :: sss ! var for Penman-Monteith soil evap
    REAL(r_1), DIMENSION(mp) :: cc1 ! var for Penman-Monteith soil evap
    REAL(r_1), DIMENSION(mp) :: cc2 ! var for Penman-Monteith soil evap
    REAL(r_1), DIMENSION(mp) :: qstvair ! sat spec humidity at leaf temperature
    REAL(r_1), DIMENSION(mp) :: xx ! delta-type function for sparse canopy limit, p20 SCAM manual
    REAL(r_1), DIMENSION(mp,mf) :: frac42 ! 2D frac4
    REAL(r_1), DIMENSION(mp,niter):: zetar ! stability correction
    INTEGER(i_d) :: iter ! iteration #
    INTEGER(i_d) :: iterplus !

! variables for alternative method to PM
!    REAL(r_1), DIMENSION(mp) :: dq ! sat spec hum diff.
!    REAL(r_1), DIMENSION(mp) :: qstss ! sat spec humidity at soil/snow temperature
! end variables for alternative method to PM

!!%% changes by Ashok Luhar (low wind speed) 
!    REAL(r_1), PARAMETER      :: alpha1=4.0 
!    REAL(r_1), PARAMETER      :: beta1=0.5
!    REAL(r_1), PARAMETER      :: gamma1=0.3
!    REAL(r_1), DIMENSION(mp)  :: zeta1
!    REAL(r_1), DIMENSION(mp)  :: zeta2
!!%%
    REAL(r_1), DIMENSION(mp,3) :: xi
    REAL(r_1), DIMENSION(mp,3) :: ti
    REAL(r_1), DIMENSION(mp,3) :: si

!%% changes by Ian Harman (screen temperature)
    REAL(r_1), DIMENSION(mp)  :: term1
    REAL(r_1), DIMENSION(mp)  :: term2
    REAL(r_1), DIMENSION(mp)  :: term3
    REAL(r_1), DIMENSION(mp)  :: term5
    REAL(r_1), DIMENSION(mp)  :: r_sc
    REAL(r_1), DIMENSION(mp)  :: zscl
!%%

    INTEGER(i_d)   :: ii    ! BP jul2010
    REAL(r_1), DIMENSION(mp)  :: tmp1  ! BP jul2010
    REAL(r_1), DIMENSION(mp)  :: tmp2  ! BP jul2010
    REAL(r_1), DIMENSION(mp)  :: tmp3  ! BP jul2010
    REAL(r_2), DIMENSION(mp)  :: tmp8  ! BP jul2010

! ypw testing new formulation (26oct2010)
    INTEGER(i_d)   :: ns
!!!!    REAL(r_1), parameter ::rootgamma = 0.1
!!!!    REAL(r_1), parameter ::rootgamma = 0.001
!!!!    REAL(r_1), parameter ::rootgamma = 1.0e-5
    REAL(r_1), parameter ::rootgamma = 0.01
    REAL(r_1), DIMENSION(mp)  :: dummy, normFac
    REAL(r_1), DIMENSION(mp,ms):: alpha1,alpha2
! ypw (26oct2010)

    ALLOCATE(cansat(mp),ghwet(mp),  dsx(mp))
    ALLOCATE(fwsoil(mp), tlfx(mp), tlfy(mp))
    ALLOCATE(   ecy(mp),  hcy(mp),  rny(mp))
    ALLOCATE(  gbhu(mp,mf), gbhf(mp,mf))
    ALLOCATE(gswmin(mp,mf), gswx(mp,mf), csx(mp,mf))

    ! 1-oct-2002 ypw: to keep the unit consistent for resistance or conductance
    ! s/m for r; mol/m2/s for g, and always use g where appropriate
    ! replace rh, rhr, rw  with ghdry/ghwet,ghrdry, ghrwet, gwdry, gwwet

    ! Set surface water vapour pressure deficit:
    met%da = (qsatf(met%tk-tfrz,met%pmb) - met%qv) *rmair/rmh2o *met%pmb *100.0

! ypw testing new formulation (26oct2010)
    ! Soil water limitation on stomatal conductance:
!    rwater = MAX(1.0e-4, &
!         SUM(veg%froot * MIN(1.0,REAL(ssoil%wb,r_1) - &
!         SPREAD(soil%swilt, 2, ms)),2) / (soil%sfc-soil%swilt))
    ! construct function to limit stom cond for soil water availability
!!!    ! Huqiang suggested removing the upper limit (BP may2010)
!!!    fwsoil = MAX(1.0e-4, veg%vbeta * rwater)
!    fwsoil = MAX(1.0e-4,MIN(1.0, veg%vbeta * rwater))

    ! ypw 19/may/2010 soil water uptake efficiency (see Lai and Ktaul 2000)
!    fwsoil(:) = 1.0e-4
    normFac(:) = 0.0
    do ns=1,ms
      dummy(:) = rootgamma/max(1.0e-3,ssoil%wb(:,ns)-soil%swilt(:))
!  WRITE(46,*) 'ns, dummy = ', ns, dummy
!  WRITE(46,*) 'ssoil%wb(:,ns) = ', ssoil%wb(:,ns)
!  WRITE(46,*) 'soil%swilt(:)  = ', soil%swilt(:)
!  WRITE(46,*) 'soil%ssat(:)   = ', soil%ssat(:)
      alpha1(:,ns) = MAX(1.0e-4,ssoil%wb(:,ns)/(soil%ssat(:)-soil%swilt(:)))
      alpha2(:,ns) =MAX(1.0e-4,((ssoil%wb(:,ns)-soil%swilt(:))/soil%ssat(:)) &
                    ** dummy)
      normFac(:) = normFac(:) + alpha2(:,ns) * alpha1(:,ns)* veg%froot(:,ns)
    enddo
    do ns=1,ms
      canopy%rwater(:,ns) = alpha2(:,ns)*alpha1(:,ns)*veg%froot(:,ns)/normFac(:)
!  WRITE(46,*) 'rwater = ', canopy%rwater(:,ns)
!      fwsoil(:) = min(1.0,max(fwsoil(:),alpha2(:,ns)))
    enddo
! ypw (26oct2010)

    ! Soil water limitation on stomatal conductance:
    ! replace linear approximation to fwsoil (and rwater) by polynomial fitting    EAK 3/08/10
    !-------------------------------------------------------------------------
    rwater = MAX(1.0e-4_r_2, &
         SUM(veg%froot * MAX(0.024,MIN(1.0_r_2,ssoil%wb - SPREAD(soil%swilt, 2, ms))),2) &
         /(soil%sfc-soil%swilt))
    rwater = soil%swilt + rwater * (soil%sfc-soil%swilt)
    ! construct function to limit stom cond for soil water availability
!    fwsoil = MAX(1.0e-4,MIN(1.0, veg%vbeta * rwater))
      xi(:,1) = soil%swilt
      xi(:,2) = soil%swilt + (soil%sfc - soil%swilt)/2.0
      xi(:,3) = soil%sfc
      ti(:,1) = 0.
      ti(:,2) = 0.9
      ti(:,3) = 1.0
      si(:,1) = (rwater - xi(:,2)) / ( xi(:,1) - xi(:,2)) * ( rwater  - xi(:,3)) / (xi(:,1) - xi(:,3))
      si(:,2) = (rwater - xi(:,1)) / ( xi(:,2) - xi(:,1)) * ( rwater  - xi(:,3)) / (xi(:,2) - xi(:,3))
      si(:,3) = (rwater - xi(:,1)) / ( xi(:,3) - xi(:,1)) * ( rwater  - xi(:,2)) / (xi(:,3) - xi(:,2))
      fwsoil = 1.
      where (rwater < soil%sfc - 0.02)
        fwsoil = max(0.,min(1., ti(:,1)*si(:,1) + ti(:,2)*si(:,2) + ti(:,3)*si(:,3)))
      endwhere
    !-------------------------------------------------------------------------

    ! BATS-type canopy saturation proportional to LAI:
    cansat = veg%canst1 * canopy%vlaiw
    ! Set previous time step canopy water storage:
    oldcansto=canopy%cansto
    ! Rainfall variable is limited so canopy interception is limited,
    ! used to stabilise latent fluxes.
    ! to avoid excessive direct canopy evaporation (EK nov2007, snow scheme)
    !      cc =MIN(met%precip-met%precip_s, 4./(1440./(dels/60.)))
    ! modified further by timestep requirement (EAK aug08)
    cc =MIN(met%precip-met%precip_s, 4.0 * MIN(dels,1800.0) / 60.0 /1440.0 )
    !    ! to avoid canopy temp. oscillations
    !    cc =MIN(met%precip, 4./(1440./(dels/60.)))
    ! Calculate canopy intercepted rainfall, equal to zero if temp < 0C:
    canopy%wcint = MERGE(MIN(MAX(cansat - canopy%cansto,0.0), cc), 0.0, &
         cc > 0.0  )  ! EK nov2007, snow scheme
    !         cc > 0.0  .AND. met%tk > tfrz)
    ! Define canopy throughfall (100% of precip if temp < 0C, see above):
    canopy%through = met%precip_s + MIN( met%precip - met%precip_s , &
         MAX(0.0, met%precip - met%precip_s - canopy%wcint) )  ! EK nov2007
    ! Delete line below in case it affects snow sites (precip_s) (EK Jul08)
    !      canopy%through = MIN(met%precip,MAX(0.0, met%precip - canopy%wcint))
    ! Add canopy interception to canopy storage term:
    canopy%cansto = canopy%cansto + canopy%wcint
    ! Calculate fraction of canopy which is wet:
    canopy%fwet   = MAX(0.0,MIN(1.0,0.8*canopy%cansto/MAX(cansat,0.01)))

!    ssoil%wetfac = 0.4 * MAX(0.0, MIN(1.0, &
!         (REAL(ssoil%wb(:,1),r_1) - soil%swilt) / (soil%sfc - soil%swilt)))
    ssoil%wetfac = MAX(0.0, MIN(1.0, &
         (REAL(ssoil%wb(:,1),r_1) - soil%swilt) / (soil%sfc - soil%swilt)))

    ! owetfac introduced to reduce sharp changes in dry regions,
    ! especially in offline runs where there may be discrepancies between
    ! timing of precip and temperature change (EAK apr2009)
    ssoil%wetfac = 0.5*(ssoil%wetfac + ssoil%owetfac)
    ! Temporay fixer for accounting the reduction of soil evap due to freezing
    WHERE ( ssoil%wbice(:,1) > 0.0 ) ! Prevents divide by zero at glaciated
                                    ! points where wb and wbice=0.
      ssoil%wetfac = ssoil%wetfac &
                   * REAL(((1.0-ssoil%wbice(:,1)/ssoil%wb(:,1)))**2,r_1)
    END WHERE
    zetar(:,:) = zeta0 ! initialize all (BP jul2010)
!    zetar(:,1) = zeta0 ! stability correction terms
    zetar(:,2) = zetpos + 1
    frac42 = SPREAD(veg%frac4, 2, mf) ! frac C4 plants
    gswx = 1e-3     ! default stomatal conuctance 
    gbhf = 1e-3     ! default free convection boundary layer conductance
    gbhu = 1e-3     ! default forced convection boundary layer conductance
    ghwet = 1.0e-3
    ! Initialise in-canopy temperatures and humidity:
    csx = SPREAD(met%ca, 2, mf) ! initialise leaf surface CO2 concentration
    met%tvair = met%tk
    met%tvrad = met%tk
    met%qvair = met%qv
    CALL define_air(met, air)    ! function of tvair & pmb to update dsatdk & psyc
    qstvair = qsatf((met%tvair-tfrz),met%pmb)
    met%dva = (qstvair - met%qvair) *  rmair/rmh2o * met%pmb * 100.0
    dsx = met%dva     ! init. leaf surface vpd
    tlfx = met%tvair  ! initialise leaf temp iteration memory variable (K)
    tlfy = met%tvair  ! initialise current leaf temp (K)
    ortsoil = ssoil%rtsoil
    ssoil%tss = (1-ssoil%isflag)*ssoil%tgg(:,1)+ssoil%isflag*ssoil%tggsn(:,1)
    tss4 = ssoil%tss**4

    DO iter = 1, niter
       CALL radiation(bal, soil, ssoil, veg, air, met, rad, canopy)
       gswmin = rad%scalex * (gsw03 * (1.0 - frac42) + gsw04 * frac42)
       ! AERODYNAMIC PROPERTIES: friction velocity us, thence turbulent
       ! resistances rt0, rt1 (elements of dispersion matrix):
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.46:
       canopy%us = MAX( 1.0e-6, vonk * MAX(met%ua,umin) &
            / ( LOG(rough%zref_uv/rough%z0m) - psim(zetar(:,iter)) &
            + psim(zetar(:,iter)*rough%z0m/rough%zref_uv) ) )
!!%%change by Ashok Luhar - low wind formulation
!       WHERE (zetar(:,iter) > 0.7)
!         zeta1=zetar(:,iter) * rough%z0m / rough%zref_uv
!         canopy%us = MAX(1.0e-6, &
!            vonk * MAX(met%ua,umin) / ( &
!            alpha1* ((zetar(:,iter)**beta1*  &
!               (1.0+gamma1*zetar(:,iter)**(1.0-beta1)))  &
!             - (zeta1**beta1*(1.0+gamma1*zeta1**(1.0-beta1))))))
!       ENDWHERE
!!%%
       ! Turbulent aerodynamic resistance from roughness sublayer depth to
       ! reference height,
       ! x=1 if zref+disp>zruffs, 0 otherwise:
       ! thus rt1usc = 0 if zref+disp<zruffs
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.53:
       xx = 0.5 + SIGN(0.5,rough%zref_tq+rough%disp-rough%zruffs)
       !              correction  by Ian Harman to the 2nd psis term
       rt1usc = xx * ( LOG( rough%zref_tq/MAX(rough%zruffs-rough%disp, &
            rough%z0soilsn) ) - psis(zetar(:,iter)) &
            + psis( zetar(:,iter)*MAX(rough%zruffs-rough%disp, &
            rough%z0soilsn)/rough%zref_tq ) )/vonk
       ! rt0 = turbulent resistance from soil to canopy:
!!$     ! correction  by Ian Harman to rough%rt0us = f( zetar )
!!$     WHERE (canopy%vlaiw.LT.0.01 .OR. rough%hruff.LT. rough%z0soilsn)
!!$       rough%rt0us  = 0.0
!!$       rt0old  = 0.0
!!$     ELSEWHERE
!!$!      rough%term6=exp(3.*rough%coexp*(rough%disp/rough%hruff -1.))
!!$       rt0old  = rough%term5*(zdlin*LOG(zdlin*rough%disp/rough%z0soilsn) &
!!$            + (1-zdlin))*(EXP(2*csw*canopy%vlaiw) - rough%term2)/rough%term3
!!$       rough%rt0us=rough%term5*(zdlin*LOG(zdlin*rough%disp/rough%z0soilsn) &
!!$!          - psis( zetar(:,iter) * rough%disp/rough%zref/rough%term6)  &
!!$!          + psis( zetar(:,iter) * rough%z0soilsn/rough%zref/rough%term6) &
!!$           + (1-zdlin))*(EXP(2*csw*canopy%vlaiw) - rough%term2)/rough%term3 &
!!$              / rough%term6
!!$     ENDWHERE
!!$     rt0old = rt0old / canopy%us
!!$     rt0 = MAX(5.,rough%rt0us / canopy%us)
       rt0 = rough%rt0us / canopy%us
!       rt0 = MIN(500.0, rt0)      ! for testing (BP aug2010)
       ! Aerodynamic resistance (sum 3 height integrals)/us
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.50:
       rough%rt1 = MAX(5.,(rough%rt1usa + rough%rt1usb + rt1usc) / canopy%us)
       WHERE (ssoil%snowd > 0.1)
          ssoil%wetfac = 1.0
       END WHERE
       ! change canopy%vlaiw requirement to 0.01 for conformity (BP may 2009)
!       ssoil%rtsoil = rt0 + rough%rt1*(0.5 + SIGN(0.5,0.02-canopy%vlaiw)) 
       WHERE (canopy%vlaiw > 0.01)
         ssoil%rtsoil = rt0
       ELSEWHERE
         ssoil%rtsoil = rt0 + rough%rt1
       END WHERE
       ssoil%rtsoil = MAX(25.,ssoil%rtsoil)   
       WHERE (ssoil%rtsoil > 2.*ortsoil .OR. ssoil%rtsoil < 0.5*ortsoil)
          ssoil%rtsoil = MAX(25.,0.5*(ssoil%rtsoil + ortsoil))
       END WHERE
       ! Vegetation boundary-layer conductance (mol/m2/s)
       ! prandt = kinematic viscosity/molecular diffusivity
       ! See CSIRO SCAM, Raupach et al 1997, eq. 3.12. Top leaf:
       gbvtop = air%cmolar*apol * air%visc / prandt / veg%dleaf &
            * (canopy%us / MIN(rough%usuh, 0.2) &
            * veg%dleaf / air%visc)**0.5 * prandt**(1.0/3.0) / veg%shelrb
!       gbvtop = MAX(0.05,gbvtop)      ! for testing (BP aug2010)
       ! Forced convection boundary layer conductance
       ! (see Wang & Leuning 1998, AFM):
       ! gbhu corrected by F.Cruz & A.Pitman on 13/03/07
       tmp8(:) = -canopy%vlaiw*(0.5*rough%coexp+rad%extkb)
       gbhu(:,1) = gbvtop*(1.0-EXP(tmp8)) / (rad%extkb+0.5*rough%coexp)
       tmp8(:) = -0.5*rough%coexp*canopy%vlaiw
       gbhu(:,2) = (2.0/rough%coexp)*gbvtop*(1.0-EXP(tmp8))-gbhu(:,1)
!       gbhu(:,1) = gbvtop*(1.0-EXP(-canopy%vlaiw*(0.5*rough%coexp+rad%extkb)))&
!            / (rad%extkb+0.5*rough%coexp)
!       gbhu(:,2) = (2.0/rough%coexp)*gbvtop*  &
!            (1.0 - EXP(-0.5*rough%coexp*canopy%vlaiw))-gbhu(:,1)

       hcy = 0.0              ! init current estimate lat heat
       rny = SUM(rad%rniso,2) ! init current estimate net rad
       ecy = rny - hcy        ! init current estimate lat heat
       
       CALL dryLeaf(dels,rad,air,met,ssoil,soil,veg,canopy,ktau)

       CALL wetLeaf(dels,rad,air,met,canopy)

      ! Calculate latent heat from vegetation:
      canopy%fev = REAL(canopy%fevc + canopy%fevw,r_1)
      ! Calculate sensible heat from vegetation:
      canopy%fhv = (1.0 - canopy%fwet) * REAL(hcy, r_1) + REAL(canopy%fhvw, r_1)
      ! Calculate net rad absorbed by canopy:
      canopy%fnv = (1.0 - canopy%fwet) * REAL(rny, r_1) &
           + REAL(canopy%fevw + canopy%fhvw , r_1)
      ! canopy radiative temperature is calculated based on long-wave
      ! radiation balance
      !   Q_lw=Q_lw_iso - (1.0-fwet)*SUM(capp*rmair*(tlfy-tair)*gri &
      !                 - canopy%fhvw*gr/ghw
      !   Q_lw=(1-transd)*(L_soil+L_sky-2.0*L_can)
      ! therefore
      !   Q_lw_iso-Q_lw=2(1-transd)*emleaf*(Tv^4-Tc^4)
      
      !rad%lwabv = (1.0-canopy%fwet)*(capp*rmair*(tlfy(:,1) &
      !          - (met%tk-tfrz))*rad%gradis(:,1) &
      !          + capp*rmair*(tlfy(:,2) - (met%tk-tfrz))*rad%gradis(:,2)) &
      !          + canopy%fhvw*SUM(rad%gradis,2)/ghwet
      rad%lwabv = capp * rmair * (tlfy-met%tk) * SUM(rad%gradis,2) ! YP nov2009
!      rad%lwabv = (1.0-canopy%fwet) * capp * rmair * (tlfy-met%tvair) &
!        * SUM(rad%gradis,2) + canopy%fhvw * SUM(rad%gradis,2) / MAX(ghwet,0.001)

      ! add if condition here to avoid dividing by zero
      ! ie when rad%transd=1.0 Ypw:24-02-2003
      WHERE (canopy%vlaiw > 0.01 .AND. rough%hruff > rough%z0soilsn)
         canopy%tv = (rad%lwabv / (2.0*(1.0-rad%transd)*sboltz*emleaf) &
              + met%tk**4)**0.25   ! revert back to met%tk (BP nov2009)
!         canopy%tv = (rad%lwabv / (2.0*(1.0-rad%transd)*sboltz*emleaf) &
!              + met%tvair**4)**0.25
         !                 YP & Mao (jun08) replaced met%tk with tvair
         !                  & + met%tk**4)**0.25
      ELSEWHERE ! sparse canopy
         canopy%tv = met%tvair
         !         YP & Mao (jun08) replaced met%tk with tvair
         !          canopy%tv = met%tk
      END WHERE
      ! Calculate ground heat flux:
      canopy%ghflux = (1-ssoil%isflag)*canopy%ghflux &
           + ssoil%isflag*canopy%sghflux
      ! Calculate net rad to soil:
      canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd) &
           *emleaf*sboltz*canopy%tv**4 - emsoil*sboltz* tss4
      ! Penman-Monteith formula
      sss=air%dsatdk
      cc1=sss/(sss+air%psyc )
      cc2=air%psyc /(sss+air%psyc )
      ssoil%potev = cc1 * (canopy%fns - canopy%ghflux) + cc2 * air%rho &
                 * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil
!      ! method alternative to P-M formula above
!      ! Saturation specific humidity at soil/snow surface temperature:
!      qstss = qsatf((ssoil%tss - tfrz),met%pmb)
!      ! Spec hum deficit at soil/snow surface:
!      dq = qstss - met%qv
!      ! excessive dew over snow area
!      WHERE (ssoil%snowd > 0.1)
!         dq = MAX( -0.1e-3, dq)
!      END WHERE
!      ssoil%potev =air%rho * air%rlam * dq /ssoil%rtsoil
!      ! end method alternative to P-M formula above
      ! Soil latent heat:
      canopy%fes= ssoil%wetfac * ssoil%potev
      WHERE (ssoil%snowd < 0.1 .AND. canopy%fes > 0.0)
         ! Reduce for wilting point limitation:
         canopy%fes= MIN( canopy%fes, &
              MAX( 0.0, (REAL(ssoil%wb(:,1),r_1)-soil%swilt) ) &
              *soil%zse(1)*1000.0*air%rlam/dels )
         ! Reduce for soil ice limitation:
         canopy%fes = MIN(canopy%fes,REAL(ssoil%wb(:,1)-ssoil%wbice(:,1),r_1) &
              * soil%zse(1) * 1000.0 * air%rlam / dels)
      END WHERE
      ssoil%cls=1.
      WHERE (ssoil%snowd >= 0.1)
         ssoil%cls =1.1335
         canopy%fes=MIN(ssoil%wetfac*ssoil%potev, &
                        ssoil%snowd/dels*air%rlam*ssoil%cls)
      END WHERE
      ! Calculate soil sensible heat:
      canopy%fhs = air%rho*capp*(ssoil%tss - met%tvair) /ssoil%rtsoil
      ! Calculate ground heat flux:
      canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssoil%cls
      ! Calculate total latent heat:
      canopy%fe = canopy%fev + canopy%fes
      ! Calculate total sensible heat:
      canopy%fh = canopy%fhv + canopy%fhs
     
      WHERE ( veg%meth > 0 .AND. canopy%vlaiw > 0.01 .AND. &
           rough%hruff > rough%z0soilsn ) 
         ! use the dispersion matrix (DM) to find the air temperature 
         ! and specific humidity (Raupach, Finkele and Zhang 1997, pp 17)
         ! leaf boundary layer resistance for water
         rbw = air%cmolar/REAL(SUM(gbhu+gbhf,2),r_1)
         ! leaf stomatal resistance for water
         rsw = air%cmolar/REAL(SUM(gswx,2),r_1)
         ! A_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmah = (rt0+rough%rt1) * ((1.+air%epsi)/rsw + 1.0/rbw) &
              + air%epsi * (rt0*rough%rt1) / (rbw*rsw)
         ! B_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmbh = (-air%rlam/capp)*(rt0*rough%rt1)/(rbw*rsw)
         ! C_{H} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmch = ((1.+air%epsi)/rsw + 1.0/rbw) * rt0 * rough%rt1 &
              * (canopy%fhv+canopy%fhs) / (air%rho*capp)
         ! A_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmae = (-air%epsi*capp/air%rlam)*(rt0*rough%rt1)/(rbw*rsw)
         ! B_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmbe = (rt0 + ssoil%wetfac*rough%rt1) * ((1.+air%epsi)/rsw + 1.0/rbw) &
              + (rt0*rough%rt1) / (rbw*rsw)
         ! C_{E} in eq. 3.41, SCAM manual, CSIRO tech doc 132
         dmce = ((1.+air%epsi)/rsw + 1.0/rbw) * rt0 * rough%rt1 &
              * (canopy%fev + canopy%fes) / (air%rho*air%rlam)
         ! Within canopy air temperature:
         met%tvair = met%tk &
              + (dmbe*dmch-dmbh*dmce)/(dmah*dmbe-dmae*dmbh+1.0e-12)
         ! Within canopy specific humidity:
         met%qvair = met%qv &
              + (dmah*dmce-dmae*dmch)/(dmah*dmbe-dmae*dmbh+1.0e-12)
         met%qvair = MAX(0.0,met%qvair)
      END WHERE
!      WRITE(46,*) ktau, iter, mp,'rt0', rt0, 'vlaiw', canopy%vlaiw, 'hc', veg%hc

      ! Saturated specific humidity in canopy:
      qstvair = qsatf((met%tvair-tfrz),met%pmb)
      ! Saturated vapour pressure deficit in canopy:
      met%dva = (qstvair - met%qvair) *  rmair/rmh2o * met%pmb * 100.
      ! Set radiative temperature as within canopy air temp:
      met%tvrad = met%tvair
      ! recalculate air%dsatdk (EAK aug08)
      CALL define_air (met, air)

      ! Ground heat flux:
      canopy%ghflux = (1-ssoil%isflag)*canopy%ghflux &
           + ssoil%isflag*canopy%sghflux
      ! Net radiation absorbed by soil: 
      canopy%fns = rad%qssabs + rad%transd*met%fld + (1.0-rad%transd) &
           *emleaf*sboltz*canopy%tv**4 - emsoil*sboltz* tss4
      ! YP & Mao (jun08) reinstating Penman-Monteith formula
      ! Penman-Monteith formula
      sss=air%dsatdk
      cc1=sss/(sss+air%psyc )
      cc2=air%psyc /(sss+air%psyc )
      ssoil%potev = cc1 * (canopy%fns - canopy%ghflux) + cc2 * air%rho &
              * air%rlam*(qsatf((met%tk-tfrz),met%pmb) - met%qv)/ssoil%rtsoil
!      ! method alternative to P-M formula above
!      dq = qstss - met%qvair
!      WHERE (ssoil%snowd > 0.1)
!         dq = MAX( -0.1e-3, dq)
!      END WHERE
!      ssoil%potev =air%rho * air%rlam * dq /ssoil%rtsoil
      ! Soil evaporation:
      canopy%fes= ssoil%wetfac * ssoil%potev
      WHERE (ssoil%snowd < 0.1 .AND. canopy%fes > 0.0 )
         ! Reduce for wilting point limitation:
         canopy%fes = MIN( canopy%fes , &
              MAX( 0., (REAL(ssoil%wb(:,1),r_1)-soil%swilt)) &
              * soil%zse(1) * 1000.0 * air%rlam / dels )
         ! Reduce for soil ice limitation:
         canopy%fes = MIN( canopy%fes , &
              REAL( ssoil%wb(:,1)-ssoil%wbice(:,1) , r_1 ) &
              * soil%zse(1) * 1000.0 * air%rlam / dels )
      END WHERE
      ssoil%cls=1.
      WHERE (ssoil%snowd >= 0.1)
         ssoil%cls  = 1.1335
         canopy%fes = MIN( ssoil%wetfac*ssoil%potev , &
              ssoil%snowd/dels*air%rlam*ssoil%cls)
      END WHERE
      ! Soil sensible heat:
      canopy%fhs = air%rho*capp*(ssoil%tss - met%tvair) /ssoil%rtsoil
      canopy%ga = canopy%fns-canopy%fhs-canopy%fes*ssoil%cls
      ! Set total latent heat:
      canopy%fe = canopy%fev + canopy%fes
      ! Set total sensible heat:
      canopy%fh = canopy%fhv + canopy%fhs

      ! monin-obukhov stability parameter zetar=zref/l
      !	recompute zetar for the next iteration, except on last iteration
      IF (iter < niter) THEN ! dont compute zetar on the last iter
         iterplus = MAX(iter+1,2)
         zetar(:,iterplus) = -( vonk*grav*rough%zref_tq &
              * (canopy%fh+0.07*canopy%fe) ) &
              / (air%rho*capp*met%tk*canopy%us**3)
         ! case niter=2: final zetar=zetmul*zetar(2) (compute only when iter=1)
         IF (niter == 2) THEN
            zetar(:,2) = zetmul * zetar(:,2)
            WHERE (met%fsd ==  0.0)
               zetar(:,2) = 0.5 * zetar(:,2)
            END WHERE
         END IF
         ! constrain zeta to zetpos and zetneg (set in param0)
         zetar(:,iterplus) = MIN(zetpos,zetar(:,iterplus))  ! zetar too +
         zetar(:,iterplus) = MAX(zetneg,zetar(:,iterplus))  ! zetar too -
      END IF
   END DO ! DO iter = 1, niter

   ! screen temp., windspeed and relative humidity at 1.8m
   tstar = - canopy%fh / ( air%rho*capp*canopy%us)
   qstar = - canopy%fe / ( air%rho*air%rlam *canopy%us)
   zscrn = MAX(rough%z0m,2.0-rough%disp)    ! change screen temp to be measured at 2 m (BP nov2009)
!   zscrn = MAX(rough%z0m,1.8-rough%disp)
   denom = ( LOG(rough%zref_tq/zscrn)- psim(zetar(:,iterplus)) &
        + psim(zetar(:,iterplus) * zscrn / rough%zref_tq) ) /vonk

!!%% change by Ashok Luhar !! need to comment out the above line with MIN(zetpos,zetar)
!    WHERE (zetar(:,iterplus) > 0.7)
!      zeta2=zetar(:,iterplus) * zscrn / rough%zref_tq
!      denom =alpha1* ((zetar(:,iterplus)**beta1* &
!               (1.0+gamma1*zetar(:,iterplus)**(1.0-beta1)))  &
!             - (zeta2*beta1*(1.0+gamma1*zeta2**(1.0-beta1)))) /vonk
!    ENDWHERE
!!%%

   ! Calculate screen temperature for bare ground:
   canopy%tscrn = met%tk-tfrz - tstar * denom
   rsts = qsatf(canopy%tscrn, met%pmb)
   qtgnet = rsts * ssoil%wetfac - met%qv
   canopy%cduv = canopy%us * canopy%us / (MAX(met%ua,umin))**2 ! EK jun08
   !      canopy%cduv = canopy%us * canopy%us / MAX(met%ua,umin)
   WHERE (qtgnet > 0.0)
      qsurf = rsts * ssoil%wetfac
   ELSEWHERE
      qsurf = 0.1*rsts*ssoil%wetfac + 0.9*met%qv
   END WHERE
   canopy%qscrn = qsurf + qstar * denom

   ! calculation of screen temepratures for LAI > 0.1 . Method by Ian Harman
   term1=0.
   term2=0.
   term5=0.
   r_sc =0.
   zscl = MAX(rough%z0soilsn,2.0)

   WHERE ( canopy%vlaiw > 0.01)    ! BP jul2010
!   WHERE ( canopy%vlaiw > 0.1 )
     WHERE ( rough%hruff  > 0.0 .and. rough%disp  > 0.0 )
       term1 = EXP(2*csw*canopy%vlaiw*(1-zscl/rough%hruff))
       term2 = EXP(2*csw*canopy%vlaiw*(1-rough%disp/rough%hruff))
       term5 = MAX(2./3.*rough%hruff/rough%disp, 1.)
     ENDWHERE
     term3 = a33**2*ctl*2*csw*canopy%vlaiw
     WHERE( zscl < rough%disp )
       r_sc = term5 * LOG(zscl/rough%z0soilsn) &
            * (EXP(2*csw*canopy%vlaiw) - term1) / term3
     ELSEWHERE ( rough%disp <= zscl .and. zscl < rough%hruff )
       r_sc = rough%rt0us + term5 * ( term2 - term1 ) / term3
     ELSEWHERE ( rough%hruff <= zscl .and. zscl <  rough%zruffs )
       r_sc = rough%rt0us + rough%rt1usa + term5 * ( zscl - rough%hruff )  &
                                                 / (a33**2*ctl*rough%hruff)
     ELSEWHERE (zscl >= rough%zruffs )
       r_sc = rough%rt0us + rough%rt1usa + rough%rt1usb  &
            + ( LOG( (zscl - rough%disp) &
                   / MAX(rough%zruffs - rough%disp, rough%z0soilsn) ) &
              - psis((zscl - rough%disp) * zetar(:,iterplus) / rough%zref_tq) &
              + psis((rough%zruffs - rough%disp) * zetar(:,iterplus) &
                    / rough%zref_tq )  &
            ) / vonk
     ENDWHERE   
     canopy%tscrn = ssoil%tss + (met%tk - ssoil%tss) * r_sc  &
            / (rough%rt0us + rough%rt1usa + rough%rt1usb + rt1usc) - tfrz  
   ENDWHERE    ! end tscrn calculation from Ian Harman

   canopy%uscrn = MAX(0.0, MAX(met%ua,umin) - canopy%us*denom ) ! at present
   ! incorrect
   poolcoef1=(SUM(spread(bgc%ratecp,1,mp)*bgc%cplant,2) - &
        bgc%ratecp(1)*bgc%cplant(:,1))
   poolcoef1w=(SUM(spread(bgc%ratecp,1,mp)*bgc%cplant,2) -  &
        bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(3)*bgc%cplant(:,3))
   poolcoef1r=(SUM(spread(bgc%ratecp,1,mp)*bgc%cplant,2) -  &
        bgc%ratecp(1)*bgc%cplant(:,1) - bgc%ratecp(2)*bgc%cplant(:,2))

   ! Carbon uptake from photosynthesis: 
   tmp1(:) = 3.22 - 0.046 * (met%tk(:)-tfrz)
   IF(ANY(tmp1(:) < 0.001)) THEN
     PRINT *, 'Warning: extreme air temperature and ktau =',MAXVAL(met%tk),ktau
     PRINT *, 'Warning: latitude and landpoint # = ', &
                     rad%latitude(MAXLOC(met%tk)), MAXLOC(met%tk)
   ENDIF
   tmp1(:) = MAX(tmp1(:), 0.001)
   tmp2(:) = 0.1 * (met%tk(:)-tfrz-20.0)
   tmp3(:) = tmp1(:) ** tmp2(:)
   canopy%frp  = veg%rp20 * tmp3 * poolcoef1  / (365.0*24.0*3600.0)
   canopy%frpw = veg%rp20 * tmp3 * poolcoef1w / (365.0*24.0*3600.0)
   canopy%frpr = veg%rp20 * tmp3 * poolcoef1r / (365.0*24.0*3600.0)
!   canopy%frp = veg%rp20*((3.22-0.046*(met%tk-tfrz))**(0.1*(met%tk-tfrz-20.))) &
!        * poolcoef1 /(365.0*24.0*3600.0)   ! 24/05
!   canopy%frpw = veg%rp20*((3.22-0.046*(met%tk-tfrz))**(0.1*(met%tk-tfrz-20.))) &
!        * poolcoef1w /(365.0*24.0*3600.0)    ! 24/05
!   canopy%frpr = veg%rp20*((3.22-0.046*(met%tk-tfrz))**(0.1*(met%tk-tfrz-20.))) &
!        * poolcoef1r /(365.0*24.0*3600.0)    ! 24/05

    ! Calculate dewfall: from negative lh wet canopy + neg. lh dry canopy:
    canopy%dewmm = - REAL(MIN(0.0_r_2,canopy%fevw) + &
        MIN(0.0_r_2,canopy%fevc),r_1) * dels * 1.0e3 / (rhow*air%rlam)
    ! Add dewfall to canopy water storage:
    canopy%cansto = canopy%cansto + canopy%dewmm
    ! Modify canopy water storage for evaporation:
    canopy%cansto = MAX(canopy%cansto-MAX(0.0,REAL(canopy%fevw,r_1))*dels*1.0e3/ &
        (rhow*air%rlam), 0.0)
    ! Calculate canopy water storage excess:
    !canopy%spill=MAX(0.,MIN(0.2*canopy%cansto,MAX(0.0, canopy%cansto-cansat)))
    canopy%spill=MAX(0.0, canopy%cansto-cansat)
    ! Move excess canopy water to throughfall:
    canopy%through = canopy%through + canopy%spill
    ! Initialise 'throughfall to soil' as 'throughfall from canopy';
    ! snow may absorb
    canopy%precis = canopy%through
    ! Update canopy storage term:
    canopy%cansto=canopy%cansto - canopy%spill
    ! Modify canopy water storage for evaporation:
!    canopy%cansto = MAX(canopy%cansto-MAX(0.0,REAL(canopy%fevw,r_1))*dels*1.0e3/ &
!        (rhow*air%rlam), 0.0)
    ! Calculate the total change in canopy water store (mm/dels):
    canopy%delwc = canopy%cansto-oldcansto
    ! calculate dgdtg, derivative of ghflux
    ssoil%dfn_dtg = (-1.)*4.*emsoil*sboltz*tss4/ssoil%tss ! d(canopy%fns)
    ! /d(ssoil%tgg)
    ssoil%dfh_dtg = air%rho*capp/ssoil%rtsoil  ! d(canopy%fhs)/d(ssoil%tgg)
    ssoil%dfe_ddq = ssoil%wetfac*air%rho*air%rlam/ssoil%rtsoil !d(canopy%fes)/d(dq)
    ssoil%ddq_dtg = (rmh2o/rmair)/met%pmb *tetena*tetenb*tetenc &
        / ((tetenc+ssoil%tss-tfrz)**2) &
        * EXP(tetenb*(ssoil%tss-tfrz)/(tetenc+ssoil%tss-tfrz))
    canopy%dgdtg = ssoil%dfn_dtg - ssoil%dfh_dtg &
        - ssoil%cls * ssoil%dfe_ddq * ssoil%ddq_dtg
    !ypw: energy balance of the dry canopy
    !    bal%drybal=ecy(:,1)+ecy(:,2)+hcy(:,1)+hcy(:,2) &
    !         -rad%rniso(:,1)-rad%rniso(:,2) &
    !         +capp*rmair*(tlfy(:,1)-(met%tvair-tfrz))*rad%gradis(:,1) &
    !         +capp*rmair*(tlfy(:,2)-(met%tvair-tfrz))*rad%gradis(:,2)
    bal%drybal=REAL(ecy+hcy,r_1)-SUM(rad%rniso,2) &
         +capp*rmair*(tlfy-met%tk)*SUM(rad%gradis,2)  ! YP nov2009
!         +capp*rmair*(tlfy-met%tvair)*SUM(rad%gradis,2)

    !ypw: energy balance of the wet canopy
    bal%wetbal=canopy%fevw+canopy%fhvw-SUM(rad%rniso,2)*canopy%fwet &
         +capp*rmair*(tlfy-met%tk)*SUM(rad%gradis,2)*canopy%fwet  ! YP nov2009
!         +canopy%fhvw*SUM(rad%gradis,2)/MAX(0.001,ghwet)

    DEALLOCATE(cansat, ghwet, dsx)
    DEALLOCATE(fwsoil, tlfx, tlfy)
    DEALLOCATE(ecy, hcy, rny)
    DEALLOCATE(gbhu, gbhf, gswmin, gswx, csx)

  END SUBROUTINE define_canopy

  !--------------------------------------------------------------------------
  ELEMENTAL FUNCTION qsatf(tair,pmb) RESULT(r)
    ! MRR, 1987
    ! AT TEMP tair (DEG C) AND PRESSURE pmb (MB), GET SATURATION SPECIFIC
    ! HUMIDITY (KG/KG) FROM TETEN FORMULA
    REAL(r_1), INTENT(IN) :: tair ! air temperature (C)
    REAL(r_1), INTENT(IN) :: pmb  ! pressure PMB (mb)
    REAL(r_1)             :: r    ! result; sat sp humidity
    r = (rmh2o/rmair) * (tetena*EXP(tetenb*tair/(tetenc+tair))) / pmb
  END FUNCTION qsatf
  !---------------------------------------------------------
  ELEMENTAL FUNCTION ej3x(parx,x) result(z)
    REAL(r_1), INTENT(IN) :: parx
    REAL(r_1), INTENT(IN) :: x
    REAL(r_1)             :: z
    z = MAX( 0.0 , &
         0.25*((alpha3*parx+x-SQRT((alpha3*parx+x)**2 - &
         4.0*convx3*alpha3*parx*x)) /(2.0*convx3)) )
  END FUNCTION ej3x
  !---------------------------------------------------------
  ELEMENTAL FUNCTION ej4x(parx,x) result(z)
    REAL(r_1), INTENT(IN) :: parx
    REAL(r_1), INTENT(IN) :: x
    REAL(r_1)             :: z
    z = MAX( 0.0 , &
         (alpha4*parx+x-SQRT((alpha4*parx+x)**2 - &
         4.0*convx4*alpha4*parx*x))/(2.0*convx4))
  END FUNCTION ej4x
  !---------------------------------------------------------
  ! Explicit array dimension as temporary work around for NEC inlining problem
  FUNCTION xvcmxt4(x) result(z)
    REAL(r_1), PARAMETER :: q10c4 = 2.0
    REAL(r_1), INTENT(IN):: x
    REAL(r_1)            :: z
    z = q10c4 ** (0.1 * x - 2.5) / &
         ((1.0 + EXP(0.3 * (13.0 - x))) * (1.0 + EXP(0.3 * (x - 36.0))))
  END FUNCTION xvcmxt4
  !---------------------------------------------------------
  FUNCTION xvcmxt3(x) result(z)
    !  leuning 2002 (p c & e) equation for temperature response
    !  used for vcmax for c3 plants
    REAL(r_1), INTENT(IN) :: x
    REAL(r_1), PARAMETER  :: EHaVc  = 73637.0  ! J/mol (Leuning 2002)
    REAL(r_1), PARAMETER  :: EHdVc  = 149252.0 ! J/mol (Leuning 2002)
    REAL(r_1), PARAMETER  :: EntropVc = 486.0  ! J/mol/K (Leuning 2002)
    REAL(r_1), PARAMETER  :: xVccoef = 1.17461 ! derived parameter
    REAL(r_1)             :: xvcnum
    REAL(r_1)             :: xvcden
    REAL(r_1)             :: z
    xvcnum = xVccoef * EXP((EHaVc/(rgas*trefk))*(1.-trefk/x))
    xvcden = 1.0 + EXP((EntropVc*x-EHdVc)/(rgas*x))
    z = MAX(0.0,xvcnum/xvcden)
  END FUNCTION xvcmxt3
  !---------------------------------------------------------
  FUNCTION xejmxt3(x) result(z)
    !  leuning 2002 (p c & e) equation for temperature response
    !  used for jmax for c3 plants
    REAL(r_1), INTENT(IN) :: x
    REAL(r_1), PARAMETER  :: EHaJx  = 50300.0  ! J/mol (Leuning 2002)
    REAL(r_1), PARAMETER  :: EHdJx  = 152044.0 ! J/mol (Leuning 2002)
    REAL(r_1), PARAMETER  :: EntropJx = 495.0  ! J/mol/K (Leuning 2002)
    REAL(r_1), PARAMETER  :: xjxcoef = 1.16715 ! derived parameter
    REAL(r_1)             :: xjxnum
    REAL(r_1)             :: xjxden
    REAL(r_1)             :: z
    xjxnum = xjxcoef * EXP((ehajx/(rgas*trefk))*(1.-trefk/x))
    xjxden = 1.0 + EXP((entropjx*x-ehdjx)/(rgas*x))
    z = MAX(0.0, xjxnum/xjxden)
  END FUNCTION xejmxt3
  !---------------------------------------------------------
  ELEMENTAL FUNCTION psim(zeta) result(r)
    ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
    ! computes integrated stability function psim(z/l) (z/l=zeta)
    ! for momentum, using the businger-dyer form for unstable cases
    ! and the webb form for stable cases. see paulson (1970).
    USE math_constants
    REAL(r_1), INTENT(IN) :: zeta
    REAL(r_1)             :: r
    REAL(r_1)             :: x
    REAL(r_1), PARAMETER  :: gu = 16.0
    REAL(r_1), PARAMETER  :: gs = 5.0
    x = (1.0 + gu * ABS(zeta))**0.25
    r = MERGE(LOG((1.0+x*x)*(1.0+x)**2/8.0) - 2.0*ATAN(x) &
         + pi_c *0.5, -1.0*gs*zeta, zeta < 0.0)
  END FUNCTION psim
  !---------------------------------------------------------
  ELEMENTAL FUNCTION psis(zeta) result(r)
    ! mrr, 16-sep-92 (from function psi: mrr, edinburgh 1977)
    ! computes integrated stability function psis(z/l) (z/l=zeta)
    ! for scalars, using the businger-dyer form for unstable cases
    ! and the webb form for stable cases. see paulson (1970).
    REAL(r_1), INTENT(IN) :: zeta
    REAL(r_1)             :: r
    REAL(r_1), PARAMETER :: gu = 16.0
    REAL(r_1), PARAMETER :: gs = 5.0
    r = MERGE(2.0 * LOG((1.0 + SQRT(1.0 + gu * ABS(zeta))) * 0.5), &
         -1.0 * gs * zeta, zeta < 0.0)
  END FUNCTION psis
  !---------------------------------------------------------

  SUBROUTINE dryLeaf(dels,rad,air,met,ssoil,soil,veg,canopy,ktau)
    TYPE (radiation_type),      INTENT(INOUT) :: rad
    TYPE (air_type),            INTENT(INOUT) :: air
    TYPE (met_type),            INTENT(INOUT) :: met
    TYPE (soil_snow_type),      INTENT(INOUT) :: ssoil
    TYPE (soil_parameter_type), INTENT(INOUT) :: soil
    TYPE (veg_parameter_type),  INTENT(INOUT) :: veg
    TYPE (canopy_type),         INTENT(INOUT) :: canopy
    INTEGER, INTENT(IN) :: ktau
    REAL(r_1), INTENT(IN)     :: dels ! integration time step (s)
    REAL(r_1), PARAMETER  :: co2cp3 = 0.0 ! CO2 compensation pt C3
    REAL(r_1), PARAMETER  :: jtomol = 4.6e-6 ! Convert from J to Mol for light
    REAL(r_1), DIMENSION(mp)  :: conkct ! Michaelis Menton const.
    REAL(r_1), DIMENSION(mp)  :: conkot ! Michaelis Menton const.
    REAL(r_2), DIMENSION(mp)  :: cx1  ! "d_{3}" in Wang and Leuning,
    REAL(r_2), DIMENSION(mp)  :: cx2  !     1998, appendix E
    REAL(r_1), DIMENSION(mp)  :: tdiff ! leaf air temp diff.
    REAL(r_1), DIMENSION(mp)  :: tlfxx ! leaf temp of current iteration (K)
    REAL(r_1), DIMENSION(mp)  :: abs_deltlf ! ABS(deltlf)
    REAL(r_1), DIMENSION(mp)  :: deltlf ! deltlfy of prev iter.
    REAL(r_1), DIMENSION(mp)  :: deltlfy ! del temp successive iter.
    REAL(r_1), DIMENSION(mp)  :: gras ! Grashof coeff
    REAL(r_1), DIMENSION(mp)  :: evapfb !
    REAL(r_1), DIMENSION(mp)  :: temp
    REAL(r_2), DIMENSION(mp)  :: ecx ! lat. hflux big leaf
    REAL(r_2), DIMENSION(mp)  :: hcx ! sens heat fl big leaf prev iteration
    REAL(r_2), DIMENSION(mp)  :: rnx ! net rad prev timestep
! ypw testing new formulation (26oct2010), move evapfbl to canopy_type
!    REAL(r_1), DIMENSION(mp,ms)  :: evapfbl !
    REAL(r_2), DIMENSION(mp,mf)  :: gw  ! cond for water for a dry canopy
    REAL(r_2), DIMENSION(mp,mf)  :: gh  ! cond for heat for a dry canopy
    REAL(r_2), DIMENSION(mp,mf)  :: ghr ! dry canopy cond for heat & thermal rad
    REAL(r_2), DIMENSION(mp,mf)  :: anx ! net photos. prev iteration
    REAL(r_2), DIMENSION(mp,mf)  :: an_y ! net photosynthesis soln
    REAL(r_1), DIMENSION(mp,mf)  :: rdx ! daytime leaf resp rate, prev iteration
    REAL(r_1), DIMENSION(mp,mf)  :: rdy ! daytime leaf resp rate
    REAL(r_1), DIMENSION(mp,mf)  :: ejmax2 ! jmax of big leaf
    REAL(r_1), DIMENSION(mp,mf)  :: ejmxt3 ! jmax big leaf C3 plants
    REAL(r_1), DIMENSION(mp,mf)  :: vcmxt3 ! vcmax big leaf C3
    REAL(r_1), DIMENSION(mp,mf)  :: vcmxt4 ! vcmax big leaf C4
    REAL(r_2), DIMENSION(mp,mf)  :: vx3 ! carboxylation C3 plants
    REAL(r_2), DIMENSION(mp,mf)  :: vx4 ! carboxylation C4 plants
    REAL(r_2), DIMENSION(mp,mf)  :: xleuning ! leuning stomatal coeff
    REAL(r_1), DIMENSION(mp,mf)  :: psycst ! modified pych. constant
    REAL(r_1), DIMENSION(mp,mf)  :: temp2
    INTEGER(i_d)   :: k  ! iteration count
    INTEGER(i_d)   :: i, kk  ! iteration count
    REAL, PARAMETER :: min_lai = 0.01
    REAL(r_2) :: oldevapfbl(mp,ms)

    gw  = 1.0e-3 ! default values of conductance
    gh  = 1.0e-3
    ghr = 1.0e-3
    gras = 1.0e-6
    anx = 0.0
    an_y= 0.0
    hcx = 0.0              ! init sens heat iteration memory variable
    hcy = 0.0
    rdx = 0.0
    rdy = 0.0
    ecx = SUM(rad%rniso,2) ! init lat heat iteration memory variable
    rnx = SUM(rad%rniso,2)
    tlfxx = tlfx
    abs_deltlf = 999.0
    psycst(:,:) = SPREAD(air%psyc,2,mf)
    canopy%fevc    = 0.0
    canopy%evapfbl = 0.0
    DO kk = 1, mp
      IF(canopy%vlaiw(kk) <=1.0e-2) THEN
        abs_deltlf(kk)=0.0
        ecx(kk) = 0.0 ! intialise
        rnx(kk) = 0.0 ! intialise
        rny(kk) = rnx(kk) ! store initial values
        ecy(kk) = ecx(kk) ! store initial values
      END IF
    ENDDO
    deltlfy = abs_deltlf
    k = 0
    DO WHILE (ANY(abs_deltlf > 0.1)  .AND.  k < maxiter)
      k = k + 1
      DO i = 1, mp
      IF (canopy%vlaiw(i) > min_lai .AND. abs_deltlf(i) > 0.1) THEN
        ! Grashof number (Leuning et al, 1995) eq E4:
        gras(i) = MAX(1.0e-6, &
                  1.595E8*ABS(tlfx(i)-met%tvair(i))*(veg%dleaf(i)**3.0))
        ! See Appendix E in (Leuning et al, 1995):
        gbhf(i,:) = rad%fvlai(i,:) * air%cmolar(i) * 0.5 * dheat &
             * (gras(i)**0.25) / veg%dleaf(i)
        ! Conductance for heat:
        gh(i,:) = 2.0 * (gbhu(i,:) + gbhf(i,:))
        ! Conductance for heat and longwave radiation:
        ghr(i,:) = rad%gradis(i,:)+gh(i,:)
        !  Leuning 2002 (P C & E) equation for temperature response
        !  used for Vcmax for C3 plants:
        temp(i) =  xvcmxt3(tlfx(i)) * veg%vcmax(i) * (1.0-veg%frac4(i))
        vcmxt3(i,:) = rad%scalex(i,:) * temp(i)
        ! Temperature of Vcmax for C4 plants (Collatz et al 1989):
        temp(i) =  xvcmxt4(tlfx(i)-tfrz) * veg%vcmax(i) * veg%frac4(i)
        vcmxt4(i,:) = rad%scalex(i,:) * temp(i)
        !  Leuning 2002 (P C & E) equation for temperature response
        !  used for Jmax for C3 plants:
        temp(i) = xejmxt3(tlfx(i)) * veg%ejmax(i) * (1.0-veg%frac4(i))
        ejmxt3(i,:) = rad%scalex(i,:) * temp(i)
        ! Difference between leaf temperature and reference temperature:
        tdiff(i) =  tlfx(i) - trefk
        ! Michaelis menten constant of Rubisco for CO2:
        conkct(i) = conkc0 * EXP((ekc/(rgas* trefk)) * (1.0-trefk/tlfx(i)))
        ! Michaelis menten constant of Rubisco for oxygen:
        conkot(i) = conko0 * EXP((eko/(rgas* trefk)) * (1.0-trefk/tlfx(i)))
        ! Store leaf temperature:
        tlfxx(i) = tlfx(i)
        ! "d_{3}" in Wang and Leuning, 1998, appendix E:
        cx1(i) = conkct(i) * (1.0+0.21/conkot(i))
        cx2(i) = 2.0 * gam0 * (1.0+gam1*tdiff(i) + gam2*tdiff(i)*tdiff(i))

        ! All equations below in appendix E in Wang and Leuning 1998 are
        ! for calculating anx, csx and gswx for Rubisco limited,
        ! RuBP limited, sink limited
!! ypw testing new formulation (26oct2010)
!        temp2(i,:) =rad%qcan(i,:,1)*jtomol*(1.0-veg%frac4(i))
!        vx3(i,:) = ej3x(temp2(i,:), ejmxt3(i,:)) * fwsoil(i)
!        temp2(i,:) = rad%qcan(i,:,1)*jtomol*veg%frac4(i)
!        vx4(i,:) = ej4x(temp2(i,:), vcmxt4(i,:)) * fwsoil(i)
!        rdx(i,:) = (cfrd3 * vcmxt3(i,:) + cfrd4 * vcmxt4(i,:)) * fwsoil(i)
!        xleuning(i,:) = (1.0-veg%frac4(i)) * a1c3 / (1.0 + dsx(i)/d0c3) &
!                       + veg%frac4(i) * a1c4 / (1.0 + dsx(i)/d0c4)
!        xleuning(i,:) = xleuning(i,:) / (csx(i,:)-co2cp3)

        temp2(i,:) =rad%qcan(i,:,1)*jtomol*(1.0-veg%frac4(i))
        vx3(i,:) = ej3x(temp2(i,:), ejmxt3(i,:))
        temp2(i,:) = rad%qcan(i,:,1)*jtomol*veg%frac4(i)
        vx4(i,:) = ej4x(temp2(i,:), vcmxt4(i,:))
        rdx(i,:) = (cfrd3 * vcmxt3(i,:) + cfrd4 * vcmxt4(i,:)) * fwsoil(i)
        xleuning(i,:) = fwsoil(i) &
                   * ( (1.0-veg%frac4(i)) * a1c3 / (1.0 + dsx(i)/d0c3) &
                       + veg%frac4(i) * a1c4 / (1.0 + dsx(i)/d0c4) )
        xleuning(i,:) = xleuning(i,:) / (csx(i,:)-co2cp3)
!! ypw (26oct2010)
      END IF
      ENDDO !i=1,mp

      CALL photosynthesis(csx(:,:),SPREAD(cx1(:),2,mf), &
                 SPREAD(cx2(:),2,mf),gswmin(:,:),rdx(:,:), &
                 vcmxt3(:,:),vcmxt4(:,:),vx3(:,:),vx4(:,:),xleuning(:,:), &
                 rad%fvlai(:,:),SPREAD(abs_deltlf,2,mf),anx(:,:))

      DO i = 1, mp
      IF (canopy%vlaiw(i) > min_lai .AND. abs_deltlf(i) > 0.1) THEN
        ecx(i) = 0.0
        DO kk = 1, mf
        IF (rad%fvlai(i,kk) > 0.0) THEN
          csx(i,kk) = met%ca(i) - anx(i,kk) * rgbwc / (gbhu(i,kk) + gbhf(i,kk))
          csx(i,kk) = MAX(1.0e-4,csx(i,kk))   ! BP (aug2010)
          gswx(i,kk) = gswmin(i,kk)+MAX(0.0_r_2,rgswc*xleuning(i,kk) *anx(i,kk))
          gswx(i,kk) = MAX(1.0e-3,gswx(i,kk)) ! Q.Zhang (28/03/2011)
          ! Recalculate conductance for water:
          gw(i,kk) = 1.0/(1.0/gswx(i,kk) + 1.0/(1.075*(gbhu(i,kk)+gbhf(i,kk))))
          ! Modified psychrometric constant (Monteith and Unsworth, 1990)
          psycst(i,kk) = air%psyc(i) * REAL( ghr(i,kk)/gw(i,kk) ,r_1)
          ! Update canopy latent heat flux:
          ecx(i) = ecx(i) + ( air%dsatdk(i) * (rad%rniso(i,kk) &
                 - capp*rmair*(met%tvair(i)-met%tk(i))*rad%gradis(i,kk)) &
                 + capp * rmair * met%dva(i) * ghr(i,kk)) &
                 / (air%dsatdk(i) + psycst(i,kk))
        END IF
        ENDDO ! kk = 1, mf

!        ! Update canopy latent heat flux:
!        ecx(i) = ( air%dsatdk(i) * (rad%rniso(i,1) &
!                 - capp*rmair*(met%tvair(i)-met%tk(i))*rad%gradis(i,1)) &
!                 + capp * rmair * met%dva(i) * ghr(i,1)) &
!                 / (air%dsatdk(i) + psycst(i,1)) &
!               + ( air%dsatdk(i) * (rad%rniso(i,2) &
!                 - capp*rmair*(met%tvair(i)-met%tk(i))*rad%gradis(i,2)) &
!                 + capp * rmair * met%dva(i) * ghr(i,2)) &
!                 / (air%dsatdk(i) + psycst(i,2))  ! YP nov2009

!        canopy%fevc    = 0.0
!        canopy%evapfbl = 0.0
        ! convert W/m2 to mm/dt using                dels/air%rlam
        evapfb(i) = (1.0-canopy%fwet(i)) * ecx(i) * dels / air%rlam(i)
        IF (ecx(i) > 0.0 .AND. canopy%fwet(i) < 1.0) THEN
          ! Calcualte contribution by each soil layer to canopy transpiration:
          DO kk = 1, ms
! ypw testing new formulation (26oct2010)
            canopy%evapfbl(i,kk) = MIN(evapfb(i)*veg%froot(i,kk), &
              MAX(0.0, REAL(MIN(ssoil%wb(i,kk)-soil%swilt(i),ssoil%wb(i,kk) &
              -1.05*ssoil%wbice(i,kk)), r_1))*soil%zse(kk)*1000.0)
!            canopy%evapfbl(i,kk) =MIN(evapfb(i)*canopy%rwater(i,kk), &
!              MAX(0.0, REAL(MIN(ssoil%wb(i,kk)-soil%swilt(i),ssoil%wb(i,kk) &
!              -1.05*ssoil%wbice(i,kk)), r_1))*soil%zse(kk)*1000.0)
          ENDDO
!          canopy%fevc(i) = SUM(evapfbl(i,:)) * air%rlam(i) / dels
          canopy%fevc(i) = SUM(canopy%evapfbl(i,:)) * air%rlam(i) / dels
! ypw (26oct2010)
          ecx(i)         = canopy%fevc(i) / (1.0-canopy%fwet(i))
        END IF

        ! Update canopy sensible heat flux:
        hcx(i)=(SUM(rad%rniso(i,:))-ecx(i)-capp*rmair*(met%tvair(i)-met%tk(i)) &
               *SUM(rad%gradis(i,:))) * SUM(gh(i,:))/SUM(ghr(i,:)) ! YP nov2009
        ! Update leaf temperature:
        tlfx(i)=met%tvair(i)+REAL(hcx(i),r_1)/(capp*rmair*SUM(gh(i,:)))

        ! Update net radiation for canopy:
        rnx(i) = SUM(rad%rniso(i,:)) - capp * rmair * (tlfx(i)-met%tk(i)) &
                                         * SUM(rad%gradis(i,:)) ! YP nov2009
        ! Update leaf surface vapour pressure deficit:
        dsx(i) = met%dva(i) + air%dsatdk(i) * (tlfx(i)-met%tvair(i))

        ! Store change in leaf temperature between successive iterations:
        deltlf(i) = tlfxx(i) - tlfx(i)
        abs_deltlf(i) = ABS(deltlf(i))

      END IF
      ENDDO !second i=1,mp

      ! Where leaf temp change b/w iterations is significant, and difference
      ! is smaller than the previous iteration, store results:
      WHERE (abs_deltlf > 0.1 .AND. abs_deltlf < ABS(deltlfy) )
        deltlfy = deltlf
        tlfy = tlfx
        rny = rnx
        hcy = hcx
        ecy = ecx
        rdy(:,1) = rdx(:,1)
        rdy(:,2) = rdx(:,2)
        an_y(:,1) = anx(:,1)
        an_y(:,2) = anx(:,2)
        ! save last values calculated for evapfbl 
        oldevapfbl(:,1) = canopy%evapfbl(:,1)
        oldevapfbl(:,2) = canopy%evapfbl(:,2)
        oldevapfbl(:,3) = canopy%evapfbl(:,3)
        oldevapfbl(:,4) = canopy%evapfbl(:,4)
        oldevapfbl(:,5) = canopy%evapfbl(:,5)
        oldevapfbl(:,6) = canopy%evapfbl(:,6)
      END WHERE
      WHERE (abs_deltlf > 0.1)
        ! after four iterations, take the mean value of current and previous
        ! estimates as the next estimate of leaf temperature, to avoid
        ! oscillation
        tlfx = (0.5*(MAX(0,k-5)/(k-4.9999))) *tlfxx &
             & + (1.0- (0.5*(MAX(0,k-5)/(k-4.9999))))*tlfx
      END WHERE
      IF(k==1) THEN
        !        taken the first iterated estimates as the defaults
        tlfy = tlfx
        rny = rnx
        hcy = hcx
        ecy = ecx
        rdy = rdx
        an_y = anx
        ! save last values calculated for evapfbl
        oldevapfbl(:,:) = canopy%evapfbl(:,:)
      END IF
    END DO  ! DO WHILE (ANY(abs_deltlf > 0.1) .AND. k < maxiter)
    ! dry canopy flux
    canopy%fevc = (1.0-canopy%fwet) * ecy
    ! Recalculate evapfbl as ecy may not be updated with the ecx
    ! calculated in the last iteration.
    ! DO NOT use simple scaling as there are times that canopy%evapfbl is zero.
    ! ** canopy%evapfbl(i,:) = canopy%evapfbl(i,:) * ecy(i) / ecx(i) **
    DO i = 1, mp
      IF (ecy(i) > 0.0 .AND. canopy%fwet(i) < 1.0) THEN
      IF (ABS(ecy(i)-ecx(i)) > 1.0e-6 ) THEN
        IF (ABS(canopy%fevc(i) - (SUM(oldevapfbl(i,:))*air%rlam(i)/dels)) &
              > 1.0e-4 ) THEN
          PRINT *, 'Error! oldevapfbl not right.', ktau, i
          PRINT *, 'ecx, ecy = ', ecx(i), ecy(i)
          PRINT *, 'or in mm = ',ecx(i)*(1.0-canopy%fwet(i))/air%rlam(i)*dels, &
                                 ecy(i)*(1.0-canopy%fwet(i))/air%rlam(i)*dels
          PRINT *,'fevc = ',canopy%fevc(i),SUM(oldevapfbl(i,:))*air%rlam(i)/dels
          PRINT *, 'fwet = ', canopy%fwet(i)
          PRINT *, 'oldevapfbl = ', oldevapfbl(i,:)
          PRINT *, 'evapfbl before rescaling: ', canopy%evapfbl(i,:)  
!          STOP
        ELSE
          canopy%evapfbl(i,:) = oldevapfbl(i,:)
        END IF
      END IF
      END IF
    END DO
    canopy%frday = 12.0 * SUM(rdy, 2)
    canopy%fpn = -12.0 * REAL(SUM(an_y, 2),r_1)

  END SUBROUTINE dryLeaf

  SUBROUTINE photosynthesis(csxz,cx1z,cx2z,gswminz,rdxz,vcmxt3z,vcmxt4z, &
                    vx3z,vx4z,xleuningz,vlaiz,deltlfz,anxz)
    ! inputs:
    REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: csxz
    REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: cx1z
    REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: cx2z
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: gswminz
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: rdxz
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: vcmxt3z
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: vcmxt4z
    REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: vx4z
    REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: vx3z
    REAL(r_2), DIMENSION(mp,mf), INTENT(IN) :: xleuningz
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: vlaiz
    REAL(r_1), DIMENSION(mp,mf), INTENT(IN) :: deltlfz
    REAL(r_2), DIMENSION(mp,mf), INTENT(INOUT) :: anxz
    !local variables
    REAL(r_2), DIMENSION(mp,mf) :: coef0z,coef1z,coef2z
    REAL(r_2), DIMENSION(mp,mf) :: ciz,delcxz
    REAL(r_2), DIMENSION(mp,mf) :: anrubiscoz,anrubpz,ansinkz
    REAL(r_1), PARAMETER  :: effc4 = 4000.0  ! Vc=effc4*Ci*Vcmax (see
                                             ! Bonan,LSM version 1.0, p106)
    INTEGER(i_d) :: i, j
    REAL(r_1), PARAMETER  :: min_lai = 0.01

    ! rgswc    - inherited from canopy_module's USE photosynthetic_constants
    ! mp       - inherited from canopy_module
    ! mf       - inherited from canopy_module
   
   DO i = 1, mp
   IF (SUM(vlaiz(i,:)) > min_lai) THEN
   DO j = 1, mf
   IF (vlaiz(i,j) > 0.0 .AND. deltlfz(i,j) > 0.1) THEN
     ! Rubisco limited:
     coef2z(i,j) = gswminz(i,j)/rgswc+xleuningz(i,j) &
                 * (vcmxt3z(i,j)-(rdxz(i,j)-vcmxt4z(i,j)))
     coef1z(i,j) = (1.0-csxz(i,j)*xleuningz(i,j)) &
                 * (vcmxt3z(i,j)+vcmxt4z(i,j)-rdxz(i,j))      &
                 + (gswminz(i,j)/rgswc)*(cx1z(i,j)-csxz(i,j)) &
                 - xleuningz(i,j)*(vcmxt3z(i,j)*cx2z(i,j)/2.0 &
                 + cx1z(i,j)*(rdxz(i,j)-vcmxt4z(i,j)))
     coef0z(i,j) = -(1.0-csxz(i,j)*xleuningz(i,j)) &
                 * (vcmxt3z(i,j)*cx2z(i,j)/2.0  &
                 + cx1z(i,j)*(rdxz(i,j)-vcmxt4z(i,j))) &
                 - (gswminz(i,j)/rgswc)*cx1z(i,j)*csxz(i,j)
     IF (ABS(coef2z(i,j)) < 1.0e-9 .AND. ABS(coef1z(i,j)) < 1.0e-9) THEN
       ! no solution, give it a huge number
       ciz(i,j) = 99999.0   ! quadratic below cannot handle zero denominator
       anrubiscoz(i,j) = 99999.0    ! OR should ciz=0 and anrubiscoz calculated?
     END IF 
     IF  (ABS(coef2z(i,j)) < 1.e-9 .AND. ABS(coef1z(i,j)) >= 1.e-9) THEN
       ! solve linearly
       ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)   ! same reason as above
       ciz(i,j) = MAX(0.0_r_2,ciz(i,j))
       anrubiscoz(i,j) = vcmxt3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) &
                       / (ciz(i,j) + cx1z(i,j)) + vcmxt4z(i,j) - rdxz(i,j)
     END IF 
     IF  (ABS(coef2z(i,j)) >= 1.e-9) THEN
       ! solve quadratic (only take the more positive solution)
       delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)
       ciz(i,j)    = (-coef1z(i,j)+SQRT(MAX(0.0_r_2,delcxz(i,j)))) &
                   / (2.0*coef2z(i,j))
       ciz(i,j)    = MAX(0.0_r_2,ciz(i,j))   ! must be positive, why?
       anrubiscoz(i,j) = vcmxt3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) &
                       / (ciz(i,j) + cx1z(i,j)) + vcmxt4z(i,j) - rdxz(i,j)
     END IF 

   ! RuBP limited:
     coef2z(i,j) = gswminz(i,j)/rgswc+xleuningz(i,j) &
                 * (vx3z(i,j)-(rdxz(i,j)-vx4z(i,j)))
     coef1z(i,j) = (1.0-csxz(i,j)*xleuningz(i,j)) &
                 * (vx3z(i,j)+vx4z(i,j)-rdxz(i,j))    &
                 + (gswminz(i,j)/rgswc)*(cx2z(i,j)-csxz(i,j)) &
                 - xleuningz(i,j)*(vx3z(i,j)*cx2z(i,j)/2.0 &
                 + cx2z(i,j)*(rdxz(i,j)-vx4z(i,j)))
     coef0z(i,j) = -(1.0-csxz(i,j)*xleuningz(i,j)) *(vx3z(i,j)*cx2z(i,j)/2.0  &
                 + cx2z(i,j)*(rdxz(i,j)-vx4z(i,j))) &
                 - (gswminz(i,j)/rgswc)*cx2z(i,j)*csxz(i,j)
     IF (ABS(coef2z(i,j)) < 1.0e-9 .AND. ABS(coef1z(i,j)) < 1.0e-9) THEN
       ! no solution, give it a huge number
       ciz(i,j) = 99999.0
       anrubpz(i,j)  = 99999.0
     END IF 
     IF  (ABS(coef2z(i,j)) < 1.e-9 .AND. ABS(coef1z(i,j)) >= 1.e-9) THEN
       ! solve linearly
       ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)
       ciz(i,j)    = MAX(0.0_r_2,ciz(i,j))
       anrubpz(i,j) = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) &
                    / (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)
     END IF 
     IF  (ABS(coef2z(i,j)) >= 1.e-9) THEN
       ! solve quadratic (only take the more positive solution)
       delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)
       ciz(i,j)    = (-coef1z(i,j)+SQRT(MAX(0.0_r_2,delcxz(i,j)))) &
                   / (2.0*coef2z(i,j))
       ciz(i,j)    = MAX(0.0_r_2,ciz(i,j)) 
       anrubpz(i,j) = vx3z(i,j)*(ciz(i,j)-cx2z(i,j)/2.0) &
                    / (ciz(i,j)+cx2z(i,j)) +vx4z(i,j)-rdxz(i,j)
     END IF 

   ! Sink limited:
     coef2z(i,j) = xleuningz(i,j)
     coef1z(i,j) = gswminz(i,j)/rgswc + xleuningz(i,j) &
                 * (rdxz(i,j) - 0.5*vcmxt3z(i,j)) &
                 + effc4 * vcmxt4z(i,j) - xleuningz(i,j) &
                 * csxz(i,j) * effc4 *vcmxt4z(i,j)
     coef0z(i,j) = -(gswminz(i,j)/rgswc)*csxz(i,j) *effc4*vcmxt4z(i,j) +    &
                    (rdxz(i,j) -0.5*vcmxt3z(i,j))*gswminz(i,j)/rgswc
     IF (ABS(coef2z(i,j)) < 1.0e-9 .AND. ABS(coef1z(i,j)) < 1.0e-9) THEN
       ! no solution, give it a huge number
       ciz(i,j) = 99999.0
       ansinkz(i,j)  = 99999.0
     END IF 
     IF  (ABS(coef2z(i,j)) < 1.e-9 .AND. ABS(coef1z(i,j)) >= 1.e-9) THEN
       ! solve linearly
       ciz(i,j) = -1.0 * coef0z(i,j) / coef1z(i,j)
       ansinkz(i,j)  = ciz(i,j)
     END IF 
     IF  (ABS(coef2z(i,j)) >= 1.e-9) THEN
       ! solve quadratic (only take the more positive solution)
       delcxz(i,j) = coef1z(i,j)**2 -4.0*coef0z(i,j)*coef2z(i,j)
       ciz(i,j)    = (-coef1z(i,j)+SQRT(MAX(0.0_r_2,delcxz(i,j)))) &
                   / (2.0*coef2z(i,j))
       ansinkz(i,j) = ciz(i,j)
     END IF 

   ! minimal of three limited rates
     anxz(i,j) = MIN(anrubiscoz(i,j),anrubpz(i,j),ansinkz(i,j))
   END IF
   END DO
   END IF
   END DO

  END SUBROUTINE photosynthesis

  SUBROUTINE wetLeaf(dels,rad,air,met,canopy)
   ! assuming the temperature of wet leaf is equal that of dry leaf ="tlfy"
    TYPE (radiation_type),      INTENT(INOUT) :: rad
    TYPE (air_type),            INTENT(INOUT) :: air
    TYPE (met_type),            INTENT(INOUT) :: met
    TYPE (canopy_type),         INTENT(INOUT) :: canopy
    REAL(r_1), INTENT(IN)           :: dels ! integration time step (s)
    REAL(r_2), DIMENSION(mp)  :: ccfevw ! limitation term for
                                              ! wet canopy evaporation rate  
    REAL(r_2), DIMENSION(mp)  :: gwwet  ! cond for water for a wet canopy
    REAL(r_2), DIMENSION(mp)  :: ghrwet ! wet canopy cond: heat & thermal rad
    gwwet = 1.0e-3
    ghrwet= 1.0e-3
    canopy%fevw = 0.0
    canopy%fhvw = 0.0
    canopy%potev_c = 0.0  ! Added so that bareground values won't go haywire
    WHERE (canopy%vlaiw > 0.01)
    ! VEG SENSIBLE AND LATENT HEAT FLUXES fevw, fhvw (W/m2) for a wet canopy
    ! calculate total thermal resistance, rthv in s/m
      ghwet = 2.0   * SUM((gbhu+gbhf),2)
      gwwet = 1.075 * SUM((gbhu+gbhf),2)
      ghrwet = SUM(rad%gradis,2) + ghwet
      ! Calculate fraction of canopy which is wet:
      canopy%fwet = MAX(0.0,MIN(1.0,0.8*canopy%cansto/MAX(cansat,0.01)))
      ! Calculate lat heat from wet canopy, may be neg. if dew on wet canopy
      ! to avoid excessive evaporation:
      ccfevw = MIN(canopy%cansto * air%rlam / dels, &
                   2.0 / (1440.0 / (dels/60.0)) * air%rlam)
      canopy%fevw = MIN(canopy%fwet * (air%dsatdk * (SUM(rad%rniso,2) &
                  - capp*rmair*(met%tvair(:)-met%tk(:))*sum(rad%gradis,2)) &
                  + capp*rmair*met%dva*ghrwet) &
                  / (air%dsatdk+air%psyc*ghrwet/gwwet), ccfevw)
      ! canopy potential evapotranspiratn for output purposes (YP & Mao jun08)
      canopy%potev_c = MAX(0.0_r_2,(air%dsatdk*SUM(rad%rniso,2) &
           + capp*rmair*met%da*ghrwet)/(air%dsatdk+air%psyc*ghrwet/gwwet))
      ! Calculate sens heat from wet canopy:
      canopy%fhvw(:) = canopy%fwet(:) * (SUM(rad%rniso,2) - capp * rmair &
                  * (tlfy(:) - met%tk(:)) * SUM(rad%gradis,2)) - canopy%fevw(:)
    END WHERE
  END SUBROUTINE wetLeaf

END MODULE canopy_module
