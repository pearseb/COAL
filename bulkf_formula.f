*   The parameterisations for bulk formula of shortwave, long-
*   wave, sensible and latent heat fluxes, and freshwater 
*   fluxes are taken from Stephanie Dutkiewicz implementation
*   of bulk formula in the MITgcm, which is open source.
*   
*   This option is a recommended alternative to restoring
*   towards a prescribed SST and SSS climatology (Haney, 1971).
*   Spurious multiple equilibria in oceanic circulation states
*   have been found from using restoring boundary conditions.
*   Thus, we highly recommend enabling bulk transfer formula
*   for heat and freshwater if the appropriate boundary
*   conditions are known. For more information regarding the
*   bulk formula presented here and their advantage over a
*   restoring approach see:
*      Large et al (1997) Journal of Physical Oceanography
*      Herold (2009) Paleocean modelling, in "Encyclopedia of
*          Paleoclimatology and Ancient Environments
*      Rahmstorf & Willebrand (1995) Journal of Physical
*          Oceanography
*   
*   The following conditions are required for heat fluxes:
*      1. Surface wind speeds (m/s)
*      2. Sea surface temperature (deg C) --> converted to K
*      3. Surface specific humidity (%)
*      4. Surface air temperature (K)
*      5. Net downwelling shortwave radiation at surface
*         assuming clear sky (W/m2)
*      6. Fractional cloud cover
*      7. Fractional sea ice cover
*      8. Downward longwave radiation (W/m2)
*
*  [ These inputs are in order or appearance in subroutine call ]
*   
*   Net heat from shortwave radiation     --> 5,6,7
*   Net heat from longwave radiation      --> 2,8
*   Net heat from sensible flux           --> 1,2,3,4,7
*   Net heat from latent flux             --> 1,2,3,4,7
*   
*   The following climatological boundary conditions are
*   required for freshwater transfer:
*      1. Evaporation (mg/m2/s)
*      2. Precipitation (mg/m2/s)
*      3. Runoff (mg/m2/s)
*      4. Sea surface salinity (psu)
*
*   Sea surface salinity is required because we retain a small
*   restorative term for surface salinity to ensure stability
*   in the global salt budget, such that if the intergral of
*   the global freshwater flux climatology is not equal to 0
*   then salinity will not drift indefinitely.
*   
*   IMPLEMENTED BY PEARSE J BUCHANAN, July 2018, while at Princeton
*      pearse.buchanan@utas.edu.au
*      pearseb@princeton.edu
*      pearse.buchanan@liverpool.ac.uk

      SUBROUTINE bulkf_formula ( 
     I             xwnd, ywnd, tos, hus, tas, swdn, fcloud, sico, lwdn,
     I             rai, sno, rof,
     O             swflxx, lwupx, lwflxx, senflxx, latflxx, 
     O             hflxx, fflxx, ustr, vstr,
     I             i,j, ncar, mit)

C ***    *** Adapted from MITgcm routine ***   ***
C
C     !DESCRIPTION: 
C     *==========================================================*
C     | SUBROUTINE BULKF_FORMULA_LAY
C     | o Calculate bulk formula fluxes over open ocean or seaice
C     |   Large and Yeager, 2004, NCAR/TN-460+STR.
C     *==========================================================*
C
C ***    *** Adapted from MITgcm routine ***   ***
C
C
C === Turbulent Fluxes :
C  * use the approach "B": shift coeff to height & stability of the
C    atmosphere state (instead of "C": shift temp & humid to the height
C    of wind, then shift the coeff to this height & stability of the atmos).
C  * similar to EXF (except over sea-ice) ; default parameter values
C    taken from Large & Yeager.
C  * assume that Qair & Tair inputs are from the same height (zq=zt)
C  * formulae in short:
C     wind stress = (ust,vst) = rhoA * Cd * Ws * (del.u,del.v)
C     Sensib Heat flux = fsha = rhoA * Ch * Ws * del.T * CpAir
C     Latent Heat flux = flha = rhoA * Ce * Ws * del.Q * Lvap
C                      = -Evap * Lvap
C   with Ws = wind speed = sqrt(del.u^2 +del.v^2) ;
C        del.T = Tair - Tsurf ; del.Q = Qair - Qsurf ;
C        Cd,Ch,Ce = drag coefficient, Stanton number and Dalton number
C              respectively [no-units], function of height & stability

C     !USES:
      IMPLICIT NONE

C     === Local constants ===
      include 'bulkf_const.f'

C     === Local input to subroutine (from tracer.F) ===
      real :: xwnd,     ! zonal wind velocity at 10 m [m/s]
     +        ywnd,     ! meridional wind velocity at 10 m [m/s]
     +        tos,      ! sea surface temperature [deg C]
     +        hus,      ! specific humidity at surface [kg/kg]
     +        tas,      ! air surface temperature [K]
     +        swdn,     ! shortwave radiation down [W/m2]
     +        fcloud,   ! fractional cloud coverage 
     +        sico,     ! fractional sea ice coverage
     +        lwdn,     ! longwave radiation down [W/m2]
     +        rai,      ! liquid precipitation [kg/(m2*s)]
     +        sno,      ! solid precipitation [kg/(m2*s)]
     +        rof       ! continental runoff [kg/(m2*s)]

C     === Local output of subroutine (to tracer.F) ===
      real :: swflxx,   ! net shortwave heat flux [W/m2]
     +        lwupx,    ! upwards longwave heat flux [W/m2]
     +        lwflxx,   ! net longwave heat flux [W/m2]
     +        senflxx,  ! sensible heat flux [W/m2]
     +        latflxx,  ! latent heat flux [W/m2]
     +        hflxx,    ! net heat flux [W/m2]
     +        fflxx,    ! freshwater flux [kg/m2/s]
     +        ustr,     ! wind stress in x direction
     +        vstr      ! wind stress in y direction

C     === Local variables ===
      real :: zwln,     ! log(measurement height/reference height)
     +        ztln,     ! log(measurement height/reference height)
     +        tsf,      ! sea surface temperature [K]
     +        czol,     ! height * karman constant * gravity
     +        alb,      ! albedo
     +        lath,     ! latent heat of vaporization or sublimation [J/kg]
     +        wnd,      ! wind velocity at 10 m height [m/s]
     +        wndm,     ! wind velocity at 10 m height [m/s] >= umin
     +        emiss,    ! emissivity of surface given ice concentration
     +        t0,       ! virtual air temperature [K]
     +        ssq,      ! saturated specific humidity at surface [kg/kg]
     +        dEvdT,    ! derivative of evap. with respect to tsf [kg/m2/s/K]
     +        deltap,   ! temperature difference between sea and air [K]
     +        delq,     ! specific humidity difference between sea and air [kg/kg]
     +        stable,   ! = 1 if stable ; = 0 if unstable
     +        tmpblk,   ! temporary value for use in multiple places
     +        rdn,      ! neutral transfer coefficient (momentum flux) [-]
     +        rhn,      ! neutral transfer coefficient (heat flux) [W/(m2*K)]
     +        ren,      ! neutral transfer coefficient (vapor flux) []
     +        ustar,    ! friction velocity [m/s]
     +        tstar,    ! temperature scale [K]
     +        qstar,    ! humidity scale [kg/kg]
     +        huol,     ! stability parameter at zwd (z/Monin-Obuklov length) [-]
     +        htol,     ! stability parameter at zth (z/Monin-Obuklov length) [-]
     +        xchi,     ! stability function [-] 
     +        xchi2,    ! xchi squared [-]
     +        psimh,    ! momentum stability function
     +        psixh,    ! latent and sensible heat stability function
     +        wnd_n,    ! neutral, zref (=10m) wind speed [m/s]
     +        wndm_n,   ! wnd_n but limited to minimum speed >= umin [m/s]
     +        rd,       ! = sqrt(Cd)          [-]
     +        rh,       ! = Ce / sqrt(Cd)     [-]
     +        re,       ! = Ch / sqrt(Cd)     [-]
     +        tau,      ! surface stress  coef = rhoA * Ws * Cd
     +        evp       ! freshwater flux associated with latent heating [kg/m2/s]

      integer :: i,j,iter

      logical :: ncar, mit



         ! Declare some useful variables
         zwln = log(zwd/zref)
         ztln = log(zth/zref)

         tsf = tos + 273.15
         czol = zref*karman*gravity  ! 39.20

         alb = alb_oce*(1.-sico) + alb_ice*sico
         
         lath = Lvap
         if (sico.gt.0.0) lath = Lvap+Lfresh
         
         wnd = sqrt(xwnd**2 + ywnd**2)
         wndm = max(wnd,umin)


c... 1. Calculate shortwave incoming radiation at surface (W/m2)
c           swflxx = swdown * ( 1 - albedo ) * ( 1 - cloud_fraction )
c           (potential to add zenith angle and this make albedo variable)
        swflxx = swdn*(1.-alb)*(1.-fcloud)

c... 2. Calculate longwave radiation exchange at surface (W/m2)
c           lwflxx = lwdown + lwup 
c               lwdown is provided as input to model in ocdatro.f 
c               lwup = ocean_emissivity * stefan_boltzmann * SST**4
        emiss = emi_oce*(1.-sico) + emi_ice*sico
        lwupx = (emiss * sboltz * tsf**4.) * (1.-sico)
        lwflxx = lwdn*emiss*(1.-sico) - lwupx

c... 3. Calculate sensible and latent radiation fluxes at surface (W/m2)
c           i) solve for virtual temperature, saturated specific humidity,
c              and differences in these variables between air and sea
c              surface
        t0 = tas * (1. + fhumid * hus)

        ssq = (saltQs * qs1w *exp(qs2w/tsf) / rhoA)*(1.-sico) + 
     +        (qs1i *exp(qs2i/tsf) / rhoA)*sico
        dEvdT = qs2w*(1.-sico) + qs2i*sico
        !ssq = ssq0*exp(lath*(ssq1-ssq2/tsf))/p0
        deltap = tas + gamma_blk*zth - tsf
        delq = hus - ssq


c        ii) initialise the exchange coefficients for sensible and
c            latent heat.
        stable = 0.5 + sign(0.5, deltap) ! stable if air > sea
        tmpblk = cdrag_1/wndm + cdrag_2 + cdrag_3*wndm
        ! depending on ice state, initialise neutral transfer coefficients
        rdn = (1.-sico)*sqrt(tmpblk) + sico*rdn_i
        rhn = (stable*cStantonS + (1.0-stable)*cStantonU)*(1.-sico) +
     +        sico*rdn_i
        ren = cDalton*(1.-sico) + sico*rdn_i

c        iii) calculate turbulent scales
        ustar = rdn*wndm
        tstar = rhn*deltap
        qstar = ren*delq

c        iv) iterate and find the transfer coefficients by finding if
c            conditions are stable or unstable, estimating the
c            integrated flux profile (psimh & psixh), new wind speed
c            and update the turbulent scales of momentum, heat and
c            water 
        do iter=1,niter_bulk
           
           huol = czol/ustar**2 * (tstar/t0 + 
     +            qstar/(1./fhumid + hus))
           huol = sign( min(abs(huol), 10.0), huol) ! -10 < huol < 10 
           xchi2 = sqrt( abs(1. - 16.*huol))
           xchi2 = max(xchi2, 1.0)
           xchi = sqrt(xchi2)
           
           ! NCAR method (suggested for COREv2 forcing) based on Large
           ! and Yeager 2004 
           !    see ncar_ocean_fluxes.f90 
           if (ncar) then
           if (huol > 0) then
             psimh = -5.0*huol
             psixh = -5.0*huol
           else
             psimh = log( (1.+2.*xchi+xchi2)*(1.+xchi2)*0.125 ) -
     +               2.0*(atan(xchi) - atan(1.0))
             psixh = 2.0*log(0.5*(1.+xchi2))
           endif
           endif

           stable = 5e-1 + sign(5e-1, huol) ! stable = 1, unstable = 0
           
           ! MITgcm method, written by Stephanie Dutkiewicz, based on
           ! Large and Yeager 2004
           !    However, DIFFERENT in calculation of stability
           !    parameters (psimh and psixh)
           if (mit) then
           psimh = -5.0*huol*stable + (1.-stable) *
     +             ( log( (1.+2.*xchi+xchi2)*(1.+xchi2)*0.125 ) -
     +               2.0*atan(xchi) + 3.14159265*0.5 )
           htol = huol*zth/zref
           xchi2 = sqrt( abs(1. - 16.*htol))
           psixh = -5.0*huol*stable +
     +             (1.-stable) * ( 2.0*log(0.5*(1.+xchi2))) 
           endif

           ! shift wind speed to 10 m height and neutral stability
           wnd_n = wnd/(1.0 + rdn/karman*(zwln-psimh))
           wndm_n = max(wnd_n,umin)

           ! update 10m , neutral stability transfer coefficients
           tmpblk = cdrag_1/wndm_n + cdrag_2 + cdrag_3*wndm_n
           rdn = sqrt(tmpblk)
           rhn = (stable*cStantonS + (1.0-stable)*cStantonU)*rdn
           ren = cDalton*rdn

           ! shift coefficients to the measurement height and stability
           rd = tmpblk/(1.0 + rdn*(zwln-psimh)/karman)**2
           rh = rhn/(1.0 + rhn*(ztln-psixh)/karman/rdn)*sqrt(rd/tmpblk)
           re = ren/(1.0 + rhn*(ztln-psixh)/karman/rdn)*sqrt(rd/tmpblk)

           ! for next iteration
           ustar = rd*wndm
           tstar = rh*deltap
           qstar = re*delq

        enddo

c       vi) calculate the sensible and latent heat fluxes using the
c           newly calculated transfer coefficients according to:
        tau = rhoA*rd*wnd
        ustr = rhoA*rd*abs(xwnd)*xwnd
        vstr = rhoA*rd*abs(ywnd)*ywnd
        senflxx = rhoA * cpair * rh * deltap * abs(wnd)
        latflxx = rhoA * lath * re * delq * abs(wnd)

c... 5. Solve for the total heat flux
        hflxx = swflxx+lwflxx+senflxx+latflxx

c... 6. Get the (mostly) evaporative flux from latent heat flux 
        evp = -latflxx/lath !kg/m2/s 
        ! and solve for freshwater flux in m/s
        fflxx = (evp - rai - sno - rof)/1025.0 ! density of freshwater


      RETURN
      END

















