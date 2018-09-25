*   The following routine takes into account changes in surface
*   heat fluxes in high latitudes where melting and freezing of ice are
*   important for maintaining realistic surface temperatures.
*
*   Two scenarios are considered:
*       1. If cooling of surface waters is sufficient to lower 
*          temperatures below the freezing point, here set in 
*          bulkf_const.f as ocefrz = -1.85 deg C, then we solve
*          for the heating due to freezing (latent heat of fusion)
*          that is required to conserve surface temperatures == ocefrz
*       2. If warming of surface waters occurs when sea ice is
*          present, then we solve for the cooling due to melting 
*          that is required to conserve surface temperatures ~ ocefrz
*   
*   The following conditions are required as inputs:
*      1. Net heat fluxes [W/m2] (from bulkf_formula.f)
*      2. CTOP, which = RHO * Cw * Z(surface layer)
*      3. Change in time over one timestep for tracers [s]
*      4. Sea ice concentration 
*      5. Sea surface temperature [deg C]
*
*   IMPLEMENTED BY PEARSE J BUCHANAN, July 2018, while at Princeton
*      pearse.buchanan@utas.edu.au
*      pearseb@princeton.edu
*      pearse.buchanan@liverpool.ac.uk

      SUBROUTINE bulkf_seaice ( 
     I             hflxx, ctop2, dtime, sico, tos, tos_b,
     O             frzflxx, fwiflxx,
     I             pnux,i,j )


C     !USES:
      IMPLICIT NONE

C     === Local constants ===
      include 'bulkf_const.f'

C     === Local input to subroutine (from tracer.F) ===
      real :: hflxx,    ! net heat flux [W/m2]
     +        tos,      ! sea surface temperature [deg C] current timestep
     +        tos_b,    ! sea surface temperature [deg C] previous timestep
     +        sico,     ! fractional sea ice coverage
     +        ctop2,    ! RHO * Cw * Z(surface layer)
     +        dtime,    ! change in time [s]
     +        pnux      ! robert time filter coefficient "pnu" (default = 0.005)

C     === Local output to subroutine (to tracer.F) ===
      real :: frzflxx,  ! latent heat fluxes due to freezing and melting [W/m2]
     +        fwiflxx   ! latent heat fluxes due to freezing and melting [W/m2]

C     === Local variables ===
      real :: tos_a,   ! sea surface temperature [deg C] next timestep
     +        pnu2mx,   ! robert time filter coefficient "pnu" (default = 0.99)
     +        temp,    ! adjusted temperature using heat fluxes [degC]
     +        temp_new ! freezing point - surface temperature [degC]

      integer :: i,j

c... 1. Adjust temperature using heat fluxes calculated previously 
      ! Note that the adjusted temperature is an ESTIMATE of change for
      ! the next timestep by emulating the robert time filter. Because
      ! this estimate occurs prior to vertical diffusion, the
      ! temperature change is imperfect.
      tos_a = tos+hflxx/ctop2*dtime
      pnu2mx = 1.0 - 2.0*pnux
      temp = tos*pnu2mx + pnux*(tos_b + tos_a) ! see step.f (Robert time filtering)

!      if (i.eq.124 .and. j.eq.105) then
!      print*, "Call from bulkf_seaice.f"
!      print*, tos_b,tos,tos_a
!      print*, temp, ocefrz
!      print*, " "
!      endif

      if (temp.le.ocefrz .and. temp.lt.tos) then !freezing and cooling
c... 2. Calculate latent heating due to freezing
         temp_new = tos-temp ! positive delta T for warming
         frzflxx = temp_new * ctop2/dtime 
      ! Second, melting when warming occurs... 
      !   else 
      !      temp_new = tos-temp ! negative delta T for cooling
      !      frzflxx = temp_new * ctop2/dtime * sico
      else
         frzflxx = 0.0
      endif

c... 3. Account for freshwater fluxes from latent heat change
      fwiflxx = frzflxx/Lfresh ![kg/m2/s] 

      RETURN
      END
