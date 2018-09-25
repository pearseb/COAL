! $Log: bvnfcalc.f90,v $
! Revision 1.4  2001/10/12 02:06:57  rot032
! HBG changes, chiefly for new river-routing scheme and conservation
!
! Revision 1.1  2001/02/18 23:35:53  rot032
! Initial revision
!
! Revision 1.3  1999/05/20 06:23:53  rot032
! HBG changes to V5-2
!
! Revision 1.2  1997/12/17  23:22:52  ldr
! Changes from MRD for parallel execution on NEC.
!
! Revision 1.1  1996/12/23  03:53:08  mrd
! Initial revision
!

      subroutine bvnfcalc(tg,ttg,prf,dprf,pdpsk,bvnf)

! This routine calculates the Brunt-Vaisalla frequency at model full levels.
! It's needed to replace the calculation in gwdrag when the new gravity wave 
! drag scheme is used.

! The routine is written in fixed format f90

      implicit none

      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

      include 'CNSTA.f'
      include 'FEWFLAGS.f'

      real, intent(in),  dimension(ln2)    :: tg   ! Surface temperature
      real, intent(in),  dimension(ln2,nl) :: ttg  ! Atmospheric temperature
      real, intent(in),  dimension(ln2,nl) :: prf
      real, intent(in),  dimension(ln2,nl) :: dprf
      real, intent(in),  dimension(ln2,nl) :: pdpsk
      real, intent(out), dimension(ln2,nl) :: bvnf ! BV freq.

      real, parameter :: dzx=grav/rdry
      real, dimension(ln2,nl) :: dzi   ! Level thickness
      real, dimension(ln2,nl) :: thf   ! Potential temperature
      real, dimension(ln2,nl) :: dthdz ! d(theta) / dz
      integer :: k

      dzi = prf/dprf * dzx/ttg

! Calculate potential temperature
      thf = ttg/pdpsk

! Ensure real BV freq by forcing theta to increase with height
      thf(:,1)  = max ( thf(:,1), tg + 0.1 )
      do k=2,nl
	 thf(:,k)  = max ( thf(:,k), thf(:,k-1) + 0.1 )
      end do

! Calculate d(theta)/dz
      dthdz(:,1) = ( thf(:,1) - tg ) * dzi(:,1)
      do k=2,nl
	 dthdz(:,k) = ( thf(:,k) - thf(:,k-1) ) * dzi(:,k)
      end do

! Calculate full level BV freq by averaging the half level derivatives
      do k=1,nl-1
	 bvnf(:,k) = sqrt ( (0.5*grav)*(dthdz(:,k)+dthdz(:,k+1))/
     &                       thf(:,k) )
      end do
      bvnf(:,nl) = sqrt ( grav*dthdz(:,nl)/thf(:,nl) )

      return
      end subroutine bvnfcalc
