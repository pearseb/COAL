c New version of OCEND, containing the end-of-timestep work carried out within
c the coupled model only. The end-of-month work has been transferred to a new
c subroutine, OCEMON.
c SJP 2009/05/11

      subroutine ocend(nato)
C
C=======================================================================
C   Routine to do end-of-timestep work within the coupled model      ===
C=======================================================================
C
      implicit none

C Global parameters
      include 'OPARAMS.f'

C Argument list
      integer nato

C Global data blocks
      include 'A2O.f'    
      include 'CRELAX.f'

C Local work arrays and variables
      integer i, j

c...  Calculate the average sea surface fluxes for this timestep.
      do j = 1, jmt-2
        do i = 1, imt-2
          otx(i, j) = otx(i, j) / float(nato)
          oty(i, j) = oty(i, j) / float(nato)
          osalf(i, j) = osalf(i, j) / float(nato)
          osurh(i, j) = osurh(i, j) / float(nato)
          osrad(i, j) = osrad(i, j) / float(nato)
          oswnd(i, j) = oswnd(i, j) / float(nato)
          osice(i, j) = osice(i, j) / float(nato)

        end do
      end do

c...  If CRELAX_FLAG is true, add the surface fluxes to the running total for
c...  this month.
      if (crelax_flag) then
        do j = 2, jmtm1
          do i = 2, imtm1
            crelax_stfht(i, j) = crelax_stfht(i, j) +
     &                           crelax_qnf(i, j)
            crelax_stfsal(i, j) = crelax_stfsal(i, j) +
     &                            crelax_salf(i, j)
          end do
        end do
      end if

      RETURN
      END
