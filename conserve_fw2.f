c Purpose
c -------
c Corrects for the conservation error that arises when fluxes of water into the
c ocean are converted into surface salinity tendencies.
c
c The three components of the surface water flux (DSIS, DSISB and DSFW) are
c each integrated over the surface of the ocean, and then converted into the
c equivalent rate of change of salinity of the entire ocean.
c
c A uniform correction is then applied to the values of ATSF, in order to
c ensure that they give the correct rate of change of salinity of the entire
c ocean, thus conserving freshwater.
c
c Inputs
c ------
c DSIS		Flux of water arising from ice growth/melt
c DSISB		Flux of water arising from ice sublimation
c DSFW		Flux of water arising from (P-E+RUNOFF)
c ATSF		Initial values for the total surface salinity tendency
c
c Outputs
c -------
c ATSF		Corrected values for the total surface salinity tendency
c
c Notes
c -----
c (1) We assume that the average salinity of the ocean is 34.70 psu when
c     calculating the rate of change of salinity of the entire ocean, this
c     being the value according to the World Ocean Atlas 1998. In theory, we
c     should use the actual average salinity of the ocean when we are running
c     in coupled mode. However, this should never differ significantly from
c     34.70 psu, and so any error introduced should be utterly negligible.
c     Furthermore, this approach ensures consistency between the stand-alone
c     AGCM and the coupled model.
c
c History
c -------
c 2007 Dec 19	Steven Phipps	Original version
c 2009 Apr 21	Steven Phipps	Modified to get the volume of the ocean via the
c				COMMON block /COUPLING/

      subroutine conserve_fw2(dsis, dsisb, dsfw, atsf)

      implicit none

C Global parameters
      include 'PARAMS.f' 
      include 'PHYSPARAMS.f'
      include 'AOGCM.f'

C Argument list
      real dsis(ln2, lat), dsisb(ln2, lat), dsfw(ln2, lat),
     &     atsf(ln2, lat)

C Global data blocks
      include 'GAUSL.f'
      include 'LSMI.f'
      include 'TIMEX.f'

C Local work arrays and variables
      integer i, j
      real corr, dsdt1, dsdt2, total, total1, total2, total3

C Local data, functions etc
      logical first, masko(ln2, lat)
      real aocean, area(lat), dz, salice, saloce
      data first /.true./
      parameter (dz = 25.0)
      parameter (salice = 0.01)
      parameter (saloce = 0.0347)
      save aocean, area, first, masko

C Start code : ------------------------------------------------------------

c...  On the first call to this routine, initialise various arrays
      if (first) then

c...  Calculate the areas of the gridboxes at each latitude
        do j = 1, lat
          area(j) = 2.0 * pi * w(j) * eradsq / float(lon)
        end do

c...  Create a land/sea mask
        do j = 1, lat
          do i = 1, ln2
            if (imsl(i, j) .eq. 4) then
              masko(i, j) = .false.
            else
              masko(i, j) = .true.
            end if
          end do
        end do

c...  Calculate the surface area of the ocean
        aocean = 0.0
        do j = 1, lat
          do i = 1, ln2
            if (masko(i, j)) aocean = aocean + area(j)
          end do
        end do

c...  Set FIRST to FALSE so that this section of code is only executed once
        first = .false.

      end if

c...  Integrate the water fluxes, and the initial values for ATSF, over the
c...  surface of the ocean
      total = 0.0
      total1 = 0.0
      total2 = 0.0
      total3 = 0.0
      do j = 1, lat
        do i = 1, ln2
          if (masko(i, j)) then
            total = total + atsf(i, j) * area(j)
            total1 = total1 + dsis(i, j) * area(j)
            total2 = total2 + dsisb(i, j) * area(j)
            total3 = total3 + dsfw(i, j) * area(j)
          end if
        end do
      end do

c...  Calculate the equivalent rates of change in the salinity of the ocean
      dsdt1 = (total1 * (salice - saloce) + total2 * salice -
     &         total3 * saloce) / (dt * ovolume)
      dsdt2 = total * dz / ovolume

c...  Calculate the correction to apply to the initial values of ATSF
      corr = (dsdt1 - dsdt2) * ovolume / (aocean * dz)

c...  Correct the initial values of ATSF
      do j = 1, lat
        do i = 1, ln2
          if (masko(i, j)) then
            atsf(i, j) = atsf(i, j) + corr
          else
            atsf(i, j) = 0.0
          end if
        end do
      end do

c...  Write a summary to standard output
C      write (*, *)
C      write (*, *) "TOTAL1 = ", total1/dt, " m^3/s^-1"
C      write (*, *) "TOTAL2 = ", total2/dt, " m^3/s^-1"
C      write (*, *) "TOTAL3 = ", total3/dt, " m^3/s^-1"
C      write (*, *)
C      write (*, *) "TOTAL = ", total*dz, " m^3/s^-1"
C      write (*, *)
C      write (*, *) "DSDT1 = ", dsdt1, " s^-1"
C      write (*, *) "DSDT2 = ", dsdt2, " s^-1"
C      write (*, *)
C      write (*, *) "CORR = ", corr, " s^-1"
C      write (*, *)

      return
      end subroutine
