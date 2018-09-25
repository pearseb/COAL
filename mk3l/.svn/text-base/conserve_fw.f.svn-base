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
c (1) CONSERVE must be called before the first call to this routine, in order
c     for AOCEANO, AREA and MASKO to be initialised.
c
c (2) We assume that the average salinity of the ocean is 34.70 psu when
c     calculating the rate of change of salinity of the entire ocean, this
c     being the value according to Levitus 1998. In theory, we should use the
c     actual average salinity of the ocean when we are running in coupled mode.
c     However, this should never differ significantly from 34.70 psu, and so
c     any error introduced should be utterly negligible. Furthermore, this
c     approach ensures consistency between the stand-alone AGCM and the coupled
c     model.
c
c History
c -------
c 2004 Sep 13	Steven Phipps	Original version
c 2009 Apr 21	Steven Phipps	Modified to get the volume of the ocean via the
c				COMMON block /COUPLING/

      subroutine conserve_fw(dsis, dsisb, dsfw, atsf)

      implicit none

C Global parameters
      include 'PARAMS.f' 
      include 'AOGCM.f'

C Argument list
      real dsis(ln2, lat), dsisb(ln2, lat), dsfw(ln2, lat),
     &     atsf(ln2, lat)

C Local shared common blocks
      logical masko
      real aoceano, area
      common /conservation/ aoceano, area(lat), masko(ln2, lat)

C Global data blocks
      include 'TIMEX.f'

C Local work arrays and variables
      integer i, j
      real corr, dsdt1, dsdt2, total, total1, total2, total3

C Local data, functions etc
      real dz
      real salice
      real saloce
      parameter (dz = 25.0)
      parameter (salice = 0.01)
      parameter (saloce = 0.0347)

C Start code : ------------------------------------------------------------

c Integrate the water fluxes, and the initial values for ATSF, over the surface
c of the ocean
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

c Calculate the equivalent rates of change in the salinity of the ocean
      dsdt1 = (total1 * (salice - saloce) + total2 * salice -
     &         total3 * saloce) / (dt * ovolume)
      dsdt2 = total * dz / ovolume

c Calculate the correction to apply to the initial values of ATSF
      corr = (dsdt1 - dsdt2) * ovolume / (aoceano * dz)

c Correct the initial values of ATSF
      do j = 1, lat
        do i = 1, ln2
          if (masko(i, j)) then
            atsf(i, j) = atsf(i, j) + corr
          else
            atsf(i, j) = 0.0
          end if
        end do
      end do

c Write a summary to standard output
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
