c Purpose
c -------
c Applies freshwater hosing within the coupled model.
c
c Inputs
c ------
c RIMSL		Original flux of runoff into the ocean
c
c Outputs
c -------
c RIMSL		Modified flux of runoff into the ocean
c
c Notes
c -----
c The flux of water arising from the freshwater hosing is added to the
c atmosphere model array RIMSL, which contains the flux of runoff into the
c ocean. This approach allow for a neat and simple implementation.
c
c History
c -------
c 2009 Aug 5	Jess Trevena	Original version
c 2009 Aug 6	Steven Phipps	Integrated into the core model source code

      subroutine hose(rimsl)

      implicit none

C Global parameters
      include 'PARAMS.f'
      include 'OPARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      real rimsl(lon, 0:lat2)

C Global data blocks
      include 'GAUSL.f'
      include 'HOSING.f'
      include 'LSMI.f'
      include 'TIMEX.f'

C Local work arrays and variables
      integer i, ierr, j, lg, ma
      real area(lat2), area_fw

C Local data, functions etc
      logical first
      integer hosemask(lon, lat2)
      real conv_extra_fw
      data first /.true./
      save conv_extra_fw, first, hosemask

C Start code : ------------------------------------------------------------

c...  On the first call to this routine, initialise various arrays
      if (first) then

c...  Calculate the areas of the gridboxes at each latitude
        do j = 1, lat
          area(j) = 2.0 * pi * w(j) * eradsq / real(lon)
          area(lat2+1-j) = 2.0 * pi * w(j) * eradsq / real(lon)
        end do

c...  Read in the mask
        write (*, *)
        write (*, *) "Reading freshwater hosing mask :"
        write (*, *)
        write (*, *) "File = hosemask"
        write (*, *)
        open (unit=99, file="hosemask", form="formatted", status="old",
     &        iostat=ierr)
        if (ierr .ne. 0) then
          write (*, *)
          write (*, *) "***  Error: Could not open file"
          write (*, *)
          stop
        end if
        do j = lat2, 1, -1
          read (99, '(64i1)') (hosemask(i, j), i=1,lon)
        end do
        close (99)

c...  Echo the mask to standard output
        write (*, *)
        write (*, *) "Freshwater hosing mask :"
        write (*, *)
        do j = lat2, 1, -1
          write (*, '(64i1)') (hosemask(i, j), i=1,lon)
        end do
        write (*, *)

c...  Check that the mask is valid i.e. that no freshwater hosing is specified
c...  at any land gridpoints
        do j = 1, lat2
          if (j .gt. lat) then
            lg = lat2 + 1 - j
          else
            lg = j
          end if
          do i = 1, lon
            if (j .gt. lat) then
              ma = i
            else
              ma = i + lon
            end if
            if ((hosemask(i, j) .eq. 1) .and. (imsl(ma, lg) .eq. 4))
     &        then
              write (*, *)
              write (*, *) "***  Error in freshwater hosing mask :"
              write (*, *) "***"
              write (*, *) "***  The gridpoint"
              write (*, *) "***"
              write (*, *) "***  i = ", i
              write (*, *) "***  j = ", j
              write (*, *) "***"
              write (*, *)
     &          "***  is land, but freshwater hosing is specified."
              write (*, *)
              stop
            end if 
          end do
        end do

c...  Calculate the area over which freshwater hosing is specified
        area_fw = 0.0
        do j = 1, lat2
          do i = 1, lon
            if (hosemask(i, j) .eq. 1) area_fw = area_fw + area(j)
          end do
        end do

c...  Convert the specified hosing rate from Sv to m per timestep
        conv_extra_fw = 6.0e7 * real(mstep) * hosing_rate / area_fw

c...  Display some diagnostics
        write (*, *)
        write (*, *) "Freshwater hosing statistics :"
        write (*, *)
        write (*, *) "Area = ", area_fw/1.0e6, " km^2"
        write (*, *)
        write (*, *) "Rate = ", hosing_rate, " Sv"
        write (*, *) "     = ", conv_extra_fw, " m/timestep"
        write (*, *) "     = ", 525600.0*conv_extra_fw/real(mstep),
     &                          " m/year"
        write (*, *)

c...  Set FIRST to false, so that this section of code is only executed once
        first = .false.

      end if

c...  Apply the freshwater hosing
      do j = 1, lat2
        do i = 1, lon
          if (hosemask(i, j) .eq. 1) rimsl(i, j) = rimsl(i, j) +
     &                                             conv_extra_fw
        end do
      end do

      return
      end subroutine
