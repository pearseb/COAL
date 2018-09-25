c Purpose
c -------
c Ensures that quantities are conserved when fields are passed between the AGCM
c and the OGCM.
c
c Values on the AGCM and OGCM grids are integrated either over land or over the
c ocean, and a correction applied to the values on the target grid in order to
c ensure that either the global integral or mean (as appropriate) is conserved.
c
c Inputs
c ------
c source	Values on the source grid
c target	Initial values on the target grid
c type		Type of conservation operation to perform:
c
c		1 = Conserve global integral AGCM ocean -> OGCM ocean
c		2 = Conserve global integral AGCM land -> OGCM land
c		3 = Conserve global mean OGCM ocean -> AGCM ocean
c
c Outputs
c -------
c source	Corrected values on the target grid
c
c Notes
c -----
c (1) We use the array SOURCE for both input and output, in order to save
c     memory.
c
c (2) Over the ocean, a uniform correction is added to/subtracted from the data
c     in order to achieve conservation. Over land, however, the data is
c     multiplied by a correction factor. These methods are the most suitable
c     for the fields which are currently being conserved - heat/freshwater
c     fluxes and SST over the ocean, and runoff over land.
c
c History
c -------
c 2003 Apr 28	Steve Phipps	Original version
c 2003 May 29	Steve Phipps	Slight modifications, including an error fix
c 2003 Jun 9	Steve Phipps	Another error fix - subroutine has now been
c                               tested and definitely works OK!
c 2003 Jun 10	Steve Phipps	Added the ability to conserve over land, as
c                               well as over the ocean
c 2004 Sep 11	Steve Phipps	Moved AOCEANO, AREA and MASKO to /CONSERVATION/
c				so that the values can be used by CONSERVE_FW
c 2004 Nov 18	Steve Phipps	Generalised in order to support conservation of
c				global-mean SST

      subroutine conserve(source, target, type)
      implicit none

C Global parameters
      include 'PARAMS.f' 
      include 'PHYSPARAMS.f'

C Argument list
      integer type
      real source(ln2, lat), target(ln2, lat)

C Local shared common blocks
      logical masko
      real aoceano, area
      common /conservation/ aoceano, area(lat), masko(ln2, lat)

C Global data blocks
      include 'GAUSL.f'
      include 'LSMI.f'

C Local work arrays and variables
      integer i, j, nland, nocean
      real corr, totala, totalo

C Local data, functions etc
      logical first, maska(ln2, lat)
      real alanda, alando, aoceana
      data first /.true./
      save alanda, alando, aoceana, first, maska

C Start code : ------------------------------------------------------------

c On the first call to this routine, set up the AGCM and OGCM masks. These are
c logical arrays, with ocean points set to .TRUE. and land points to .FALSE.
      if (first) then

c Calculate the area of gridboxes at each latitude
        do j = 1, lat
          area(j) = 2.0 * pi * w(j) * eradsq / float(lon)
        end do

c AGCM mask - get data from array IMSL. We also set the OGCM mask equal to the
c             AGCM mask at this point
        do j = 1, lat
          do i = 1, ln2
            if (imsl(i, j) .eq. 4) then
              maska(i, j) = .false.
              masko(i, j) = .false.
            else
              maska(i, j) = .true.
              masko(i, j) = .true.
            end if
          end do
        end do

c OGCM mask - derive this by adjusting the AGCM mask accordingly
c
c (1) 11 points are treated as land by the AGCM, but as ocean by the OGCM
        masko(3, 4) = .true.
        masko(62, 8) = .true.
        masko(23, 26) = .true.
        masko(84, 27) = .true.
        masko(116, 12) = .true.
        masko(118, 8) = .true.
        masko(86, 28) = .true.
        masko(87, 28) = .true.
        masko(90, 28) = .true.
        masko(91, 28) = .true.
        masko(91, 27) = .true.

c (2) 18 points are treated as ocean by the AGCM, but as land by the OGCM
        masko(91, 16) = .false.
        masko(72, 21) = .false.
        masko(72, 22) = .false.
        masko(72, 23) = .false.
        masko(91, 25) = .false.
        masko(67, 27) = .false.
        masko(85, 28) = .false.
        masko(20, 28) = .false.
        masko(21, 28) = .false.
        masko(19, 27) = .false.
        masko(19, 26) = .false.
        masko(19, 25) = .false.
        masko(23, 18) = .false.
        masko(64, 17) = .false.
        masko(52, 9) = .false.
        masko(53, 9) = .false.
        masko(35, 8) = .false.
        masko(52, 5) = .false.

c Derive and display a summary of each mask - keep the surface areas
c according to each grid, as we will need these values later
        write (*, *)
        write (*, *) "AGCM land/sea mask"
        write (*, *) "------------------"
        nland = 0
        nocean = 0
        alanda = 0.0
        aoceana = 0.0
        do j = 1, lat
          do i = 1, ln2
            if (maska(i, j)) then
              nocean = nocean + 1
              aoceana = aoceana + area(j)
            else
              nland = nland + 1
              alanda = alanda + area(j)
            end if
          end do
        end do
        write (*, *) "Number of land grid points  = ", nland
        write (*, *) "Number of ocean grid points = ", nocean
        write (*, *) "Land area  = ", alanda/1.0e6, " km^2"
        write (*, *) "Ocean area = ", aoceana/1.0e6, " km^2"
        write (*, *) "Ocean fraction = ", aoceana / (alanda + aoceana)
        write (*, *)
        write (*, *) "OGCM land/sea mask"
        write (*, *) "------------------"
        nland = 0
        nocean = 0
        alando = 0.0
        aoceano = 0.0
        do j = 1, lat
          do i = 1, ln2
            if (masko(i, j)) then
              nocean = nocean + 1
              aoceano = aoceano + area(j)
            else
              nland = nland + 1
              alando = alando + area(j)
            end if
          end do
        end do
        write (*, *) "Number of land grid points  = ", nland
        write (*, *) "Number of ocean grid points = ", nocean
        write (*, *) "Land area  = ", alando/1.0e6, " km^2"
        write (*, *) "Ocean area = ", aoceano/1.0e6, " km^2"
        write (*, *) "Ocean fraction = ", aoceano / (alando + aoceano)
        write (*, *)

c Set FIRST to .FALSE. so that this section of code is only executed once
        first = .false.

      end if

c Perform the desired conservation operation
      if (type .eq. 1) then

c (1) Conserve the global integral over the surface of the ocean, for fields
c     passed from the AGCM to the OGCM.
c
c (1.1) Integrate the total fluxes on both the AGCM and OGCM grids
        totala = 0.0
        totalo = 0.0
        do j = 1, lat
          do i = 1, ln2
            if (maska(i, j)) totala = totala + source(i, j) * area(j)
            if (masko(i, j)) totalo = totalo + target(i, j) * area(j)
          end do
        end do
C        write (*, *)
C        write (*, *) "Total value on AGCM grid = ", totala
C        write (*, *) "Total value on OGCM grid = ", totalo
C        write (*, *)
C        write (*, *) "Mean value on AGCM grid = ", totala / aoceana
C        write (*, *) "Mean value on OGCM grid = ", totalo / aoceano
C        write (*, *)

c (1.2) Derive the correction to apply to the fluxes on the OGCM grid
        corr = (totala - totalo) / aoceano
C        write (*, *) "Correction = ", corr
C        write (*, *)

c (1.3) Correct the values on the OGCM grid, putting the corrected values into
c       the array SOURCE
        do j = 1, lat
          do i = 1, ln2
            if (masko(i, j)) then
              source(i, j) = target(i, j) + corr
            else
              source(i, j) = 0.0
            end if
          end do
        end do

      else if (type .eq. 2) then

c (2) Conserve the global integral over land, for fields passed from the AGCM
c     to the OGCM.
c
c (2.1) Integrate the total fluxes on both the AGCM and OGCM grids
        totala = 0.0
        totalo = 0.0
        do j = 1, lat
          do i = 1, ln2
            if (.not. maska(i, j)) totala = totala +
     $                                        source(i, j) * area(j)
            if (.not. masko(i, j)) totalo = totalo +
     $                                        target(i, j) * area(j)
          end do
        end do
C        write (*, *)
C        write (*, *) "Total value on AGCM grid = ", totala
C        write (*, *) "Total value on OGCM grid = ", totalo
C        write (*, *)
C        write (*, *) "Mean value on AGCM grid = ", totala / alanda
C        write (*, *) "Mean value on OGCM grid = ", totalo / alando
C        write (*, *)

c (2.2) Derive the correction to apply to the fluxes on the OGCM grid
        corr = totala / totalo
C        write (*, *) "Correction factor = ", corr
C        write (*, *)

c (2.3) Correct the values on the OGCM grid, putting the corrected values into
c       the array SOURCE
        do j = 1, lat
          do i = 1, ln2
            if (.not. masko(i, j)) then
              source(i, j) = corr * target(i, j)
            else
              source(i, j) = 0.0
            end if
          end do
        end do

      else if (type .eq. 3) then

c (3) Conserve the global mean over the surface of the ocean, for fields
c     passed from the OGCM to the AGCM.
c
c (3.1) Integrate the total fluxes on both the AGCM and OGCM grids
        totalo = 0.0
        totala = 0.0
        do j = 1, lat
          do i = 1, ln2
            if (masko(i, j)) totalo = totalo + source(i, j) * area(j)
            if (maska(i, j)) totala = totala + target(i, j) * area(j)
          end do
        end do
C        write (*, *)
C        write (*, *) "Total value on OGCM grid = ", totalo
C        write (*, *) "Total value on AGCM grid = ", totala
C        write (*, *)
C        write (*, *) "Mean value on OGCM grid = ", totalo / aoceano
C        write (*, *) "Mean value on AGCM grid = ", totala / aoceana
C        write (*, *)

c (3.2) Derive the correction to apply to the values on the AGCM grid
        corr = totalo / aoceano - totala / aoceana
C        write (*, *) "Correction = ", corr
C        write (*, *)

c (3.3) Correct the values on the AGCM grid, putting the corrected values into
c       the array SOURCE
        do j = 1, lat
          do i = 1, ln2
            if (maska(i, j)) then
              source(i, j) = target(i, j) + corr
            else
              source(i, j) = 0.0
            end if
          end do
        end do

      else

c The value of TYPE that has been passed is illegal

        write (*, *)
        write (*, *) "Type = ", type
        write (*, *)
        write (*, *) "STOP: Illegal value for type"
        write (*, *)
        stop

      end if

      return
      end subroutine
