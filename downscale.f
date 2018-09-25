c Purpose
c -------
c Uses bilinear interpolation to downscale surface flux fields from the AGCM
c grid to the new, high-resolution OGCM grid.
c
c A globally-uniform correction is applied to the surface fluxes in order to
c ensure that the global integral is conserved.
c
c Inputs
c ------
c atmos		Values on the AGCM grid
c
c Outputs
c -------
c ocean		Values on the OGCM grid
c
c History
c -------
c 2007 Dec 18	Steven Phipps	Original version
c 2007 Dec 21	Steven Phipps	Modified to ensure conservation on a global
c                               basis, rather than a gridbox-by-gridbox basis

      subroutine downscale(atmos, ocean)

      implicit none

C Global parameters
      include 'PARAMS.f' 
      include 'OPARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      real atmos(ln2, lat), ocean(imt-2, jmt-2)

C Global data blocks
      include 'GAUSL.f'
      include 'LSMI.f'
      include 'ONEDIM.f'

C Local work arrays and variables
      logical masko(lon, lat2)
      real aflux(lon, lat2), areaa(lat2), areao(jmt-2)
      integer i, ia, im1, ip1, j, ja, jm1, jp1
      logical first
      real aoceana, aoceano, corr, error, sum, total1, total2, weight
      data first /.true./
      save aoceana, aoceano, areaa, areao, first, masko

C Start code : ------------------------------------------------------------

c...  On the first call to this routine, initialise various arrays
      if (first) then

c...  Calculate the areas of the AGCM gridboxes
        do j = 1, lat
          areaa(j) = 2.0 * pi * w(j) * eradsq / float(lon)
          areaa(lat2+1-j) = 2.0 * pi * w(j) * eradsq / float(lon)
        end do

c...  Calculate the areas of the OGCM gridboxes 
        do j = 1, jmt-2
          areao(j) = 1.0e-4 * cst(j+1) * dxt(1) * dyt(j+1)
        end do

c...  Create a land/sea mask
        do j = 1, lat2
          if (j .gt. lat) then
            ja = lat2 + 1 - j
          else
            ja = j
          endif

          do i = 1, lon
            if (j .gt. lat) then
              ia = i
            else
              ia = i + lon
            endif

            if (imsl(ia, ja) .eq. 4) then
              masko(i, j) = .false.
            else
              masko(i, j) = .true.
            end if

          end do
        end do

c...  Calculate the surface area of the ocean according to the AGCM
        aoceana = 0.0
        do j = 1, lat2
          do i = 1, lon
            if (masko(i, j)) aoceana = aoceana + areaa(j)
          end do
        end do

c...  Calculate the surface area of the ocean according to the OGCM
        aoceano = 0.0
        do j = 1, jmt-2
          ja = (j + 1) / 2
          do i = 1, imt-2
            ia = (i + 1) / 2
            if (masko(ia, ja)) aoceano = aoceano + areao(j)
          end do
        end do

c...  Display a summary
      write (*, *)
      write (*, *) "Ocean area according to AGCM = ",
     &             aoceana/1.0e6, " km^2"
      write (*, *) "Ocean area according to OGCM = ",
     &             aoceano/1.0e6, " km^2"
      write (*, *)

c...  Set FIRST to FALSE so that this section of code is only executed once
        first = .false.

      end if

c...  Transfer the AGCM surface fluxes to an array with a sensible shape
      do j = 1, lat2
        if (j .gt. lat) then
          ja = lat2 + 1 - j
        else
          ja = j
        endif

        do i = 1, lon
          if (j .gt. lat) then
            ia = i
          else
            ia = i + lon
          endif

          aflux(i, j) = atmos(ia, ja)

        end do
      end do

c...  Calculate the global integral of the AGCM surface fluxes
      total1 = 0.0
      do j = 1, lat2
        do i = 1, lon
          if (masko(i, j)) total1 = total1 + areaa(j) * aflux(i, j)
        end do
      end do

c...  Loop over all AGCM gridboxes
      do j = 1, lat2
        jm1 = j - 1
        jp1 = j + 1

        do i = 1, lon
          im1 = i - 1
          ip1 = i + 1
          if (i .eq. 1) im1 = lon
          if (i .eq. lon) ip1 = 1

          if (masko(i, j)) then

c...  This is an ocean gridbox, so use bilinear interpolation to calculate the
c...  surface flux on the OGCM grid
c...
c...  (1) The OGCM gridbox at (2*i-1, 2*j-1)
            sum = 9.0 * areaa(j) * aflux(i, j)
            weight = 9.0 * areaa(j)

            if (masko(im1, j)) then
              sum = sum + 3.0 * areaa(j) * aflux(im1, j)
              weight = weight + 3.0 * areaa(j)
            end if

            if (j .gt. 1) then

              if (masko(i, jm1)) then
                sum = sum + 3.0 * areaa(jm1) * aflux(i, jm1)
                weight = weight + 3.0 * areaa(jm1)
              end if

              if (masko(im1, jm1)) then
                sum = sum + areaa(jm1) * aflux(im1, jm1)
                weight = weight + areaa(jm1)
              end if

            end if

            ocean(2*i-1, 2*j-1) = sum / weight

c...  (2) The OGCM gridbox at (2*i, 2*j-1)
            sum = 9.0 * areaa(j) * aflux(i, j)
            weight = 9.0 * areaa(j)

            if (masko(ip1, j)) then
              sum = sum + 3.0 * areaa(j) * aflux(ip1, j)
              weight = weight + 3.0 * areaa(j)
            end if

            if (j .gt. 1) then

              if (masko(i, jm1)) then
                sum = sum + 3.0 * areaa(jm1) * aflux(i, jm1)
                weight = weight + 3.0 * areaa(jm1)
              end if

              if (masko(ip1, jm1)) then
                sum = sum + areaa(jm1) * aflux(ip1, jm1)
                weight = weight + areaa(jm1)
              end if

            end if

            ocean(2*i, 2*j-1) = sum / weight

c...  (3) The OGCM gridbox at (2*i-1, 2*j)
            sum = 9.0 * areaa(j) * aflux(i, j)
            weight = 9.0 * areaa(j)

            if (masko(im1, j)) then
              sum = sum + 3.0 * areaa(j) * aflux(im1, j)
              weight = weight + 3.0 * areaa(j)
            end if

            if (j .lt. lat2) then

              if (masko(i, jp1)) then
                sum = sum + 3.0 * areaa(jp1) * aflux(i, jp1)
                weight = weight + 3.0 * areaa(jp1)
              end if

              if (masko(im1, jp1)) then
                sum = sum + areaa(jp1) * aflux(im1, jp1)
                weight = weight + areaa(jp1)
              end if

            end if

            ocean(2*i-1, 2*j) = sum / weight

c...  (4) The OGCM gridbox at (2*i, 2*j)
            sum = 9.0 * areaa(j) * aflux(i, j)
            weight = 9.0 * areaa(j)

            if (masko(ip1, j)) then
              sum = sum + 3.0 * areaa(j) * aflux(ip1, j)
              weight = weight + 3.0 * areaa(j)
            end if

            if (j .lt. lat2) then

              if (masko(i, jp1)) then
                sum = sum + 3.0 * areaa(jp1) * aflux(i, jp1)
                weight = weight + 3.0 * areaa(jp1)
              end if

              if (masko(ip1, jp1)) then
                sum = sum + areaa(jp1) * aflux(ip1, jp1)
                weight = weight + areaa(jp1)
              end if

            end if

            ocean(2*i, 2*j) = sum / weight

          else

c...  This is a land gridbox, so set the surface flux on the OGCM grid to zero
            ocean(2*i-1, 2*j-1) = 0.0
            ocean(2*i, 2*j-1) = 0.0
            ocean(2*i-1, 2*j) = 0.0
            ocean(2*i, 2*j) = 0.0

          end if

        end do
      end do

c...  Calculate the global integral of the OGCM surface fluxes
      total2 = 0.0
      do j = 1, jmt-2
        ja = (j + 1) / 2
        do i = 1, imt-2
          ia = (i + 1) / 2
          if (masko(ia, ja)) total2 = total2 + areao(j) * ocean(i, j)
        end do
      end do

c...  Calculate the conservation error, and correct the OGCM surface fluxes
      error = total2 - total1
      corr = error / aoceano
      do j = 1, jmt-2
        ja = (j + 1) / 2
        do i = 1, imt-2
          ia = (i + 1) / 2
          if (masko(ia, ja)) ocean(i, j) = ocean(i, j) - corr
        end do
      end do

c...  Display a summary
c      write (*, *)
c      write (*, *) "TOTAL1 = ", total1
c      write (*, *) "TOTAL2 = ", total2
c      write (*, *)
c      write (*, *) "ERROR = ", error
c      write (*, *)
c      write (*, *) "CORR = ", 0.0-corr
c      write (*, *)

      return
      end subroutine
