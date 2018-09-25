c Modified for the conversion of the coupled model auxiliary files to netCDF.
c In the process, these files have been renamed as follows:
c   dtm1av      -> dtm.nc
c   hfcor.dat12 -> hfcor.nc
c   s1coravth   -> ssscor.nc
c   sfcor.dat12 -> sfcor.nc
c   t1coravth   -> sstcor.nc
c   txcor.dat12 -> txcor.nc
c   tycor.dat12 -> tycor.nc
c SJP 2008/02/05
c
c Added IMPLICIT NONE statement, and replaced declarations of IMT and JMT with
c the header file OPARAMS.f.
c SJP 2007/05/30
c
c Re-inserted wind stress adjustment files txcor.dat12 and tycor.dat12.
c SJP 2004/11/18
c
c Split off from step.f; redundant variable declarations deleted.
c SJP 2003/12/30

      subroutine ocdatra

c---- Routine to read 2 months of data for coupled model

      implicit none

      include 'OPARAMS.f'
      include 'LEVD.f'
      include 'FCOR.f'

      character title*60
      integer i, j
      logical got_data
      real hfcor_save(imt-2, jmt-2, 12), sfcor_save(imt-2, jmt-2, 12),
     &     txcor_save(imt-2, jmt-2, 12), tycor_save(imt-2, jmt-2, 12)

      data got_data /.false./
      save got_data, hfcor_save, sfcor_save, txcor_save, tycor_save

c---- Read in data from files holding 12 months of flux corrections.
c---- Put data for months lmon1,lmon2 into global arrays
c---- with last index locations of 1,2 respectively.

c...  ***********************************
c...  ***  READ FLUX ADJUSTMENT DATA  ***
c...  ***********************************

c...  Read the flux adjustment data, if this has not already been done
      if (.not. got_data) then

c...  Surface heat flux adjustments
        write (*, *)
        write (*, *) "Reading surface heat flux adjustment data :"
        write (*, *)
        write (*, *) "File     =  hfcor.nc"
        write (*, *) "Variable =  hfcor"

        call read_real_3d("hfcor.nc", "hfcor", imt-2, jmt-2, 12,
     &                    hfcor_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

c...  Surface salinity tendency adjustments
        write (*, *)
        write (*, *)
     &    "Reading surface salinity tendency adjustment data :"
        write (*, *)
        write (*, *) "File     =  sfcor.nc"
        write (*, *) "Variable =  sfcor"

        call read_real_3d("sfcor.nc", "sfcor", imt-2, jmt-2, 12,
     &                    sfcor_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

c...  Zonal wind stress adjustments
        write (*, *)
        write (*, *) "Reading zonal wind stress adjustment data :"
        write (*, *)
        write (*, *) "File     =  txcor.nc"
        write (*, *) "Variable =  txcor"

        call read_real_3d("txcor.nc", "txcor", imt-2, jmt-2, 12,
     &                    txcor_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

c...  Meridional wind stress adjustments
        write (*, *)
        write (*, *) "Reading meridional wind stress adjustment data :"
        write (*, *)
        write (*, *) "File     =  tycor.nc"
        write (*, *) "Variable =  tycor"

        call read_real_3d("tycor.nc", "tycor", imt-2, jmt-2, 12,
     &                    tycor_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_data = .true.

      end if

c...  Transfer the data to the arrays HFCOR, SFCOR, TXCOR and TYCOR. It is not
c...  necessary to pad in the x-direction, as the values for I=1 and I=IMT are
c...  not used.
      do j = 2, jmtm1
        do i = 2, imtm1
          hfcor(i, j, 1) = hfcor_save(i-1, j-1, lmon1)
          hfcor(i, j, 2) = hfcor_save(i-1, j-1, lmon2)
          sfcor(i, j, 1) = sfcor_save(i-1, j-1, lmon1)
          sfcor(i, j, 2) = sfcor_save(i-1, j-1, lmon2)
          txcor(i, j, 1) = txcor_save(i-1, j-1, lmon1)
          txcor(i, j, 2) = txcor_save(i-1, j-1, lmon2)
          tycor(i, j, 1) = tycor_save(i-1, j-1, lmon1)
          tycor(i, j, 2) = tycor_save(i-1, j-1, lmon2)
        end do
      end do

      return
      end
