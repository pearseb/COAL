c Added bulk formula functionality to ocean GCM
c Pearse J Buchanan 16/07/2018
c
c Modified for the conversion of the OGCM auxiliary files to netCDF. In the
c process, these files have been renamed as follows:
c   sss.dat    -> sss.nc
c   sst.dat    -> sst.nc
c   stress.dat -> stress.nc
c SJP 2008/02/03
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/17
c
c Added IMPLICIT NONE statement, plus variable declarations as necessary.
c SJP 2007/05/30
c
c Transferred COMMON block declarations to separate header files.
c SJP 2007/05/29
c
c Modified to enable relaxation of the coupled model SSTs and SSSs towards
c prescribed values.
c SJP 2006/01/05
c
c Split off from step.f; redundant code deleted.
c SJP 2003/12/30

      subroutine ocdatro

c---- Routine to read 12 months of data for ocean model
c---- Select 2 relevant months for interpolation.

      implicit none

      include 'OPARAMS.f'
      include 'SFC1.f'
      include 'LEVD.f'
      include 'SSTSAL.f'
      include 'FEWFLAGS.f'
      include 'OCEAN_NML.f'

      character title*60
      integer i, j, k
      logical got_sss, got_sst, got_stress, got_spechum, got_airtem,
     &        got_slp, got_swdown, got_lwdown, got_fclt, got_rain,
     &        got_snow, got_roff, got_fice, got_uwnd, got_vwnd 
      real sss_save(imt-2, jmt-2, 12), sst_save(imt-2, jmt-2, 12),
     &     strx_save(imt-2, jmt-2, 12), stry_save(imt-2, jmt-2, 12),
     &     spechum_save(imt-2, jmt-2, 12),airtem_save(imt-2, jmt-2, 12),
     &     slp_save(imt-2, jmt-2, 12), swdown_save(imt-2, jmt-2, 12), 
     &     lwdown_save(imt-2, jmt-2, 12), fclt_save(imt-2, jmt-2, 12), 
     &     rain_save(imt-2, jmt-2, 12), snow_save(imt-2, jmt-2, 12), 
     &     roff_save(imt-2, jmt-2, 12), fice_save(imt-2, jmt-2, 12), 
     &     uwnd_save(imt-2, jmt-2, 12), vwnd_save(imt-2, jmt-2, 12)

      data got_sss /.false./
      data got_sst /.false./
      data got_stress /.false./
      data got_spechum /.false./
      data got_airtem /.false./
      data got_slp /.false./
      data got_swdown /.false./
      data got_lwdown /.false./
      data got_fclt /.false./
      data got_rain /.false./
      data got_snow /.false./
      data got_roff /.false./
      data got_fice /.false./
      data got_uwnd /.false./
      data got_vwnd /.false./

      save got_sss, got_sst, got_stress, got_spechum, got_airtem,
     &     got_slp, got_swdown, got_lwdown, got_fclt, got_rain,
     &     got_snow, got_roff, got_fice, got_uwnd, got_vwnd,
     &     sss_save, sst_save, strx_save, stry_save, spechum_save,
     &     airtem_save, slp_save, swdown_save, lwdown_save, fclt_save,
     &     rain_save, snow_save, roff_save, fice_save, uwnd_save, 
     &     vwnd_save

c---- Read in data from files holding 12 months of data.
c---- Put data for months lmon1,lmon2 into global arrays
c---- with last index locations of 1,2 respectively.

c...  ***************************************
c...  ***  READ SURFACE WIND STRESS DATA  ***
c...  ***************************************
      if (.not. lcouple) then

c...  Read the surface wind stress data, if this has not already been done
      if (.not. bulk_force) then
      if (.not. got_stress) then

        write (*, *)
        write (*, *) "Reading zonal wind stress data :"
        write (*, *)
        write (*, *) "File     =  stress.nc"
        write (*, *) "Variable =  strx"

        call read_real_3d("stress.nc", "strx", imt-2, jmt-2, 12,
     &                    strx_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        write (*, *)
        write (*, *) "Reading meridional wind stress data :"
        write (*, *)
        write (*, *) "File     =  stress.nc"
        write (*, *) "Variable =  stry"

        call read_real_3d("stress.nc", "stry", imt-2, jmt-2, 12,
     &                    stry_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_stress = .true.

      end if

c...  Transfer the surface wind stress data to the arrays STRX and STRY, and
c...  pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          strx(i, j, 1) = strx_save(i-1, j-1, lmon1)
          strx(i, j, 2) = strx_save(i-1, j-1, lmon2)
          stry(i, j, 1) = stry_save(i-1, j-1, lmon1)
          stry(i, j, 2) = stry_save(i-1, j-1, lmon2)
        end do
        strx(1, j, 1) = strx(imtm1, j, 1)
        strx(1, j, 2) = strx(imtm1, j, 2)
        strx(imt, j, 1) = strx(2, j, 1)
        strx(imt, j, 2) = strx(2, j, 2)
        stry(1, j, 1) = stry(imtm1, j, 1)
        stry(1, j, 2) = stry(imtm1, j, 2)
        stry(imt, j, 1) = stry(2, j, 1)
        stry(imt, j, 2) = stry(2, j, 2)
      end do

      endif    ! .not. bulk_force
      end if   ! .not. lcouple

c...  *******************************************
c...  ***  READ SEA SURFACE TEMPERATURE DATA  ***
c...  *******************************************

c...  Read the SST data, if this has not already been done
      if (.not. got_sst) then

        write (*, *)
        write (*, *) "Reading sea surface temperature data :"
        write (*, *)
        write (*, *) "File     =  sst.nc"
        write (*, *) "Variable =  sst"

        call read_real_3d("sst.nc", "sst", imt-2, jmt-2, 12, sst_save,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_sst = .true.

      end if

c...  Transfer the SST data to the array SSTM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          sstm(i, j, 1) = sst_save(i-1, j-1, lmon1)
          sstm(i, j, 2) = sst_save(i-1, j-1, lmon2)
        end do
        sstm(1, j, 1) = sstm(imtm1, j, 1)
        sstm(1, j, 2) = sstm(imtm1, j, 2)
        sstm(imt, j, 1) = sstm(2, j, 1)
        sstm(imt, j, 2) = sstm(2, j, 2)
      end do

c...  ****************************************
c...  ***  READ SEA SURFACE SALINITY DATA  ***
c...  ****************************************

c...  Read the SSS data, if this has not already been done, subtract 35 psu and
c...  convert to absolute units
      if (.not. got_sss) then

        write (*, *)
        write (*, *) "Reading sea surface salinity data :"
        write (*, *)
        write (*, *) "File     =  sss.nc"
        write (*, *) "Variable =  sss"

        call read_real_3d("sss.nc", "sss", imt-2, jmt-2, 12, sss_save,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_sss = .true.

        do k = 1, 12
          do j = 1, jmt-2
            do i = 1, imt-2
              sss_save(i, j, k) = (sss_save(i, j, k) - 35.0) / 1000.0
            end do
          end do
        end do

      end if

c...  Transfer the SSS data to the array SALM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          salm(i, j, 1) = sss_save(i-1, j-1, lmon1)
          salm(i, j, 2) = sss_save(i-1, j-1, lmon2)
        end do
        salm(1, j, 1) = salm(imtm1, j, 1)
        salm(1, j, 2) = salm(imtm1, j, 2)
        salm(imt, j, 1) = salm(2, j, 1)
        salm(imt, j, 2) = salm(2, j, 2)
      end do

c...  ****************************************
c...  ***  READ 10 m WIND VELOCITY DATA  ***
c...  ****************************************

c...  Read the WIND data, if this has not already been done
      if (.not. got_uwnd) then

        write (*, *)
        write (*, *) "Reading 10 m zonal wind data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  uwnd"

        call read_real_3d("forcings.nc","uwnd",imt-2,jmt-2,12,uwnd_save,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_uwnd = .true.

        do k = 1, 12
          do j = 1, jmt-2
            do i = 1, imt-2
              uwnd_save(i, j, k) = uwnd_save(i, j, k)
            end do
          end do
        end do

      end if

c...  Transfer the UWND data to the array UWNDM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          uwndm(i, j, 1) = uwnd_save(i-1, j-1, lmon1)
          uwndm(i, j, 2) = uwnd_save(i-1, j-1, lmon2)
        end do
        uwndm(1, j, 1) = uwndm(imtm1, j, 1)
        uwndm(1, j, 2) = uwndm(imtm1, j, 2)
        uwndm(imt, j, 1) = uwndm(2, j, 1)
        uwndm(imt, j, 2) = uwndm(2, j, 2)
      end do

      if (.not. got_vwnd) then

        write (*, *)
        write (*, *) "Reading 10 m zonal wind data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  vwnd"

        call read_real_3d("forcings.nc","vwnd",imt-2,jmt-2,12,vwnd_save,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_vwnd = .true.

        do k = 1, 12
          do j = 1, jmt-2
            do i = 1, imt-2
              vwnd_save(i, j, k) = vwnd_save(i, j, k)
            end do
          end do
        end do

      end if

c...  Transfer the vwnd data to the array vwndM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          vwndm(i, j, 1) = vwnd_save(i-1, j-1, lmon1)
          vwndm(i, j, 2) = vwnd_save(i-1, j-1, lmon2)
        end do
        vwndm(1, j, 1) = vwndm(imtm1, j, 1)
        vwndm(1, j, 2) = vwndm(imtm1, j, 2)
        vwndm(imt, j, 1) = vwndm(2, j, 1)
        vwndm(imt, j, 2) = vwndm(2, j, 2)
      end do

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if (bulk_force .or. nt.gt.2) then

C... If bulk formula for heat and freshwater fluxes are enabled, read in
C      the required variables from file forcings.nc
      
      if (.not. got_spechum) then

        write (*, *)
        write (*, *) "Reading surface specific humidity data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  spechum"

        call read_real_3d("forcings.nc","spechum",imt-2,jmt-2,12,
     &                    spechum_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_spechum = .true.

      end if

c...  Transfer the spechum data to the array SPECHUMM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          spechumm(i, j, 1) = spechum_save(i-1, j-1, lmon1)
          spechumm(i, j, 2) = spechum_save(i-1, j-1, lmon2)
        end do
        spechumm(1, j, 1) = spechumm(imtm1, j, 1)
        spechumm(1, j, 2) = spechumm(imtm1, j, 2)
        spechumm(imt, j, 1) = spechumm(2, j, 1)
        spechumm(imt, j, 2) = spechumm(2, j, 2)
      end do
      
      
      if (.not. got_airtem) then

        write (*, *)
        write (*, *) "Reading surface air temperature data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  airtem"

        call read_real_3d("forcings.nc","airtem",imt-2,jmt-2,12,
     &                    airtem_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_airtem = .true.

      end if

c...  Transfer the airtem data to the array AIRTEMM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          airtemm(i, j, 1) = airtem_save(i-1, j-1, lmon1)
          airtemm(i, j, 2) = airtem_save(i-1, j-1, lmon2)
        end do
        airtemm(1, j, 1) = airtemm(imtm1, j, 1)
        airtemm(1, j, 2) = airtemm(imtm1, j, 2)
        airtemm(imt, j, 1) = airtemm(2, j, 1)
        airtemm(imt, j, 2) = airtemm(2, j, 2)
      end do
      
      if (.not. got_slp) then

        write (*, *)
        write (*, *) "Reading surface air temperature data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  slp"

        call read_real_3d("forcings.nc","slp",imt-2,jmt-2,12,
     &                    slp_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_slp = .true.

      end if

c...  Transfer the slp data to the array slpM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          slpm(i, j, 1) = slp_save(i-1, j-1, lmon1)
          slpm(i, j, 2) = slp_save(i-1, j-1, lmon2)
        end do
        slpm(1, j, 1) = slpm(imtm1, j, 1)
        slpm(1, j, 2) = slpm(imtm1, j, 2)
        slpm(imt, j, 1) = slpm(2, j, 1)
        slpm(imt, j, 2) = slpm(2, j, 2)
      end do
      
      if (.not. got_swdown) then

        write (*, *)
        write (*, *) "Reading downward shortwave radiation data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  swdown"

        call read_real_3d("forcings.nc","swdown",imt-2,jmt-2,12,
     &                    swdown_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_swdown = .true.

      end if

c...  Transfer the swdown data to the array SWDOWNM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          swdownm(i, j, 1) = swdown_save(i-1, j-1, lmon1)
          swdownm(i, j, 2) = swdown_save(i-1, j-1, lmon2)
        end do
        swdownm(1, j, 1) = swdownm(imtm1, j, 1)
        swdownm(1, j, 2) = swdownm(imtm1, j, 2)
        swdownm(imt, j, 1) = swdownm(2, j, 1)
        swdownm(imt, j, 2) = swdownm(2, j, 2)
      end do
      
      
      if (.not. got_lwdown) then

        write (*, *)
        write (*, *) "Reading downward longwave radiation data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  lwdown"

        call read_real_3d("forcings.nc","lwdown",imt-2,jmt-2,12,
     &                    lwdown_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_lwdown = .true.

      end if

c...  Transfer the lwdown data to the array LWDOWNM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          lwdownm(i, j, 1) = lwdown_save(i-1, j-1, lmon1)
          lwdownm(i, j, 2) = lwdown_save(i-1, j-1, lmon2)
        end do
        lwdownm(1, j, 1) = lwdownm(imtm1, j, 1)
        lwdownm(1, j, 2) = lwdownm(imtm1, j, 2)
        lwdownm(imt, j, 1) = lwdownm(2, j, 1)
        lwdownm(imt, j, 2) = lwdownm(2, j, 2)
      end do
      
      
      if (.not. got_fclt) then

        write (*, *)
        write (*, *) "Reading fractional cloud cover data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  fclt"

        call read_real_3d("forcings.nc","fclt",imt-2,jmt-2,12,
     &                    fclt_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_fclt = .true.

      end if

c...  Transfer the fclt data to the array FCLTM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          fcltm(i, j, 1) = fclt_save(i-1, j-1, lmon1)
          fcltm(i, j, 2) = fclt_save(i-1, j-1, lmon2)
        end do
        fcltm(1, j, 1) = fcltm(imtm1, j, 1)
        fcltm(1, j, 2) = fcltm(imtm1, j, 2)
        fcltm(imt, j, 1) = fcltm(2, j, 1)
        fcltm(imt, j, 2) = fcltm(2, j, 2)
      end do
      
      
      if (.not. got_rain) then

        write (*, *)
        write (*, *) "Reading precipitation data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  rain"

        call read_real_3d("forcings.nc","rain",imt-2,jmt-2,12,
     &                    rain_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_rain = .true.

      end if

c...  Transfer the rain data to the array RAINM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          rainm(i, j, 1) = rain_save(i-1, j-1, lmon1)
          rainm(i, j, 2) = rain_save(i-1, j-1, lmon2)
        end do
        rainm(1, j, 1) = rainm(imtm1, j, 1)
        rainm(1, j, 2) = rainm(imtm1, j, 2)
        rainm(imt, j, 1) = rainm(2, j, 1)
        rainm(imt, j, 2) = rainm(2, j, 2)
      end do
      
      
      if (.not. got_snow) then

        write (*, *)
        write (*, *) "Reading runoff data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  snow"

        call read_real_3d("forcings.nc","snow",imt-2,jmt-2,12,
     &                    snow_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_snow = .true.

      end if

c...  Transfer the snow data to the array snowM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          snowm(i, j, 1) = snow_save(i-1, j-1, lmon1)
          snowm(i, j, 2) = snow_save(i-1, j-1, lmon2)
        end do
        snowm(1, j, 1) = snowm(imtm1, j, 1)
        snowm(1, j, 2) = snowm(imtm1, j, 2)
        snowm(imt, j, 1) = snowm(2, j, 1)
        snowm(imt, j, 2) = snowm(2, j, 2)
      end do
      
      
      if (.not. got_roff) then

        write (*, *)
        write (*, *) "Reading runoff data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  roff"

        call read_real_3d("forcings.nc","roff",imt-2,jmt-2,12,
     &                    roff_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_roff = .true.

      end if

c...  Transfer the roff data to the array ROFFM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          roffm(i, j, 1) = roff_save(i-1, j-1, lmon1)
          roffm(i, j, 2) = roff_save(i-1, j-1, lmon2)
        end do
        roffm(1, j, 1) = roffm(imtm1, j, 1)
        roffm(1, j, 2) = roffm(imtm1, j, 2)
        roffm(imt, j, 1) = roffm(2, j, 1)
        roffm(imt, j, 2) = roffm(2, j, 2)
      end do
      
      
      if (.not. got_fice) then

        write (*, *)
        write (*, *) "Reading fice data :"
        write (*, *)
        write (*, *) "File     =  forcings.nc"
        write (*, *) "Variable =  fice"

        call read_real_3d("forcings.nc","fice",imt-2,jmt-2,12,
     &                    fice_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_fice = .true.

      end if

c...  Transfer the fice data to the array ficeM, and pad in the x-direction
      do j = 2, jmtm1
        do i = 2, imtm1
          ficem(i, j, 1) = min(1.0, max(0.0, fice_save(i-1, j-1,lmon1)))
          ficem(i, j, 2) = min(1.0, max(0.0, fice_save(i-1, j-1,lmon2)))
        end do
        ficem(1, j, 1) = ficem(imtm1, j, 1)
        ficem(1, j, 2) = ficem(imtm1, j, 2)
        ficem(imt, j, 1) = ficem(2, j, 1)
        ficem(imt, j, 2) = ficem(2, j, 2)
      end do
      
      endif ! --> bulk_force
      
      RETURN
      END
