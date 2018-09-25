c Modified for the conversion of the atmosphere model auxiliary files to
c netCDF. In the process, these files have been renamed as follows:
c
c   albnew21f/albnew42f.dat2/sibalbT63.dat4     ->  albedo.nc
c   clim3f.sss/clim242f.sss/clim3T63f.sss       ->  sssa.nc
c   clim3f.sst/clim242f.sst/clim3T63f.sst       ->  ssta.nc
c   icedivl.R21/icedivl.R42/icedivl.T63         ->  icedivl.nc
c   ocuv.3st/ocuv.242/ocuv.T63                  ->  ocuv.nc
c   psrk21f.dat/psrk42f.dat/psrkT63f.dat3       ->  psrk.nc
c   sibrs.dat/sibrs42.dat2/sibrsT63.dat3        ->  sibrs.nc
c   sibsig.dat/sibsig42.dat2/sibsigT63.dat3     ->  sibsig.nc
c   sibsoil.dat/sibsoilr42.dat/sibsoilT63.dat3  ->  sibsoil.nc
c   sibvegt.dat/sibvegt42.dat2/sibvegtT63.dat4  ->  sibvegt.nc
c   sibz0.dat/sibz042.dat2/sibz0T63.dat4        ->  sibz0.nc
c
c The units of some variables have also been changed as follows:
c
c   ocuv.nc :   cm/s     ->  m/s
c   ssta.nc :   minus K  ->  K
c
c SJP 2008/02/06
c
c Modified to read the ocean currents from a single input file, rather than
c having a separate input file for each month.
c SJP 2004/09/22
c
c $Log: flatset.f,v $
c Revision 1.20  2001/02/22 05:34:37  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.19  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.18  1998/12/10  00:55:30  ldr
c HBG changes to V5-1-21
c
c Revision 1.17  1998/05/27  02:07:35  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.16  1996/10/24  01:02:45  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.15  1996/03/21  03:18:43  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.14  1994/08/08  17:21:19  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.13  94/03/30  12:34:25  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.12  93/12/17  11:51:36  ldr
c Changes to V4-4-45l from HBG for coupled model
c 
c Revision 1.11  93/11/29  11:38:27  ldr
c Changes to V4-4-32l from HBG for coupled model
c 
c Revision 1.10  93/10/05  13:06:15  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.9  93/08/13  15:21:03  ldr
c Corrected name of second ocuv file for R42.
c 
c Revision 1.8  93/07/20  10:01:24  ldr
c New ocuv files for R42 (242 suffix).
c 
c Revision 1.7  93/07/14  17:36:33  ldr
c Corrected name of 2nd ocuv file.
c 
c Revision 1.6  93/07/12  15:26:27  ldr
c Use new ocean current files ocuv??.2st. Still using old R42 files.
c 
c Revision 1.5  93/07/06  16:37:52  ldr
c      changes for setting all values in common blocks
c
c     INPUT/OUTPUT:
c     Input:   from arguments
c                  mnth - month counter
c
c     Output:  from common/ocwind in COMDICE.f
c                  fav - average ice concentration
c                  uav - average zonal velocities
c                  vav - average meridional velocities
c                  wd,wl cumulative arrays for ice divergence
c
c              from common/red in this subroutine
c                  groice - new ice growth
c                  redice - ice redistribution
c                  hav - advective ice
c
c     In/Out:  from common/ocwind in COMDICE.f
c                  ocu - ocean currents, zonal
c                  ocv - ocean currents, meridional
c
c 
      subroutine flatset(mnth)

c To get ocean model average velocities for driving the ice model.
c The ocu,ocv arrays are now dimensioned as for the ocean model :
c (1:plon+2,1:plat+2) which is different to the ice model arrays.

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'COMDICE.f'

C Argument list
      integer mnth

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      real redice,groice,hav
      common/red/redice(ln2,lat),groice(ln2,lat),hav(ln2,lat)

C Local work arrays and variables
      integer i, j, lg, mg, mnthp1
      character title*60

C Local data, functions etc
      logical got_data
      real ocu_save(lon, lat2, 12), ocv_save(lon, lat2, 12)

      data got_data /.false./
      save got_data, ocu_save, ocv_save

C Start code : ----------------------------------------------------------

c...  Read the ocean current data, if this has not already been done
      if (.not. got_data) then

        write (*, *)
        write (*, *) "Reading zonal ocean current data :"
        write (*, *)
        write (*, *) "File     =  ocuv.nc"
        write (*, *) "Variable =  ocu"

        call read_real_3d("ocuv.nc", "ocu", lon, lat2, 12, ocu_save,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        write (*, *)
        write (*, *) "Reading meridional ocean current data :"
        write (*, *)
        write (*, *) "File     =  ocuv.nc"
        write (*, *) "Variable =  ocv"

        call read_real_3d("ocuv.nc", "ocv", lon, lat2, 12, ocv_save,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_data = .true.

      end if

c...  Transfer the data to the arrays OCU and OCV, and pad in the x-direction
      mnthp1 = mnth + 1
      if (mnthp1 .eq. 13) mnthp1 = 1
      do j = 2, plat+1
        do i = 2, plon+1
          ocu(i, j, 1) = ocu_save(i-1, j-1, mnth)
          ocu(i, j, 2) = ocu_save(i-1, j-1, mnthp1)
          ocv(i, j, 1) = ocv_save(i-1, j-1, mnth)
          ocv(i, j, 2) = ocv_save(i-1, j-1, mnthp1)
        end do
        ocu(1, j, 1) = ocu(plon+1, j, 1)
        ocu(1, j, 2) = ocu(plon+1, j, 2)
        ocv(1, j, 1) = ocv(plon+1, j, 1)
        ocv(1, j, 2) = ocv(plon+1, j, 2)
        ocu(plon+2, j, 1) = ocu(2, j, 1)
        ocu(plon+2, j, 2) = ocu(2, j, 2)
        ocv(plon+2, j, 1) = ocv(2, j, 1)
        ocv(plon+2, j, 2) = ocv(2, j, 2)
      end do

      return

      entry flatset1
        do 520 j = 1,plat    
        do 520 i = 1,plon    
        fav(i,j)=0.0
        wd(i,j)=0.0
        wl(i,j)=0.0
        uav(i,j)=0.0
  520   vav(i,j)=0.0
        do 521 lg=1,lat
        do 521 mg=1,ln2
          redice(mg,lg)=0.0
          groice(mg,lg)=0.0
  521        hav(mg,lg)=0.0
        return
        end
