c Re-implementing the "limited free slip" boundary condition at R21 and R42
c resolutions. The new implementation is simpler, and corrects the bugs that
c were present in the previous version.
c SJP 2009/08/04
c
c Tidying up the source code for the sea ice model, particularly: fixing the
c bug in the calculation of the latitudes; and removing the buggy and 
c experimental code which sought to impose "limited slip at low resolution".
c SJP 2009/04/25
c
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
c $Log: icesetup.f,v $
c Revision 1.20  2001/02/22 05:34:46  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.19  1999/05/20 06:23:58  rot032
c HBG changes to V5-2
c
c Revision 1.18  1998/12/10  00:56:03  ldr
c HBG changes to V5-1-21
c
c Revision 1.17  1997/12/19  02:03:20  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.16  1996/10/24  01:02:56  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.15  1996/03/21  03:18:51  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.14  1995/08/31  04:30:45  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.13  1995/07/26  07:28:38  ldr
c Merge Mickles speedups to ice scheme (V4-5-30mic) into V4-7-3l.
c
c Revision 1.12  1994/08/08  17:21:35  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.11.1.1  1995/07/26  07:24:08  ldr
c Speedups from Mickles replacing character mask with integers.
c
c Revision 1.11  94/03/30  12:34:34  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.10  93/12/17  15:32:54  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.9  93/10/05  13:06:25  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.8  93/07/12  14:27:52  ldr
c Merge of SPO's changes to V4-3-10l with other changes since then.
c 
c Revision 1.7  93/07/06  16:38:08  ldr
c      slat(0) and slat(plat+1) moved to correct spot (HBG).
c 
c Revision 1.6.1.1  93/07/12  14:13:58  ldr
c Minor changes from SPO for coupled model V4-4.
c 
c Revision 1.6  93/04/26  15:03:07  ldr
c Removed error message for Svalbard fix at R42 resolution.
c 
c Revision 1.5  93/03/12  10:09:32  ldr
c Minor SPO changes from V4-2.
c 
c Revision 1.4  92/12/09  14:43:36  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.3  92/09/01  12:04:59  ldr
c Special fix for points around Svalbard.
c 
c Revision 1.2  92/07/16  12:52:29  ldr
c Minor changes to dynamical seaice model to correct runs not starting in
c January, make general for R42 and remove need for 'mask' file.
c 
c Revision 1.1  92/06/16  11:54:56  ldr
c Initial revision
c
c     INPUT/OUTPUT:
c     Input:   from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c     In/Out:  from common/dicegrid in COMDICE.f
c                  alatu - latitudes on u-grid
c                  alonu - longitudes on u-grid
c                  alat - latitude (from -pi/2 to pi/2)
c                  deltxu - east-west widths of u-grid boxes
c                  deltyu - nl-south widths of u-grid boxes
c                  deltx - east-west widths of h-grid boxes
c                  delty - nl-south widths of h-grid boxes
c                  slat, slatu - sine latitude scales
c                  srfarea - surface area of grid boxes
c
c              from common/dicemask in COMDICE.f
c                  imask - mask of land surface type(imsl)
c                  umask - u-grid land-ocean mask 
c
c              from common/diceturn in COMDICE.f
c                  coswturn - water stress turning angle
c
c     Output:  from common/dicegrid in COMDICE.f
c                  alon  - longitudes on h-grid
c                  blatu - latitudes on u-grid (alatu) in degrees 
c                  blonu - longitudes on u-grid (alonu) in degrees
c                  clonu - cos longitudes on u-grid (alonu)
c                  slonu - sin longitudes on u-grid (alonu)
c
c              from common/dicemask in COMDICE.f
c                  dmask - h-grid "divergence" mask
c
c              from common/diceturn in COMDICE.f
c                  sinwturn - water stress turning angle
c 
c 
      subroutine icesetup(kk)

c Initialization for dynamical sea ice, called once at (re)start.

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer kk

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'COMDICE.f'
      include 'GRADNS.f'
      include 'LSMI.f'

C Local work arrays and variables
      integer i
      integer j
      integer lg

      real a
      real asrat
      real c
      real ca
      real s
      real sa
      real twopi

      character title*60

C Local data, functions etc

C Start code : ----------------------------------------------------------

      IF(kk.eq.1)THEN

      twopi = 2.*pi 

c Set grid metrics. deltx,delty are e-w and nl-s widths of h-grid boxes,
c and deltxu,deltyu are e-w and nl-s widths of u-grid boxes. 

      do 10 i=0,plon+1
        alonu(i) =  (i-.5)*twopi/plon
        blonu(i) = alonu(i)*180./pi
        alon(i) = (i-1.)*twopi/plon
        slonu(i) = sin(alonu(i))
        clonu(i) = cos(alonu(i))
   10 continue

      do lg = 1, lat
        alat(lg) = glats(lg)
        alat(lat2+1-lg) = -glats(lg)
      end do

      alat(0) = 2*alat(1) - alat(2) 
      alat(plat+1) = 2*alat(plat) - alat(plat-1)

      do 12 j=1,plat
        slat(j) = sin(alat(j))
        deltx(j)  = erad * cos(alat(j)) * (twopi/plon)
        delty(j)  = erad * abs(alat(j+1)-alat(j))
   12 continue

      slat(plat+1) = 2*slat(plat)-slat(plat-1)
      slat(0) =2*slat(1)-slat(2) 
      deltx(0) = deltx(1)
      delty(0) = delty(1)
      deltx(plat+1) = deltx(plat)
      delty(plat+1) = delty(plat)

      do 14 j=1,plat+1
        alatu(j) = (alat(j)+alat(j-1))*0.5   
        blatu(j)= alatu(j)*180./pi
        slatu(j) = sin(alatu(j))
        deltxu(j) = erad * cos(alatu(j)) * (twopi/plon)
  14    continue
        alatu(plat+2)=alat(plat+1)+0.5*(alat(plat+1)-alat(plat))
        do 15 j=1,plat+1
        deltyu(j) = erad * abs(alatu(j+1)-alatu(j))
   15 continue

c Set the surface area of each h-grid box

      do 100 j=1,plat
        srfarea(j) = abs(slatu(j+1)-slatu(j))*(twopi/plon)*erad*erad
  100 continue

c Set the surface area of each u-grid box, including triangular pole 
c areas

      do 120 j=2,plat
        srfareau(j) =abs(slat(j)-slat(j-1))*(twopi/plon) * erad*erad
  120 continue

      srfareau(1) = abs(slat(1)-slatu(1))*(twopi/plon) * erad*erad
      srfareau(plat+1) = abs(slatu(plat+1)-slat(plat))*(twopi/plon)
     &                   * erad*erad

c Set the water turning parameters on u-grid

      c = cos(wturn*pi/180.)
      s = sin(wturn*pi/180.)
      ca = cos(aturn*pi/180.)
      sa = sin(aturn*pi/180.)
      do 200 j=1,plat+1
        coswturn(j) = c
        sinwturn(j) = s
        cosaturn(j) = ca
        sinaturn(j) = -sa
        if (alatu(j).lt.0.) sinwturn(j) = -s
        if (alatu(j).lt.0.) sinaturn(j) = sa
  200 continue

      ELSE ! (kk=2)

c Set u-grid land-ocean mask: 1 for ocean points, 0 for any point
c on the corner of a land h-grid box.

      do 300 j=1,plat+1
      do 300 i=0,plon+1
        umask(i,j) = 1.
  300 continue

c Set up imask array which runs from South pole to north

      do j=1,lat
        do i=1,lon
          imask(i,j)=imsl(i+lon,j)
          imask(i,lat2+1-j)=imsl(i,j)
        enddo
      enddo

      do 310 j=1,plat
        do 312 i=1,plon
          if(imask(i,j).gt.3) then
            umask(i,j) = 0.
            umask(i+1,j) = 0.
            umask(i,j+1) = 0.
            umask(i+1,j+1) = 0.
          endif
  312   continue
  310 continue

c At this point, umask on the Greenwich meridian (i=1 on u-grid) has
c only been affected by the eastern h-grid box. Remedy this here.
c Also set umask wraparound.

      do 320 j=1,plat+1
        umask(1,j) = min (umask(1,j), umask(plon+1,j))
  320 continue

      do 330 j=1,plat+1
        umask(0,j) = umask(plon,j)
        umask(plon+1,j) = umask(1,j)
  330 continue

c...  Set up a logical mask of coastal gridpoints on the u-grid. CMASK is true
c...  if a gridpoint is land and at least one of the neighbouring gridpoints is
c...  ocean; it is false otherwise. This mask is used to apply the "limited
c...  free slip" boundary condition at R21 and R42 resolutions.
c...
c...  1. Initialise the mask
      do j = 1, plat+1
        do i = 0, plon+1
          cmask(i, j) = .false.
        end do
      end do

c...  2. Set the mask to true for all coastal gridpoints
      do j = 2, plat
        do i = 1, plon
          if (umask(i, j) .lt. 0.5) then
            if (umask(i, j+1) .gt. 0.5) cmask(i, j) = .true.
            if (umask(i+1, j+1) .gt. 0.5) cmask(i, j) = .true.
            if (umask(i+1, j) .gt. 0.5) cmask(i, j) = .true.
            if (umask(i+1, j-1) .gt. 0.5) cmask(i, j) = .true.
            if (umask(i, j-1) .gt. 0.5) cmask(i, j) = .true.
            if (umask(i-1, j-1) .gt. 0.5) cmask(i, j) = .true.
            if (umask(i-1, j) .gt. 0.5) cmask(i, j) = .true.
            if (umask(i-1, j+1) .gt. 0.5) cmask(i, j) = .true.
          end if
        end do
      end do

c Set h-grid "divergence" mask, for use in cavit or cavit2. If for 
c cavit (flato=.true.), it is the divergence corresponding to unit
c outward velocities at all non-masked corners of an h-grid box, except
c that the u-velocities are weighted by asrat, the mean x-to-y aspect
c ratio of the grid box (see comment in cavit). If for cavit2
c (flato=.false.), it is the divergence produced by solving icefree
c with unit pressure pres inside an h-grid box, with pres outside and
c all other forcings zero. 

      if (flato) then

        do 400 j=1,plat
          asrat = 0.5*(deltxu(j)+deltxu(j+1)) / delty(j)
          do 402 i=1,plon
            dmask(i,j) = 
     &        (   deltxu(j)  *0.5*(umask(i,j)  +umask(i+1,j))
     &          + deltxu(j+1)*0.5*(umask(i,j+1)+umask(i+1,j+1))
     &          + delty(j)   *0.5*(umask(i,j)  +umask(i,j+1))  *asrat
     &          + delty(j)   *0.5*(umask(i+1,j)+umask(i+1,j+1))*asrat
     &        ) / srfarea(j)
  402     continue
  400   continue

      else

        do 410 j=1,plat
          a = csubw * 0.5*(coswturn(j)+coswturn(j+1))
          do 412 i=1,plon
            dmask(i,j) = 
     &        (   0.25*delty(j)*umask(i,j)        /deltx(j) 
     &          + 0.25*delty(j)*umask(i,j+1)      /deltx(j) 
     &          + 0.25*delty(j)*umask(i+1,j)      /deltx(j) 
     &          + 0.25*delty(j)*umask(i+1,j+1)    /deltx(j) 

     &          + 0.25*deltxu(j+1)*umask(i+1,j+1) /deltyu(j+1)
     &          + 0.25*deltxu(j+1)*umask(i,j+1)   /deltyu(j+1)
     &          + 0.25*deltxu(j)*umask(i,j)       /deltyu(j)
     &          + 0.25*deltxu(j)*umask(i+1,j)     /deltyu(j)
     &        ) / (a*srfarea(j))
  412     continue
  410   continue

      endif

c Read in ice divergence limiter mask
      write (*, *)
      write (*, *) "Reading sea ice divergence limiter mask :"
      write (*, *)
      write (*, *) "File     =  icedivl.nc"
      write (*, *) "Variable =  kma"

      call read_int_2d("icedivl.nc", "kma", plon, plat, kma, title)

      write (*, *) "Title    =  ", trim(title)
      write (*, *)

      ENDIF

      return
      end
