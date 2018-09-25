c Removing a loop which should NOT have been present, which had the effect of
c reducing the vegetation fractions specified via the auxiliary file sibsig.nc.
c In particular, this loop had the strongly undesirable effect of capping the
c vegetation fraction at a maximum value of 0.80.
c SJP 2009/04/27
c
c Add a call to NINT to resolve a warning issued by the g95 Fortran compiler.
c SJP 2009/04/15
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified so that albedos are only read when using the NSiB land
c surface scheme; CABLE calculates albedos itself.
c SJP 2009/04/01
c
c nsib comparison is modified to lsm_type = "nsib "
c AJA 2009/01/22
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
c Modified for coupling of the new, high-resolution version of the ocean model
c to the atmosphere.
c SJP 2007/12/20
c
c Replaced the line "include 'OCEANPARS.f'" with "include 'OPARAMS.f'",
c enabling the header file OCEANPARS.f to be removed from the model source
c code. Renamed one loop variable NT -> NTLOOP to avoid a conflict.
c SJP 2007/05/31
c
c Transferred COMMON blocks to separate header files, as follows:
c /MDAY/  ->  MDAY.f
c SJP 2007/05/29
c
c Modified for fix-digit year numbers.
c SJP 2004/09/29
c
c Modified to read in climatological salinities when SAVEFCOR=T.
c SJP 2004/09/10
c
c $Log: datard.f,v $
c Revision 1.42  2001/10/12 02:28:09  rot032
c Merge HBG and LDR changes.
c
c Revision 1.41  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.40.1.1  2001/10/12 02:13:44  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.40  2001/02/28 04:36:36  rot032
c Further tidy ups from HBG
c
c Revision 1.39  2001/02/22 05:34:37  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.38  2001/02/12 05:39:43  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.37  2000/11/14 03:11:35  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.36  2000/08/16 02:59:24  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.35  2000/06/20 02:08:31  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.34  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.33  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.32  1998/12/10  00:55:30  ldr
c HBG changes to V5-1-21
c
c Revision 1.31  1998/05/27  02:07:35  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.30  1997/12/23  00:35:59  ldr
c Fix the little error recently introduced regarding isoilm.ge.9.
c
c Revision 1.29  1997/12/23  00:23:34  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.28  1997/12/03  05:35:42  ldr
c Put IGW's qflux code with variable MLO depth into main version of model.
c
c Revision 1.27  1997/11/27  05:35:00  ldr
c Flexible treatment of direct and indirect aerosol effects.
c
c Revision 1.26  1997/11/24  23:25:24  ldr
c Indirect sulfates changes to give version used for GRL paper
c
c Revision 1.25.1.1  1997/12/19  02:03:10  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.25  1997/10/06  04:15:03  ldr
c Final corrections to V5-1.
c
c Revision 1.24  1997/10/03  05:51:08  ldr
c Merge of LDR sulfate stuff with HBG/MRD NEC stuff.
c
c Revision 1.23  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.22.1.1  1997/10/03  05:45:45  ldr
c Changes for sulfates from LDR
c
c Revision 1.22  1997/07/24  05:42:49  ldr
c Tuning of ukconv version and T63 fixes from HBG
c
c Revision 1.21  1997/06/11  02:21:28  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.20  1996/12/23  05:42:40  mrd
c Add new gravity wave drag scheme as an option.
c
c Revision 1.19  1996/10/24  01:02:36  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.18  1996/02/19  04:09:46  ldr
c Generalize for 24 levels.
c
c Revision 1.17  1995/12/12  04:47:17  ldr
c Fix up the line that corrects for T63 mask discrepancies in sibsoil.
c
c Revision 1.16  1995/11/16  03:54:16  ldr
c New 9 soil types from EAK.
c
c Revision 1.15  1995/11/10  05:23:03  ldr
c EAK's changes to soil types (increase to 9 soil types).
c
c Revision 1.14  1995/05/02  07:49:15  ldr
c Merge of HBG stuff in V4-5-43h3 into latest version.
c
c Revision 1.13  1994/08/08  17:20:59  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.12  94/08/04  16:54:52  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c 
c Revision 1.12.1.1  1995/05/02  07:31:32  ldr
c Extra changes from HBG to files which already went in as V4-5-43h.
c
c Revision 1.11  94/03/30  12:34:09  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.10  94/01/17  11:29:37  ldr
c Changes from Hal interpolate qfluxes and make them implicit. Also move
c read of qfluxes to datard so it is done every month.
c 
c Revision 1.9  93/12/23  15:31:25  ldr
c Changes to V4-4-54l from HBG for coupled model of V4-5
c 
c Revision 1.8  93/12/17  15:32:01  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.7  93/12/10  16:37:44  ldr
c New albedo and SiB data sets to go with new mask in RS file.
c 
c Revision 1.6  93/12/06  16:55:23  ldr
c Changes to V4-4-45l to get T63 running on SGI (Seca).
c 
c Revision 1.5  93/10/05  13:05:44  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c
c     INPUT/OUTPUT
c     Input:   from common/fewflags in FEWFLAGS.f
c                  chenflag - flag for CHEN version of model
c                  lcouple  - flag to run model in coupled mode, default F,
c                  lsm_type - if "nsib ", use "New SIB" land surface scheme
c                  qflux    - set to T for qflux run, default F
c                  sltrace  - if T compute tracer transport by semi-
c                             Langrangian advection and by mixing and
c                             convection in the vertical
c                  ngwd     - if .ne. 1 read data files for new GWD scheme
c
c              from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c              from common/mday in this subroutine
c                  mdays - days in a month
c
c     Output:  from common/compsst in this subroutine
c                  tsst1 - sea surface temp at beginning of current month
c                  tsst2 - sea surface temp at end of current month
c
c              from common/curent in this subroutine
c                  g50dt - qflux correction array
c
c              from common/msp21 in this subroutine
c                  psrk - chen surface pressure variable
c
c              from common/qflxdat in QFLXDAT.f
c                  gm1cur)qflux interpolate data
c                  gm2cur)
c
c              from common/traceblk in TRACEBLK.f
c                  conmnth - monthly mean tracer concentration, tracer no. nt
c
c              from common/vegdatm in this subroutine
c                  vegm - global data sets for vegetation
c
c              from common/gwddat in GWDDATA.f
c                  sd_gwd        ! Standard deviation
c                  slope_gwd     ! Slope
c                  gamma_gwd     ! Anisotropy
c                  theta_gwd     ! Direction of principal axis
c                  
c
c     In/Out:  from common/masiv4 in MASIV4.f
c                  savegrid - global surface variable
c
c              from common/surf in this subroutine
c                  isoilm - soil types
c
c              from common/timex in TIMEX.f
c                  month - month counter
c
c              from common/traceblk in TRACEBLK.f
c                  source2) tracer sources
c                  source3)
c 

      subroutine datard

c Routine to replace parts of filerd that read data from non-restart files
c i.e. albedo, NSIB data and SSTs
c This allows the albedo and NSIB data to be updated correctly each month in 
c coupled runs.

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list

C Global data blocks
      include 'BCOGCM.f'
      include 'FEWFLAGS.f'
      include 'FILES.f'
      include 'GWDDATA.f'
      include 'LSMI.f'
      include 'MASIV4.f'
      include 'QFLXDAT.f'
      include 'SULFATE1.f'
      include 'TIMEX.f'
      include 'TRACEBLK.f'
      include 'MDAY.f'

      real tsst1,tsst2,tdtm1,tdtm2
      common/compsst/tsst1(ln2,lat),tsst2(ln2,lat)
     &              ,tdtm1(ln2,lat),tdtm2(ln2,lat)

      real ochf,occur,g50dt
      integer khf
      common/curent/ochf(ln2,lat),occur(ln2,lat),khf(ln2,lat)
     &,g50dt(ln2,lat)

      real psrk
      common/msp21/psrk(ln2,lat)

      real aftfh,aftfhg,ocondx,opreci,pmc
      integer isoilm
      common/surf/aftfh(ln2,lat),aftfhg(ln2,lat),ocondx(ln2,lat)
     &          ,isoilm(ln2,lat),opreci(ln2,lat)
     &          ,pmc(ln2,lat)

      real vegm
      common/vegdatm/vegm(ln2,4,lat)

C Local work arrays and variables
      real datacif(lon,lat2) ! for cif layout data input files
      character*80 header
      character*50 header2
      character*14 fname 

      integer ierr
      integer ir
      integer k
      integer kmonth
      integer lg
      integer mg
      integer ns
      integer ntloop
      integer i
      integer j

      real dtmth
      real dz1
      real hcap1
      real phpsr
      real t0ec
      real t00ec
      real t1ec

      character title*60
      integer monthp1
      integer soil_in(lon, lat2)
      real psrk_in(lon, lat2), sig_in(lon, lat2), vegt_in(lon, lat2)

C Local data, functions etc
      real hcap
      data hcap/2.095e8/

      logical got_albedo, got_nsib, got_psrk, got_sss, got_sst
      real albedo_save(lon, lat2, 12), rs_save(lon, lat2, 12),
     &     sss_save(lon, lat2, 12), sst_save(lon, lat2, 12),
     &     z0_save(lon, lat2, 12)

      data got_albedo /.false./
      data got_nsib /.false./
      data got_psrk /.false./
      data got_sss /.false./
      data got_sst /.false./
      save got_albedo, got_nsib, got_psrk, got_sss, got_sst,
     &     albedo_save, rs_save, sss_save, sst_save, z0_save

C Start code : ----------------------------------------------------------

c Read sulfate data
c First, that used for the direct effect

      if(naerosol_d.eq.1)then  !Annual mean sulfate
        open(unit=7,file=so4_direct_file,form='formatted',
     &       status='old', iostat=ierr)
        call errcheck(ierr,'so4_direct_file  ','datard    ')
        read(7,'(a)')header2
        write(6,*)'Direct sulfate file header is: ',header2(1:50)
        read(7,*) datacif
        call trandata(datacif,so4dir)
        close(unit=7)
      elseif(naerosol_d.eq.2)then  !Monthly mean sulfate
        open(unit=7,file=so4_direct_file,form='formatted',
     &       status='old', iostat=ierr)
        call errcheck(ierr,'so4_direct_file  ','datard    ')
        kmonth=0
        do while(kmonth.lt.month)
          read(7,'(a)')header2
          read(7,*) datacif
          call trandata(datacif,so4dir)
          kmonth=kmonth+1
        enddo
        write(6,*)'Direct sulfate file header is: ',header2(1:50)
        close(unit=7)
      endif

c Second, that used for the indirect effect in the radiation scheme

      if(naerosol_i(1).eq.1)then  !Annual mean sulfate
        open(unit=7,file=so4_ind_rad_file,form='formatted',
     &       status='old', iostat=ierr)
        call errcheck(ierr,'so4_ind_rad_file ','datard    ')
        read(7,'(a)')header2
        write(6,*)'Indirect sulfate (radiation) file header is: ',
     &             header2(1:50)
        read(7,*) datacif
        call trandata(datacif,so4rad)
        close(unit=7)
      elseif(naerosol_i(1).eq.2)then  !Monthly mean sulfate
        open(unit=7,file=so4_ind_rad_file,form='formatted',
     &       status='old', iostat=ierr)
        call errcheck(ierr,'so4_ind_rad_file ','datard    ')
        kmonth=0
        do while(kmonth.lt.month)
          read(7,'(a)')header2
          read(7,*) datacif
          call trandata(datacif,so4rad)
          kmonth=kmonth+1
        enddo
        write(6,*)'Indirect sulfate (radiation) file header is: ',
     &             header2(1:50)
        close(unit=7)
      endif

c Third, that used for the indirect effect in the rainfall scheme

      if(naerosol_i(2).eq.1)then  !Annual mean sulfate
        open(unit=7,file=so4_ind_rain_file,form='formatted',
     &       status='old', iostat=ierr)
        call errcheck(ierr,'so4_ind_rain_file','datard    ')
        read(7,'(a)')header2
        write(6,*)'Indirect sulfate (rainfall) file header is: ',
     &             header2(1:50)
        read(7,*) datacif
        call trandata(datacif,so4rain)
        close(unit=7)
      elseif(naerosol_i(2).eq.2)then  !Monthly mean sulfate
        open(unit=7,file=so4_ind_rain_file,form='formatted',
     &       status='old', iostat=ierr)
        call errcheck(ierr,'so4_ind_rain_file','datard    ')
        kmonth=0
        do while(kmonth.lt.month)
          read(7,'(a)')header2
          read(7,*) datacif
          call trandata(datacif,so4rain)
          kmonth=kmonth+1
        enddo
        write(6,*)'Indirect sulfate (rainfall) file header is: ',
     &             header2(1:50)
        close(unit=7)
      endif

c Read in data set needed for "Chen" temperature variable


      If(chenflag)Then
c---- Psr (Spectral form of Psg, where Psg is balanced to
c----  g.z4 where z4 is grid elevations [non-spectral] )
c---- Read in (psr/p00)**cappa

c...  Read the surface pressure data, if this has not already been done

      if (.not. got_psrk) then

        write (*, *)
        write (*, *) "Reading surface pressure data :"
        write (*, *)
        write (*, *) "File     =  psrk.nc"
        write (*, *) "Variable =  psrk"

        call read_real_2d("psrk.nc", "psrk", lon, lat2, psrk_in, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        do j = 1, lat2
          if (j .gt. lat) then
            lg = lat2 + 1 - j
          else
            lg = j
          end if
          do i = 1, lon
            if (j .gt. lat) then
              mg = i
            else
              mg = i + lon
            end if
            psrk(mg, lg) = psrk_in(i, j)
          end do
        end do

        got_psrk = .true.

      end if

c---- create equivalent z* (z4) from Psrk
          t00ec=288.0
          t1ec=0.6652*t00ec
          t0ec=t00ec-t1ec
          do 577 lg=1,lat
          do 577 mg=1,ln2
          phpsr=-cp*(t0ec*log(psrk(mg,lg))+t1ec*(psrk(mg,lg)-1.))
          savegrid(mg,1,lg)=phpsr/grav
 577      continue

      endif  ! chenflag

c Fields needed for new gravity wave drag scheme
      if ( ngwd .ne. 1 ) then
         write(fname,'(''gwd'',a,i2,''.dat'')') trunc, mw-1
         open(unit=7,file=fname,form='formatted',status='old',
     &        iostat=ierr)
         call errcheck(ierr,'gwd data file    ','datard    ')
         read(7,*) datacif
         call trandata(datacif,sd_gwd)
         read(7,*) datacif
         call trandata(datacif,slope_gwd)
         read(7,*) datacif
         call trandata(datacif,gamma_gwd)
         read(7,*) datacif
         call trandata(datacif,theta_gwd)
         close(7)
      end if
         
c Read in monthly albedos. This replaces stuff in timet - less error prone here
c First this month, then next month.
c
c CABLE calculates albedos itself. SJP 2009/04/01

      if (lsm_type .eq. "nsib ") then

c...  Read the albedo data, if this has not already been done
      if (.not. got_albedo) then

        write (*, *)
        write (*, *) "Reading albedo data :"
        write (*, *)
        write (*, *) "File     =  albedo.nc"
        write (*, *) "Variable =  albedo"

        call read_real_3d("albedo.nc", "albedo", lon, lat2, 12,
     &                    albedo_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_albedo = .true.

      end if

c...  Transfer the data for the start and end of the current month to SAVEGRID
      monthp1 = month + 1
      if (monthp1 .eq. 13) monthp1 = 1
      do ns = 1, 2
        do lg = 1, lat
          if (ns .eq. 1) then
            j = lat2 + 1 - lg
          else
            j = lg
          end if
          do i = 1, lon
            if (ns .eq. 1) then
              mg = i
            else
              mg = i + lon
            end if
            savegrid(mg, 7, lg) = albedo_save(i, j, month)
            savegrid(mg, 4, lg) = albedo_save(i, j, monthp1)
          end do
        end do
      end do

      end if     ! lsm_type = "nsib "

c Read in monthly SSTs. First this month, then next month
      if(.not.ncepagcm)then

c...  Read the sea surface temperature data, if this has not already been done
      if (.not. got_sst) then

        write (*, *)
        write (*, *) "Reading sea surface temperature data :"
        write (*, *)
        write (*, *) "File     =  ssta.nc"
        write (*, *) "Variable =  sst"

        call read_real_3d("ssta.nc", "sst", lon, lat2, 12, sst_save,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_sst = .true.

      end if

c...  Transfer the data for the start and end of the current month to SAVEGRID,
c...  if running in stand-alone atmosphere mode, and to the arrays TSST1 and
c...  TSST2, if running in coupled mode. Note that the sign of the values must
c...  be reversed, as the model requires the SSTs to be in minus Kelvin!
      monthp1 = month + 1
      if (monthp1 .eq. 13) monthp1 = 1

      if (lcouple) then

        do ns = 1, 2
          do lg = 1, lat
            if (ns .eq. 1) then
              j = lat2 + 1 - lg
            else
              j = lg
            end if
            do i = 1, lon
              if (ns .eq. 1) then
                mg = i
              else
                mg = i + lon
              end if
              tsst1(mg, lg) = 0.0 - sst_save(i, j, month)
              tsst2(mg, lg) = 0.0 - sst_save(i, j, monthp1)
            end do
          end do
        end do

      else

        do ns = 1, 2
          do lg = 1, lat
            if (ns .eq. 1) then
              j = lat2 + 1 - lg
            else
              j = lg
            end if
            do i = 1, lon
              if (ns .eq. 1) then
                mg = i
              else
                mg = i + lon
              end if
              savegrid(mg, 5, lg) = 0.0 - sst_save(i, j, month)
              savegrid(mg, 6, lg) = 0.0 - sst_save(i, j, monthp1)
            end do
          end do
        end do

      end if   ! lcouple

      endif ! ncepagcm

c Read in climatological SSSs for the start and end of the current month
      if (savefcor) then

c...  Read the sea surface salinity data, if this has not already been done
      if (.not. got_sss) then

        write (*, *)
        write (*, *) "Reading sea surface salinity data :"
        write (*, *)
        write (*, *) "File     =  sssa.nc"
        write (*, *) "Variable =  sss"

        call read_real_3d("sssa.nc", "sss", lon, lat2, 12, sss_save,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_sss = .true.

      end if

c...  Transfer the data for the start and end of the current month to ASAL
      monthp1 = month + 1
      if (monthp1 .eq. 13) monthp1 = 1
      do j = 1, lat2
        do i = 1, lon
          asal(i, j, 1) = sss_save(i, j, month)
          asal(i, j, 2) = sss_save(i, j, monthp1)
        end do
      end do

      endif ! savefcor

c******************************************************************************

c NSIB section

      if (lsm_type .eq. "nsib ") then

c...  Read the data for the New SIB land surface scheme, if this has not
c...  already been done
      if (.not. got_nsib) then

        write (*, *)
        write (*, *) "Reading stomatal resistance data :"
        write (*, *)
        write (*, *) "File     =  sibrs.nc"
        write (*, *) "Variable =  rs"

        call read_real_3d("sibrs.nc", "rs", lon, lat2, 12, rs_save,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        write (*, *)
        write (*, *) "Reading roughness length data :"
        write (*, *)
        write (*, *) "File     =  sibz0.nc"
        write (*, *) "Variable =  z0"

        call read_real_3d("sibz0.nc", "z0", lon, lat2, 12, z0_save,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        write (*, *)
        write (*, *) "Reading vegetation fraction data :"
        write (*, *)
        write (*, *) "File     =  sibsig.nc"
        write (*, *) "Variable =  sig"

        call read_real_2d("sibsig.nc", "sig", lon, lat2, sig_in, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        write (*, *)
        write (*, *) "Reading vegetation type data :"
        write (*, *)
        write (*, *) "File     =  sibvegt.nc"
        write (*, *) "Variable =  vegt"

        call read_real_2d("sibvegt.nc", "vegt", lon, lat2, vegt_in,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        write (*, *)
        write (*, *) "Reading soil type data :"
        write (*, *)
        write (*, *) "File     =  sibsoil.nc"
        write (*, *) "Variable =  soil"

        call read_int_2d("sibsoil.nc", "soil", lon, lat2, soil_in,
     &                   title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

c...  Transfer the vegetation fraction and vegetation type data to VEGM, and
c...  the soil type data to ISOILM. This only needs to be done once, as the
c...  data is time-invariant.
        do j = 1, lat2
          if (j .gt. lat) then
            lg = lat2 + 1 - j
          else
            lg = j
          end if
          do i = 1, lon
            if (j .gt. lat) then
              mg = i
            else
              mg = i + lon
            end if
            vegm(mg, 3, lg) = sig_in(i, j)
            vegm(mg, 4, lg) = vegt_in(i, j)
            isoilm(mg, lg) = soil_in(i, j)
          end do
        end do

        got_nsib = .true.

      end if

c...  Transfer the stomatal resistance and roughness length data for the
c...  current month to VEGM
      do ns = 1, 2
        do lg = 1, lat
          if (ns .eq. 1) then
            j = lat2 + 1 - lg
          else
            j = lg
          end if
          do i = 1, lon
            if (ns .eq. 1) then
              mg = i
            else
              mg = i + lon
            end if
            vegm(mg, 1, lg) = rs_save(i, j, month)
            vegm(mg, 2, lg) = z0_save(i, j, month)
          end do
        end do
      end do

c
c  put land sea mask over data and set z0 for sea to 0.001m
c
      do 240 lg=1,lat
      do 240 mg=1,ln2
      if(imsl(mg,lg).lt.4) then
         vegm(mg,1,lg)=10000.
         vegm(mg,2,lg)=0.001
      endif 
 240  continue

      call insoilveg

      end if   ! lsm_type = "nsib "

c End of NSIB section

c******************************************************************************

c QFLUX section

c.... Set up (GAMMA50 dT/dt) for the month for calculation
c.... of Q-fluxes at each step for SEA points
          dz1=100.0 ! mixed layer depth now 100m (was 50m)
          hcap1=hcap/50.*dz1
          dtmth=86400.0*mdays(month)
       do 220 lg=1,lat
       do 220 mg=1,ln2
        g50dt(mg,lg)=(hcap1/dtmth)*
     & (savegrid(mg,5,lg) -savegrid(mg,6,lg))
  220  continue

c IGW's modified qflux treatment with variable MLO depth...

      if (qflux) then
         open (unit=7,file='qflux.dat',form='unformatted',status='old',
     &    iostat=ierr)
         call errcheck(ierr,'qflux.dat        ','datard    ')
         rewind 7
           kmonth=0
           do while(kmonth.lt.month)
             read(7)header
             read(7) kmonth,ir
             read(7)((gm1cur(mg,lg),mg=1,lon),lg=1,lat)
     &              ,((gm1cur(mg+lon,lg),mg=1,lon),lg=1,lat)
           write(6,*)' on the qflux file the header reads:'
           write(6,'(a)')header
             read(7)header
             read(7) kmonth,ir
             read(7)((gm2cur(mg,lg),mg=1,lon),lg=1,lat)
     &              ,((gm2cur(mg+lon,lg),mg=1,lon),lg=1,lat)
           write(6,'(a)')header
             read(7)header
             read(7) kmonth,ir
             read(7)((gm3cur(mg,lg),mg=1,lon),lg=1,lat)
     &              ,((gm3cur(mg+lon,lg),mg=1,lon),lg=1,lat)
           write(6,'(a)')header
           enddo
         close(unit=7)
      end if

c sltrace section

      if (sltrace) then
         open (unit=7,file='source2.dat',form='formatted',status='old',
     &    iostat=ierr)
         call errcheck(ierr,'source2.dat      ','datard    ')
             read(7,'(a80)') header
           write(6,*)' on the source2 file the header reads:'
           write(6,'(a)')header
             read(7,'(5e12.5)')source2
         close(unit=7)

         open (unit=7,file='source3.dat',form='formatted',status='old',
     &    iostat=ierr)
         call errcheck(ierr,'source3.dat      ','datard    ')
           kmonth=0
           do while(kmonth.lt.month)
             read(7,'(a80)')header
             read(7,*)kmonth
             read(7,'(5e12.5)')source3
           enddo
         close(unit=7)
         write(6,*)' on the source3 file the header reads:',kmonth
         write(6,'(a)')header
c   rescale to ppmv
        do lg=1,lat2
        do mg=1,lon 
        source2(mg,lg)=source2(mg,lg)*29./12.*1.E6/3.1536E7
        source3(mg,lg)=source3(mg,lg)*29./12.*1.E-4
        enddo
        enddo
c      initialize month accumulation array
        write(6,*) 'initialize conmnth '
        do ntloop=1,ntrace
        do k=1,nl
        do lg=1,lat2
        do mg=1,lon 
          conmnth(mg,lg,k,ntloop)=0.
        enddo
        enddo
        enddo
        enddo
c      history files for sltarce opend in openfl.f

      end if

      return
      end

      subroutine insoilveg

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list

C Global data blocks
      include 'LSMI.f'
      include 'FEWFLAGS.f'
      include 'SURF1.f'
      include 'TIMEX.f'

      real prhos,pcss,pcnsd,pbch,pssat,pswilt,phsbh,psfc,psucs
     & ,phyds,pcdr3,prlaim,prlais,pslveg,pscveg
      integer ip2bp3,ipbp2
      common/psoilvegpr/prhos(ln2,lat),pcss(ln2,lat),
     &pcnsd(ln2,lat),pbch(ln2,lat),pssat(ln2,lat),
     &pswilt(ln2,lat),phsbh(ln2,lat),ip2bp3(ln2,lat),
     &ipbp2(ln2,lat),psfc(ln2,lat),psucs(ln2,lat),
     &phyds(ln2,lat),pcdr3(ln2,lat),prlaim(ln2,lat),
     &prlais(ln2,lat),pslveg(ln2,lat),pscveg(ln2,lat)

      real zs,zsh,zsdzs,zsfrac
      common/soilzs/zs(ms),zsh(ms+1),zsdzs(ms),zsfrac(ms)

      real aftfh,aftfhg,ocondx,opreci,pmc
      integer isoilm
      common/surf/aftfh(ln2,lat),aftfhg(ln2,lat),ocondx(ln2,lat)
     &          ,isoilm(ln2,lat),opreci(ln2,lat)
     &          ,pmc(ln2,lat)

      real vegm
      common/vegdatm/vegm(ln2,4,lat)

      logical lakesg,lakeind
      common/glakes/lakesg(ln2,lat),lakeind(lat)

C Local work arrays and variables
      real sand(9),clay(9),silt(9),cnsd(9)
      real bch(9),hyds(9),sucs(9),hsbh(9)
      integer ibp2(9),i2bp3(9)
      real rhos(9),css(9),ssat(9),swilt(9)
      real sfc(9),cdr3(9)
      real rlaim(13),rlais(13),slveg(13),scveg(13),rsmin(13)

      integer isoil
      integer iveg
      integer k
      integer lg
      integer mg

      real ztot

C Local data, functions etc
c     data zs/0.050,0.125,0.315,0.760,1.400,2.000/
      data zs/0.022,0.058,0.154,0.409,1.085,2.872/

      data sand/0.83,0.37,0.17,0.6,0.52,0.27,0.58,0.13,0.37/
      data clay/0.09,0.30,0.66,0.2,0.42,0.48,0.27,0.17,0.30/
      data silt/0.08,0.33,0.17,0.2,0.06,0.25,0.15,0.70,0.33/

      data bch/4.2,7.1,11.4,5.15,10.4,10.4,7.12,5.83,7.1/
      data hyds/166.e-6,4.e-6,1.e-6,21.e-6,2.e-6,1.e-6,6.e-6,
     *             800.e-6,1.e-6/
      data sucs/-0.106,-0.591,-0.405,-0.348,-0.153,-0.49,-0.299,
     &          -0.356,-0.153/

      data rhos/7*2600.,1300.,910./
      data  css/7*850.,1920.,2100./
      data  ssat/.398,.479,.482,.443,.426,.482,.420,.450,.479/
      data swilt/.072,.216,.286,.135,.219,.283,.175,.395,.216/

      data  sfc/0.143,0.301,0.367,0.218,0.310,0.370,0.255,0.450,0.301/
      data cdr3/1.255,0.334,0.138,0.521,0.231,0.199,0.375,0.623,0.334/

c     data rlaim/6.0,6.0,6.0,6.0,6.0,4.0,4.0,4.0,1.0,4.0,0.5,4.0,0.0/
      data rlaim/6.0,5.5,5.0,4.5,5.0,4.0,3.0,3.5,1.0,4.0,0.5,4.0,0.0/
      data rlais/2.0,2.0,2.0,2.0,2.0,1.5,1.5,1.5,1.0,0.5,0.5,0.5,0.0/
      data slveg/1.0,5.5,3.0,1.0,3.0,3.0,3.5,3.0,0.5,3.5,0.1,3.5,0.0/
      data scveg/.05,0.0,0.0,0.0,0.0,.05,.05,.05,0.1,0.0,0.0,0.4,0.0/
      data rsmin/300.,350.,300.,350.,400.,350.,240.,350.,150.,400.,
     &           995.,150.,995./


C Start code : ----------------------------------------------------------

      do 100 isoil=1,9
       cnsd(isoil)= sand(isoil)*0.3+clay(isoil)*0.25+silt(isoil)*0.265
       hsbh(isoil) = hyds(isoil)*abs(sucs(isoil))*bch(isoil)
       ibp2(isoil) = nint(bch(isoil))+2
       i2bp3(isoil) = 2*nint(bch(isoil))+3
100   enddo
      cnsd(9)=2.51

c.... zs  = Thickness of soil layers
c.... zsh = Distance between mid points of soil layers
      zsh(1) = .5*zs(1)
      zsh(ms+1) = 0.5*zs(ms)
      do k=2,ms
       zsh(k)= (zs(k-1)+zs(k))/2.0
      enddo
c.... zsdzs defined next
      zsdzs(1) = 0.
      do k=2,ms
        zsdzs(k) = zs(k)/(zs(k)+zs(k-1))
      enddo
c.... zsfrac = Fraction of total depth for any soil layer
      ztot=0.0
      do k=1,ms
       ztot=ztot+zs(k)
      enddo
      do k=1,ms
       zsfrac(k)=zs(k)/ztot
      enddo

c No lakes unless newriver & T63
      do lg=1,lat
        lakeind(lg)=.false.
      enddo
      do lg=1,lat
      do mg=1,ln2
        lakesg(mg,lg)=.false.
      enddo
      enddo

      if(newriver)then

c For T63, put in the Great Lakes
c Make these points veg type 11
c Treat them as swamp with a permanent very deep "puddle" on top. The 
c  puddle depth will be maintained through the new river routing scheme.
c  
c  Specify Great Lakes (L)
c      144 145 146 147 148 149 150 151 152 153
c     .........................................
c (lg).   .  Superior .   .   .   .   .   .   .(lgns)
c  23 .   . L . L . L .   .   .   .   .   .   .  74
c     .   .   .   .   .   .   .   .   .   .   .
c     .........................................
c     .   .   .   M   .   Huron   .   .   .   .
c  24 .   .   .   i L . L . L . L .   .   .   .  73
c     .   .   .   c   .   .   .   .   .   .   .
c     ............h............................
c     .   .   .   i   .   .   .   .   Ontario .
c  25 .   .   .   g L .   . L .   .   . L .   .  72
c     .   .   .   a   .   .   .   .   .   .   .
c     ............n............................
c     .   .   .   .   .   .   .   .   .   .   .
c  26 .   .   .   .   .   . L . L .   .   .   .  71
c     .   .   .   .   .   .  Erie .   .   .   .
c     .........................................
      lakesg(145,23)=.true.
      lakesg(146,23)=.true.
      lakesg(147,23)=.true.
      lakesg(147,24)=.true.
      lakesg(148,24)=.true.
      lakesg(149,24)=.true.
      lakesg(150,24)=.true.
      lakesg(147,25)=.true.
      lakesg(149,25)=.true.
      lakesg(149,26)=.true.
      lakesg(150,26)=.true.
c Niagara Falls between
      lakesg(152,25)=.true.

c  Specify Lake Winnipeg (approx by 2 grid boxes)
c      139 140 141 142 
c     .................
c (lg).   .   .   .   .(lgns)
c  19 .   .   .tr .   .  78
c     .   .   .221.   .
c     .................
c     .   Winnipeg.   .
c  20 .   . L . L .   .  77
c     .   .   .   .   .
c     .................
c       (puddle water out of trough 221)
      lakesg(140,20)=.true.
      lakesg(141,20)=.true.

c  Specify "Baker Lakes" (approx by 2 grid boxes)
c      139 140 141 142 143 
c     .....................
c (lg).   .Baker Lakes.   .(lgns)
c  14 .   . L . L .tr .   .  83
c     .   .   .   .273.   .
c     .....................
c       (puddle water out of trough 273)
      lakesg(140,14)=.true.
      lakesg(141,14)=.true.

c  Specify lakes to help regulate southern Hudson Bay rivers (3 boxes)
c      147 148 149 150 151 152 153
c     ........|___.....___|........
c (lg).   .   .   |   |   .   .   .(lgns)
c  20 .   .   . L |   | L .   .   .  77
c     .   .   .   |   |   .   .   .
c     ............|___|............
c     .   .   .   .   .   .   .   .
c  21 .   .tr .   . L .   .   .   .  76
c     .   .209.   .   .   .   .   .
c     .............................
c       (puddle water out of lake in middle)
      lakesg(149,20)=.true.
      lakesg(151,20)=.true.
      lakesg(150,21)=.true.

c
c Other points may be made into lakes later 
c  eg Lake Chad, Lake Victoria etc.
c
      do lg=1,lat
      do mg=1,ln2
        if(lakesg(mg,lg))vegm(mg,4,lg)=11
      enddo
      enddo
      do lg=1,lat
c Set indicator per double-latitude row if any lakes present
        call checkl(lakesg(1,lg),ln2,lakeind(lg))
      enddo

      endif ! newriver

c Set up soil and vegetation types
      do 1000 lg=1,lat
      do 1000 mg=1,ln2
      if( imsl(mg,lg).eq.4) then
       isoil = isoilm(mg,lg)
       prhos(mg,lg) = rhos(isoil)
       pcss(mg,lg)  = css(isoil)
       pcnsd(mg,lg) = cnsd(isoil)
       pbch(mg,lg)  = bch(isoil)
       pssat(mg,lg) = ssat(isoil)
       pswilt(mg,lg)= swilt(isoil)
       phsbh(mg,lg) = hsbh(isoil)
       ip2bp3(mg,lg)= i2bp3(isoil)
       ipbp2(mg,lg) = ibp2(isoil)
       psfc(mg,lg)  = sfc(isoil)
       psucs(mg,lg) = sucs(isoil)
       phyds(mg,lg) = hyds(isoil)
       pcdr3(mg,lg) = cdr3(isoil)

       iveg = nint(vegm(mg,4,lg))
       prlaim(mg,lg) = rlaim(iveg)
       prlais(mg,lg) = rlais(iveg)
       pslveg(mg,lg) = slveg(iveg)
       pscveg(mg,lg) = scveg(iveg)
       vegm(mg,1,lg) = rsmin(iveg) ! Non seasonally varying RS now

      else

       prhos(mg,lg) = -7777777
       pcss(mg,lg)  = -7777777
       pcnsd(mg,lg) = -7777777
       pbch(mg,lg)  = -7777777
       pssat(mg,lg) = -7777777
       pswilt(mg,lg)= -7777777
       phsbh(mg,lg) = -7777777
       ip2bp3(mg,lg)= -7777777
       ipbp2(mg,lg) = -7777777
       psfc(mg,lg)  = -7777777
       psucs(mg,lg) = -7777777
       phyds(mg,lg) = -7777777
       prlaim(mg,lg) = -7777777
       prlais(mg,lg) = -7777777
       pslveg(mg,lg) = -7777777
       pscveg(mg,lg) = -7777777
      endif
1000  continue

      if(newriver)then

c Always make sure lakes are wet, with no soil
c (wbice may have values > 0 : this implies frozen lake water)
        do lg=1,lat
        do mg=1,ln2
          if(lakesg(mg,lg))then
            pssat(mg,lg)=1.0
            pswilt(mg,lg)=0.0
            do k=1,ms
              wb(mg,k,lg)=1.0
            enddo
          endif
        enddo
        enddo

      endif ! newriver

      return
      end

      subroutine trandata(datacif,datatran)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      real datacif(lon,lat2)
      real datatran(ln2,lat)

C Local work arrays and variables
      integer lg
      integer lgns
      integer mg
      integer ma
      integer ns

C Start code : ----------------------------------------------------------

c To transform from cif type data with lgns=1 at SP
c to array with corresponding NH and SH latitude rows of
c data adjacent in the common block (for 2xlon vector length)

      do lg=1,lat
      do ns=1,2
        lgns=lg*(ns-1)+(lat2p-lg)*(2-ns)
        ma=(ns-1)*lon
        do mg=1,lon
          datatran(mg+ma,lg)=datacif(mg,lgns)
        enddo
      enddo
      enddo

      return
      end
