c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c Modified for five-character experiment names.
c SJP 2009/04/22
c
c Updated for the new version of the subroutine hstring.
c SJP 2008/11/20
c
c Modified to enable five-digit year numbers.
c SJP 2004/09/29
c
c Specify the precision of all integer arrays and real variables that are
c used as arguments to netCDF routines. This ensures portability of the code
c across 32- and 64-bit platforms.
c SJP 2001/12/13
c
c Updated from v2.4 to v3 of netCDF.
c SJP 2001/11/21
c
c $Log: openhist.f,v $
c Revision 1.27  2000/11/14 03:11:38  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.26  1999/05/20 06:23:54  rot032
c HBG changes to V5-2
c
c Revision 1.25  1997/08/13  00:03:44  mrd
c Minor changes to attributes so file can be read by GrADS.
c
c Revision 1.24  1997/07/24  22:58:55  mrd
c Fix zero length string problem.
c
c Revision 1.23  1997/07/22  05:04:40  mrd
c Get rid of rubbish characters at end of string from hstring.
c
c Revision 1.22  1996/01/28  22:40:41  mrd
c Change the valid_min and valid_max attributes from short to long to
c agree with format of monthly mean history.
c
c Revision 1.21  1995/05/02  07:25:58  ldr
c Merge HBG V4-5-43h2 stuff into latest version.
c
c Revision 1.20  94/09/13  09:51:31  mrd
c Add daily mean 10m windspeed and 2m relative humidity as history variables.
c 
c Revision 1.19  94/09/09  15:02:50  mrd
c  Increase scale factor by small fraction to ensure min and max values
c always fit within range despite rounding.
c 
c Revision 1.18  94/09/09  14:54:00  mrd
c Added daily average cloud, surface temperature and screen temperature
c as history variables.
c 
c Revision 1.17  94/09/09  14:14:48  mrd
c Added possbility of saving two history files at different frequencies.
c 
c Revision 1.16  94/08/08  17:20:47  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.15  94/06/28  14:35:19  mrd
c Add sea-level pressure as a history variable.
c 
c Revision 1.14.1.1  95/05/02  07:14:37  ldr
c Other changes from HBG that went into coupled run but not into V4-5-43h.
c 
c Revision 1.14  94/06/03  14:08:24  mrd
c Set a lower limit on qmax so stratospheric moisture is handled properly.
c 
c Revision 1.13  94/05/26  14:16:06  mrd
c Fixed error that mean scale factor and offset were not set properly
c for layer cloud amounts if total cloud was not being saved.
c 
c Revision 1.12  93/12/06  14:31:41  mrd
c Increased ranges for hfl, pev, rgd and run.
c 
c Revision 1.11  93/11/30  11:43:06  mrd
c Increase upper bound for runoff to 50 mm/day.
c 
c Revision 1.10  93/11/29  15:24:44  mrd
c Changed names of multi-level variables to lower case to be compatible
c with monthly mean files.
c 
c Revision 1.9  93/10/11  14:20:48  ldr
c Changes to get V4-4-24l running on the VP.
c 
c Revision 1.8  93/10/08  09:39:04  mrd
c Moved ncsfil call to work around netcdf bug on SGI.
c 
      subroutine openhist(ihist, iyear, month)

c     Open a netCDF file for the model history and set attributes of the
c     variables.

      implicit none

C Global parameters
      include 'PARAMS.f'

C Argument list
      integer ihist
      integer iyear
      integer month

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'
      include 'HIST.f'
      include 'netcdf.inc'
      include 'VERSION.f'
      include 'FILES.f'

C Local work arrays and variables
      integer   ierr
      integer   dims(3)
      integer   londim, latdim, recdim
      integer   latid, lonid, timeid
      integer   lngstr
      integer   id
      integer   imode
      integer   i,k
      real*4 xlon(lon/nhcomp)
      character hname*21, name*50, lname*50, hist*82
      character startstr*80

C Local data, functions etc

C Start code : ----------------------------------------------------------

      if((savehist(1).or.savehist(2)).and.(nhist.eq.0))then
         print *,'Model stopped in openhist,: savehist=T, but nhist=0'
         print *,'If savehist required, change to nhist=1 in HIST.f'
         stop
      endif

c  If there is only one file, use the name hist..., otherwise use hista...
c  and histb...
      if ( .not. savehist(2) ) then
         hname= 'hist'//str//'.nc'
      else  
         if ( ihist.eq.1 ) then
            hname= 'hista'//str//'.nc'
         else
            hname= 'histb'//str//'.nc'
         end if
      end if

c     Set up the netcdf file
      ierr = nf_create(hname, nf_clobber, histid(ihist))
      if (ierr .ne. nf_noerr) then
	print*, ' Error in creating file ', ierr
	stop
      end if

c     Index for the first set of history data.
      histset(ihist) = 1

c     Create dimensions, lon, lat and rec
c.... if nhcomp>1 save on reduced grid
c.... BUT must be at suitable resolution
c....  (T63 => mw=64, R42 => mw=43, R21 => mw=22)
      if((nhcomp.le.0).or.(nhcomp.ge.4))then
        print *,'openhist.f : nhcomp=',nhcomp
        print *,'May only be 1,2,or 3 depending on model resolution'
        stop
      endif
      if((nhcomp.eq.3).and.(mw.lt.64))then
        print *,'openhist.f : nhcomp=3'
        print *,'Can only use this at T63'
        stop
      endif
      if((nhcomp.eq.2).and.(mw.lt.43))then
        print *,'openhist.f : nhcomp=2'
        print *,'Cannot use this at R21'
        stop
      endif

      ierr = nf_def_dim(histid(ihist), 'lon', lon/nhcomp, londim)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      ierr = nf_def_dim(histid(ihist), 'lat', lat2/nhcomp, latdim)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      ierr = nf_def_dim(histid(ihist), 'time', nf_unlimited, recdim)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

c     Create the variables
      dims(1) = londim
      dims(2) = latdim
      dims(3) = recdim

c     Create global attributes

c     Model version
      ierr = nf_put_att_text(histid(ihist), nf_global, 'version',
     &                       len(trim(version)), trim(version))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

c     Experiment description
      ierr = nf_put_att_text(histid(ihist), nf_global, 'runtype', 5,
     &                       runtype)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

c     History
      call hstring(hist)
      ierr = nf_put_att_text(histid(ihist), nf_global, 'history',
     &                       len(trim(hist)), trim(hist))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

c     Year and month of history file
      ierr = nf_put_att_int(histid(ihist), nf_global, 'year', nf_int,
     &                      1, iyear)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      ierr = nf_put_att_int(histid(ihist), nf_global, 'month', nf_int,
     &                      1, month)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

c     Sigma levels
      ierr = nf_put_att_real(histid(ihist), nf_global, 'sigma_lev',
     &                        nf_float, nl, real(sig, 4))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      ierr = nf_put_att_real(histid(ihist), nf_global, 'sigma_halflev',
     &                        nf_float, nlp, real(sigh, 4))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

c     Define attributes for the dimensions
      ierr = nf_def_var(histid(ihist), 'lon', nf_float, 1, dims(1),
     &                  lonid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      ierr = nf_put_att_text(histid(ihist), lonid, 'long_name', 9,
     &                  'longitude')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      ierr = nf_put_att_text(histid(ihist), lonid, 'units', 12,
     &                       'degrees_east')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

      ierr = nf_def_var(histid(ihist), 'lat', nf_float, 1, dims(2),
     &                  latid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      ierr = nf_put_att_text(histid(ihist), latid, 'long_name', 8,
     &                  'latitude')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      ierr = nf_put_att_text(histid(ihist), latid, 'units', 13,
     &                       'degrees_north')
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

      ierr = nf_def_var(histid(ihist), 'time', nf_float, 1, dims(3),
     &                  timeid)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      write (startstr, "(a,i5.5,a,i2.2,a)")
     &         'days since ', iyear, '-', month, '-1 00:00'
      ierr = nf_put_att_text(histid(ihist), timeid, 'units',
     &                       lngstr(startstr), startstr)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

c     Check on whether hflgs are valid
      if ( ihist.eq.2 .and. all_hflg(ihist)) then
         print*,
     &     'ERROR - can''t save all fields in second history file'
         stop
      end if

c     Surface height is constant through a run and so need not be a 
c     record variable. Call attrib with ndim=2
      if ( all_hflg(ihist) .or. zht_hflg(ihist) ) then
         lname = 'Surface height'
         call attrib(histid(ihist), id, dims, 2, 'zht', lname, 'm')
      end if

c     Create the temperature variables
      if ( all_hflg(ihist) .or. t_hflg(ihist) ) then
         do k=1,nl
            write(name,1000) 't', k
            write(lname,1100) 'Temperature at sigma level',k
            call attrib(histid(ihist), tid(k, ihist), dims, 3, name,
     &                  lname, 'K')
         end do
      end if
 1000 format(a,i2.2)
 1100 format(a,1x,i2)

c     Create the wind variables
      if ( all_hflg(ihist) .or. u_hflg(ihist) ) then
         do k=1,nl
            write(name,1000) 'u', k
            write(lname,1100) 'Zonal wind at sigma level',k
            call attrib(histid(ihist), uid(k, ihist), dims, 3, name,
     &                  lname, 'm/s')
         end do
      end if
      if ( all_hflg(ihist) .or. v_hflg(ihist) ) then
         do k=1,nl
            write(name,1000) 'v', k
            write(lname,1100) 'Meridional wind at sigma level',k
            call attrib(histid(ihist), vid(k, ihist), dims, 3, name,
     &                  lname, 'm/s')
         end do
      end if

c     Specific humidity
      if ( all_hflg(ihist) .or. q_hflg(ihist) ) then
         do k=1,nl
            write(name,1000) 'q', k
            write(lname,1100) 'Specific humidity at sigma level',k
            call attrib(histid(ihist), qid(k, ihist), dims, 3, name,
     &                  lname, 'kg/kg')
         end do
      end if

      if ( all_hflg(ihist) .or. psf_hflg(ihist) ) then
         lname ='Surface pressure'
         call attrib(histid(ihist), pid(ihist), dims, 3, 'psf', lname,
     &               'mb')
      end if

      if ( all_hflg(ihist) .or. psl_hflg(ihist) ) then
         lname ='Sea-level pressure'
         call attrib(histid(ihist), plid(ihist), dims, 3, 'psl', lname,
     &               'mb')
      end if

      if ( all_hflg(ihist) .or. imsl_hflg(ihist) ) then
         lname = 'Surface type'
         ierr = nf_def_var(histid(ihist), 'imsl', nf_byte, 3, dims, id)
         if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
         ierr = nf_put_att_text(histid(ihist), id, 'long_name',
     &                          lngstr(lname), lname)
         if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      end if

      if ( all_hflg(ihist) .or. tsu_hflg(ihist) ) then
         lname = 'Surface temperature'
         call attrib(histid(ihist), id, dims, 3, 'tsu', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. tb2_hflg(ihist) ) then
         lname = 'Soil temperature at level 2'
         call attrib(histid(ihist), id, dims, 3, 'tb2', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. tb3_hflg(ihist) ) then
         lname = 'Soil temperature at level 3'
         call attrib(histid(ihist), id, dims, 3, 'tb3', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. als_hflg(ihist) ) then
         lname = 'Albedo'
         call attrib(histid(ihist), alsid(ihist), dims, 3, 'als',
     &               lname, ' ')
      end if

      if ( all_hflg(ihist) .or. wfg_hflg(ihist) ) then
         lname = 'Soil moisture upper'
         call attrib(histid(ihist), id, dims, 3, 'wfg', lname, ' ')
      end if

      if ( all_hflg(ihist) .or. wfb_hflg(ihist) ) then
         lname = 'Soil moisture lower'
         call attrib(histid(ihist), id, dims, 3, 'wfb', lname, ' ')
      end if

      if ( all_hflg(ihist) .or. snd_hflg(ihist) ) then
         lname = 'Snow depth'
         call attrib(histid(ihist), id, dims, 3, 'snd', lname, 'cm')
      end if

      if ( all_hflg(ihist) .or. sid_hflg(ihist) ) then
         lname = 'Sea-ice depth'
         call attrib(histid(ihist), id, dims, 3, 'sid', lname, 'm')
      end if

      if ( all_hflg(ihist) .or. rnd_hflg(ihist) ) then
         lname = 'Precipitation rate'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'rnd', lname, 'mm/day')
      end if

      if ( all_hflg(ihist) .or. rnc_hflg(ihist) ) then
         lname = 'Convective precipitation rate'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'rnc', lname, 'mm/day')
      end if

      if ( all_hflg(ihist) .or. evp_hflg(ihist) ) then
         lname = 'Evaporation'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'evp', lname, 'mm/day')
      end if

      if ( all_hflg(ihist) .or. hfl_hflg(ihist) ) then
         lname = 'Sensible heat flux'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'hfl', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. sgn_hflg(ihist) ) then
         lname = 'Net SW at the surface'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'sgn', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. sgc_hflg(ihist) ) then
         lname = 'Clear sky net SW at the surface'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'sgc', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. sgd_hflg(ihist) ) then
         lname = 'Downward SW at the surface'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'sgd', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. sit_hflg(ihist) ) then
         lname = 'SW insolation at TOA'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'sit', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. sot_hflg(ihist) ) then
         lname = 'Reflected SW at TOA'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'sot', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. soc_hflg(ihist) ) then
         lname = 'Clear sky reflected SW at TOA'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'soc', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. rgn_hflg(ihist) ) then
         lname = 'Net LW at the surface'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'rgn', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. rgc_hflg(ihist) ) then
         lname = 'Clear sky net LW at the surface'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'rgc', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. rgd_hflg(ihist) ) then
         lname = 'Downward LW at the surface'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'rgd', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. rtu_hflg(ihist) ) then
         lname = 'Net LW at TOA'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'rtu', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. rtc_hflg(ihist) ) then
         lname = 'Clear sky net LW at TOA'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'rtc', lname, 'W/m2')
      end if

      if ( all_hflg(ihist) .or. tax_hflg(ihist) ) then
         lname = 'Zonal surface stress'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'tax', lname, 'N/m2')
      end if

      if ( all_hflg(ihist) .or. tay_hflg(ihist) ) then
         lname = 'Meridional surface stress'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'tay', lname, 'N/m2')
      end if

      if ( all_hflg(ihist) .or. cld_hflg(ihist) ) then
         lname = 'Total cloud'
         call attrib(histid(ihist), cldid(ihist), dims, 3, 'cld',
     &               lname, ' ')
      end if

      if ( all_hflg(ihist) .or. clda_hflg(ihist) ) then
         lname = 'Average total cloud'
         call attrib(histid(ihist), id, dims, 3, 'clda', lname, ' ')
      end if

      if ( all_hflg(ihist) .or. cll_hflg(ihist) ) then
         lname = 'Low level cloud'
         call attrib(histid(ihist), cllid(ihist), dims, 3, 'cll',
     &               lname, ' ')
      end if

      if ( all_hflg(ihist) .or. clm_hflg(ihist) ) then
         lname = 'Middle level cloud'
         call attrib(histid(ihist), clmid(ihist), dims, 3, 'clm',
     &               lname, ' ')
      end if

      if ( all_hflg(ihist) .or. clh_hflg(ihist) ) then
         lname = 'High level cloud'
         call attrib(histid(ihist), clhid(ihist), dims, 3, 'clh',
     &               lname, ' ')
      end if

      if ( all_hflg(ihist) .or. ich_hflg(ihist) ) then
         lname = 'Level of high cloud'
         ierr = nf_def_var(histid(ihist), 'ich', nf_byte, 3, dims,
     &                     ichid(ihist))
         if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
         ierr = nf_put_att_text(histid(ihist), ichid(ihist),
     &                          'long_name', lngstr(lname), lname)
         if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      end if

      if ( all_hflg(ihist) .or. icm_hflg(ihist) ) then
         lname = 'Level of medium cloud'
         ierr = nf_def_var(histid(ihist), 'icm', nf_byte, 3, dims,
     &                     icmid(ihist))
         if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
         ierr = nf_put_att_text(histid(ihist), icmid(ihist),
     &                          'long_name', lngstr(lname), lname)
         if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      end if

      if ( all_hflg(ihist) .or. ict_hflg(ihist) ) then
         lname = 'Top of low cloud'
         ierr = nf_def_var(histid(ihist), 'ict', nf_byte, 3, dims,
     &                     ictid(ihist))
         if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
         ierr = nf_put_att_text(histid(ihist), ictid(ihist),
     &                          'long_name', lngstr(lname), lname)
         if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      end if

      if ( all_hflg(ihist) .or. icb_hflg(ihist) ) then
         lname = 'Bottom of low cloud'
         ierr = nf_def_var(histid(ihist), 'icb', nf_byte, 3, dims,
     &                     icbid(ihist))
         if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
         ierr = nf_put_att_text(histid(ihist), icbid(ihist),
     &                          'long_name', lngstr(lname), lname)
         if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"
      end if

      if ( all_hflg(ihist) .or. tsc_hflg(ihist) ) then
         lname = 'Screen temperature'
         call attrib(histid(ihist), id, dims, 3, 'tsc', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. tsh_hflg(ihist) ) then
         lname = 'Maximum screen temperature'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'tsh', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. tsl_hflg(ihist) ) then
         lname = 'Minimum screen temperature'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'tsl', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. tsca_hflg(ihist) ) then
         lname = 'Average screen temperature'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'tsca', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. tgh_hflg(ihist) ) then
         lname = 'Maximum surface temperature'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'tgh', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. tgl_hflg(ihist) ) then
         lname = 'Minimum surface temperature'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'tgl', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. tsua_hflg(ihist) ) then
         lname = 'Average surface temperature'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'tsua', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. rhsa_hflg(ihist) ) then
         lname = 'Screen level relative humidity'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'rhsa', lname,
     &               'percent')
      end if

      if ( all_hflg(ihist) .or. v10ma_hflg(ihist) ) then
         lname = 'Average 10m wind speed'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'v10ma', lname, 'm/s')
      end if

      if ( all_hflg(ihist) .or. run_hflg(ihist) ) then
         lname = 'Runoff'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'run', lname, 'mm/day')
      end if

      if ( all_hflg(ihist) .or. int_hflg(ihist) ) then
         lname = 'Canopy rainfall interception'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'int', lname, 'mm/day')
      end if

      if ( all_hflg(ihist) .or. tgg_hflg(ihist) ) then
         lname = 'Bare ground temperature'
         call attrib(histid(ihist), id, dims, 3, 'tgg', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. tgf_hflg(ihist) ) then
         lname = 'Canopy temperature'
         call attrib(histid(ihist), id, dims, 3, 'tgf', lname, 'K')
      end if

      if ( all_hflg(ihist) .or. pev_hflg(ihist) ) then
         lname = 'Potential evaporation'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'pev', lname, 'mm/day')
      end if

      if ( all_hflg(ihist) .or. sev_hflg(ihist) ) then
         lname = 'Scaling evaporation'
         call accum_chk(lname,ihist)
         call attrib(histid(ihist), id, dims, 3, 'sev', lname, 'mm/day')
      end if

      if ( all_hflg(ihist) .or. ico_hflg(ihist)  ) then
         lname = 'Ice concentration'
         call attrib(histid(ihist), id, dims, 3, 'ico', lname, ' ')
      end if

c     Leave define mode
      ierr = nf_enddef(histid(ihist))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

c     Turn off the data filling to save time.
      ierr = nf_set_fill(histid(ihist), nf_nofill, imode)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

c     Write the longitude data. The latitudes are not available yet.
c.... if nhcomp>1 save on reduced grid
      do i=1,lon/nhcomp
	xlon(i) = real((float(i-1) * 360./float(lon/nhcomp)), 4)
      end do
      ierr = nf_put_var_real(histid(ihist), lonid, xlon)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in openhist"

      return
      end

c-------------------------------------------------------------------------
      subroutine attrib(cdfid, vid, dims, ndims, name, lname, units)

      implicit none

C Global parameters
      include 'netcdf.inc'

C Argument list
      integer   cdfid
      integer   vid
      integer   ndims
      integer   dims(ndims)
      character name*(*)
      character lname*(*)
      character units*(*)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      integer lngstr
      integer ierr

C Local data, functions etc

C Start code : ----------------------------------------------------------

      ierr = nf_def_var(cdfid, name, nf_float, ndims, dims, vid)
      if (ierr .ne. nf_noerr) then
        print*, ' Error in variable declaration ', name
        stop
      end if

      ierr = nf_put_att_text(cdfid, vid, 'long_name', lngstr(lname),
     &                       lname)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in attrib"
      if (lngstr(units) .ne. 0) then
        ierr = nf_put_att_text(cdfid, vid, 'units', lngstr(units),
     &                         units)
        if (ierr .ne. nf_noerr) stop "***  netCDF error in attrib"
      end if

      ierr = nf_put_att_real(cdfid, vid, 'missing_value',
     &                       nf_float, 1, real(-7777777.0, 4))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in attrib"

      return
      end

c-----------------------------------------------------------------
      function lngstr( string )
      character string*(*)
      ilen = len(string)
      do 100 lngstr=ilen,1,-1
        if ( string(lngstr:lngstr) .ne. ' ' ) go to 99
  100 continue
      lngstr = 0
   99 continue
      return
      end
c-----------------------------------------------------------------
      subroutine accum_chk(lname,ihist)
c     Check whether an attempt was made to save an accumulated variable
c     in the second history file
      character lname*50

      if ( ihist.eq.2 ) then
         print*,
     &'ERROR - Can''t save accumulated variable in second history file'
         print*, lname
         stop
      end if

      return
      end
