c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c Fix the calls to NINT, resolving warnings issued by the g95 Fortran compiler.
c SJP 2009/04/15
c
c Specify the precision of all integer arrays and real variables that are
c used as arguments to netCDF routines. This ensures portability of the code
c across 32- and 64-bit platforms.
c SJP 2001/12/13
c
c Updated from v2.4 to v3 of netCDF.
c SJP 2001/11/16
c
c "integer ncvid" removed from each subroutine and "include 'netcdf.inc'"
c added to histw1, histw2 and histw3.
c SJP 2001/11/15
c
c $Log: hist_save.f,v $
c Revision 1.21  2001/02/22 05:34:36  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.20  2001/02/12 05:39:42  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.19  2000/11/14 03:11:35  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.18  1999/05/20 06:23:48  rot032
c HBG changes to V5-2
c
c Revision 1.17  1998/12/10  00:55:27  ldr
c HBG changes to V5-1-21
c
c Revision 1.16  1997/12/23  00:23:33  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.15  1997/12/17  23:22:43  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.14.1.1  1997/12/19  02:03:09  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.14  1996/10/24  01:03:26  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.13  1996/03/21  03:19:14  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.12  1996/01/02  22:30:39  mrd
c Fixed error in the flags used to control saving of atmospheric
c temperature and humidity.
c
c Revision 1.11  1995/10/04  06:38:16  mrd
c Changes for more efficient writing of history files.
c
c Revision 1.10  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.6.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.9  1994/09/13  09:51:42  mrd
c Add daily mean 10m windspeed and 2m relative humidity as history variables.
c
c Revision 1.8  94/09/09  14:53:57  mrd
c Added daily average cloud, surface temperature and screen temperature
c as history variables.
c 
c Revision 1.7  94/09/09  14:14:45  mrd
c Added possbility of saving two history files at different frequencies.
c 
c Revision 1.6  93/11/30  11:20:07  mrd
c Set proper max/min for out of range data.
c 
c Revision 1.5  93/10/08  09:14:49  mrd
c Added maximum and minimum daily surface temperature to history.
c 
c Revision 1.4  93/10/07  17:23:12  mrd
c Removed path from netcdf.inc include file.
c 
c Revision 1.3  93/10/07  17:18:24  mrd
c Added flags to control writing of individual variables to daily history.
c 
c Revision 1.2  93/08/31  17:22:04  ldr
c Put /surf1 in include file to avoid problems with declaration of mc.
c 
c Revision 1.1  93/08/18  15:59:13  mrd
c Initial revision
c 
c 
      subroutine hist_save(ihist, mstep, nrad)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'RDPARM.f'
      real radtodeg
      parameter(radtodeg=57.29577951)

C Argument list
      integer ihist
      integer mstep
      integer nrad

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'HIST.f'
      include 'LSMI.f'
      include 'GAUSL.f'
      include 'UVPGD.f'
      include 'RMGRID.f'
      include 'netcdf.inc'
      include 'SURF1.f'
      real statsice
      common/masiv5/statsice(ln2,14,lat)

C Local work arrays and variables
      real*4 ylat(lat2)
      real*4 ylatc(lat2/nhcomp)
      integer surftype(lon,lat2)
      real tmp(lon,lat2)
      real pginv(lon,lat2)
      real windfac(lon,lat2)

      character vname*3

      integer ierr
      integer id
      integer k
      integer lg
      integer lgns
      integer ma
      integer mg
      integer ns

      real csqr
      real fac
      real facr
      real facp
      real*4 fdays

C Local data, functions etc

C Start code : ----------------------------------------------------------

c Declare local copy of histset(ihist)

      if ( histset(ihist).eq.1 ) then
        ierr = nf_inq_varid(histid(ihist), 'lat', id)
        if (ierr .ne. nf_noerr) stop "***  netCDF error in hist_save"
	do lg=1,lat
	  ylat(lg) = real((-radtodeg * acos(sia(lg))), 4)
	  ylat(lat2+1-lg) = -ylat(lg)
        end do

        if(nhcomp.gt.1)then ! Save on reduced grid

         do lg=1,lat2/nhcomp
           ylatc(lg)=ylat(nhcomp*lg)
         enddo
         ierr = nf_put_var_real(histid(ihist), id, ylatc)
         if (ierr .ne. nf_noerr) stop "***  netCDF error in hist_save"

        else
     
         ierr = nf_put_var_real(histid(ihist), id, ylat)
         if (ierr .ne. nf_noerr) stop "***  netCDF error in hist_save"

        endif

c       Surface height only needs to be written once.
        if ( all_hflg(ihist) .or. zht_hflg(ihist) ) then
           call histw1('zht',1,ihist)
        end if
      end if

c     Surface type
      if ( all_hflg(ihist) .or. imsl_hflg(ihist) ) then
         do lgns=1,lat2
            ns=2-(lgns-1)/lat
            lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)
            do mg=1,lon
               ma=mg+(ns-1)*lon
               surftype(mg,lgns) = imsl(ma,lg)
            enddo
         enddo
         call histw4('imsl',surftype,ihist)
      end if

c     The sea-level pressure is already in the correct latitudinal order
      if ( all_hflg(ihist) .or. psl_hflg(ihist) ) then
        call histw3('psl',pslgrid,ihist)
      endif

c     The temperature is already in the correct latitudinal order
      if ( all_hflg(ihist) .or. t_hflg(ihist) ) then
         do k=1,nl
            write(vname,'(a,i2.2)') 't', k
            call histw3(vname,tgrid(1,1,k),ihist)
         end do
      end if

c     Write the winds. These arrays are already in the 
c     correct latitudinal order. The pressure scaling and latitude weighting
c     also has to be removed.

      if ( all_hflg(ihist) .or. u_hflg(ihist) .or. v_hflg(ihist) .or.
     &     q_hflg(ihist) ) then
         do lg=1,lat2
            do mg=1,lon
               pginv(mg,lg) = 1./pgd(mg,lg)
            end do
         end do
      end if

      if ( all_hflg(ihist) .or. u_hflg(ihist) .or. v_hflg(ihist) )
     &     then
         do lgns=1,lat2
            lg=lgns
            if(lgns.gt.lat)lg=lat2+1-lg
            csqr=sqrt(acsq(lg)*eradsq)
            do mg=1,lon
               windfac(mg,lgns) = pginv(mg,lgns)*csqr
            end do
         end do
      end if

      if ( all_hflg(ihist) .or. u_hflg(ihist) ) then
         do k=1,nl
            do lg=1,lat2
               do mg=1,lon
                  tmp(mg,lg) = ugd(mg,lg,k)*windfac(mg,lg)
               end do
            end do
            write(vname,'(a,i2.2)') 'u', k
            call histw3(vname,tmp,ihist)
         end do
      end if

      if ( all_hflg(ihist) .or. v_hflg(ihist) ) then
         do k=1,nl
            do lg=1,lat2
               do mg=1,lon
                  tmp(mg,lg) = vgd(mg,lg,k)*windfac(mg,lg)
               end do
            end do
            write(vname,'(a,i2.2)') 'v', k
            call histw3(vname,tmp,ihist)
         end do
      end if

c     The mixing ratio does not have pressure weighting
      if ( all_hflg(ihist) .or. q_hflg(ihist) ) then
         do k=1,nl
            do lgns=1,lat2
               ns=2-(lgns-1)/lat
               lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)
               do mg=1,lon
                  ma=mg+(ns-1)*lon
                  tmp(mg,lgns) = rmg(ma,k,lg)
               end do
            end do
            write(vname,'(a,i2.2)') 'q', k
            call histw3(vname,tmp,ihist)
         end do
      end if

      if ( all_hflg(ihist) .or. psf_hflg(ihist) ) then
         call histw3('psf',pgd,ihist)
      end if

c     Write the instantaneous fields from the savegrid array
      if ( all_hflg(ihist). or. tsu_hflg(ihist) )
     &     call histw1('tsu',3,ihist)
      if ( all_hflg(ihist). or. tb2_hflg(ihist) )
     &     call histw1('tb2',11,ihist)
      if ( all_hflg(ihist). or. tb3_hflg(ihist) )
     &     call histw1('tb3',9,ihist)
      if ( all_hflg(ihist). or. wfg_hflg(ihist) )
     &     call histw1('wfg',8,ihist)
      if ( all_hflg(ihist). or. wfb_hflg(ihist) )
     &     call histw1('wfb',10,ihist)
      if ( all_hflg(ihist). or. snd_hflg(ihist) )
     &     call histw1('snd',12,ihist)
      if ( all_hflg(ihist). or. sid_hflg(ihist) )
     &     call histw1('sid',13,ihist)

c     Write the fields accumulated in hist_acc
c     These have different scale factors.
      facp = 1440. / float(hist_interval(ihist))
      fac = float(mstep) / float(hist_interval(ihist))
      facr = float(mstep*nrad) / float(hist_interval(ihist))

      if ( all_hflg(ihist). or. rnd_hflg(ihist) )
     &     call histw2('rnd',rain_a,facp,ihist)
      if ( all_hflg(ihist). or. rnc_hflg(ihist) )
     &     call histw2('rnc',precc_a,facp,ihist)
      if ( all_hflg(ihist). or. evp_hflg(ihist) )
     &     call histw2('evp',evap_a,facp,ihist)
      if ( all_hflg(ihist). or. hfl_hflg(ihist) )
     &     call histw2('hfl',fg_a,fac,ihist)
      if ( all_hflg(ihist). or. sgn_hflg(ihist) )
     &     call histw2('sgn',sg_a,fac,ihist)
      if ( all_hflg(ihist). or. sgc_hflg(ihist) )
     &     call histw2('sgc',sgclr_a,facr,ihist)
      if ( all_hflg(ihist). or. sgd_hflg(ihist) )
     &     call histw2('sgd',sgdn_a,facr,ihist)
      if ( all_hflg(ihist). or. sit_hflg(ihist) )
     &     call histw2('sit',sint_a,facr,ihist)
      if ( all_hflg(ihist). or. sot_hflg(ihist) )
     &     call histw2('sot',sout_a,facr,ihist)
      if ( all_hflg(ihist). or. soc_hflg(ihist) )
     &     call histw2('soc',soutclr_a,facr,ihist)
      if ( all_hflg(ihist). or. rgn_hflg(ihist) )
     &     call histw2('rgn',rg_a,fac,ihist)
      if ( all_hflg(ihist). or. rgc_hflg(ihist) )
     &     call histw2('rgc',rgclr_a,facr,ihist)
      if ( all_hflg(ihist). or. rgd_hflg(ihist) )
     &     call histw2('rgd',rgdn_a,facr,ihist)
      if ( all_hflg(ihist). or. rtu_hflg(ihist) )
     &     call histw2('rtu',rt_a,facr,ihist)
      if ( all_hflg(ihist). or. rtc_hflg(ihist) )
     &     call histw2('rtc',rtclr_a,facr,ihist)
      if ( all_hflg(ihist). or. tax_hflg(ihist) )
     &     call histw2('tax',taux_a,fac,ihist)
      if ( all_hflg(ihist). or. tay_hflg(ihist) )
     &     call histw2('tay',tauy_a,fac,ihist)
      if ( all_hflg(ihist). or. tsc_hflg(ihist) )
     &     call histw3('tsc',tscrn_sav,ihist)
      if ( all_hflg(ihist). or. tsh_hflg(ihist) )
     &     call histw3('tsh',tscrn_max,ihist)
      if ( all_hflg(ihist). or. tsl_hflg(ihist) )
     &     call histw3('tsl',tscrn_min,ihist)
      if ( all_hflg(ihist). or. tsca_hflg(ihist) )
     &     call histw2('tsca',tscrn_a,fac,ihist)
      if ( all_hflg(ihist). or. tgl_hflg(ihist) )
     &     call histw3('tgl',tg_min,ihist)
      if ( all_hflg(ihist). or. tgh_hflg(ihist) )
     &     call histw3('tgh',tg_max,ihist)
      if ( all_hflg(ihist). or. tsua_hflg(ihist) )
     &     call histw2('tsua',tg_a,fac,ihist)
      if ( all_hflg(ihist). or. run_hflg(ihist) )
     &     call histw2('run',runoff_a,facp,ihist)
      if ( all_hflg(ihist). or. int_hflg(ihist) )
     &     call histw2('int',cint_a,facp,ihist)
      if ( all_hflg(ihist). or. clda_hflg(ihist) )
     &     call histw2('clda',cld_a,facr,ihist)
      if ( all_hflg(ihist). or. rhsa_hflg(ihist) )
     &     call histw2('rhsa',rhscrn_a,fac,ihist)
      if ( all_hflg(ihist). or. v10ma_hflg(ihist) )
     &     call histw2('v10ma',v10m_a,fac,ihist)

      if ( all_hflg(ihist). or. tgg_hflg(ihist) ) then
         do lgns=1,lat2
            ns=2-(lgns-1)/lat
            lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)
            do mg=1,lon
               ma=mg+(ns-1)*lon
               tmp(mg,lgns) = tggsl(ma,1,lg)
            end do
         end do
         call histw3('tgg',tmp,ihist)
      end if

      if ( all_hflg(ihist). or. tgf_hflg(ihist) ) then
         do lgns=1,lat2
            ns=2-(lgns-1)/lat
            lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)
            do mg=1,lon
               ma=mg+(ns-1)*lon
               tmp(mg,lgns) = tgf(ma,lg)
            end do
         end do
         call histw3('tgf',tmp,ihist)
      end if

      if ( all_hflg(ihist). or. pev_hflg(ihist) )
     &     call histw2('pev',pev_a,facp,ihist)
      if ( all_hflg(ihist). or. sev_hflg(ihist) )
     &     call histw2('sev',sev_a,facp,ihist)

      if ( all_hflg(ihist). or. ico_hflg(ihist) ) then
         do lgns=1,lat2
            ns=2-(lgns-1)/lat
            lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)
            do mg=1,lon
               ma=mg+(ns-1)*lon
               tmp(mg,lgns) = statsice(ma,8,lg)
            end do
         end do
         call histw3('ico',tmp,ihist)
      end if

c     Time in days since start of run.
      fdays = real((float(histset(ihist) * hist_interval(ihist))
     &             / 1440.0), 4)
      ierr = nf_inq_varid(histid(ihist), 'time', id)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in hist_save"
      ierr = nf_put_var1_real(histid(ihist), id, histset(ihist), fdays)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in hist_save"
      histset(ihist) = histset(ihist) + 1

      return
      end
c-----------------------------------------------------------------------
      subroutine histw1(name, index, ihist)

c     Write history fields from the savegrid array.

      implicit none

C Global parameters
      include 'PARAMS.f'

C Argument list
      character name*(*)
      integer index
      integer ihist

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'HIST.f'
      include 'MASIV4.f'
      include 'netcdf.inc'

C Local work arrays and variables
      integer id
      integer start(3), count(3)
      real*4 ipack(lon,lat2)
      real*4 ipackc(lon/nhcomp, lat2/nhcomp)
      real datagr(lon,lat2)

      integer ierr
      integer lg
      integer lgns
      integer ma
      integer mg
      integer mgcomp
      integer ns

C Local data, functions etc

C Start code : ----------------------------------------------------------

      do lgns=1,lat2
         ns=2-(lgns-1)/lat
         lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)
         do mg=1,lon
            ma=mg+(ns-1)*lon
            datagr(mg,lgns) = savegrid(ma,index,lg)
         end do
      end do
    
      if(nhcomp.gt.1)then ! Save on reduced grid

      start(1) = 1
      start(2) = 1
      start(3) = histset(ihist)
      count(1) = lon/nhcomp
      count(2) = lat2/nhcomp
      count(3) = 1
      ierr = nf_inq_varid(histid(ihist), name, id)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw1"

      do lgns=1,lat2/nhcomp
         do mg=1,lon/nhcomp
            mgcomp=nhcomp*(mg-1)+1
            ipackc(mg,lgns) = real(datagr(mgcomp, nhcomp*lgns), 4)
         end do
      end do

      ierr = nf_put_vara_real(histid(ihist), id, start, count, ipackc)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw1"

      else

      start(1) = 1
      start(2) = 1
      start(3) = histset(ihist)
      count(1) = lon
      count(2) = lat2
      count(3) = 1
      ierr = nf_inq_varid(histid(ihist), name, id)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw1"

      do lgns=1,lat2
         do mg=1,lon
            ipack(mg,lgns) = real(datagr(mg, lgns), 4)
         end do
      end do

      ierr = nf_put_vara_real(histid(ihist), id, start, count, ipack)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw1"

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine histw2(name, grid, fac, ihist)

c     Scale and write given fields to history.

      implicit none

C Global parameters
      include 'PARAMS.f'

C Argument list
      character name*(*)
      real grid(lon,lat2)
      real fac
      integer ihist

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'HIST.f'
      include 'netcdf.inc'

C Local work arrays and variables
      integer id
      integer start(3), count(3)
      real*4 ipack(lon,lat2)
      real*4 ipackc(lon/nhcomp,lat2/nhcomp)

      integer ierr
      integer lg
      integer mg
      integer mgcomp

C Local data, functions etc

C Start code : ----------------------------------------------------------
 
      if(nhcomp.gt.1)then ! Save on reduced grid

      start(1) = 1
      start(2) = 1
      start(3) = histset(ihist)
      count(1) = lon/nhcomp
      count(2) = lat2/nhcomp
      count(3) = 1
      ierr = nf_inq_varid(histid(ihist), name, id)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw2"

      do lg=1,lat2/nhcomp
         do mg=1,lon/nhcomp
            mgcomp=nhcomp*(mg-1)+1
            ipackc(mg,lg) = real(grid(mgcomp, nhcomp*lg)*fac, 4)
         end do
      end do

      ierr = nf_put_vara_real(histid(ihist), id, start, count, ipackc)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw2"

      else

      start(1) = 1
      start(2) = 1
      start(3) = histset(ihist)
      count(1) = lon
      count(2) = lat2
      count(3) = 1
      ierr = nf_inq_varid(histid(ihist), name, id)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw2"

      do lg=1,lat2
	 do mg=1,lon
	    ipack(mg,lg) = real(grid(mg, lg)*fac, 4)
	 end do
      end do

      ierr = nf_put_vara_real(histid(ihist), id, start, count, ipack)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw2"

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine histw3(name, grid, ihist)

c     Write given field to history.

      implicit none

C Global parameters
      include 'PARAMS.f'

C Argument list
      character name*(*)
      real grid(lon,lat2)
      integer ihist

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'HIST.f'
      include 'netcdf.inc'

C Local work arrays and variables
      integer id
      integer start(3), count(3)
      real*4 ipack(lon,lat2)
      real*4 ipackc(lon/nhcomp,lat2/nhcomp)

      integer ierr
      integer lg
      integer mg
      integer mgcomp

C Local data, functions etc

C Start code : ----------------------------------------------------------
 
      if(nhcomp.gt.1)then ! Save on reduced grid

      start(1) = 1
      start(2) = 1
      start(3) = histset(ihist)
      count(1) = lon/nhcomp
      count(2) = lat2/nhcomp
      count(3) = 1
      ierr = nf_inq_varid(histid(ihist), name, id)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw3"

      do lg=1,lat2/nhcomp
         do mg=1,lon/nhcomp
            mgcomp=nhcomp*(mg-1)+1
            ipackc(mg,lg) = real(grid(mgcomp, nhcomp*lg), 4)
         end do
      end do

      ierr = nf_put_vara_real(histid(ihist), id, start, count, ipackc)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw3"

      else

      start(1) = 1
      start(2) = 1
      start(3) = histset(ihist)
      count(1) = lon
      count(2) = lat2
      count(3) = 1
      ierr = nf_inq_varid(histid(ihist), name, id)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw3"

      do lg=1,lat2
	 do mg=1,lon
	    ipack(mg,lg) = real(grid(mg, lg), 4)
	 end do
      end do

      ierr = nf_put_vara_real(histid(ihist), id, start, count, ipack)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw3"

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine histw4(name,igrid,ihist)

c     Write given integer field to history.

      implicit none

C Global parameters
      include 'PARAMS.f'

C Argument list
      character name*(*)
      integer igrid(lon,lat2)
      integer ihist

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'HIST.f'
      include 'netcdf.inc'

C Local work arrays and variables
      integer id
      integer start(3), count(3)
      integer igridc(lon/nhcomp,lat2/nhcomp)

      integer ierr
      integer lg
      integer mg
      integer mgcomp

C Local data, functions etc

C Start code : ----------------------------------------------------------
 
      if(nhcomp.gt.1)then ! Save on reduced grid

      start(1) = 1
      start(2) = 1
      start(3) = histset(ihist)
      count(1) = lon/nhcomp
      count(2) = lat2/nhcomp
      count(3) = 1
      ierr = nf_inq_varid(histid(ihist), name, id)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw4"

      do lg=1,lat2/nhcomp
         do mg=1,lon/nhcomp
            mgcomp=nhcomp*(mg-1)+1
            igridc(mg,lg) = igrid(mgcomp,nhcomp*lg)
         end do
      end do

      ierr = nf_put_vara_int(histid(ihist), id, start, count,
     &                       int(igridc, 4))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw4"

      else

      start(1) = 1
      start(2) = 1
      start(3) = histset(ihist)
      count(1) = lon
      count(2) = lat2
      count(3) = 1
      ierr = nf_inq_varid(histid(ihist), name, id)
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw4"
      ierr = nf_put_vara_int(histid(ihist), id, start, count,
     &                       int(igrid, 4))
      if (ierr .ne. nf_noerr) stop "***  netCDF error in histw4"

      endif

      return
      end
