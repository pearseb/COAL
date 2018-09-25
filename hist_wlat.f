c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c $Log: hist_wlat.f,v $
c Revision 1.17  2001/02/12 05:39:43  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.16  2000/11/14 03:11:35  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.15  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.14  1998/12/10  00:55:29  ldr
c HBG changes to V5-1-21
c
c Revision 1.13  1997/12/17  23:22:44  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.12  1996/10/24  01:03:23  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.11  1996/03/21  03:19:13  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.10  1995/10/04  01:34:49  mrd
c Changes for more efficient writing of history files.
c
c Revision 1.10  1995/10/04  01:28:36  mrd
c Changes for more efficient writing of history files.
c
c Revision 1.9  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.8  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.7  1994/09/09  14:14:46  mrd
c Added possbility of saving two history files at different frequencies.
c
c Revision 1.6  94/06/28  16:21:51  mrd
c Add sea-level pressure as a history variable.
c 
c Revision 1.5  93/11/30  11:53:33  mrd
c Set proper max/min for out of range data.
c 
c Revision 1.4  93/11/25  10:16:09  mrd
c FIxed error in zonal wind scale factor.
c 
c Revision 1.3  93/10/07  17:23:14  mrd
c Removed path from netcdf.inc include file.
c 
c Revision 1.2  93/10/07  17:18:20  mrd
c Added flags to control writing of individual variables to daily history.
c 
c Revision 1.1  93/08/18  15:17:23  mrd
c Initial revision
c 
c 
      subroutine hist_wlat(lg,radstep)

      implicit none

!$OMP THREADPRIVATE ( /GIANT3/ )
!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /PHYSICAL/ )
!$OMP THREADPRIVATE ( /RADISW/ )
!$OMP THREADPRIVATE ( /RVARSG/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'RDPARM.f'

C Argument list
      integer lg
      logical radstep

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'GIANT3.f'
      include 'HYBRPR.f'
      include 'PHYSICAL.f'
      include 'RADISW.f'
      include 'RVARSG.f'

C Global data blocks
      include 'HIST.f'
      include 'TIMEX.f'

C Local work arrays and variables
      real psl(ln2),plk(ln2,2),tlk(ln2,2) ! Pmsl work arrays
      integer ich(ln2),icm(ln2),icb(ln2),ict(ln2)
      
      integer k
      integer lgn
      integer mg

      real c
      real conf
      real conr
      real pxlev
      real txlev

C Local data, functions etc

C Start code : ----------------------------------------------------------

c To save global data fields for history files

      lgn=lat2p-lg

      IF ((mod(mins+int(mstep, 8),int(hist_interval(1), 8)).eq.0_8).or.
     &    (mod(mins+int(mstep, 8),int(hist_interval(2), 8)).eq.0_8))
     &  THEN

      if(all_hflg(1).or.all_hflg(2) .or. t_hflg(1).or.t_hflg(2))then
c   Save a global array of grid temperature
       do k=1,nl
        do mg=1,lon
         tgrid(mg,lgn,k) = ttg(mg,k)
         tgrid(mg,lg,k)  = ttg(mg+lon,k)
        end do
       end do
      end if

      if(all_hflg(1).or.all_hflg(2) .or. psl_hflg(1).or.psl_hflg(2))then
c   Save a global array of psl
        c=grav/6.5e-3
        conr=c/rdry
        do mg=1,ln2
         pxlev=pg(mg)-90.0
         k=1
  999    k=k+1
         if(pxlev.lt.prf(mg,k))go to 999
         plk(mg,2)=prf(mg,k)
         plk(mg,1)=prf(mg,k-1)
         tlk(mg,2)=ttg(mg,k)
         tlk(mg,1)=ttg(mg,k-1)
        end do
        do mg=1,ln2
         pxlev=pg(mg)-90.0
         txlev=(tlk(mg,2)*alog(pxlev/plk(mg,1)) +
     &         tlk(mg,1)*alog(plk(mg,2)/pxlev))
     &              /alog(plk(mg,2)/plk(mg,1))
         conf=(pxlev/pg(mg))**(rdry/c)  /c
         psl(mg)=pg(mg)*(1.+conf*grav*z4(mg)/txlev)   **conr
        enddo

        do mg=1,lon
         pslgrid(mg,lgn) = psl(mg)
         pslgrid(mg,lg)  = psl(mg+lon)
        end do
      endif

      ENDIF

c     Save the cloud data to the history file on the radiation step
c     before the rest of the data is saved.

      IF ((mod(mins+int(nrad*mstep, 8),
     &         int(hist_interval(1), 8)).eq.0_8).or.
     &    (mod(mins+int(nrad*mstep, 8),
     &         int(hist_interval(2), 8)).eq.0_8))   THEN

      if(radstep)then

c Cloud amounts
        do mg=1,lon
         cldgrid(mg,lgn) =  cl(mg)
         cldgrid(mg,lg)  =  cl(mg+lon)
         cllgrid(mg,lgn) = cll(mg)
         cllgrid(mg,lg)  = cll(mg+lon)
         clmgrid(mg,lgn) = clm(mg)
         clmgrid(mg,lg)  = clm(mg+lon)
         clhgrid(mg,lgn) = clh(mg)
         clhgrid(mg,lg)  = clh(mg+lon)
         alsgrid(mg,lgn) = als(mg)
         alsgrid(mg,lg)  = als(mg+lon)
        enddo

c Cloud levels
        do mg=1,ln2
         ich(mg) = nlp - kbtm(mg,4)
         icm(mg) = nlp - kbtm(mg,3)
         icb(mg) = nlp - kbtm(mg,2)
         ict(mg) = nlp - ktop(mg,2)
         if ( cll(mg) .eq. 0. ) then
          ich(mg) = icm(mg)
          icm(mg) = ict(mg)
          icb(mg) = 0
          ict(mg) = 0
         end if
         if ( clm(mg) .eq. 0. ) then
          ich(mg) = icm(mg)
          icm(mg) = 0
         end if
         if ( clh(mg) .eq. 0 ) ich(mg) = 0
        end do

        do mg=1,lon
         ichgrid(mg,lgn) = ich(mg)
         ichgrid(mg,lg)  = ich(mg+lon)
         icmgrid(mg,lgn) = icm(mg)
         icmgrid(mg,lg)  = icm(mg+lon)
         icbgrid(mg,lgn) = icb(mg)
         icbgrid(mg,lg)  = icb(mg+lon)
         ictgrid(mg,lgn) = ict(mg)
         ictgrid(mg,lg)  = ict(mg+lon)
        enddo

C Albedo
        do mg=1,lon
         alsgrid(mg,lgn) = als(mg)
         alsgrid(mg,lg)  = als(mg+lon)
        enddo

      end if

      ENDIF

      return
      end
