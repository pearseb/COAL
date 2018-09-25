c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/10
c
c $Log: co2_read.f,v $
c Revision 1.12  2001/02/28 04:36:36  rot032
c Further tidy ups from HBG
c
c Revision 1.11  2001/02/22 05:34:37  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.10  1998/12/10 00:55:29  ldr
c HBG changes to V5-1-21
c
c Revision 1.9  1998/05/27  02:07:35  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.8  1997/12/17  23:22:44  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.7  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.6  93/12/17  15:31:53  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.5  93/07/20  10:08:42  ldr
c Use block data for initializing data in common, to keep the VP happy.
c 
c Revision 1.4  92/12/09  14:43:07  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.3  92/08/31  16:49:17  ldr
c Modified to use variable filename.
c 
c Revision 1.2  92/05/11  15:13:23  ldr
c Put Include PARAMS.f in source file rather than in RDPARM.f, 
c to avoid nested includes.
c 
c Revision 1.1  92/04/15  12:03:51  mrd
c Initial revision
c 
      block data co2_blk

      include 'PARAMS.f'
      include 'RDPARM.f'
      include 'CO2DTA.f'

c   The following coeffiecients don't depend on resolution or CO2 conc.
      data b0,b1,b2,b3/-.51926410e-4,-.18113332e-3,
     & -.10680132e-5,-.67303519e-7/

      end

c******************************************************************************

      subroutine co2_read(sigma)

      implicit none

!$OMP THREADPRIVATE ( /RADISW/ )

C Global parameters
      include 'PARAMS.f'
      include 'RDPARM.f'
      real sigtol
      parameter(sigtol=1e-3)

C Argument list
      real sigma(nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'RADISW.f'

C Global data blocks
      include 'CO2DTA.f'
      include 'FILES.f'

C Local work arrays and variables
      real sigin(nl)

      integer ierr
      integer k
      integer nlev

C Local data, functions etc

C Start code : ----------------------------------------------------------

c  This routine reads the CO2 transmission coefficients from the
c  co2_datafile (filename set in namelist)

      open(unit=7,file=co2_datafile,status='old',iostat=ierr)
      call errcheck(ierr,'co2_datafile     ','co2_read  ')

      read(7,*) nlev
c     Check that the number of levels is the same
      if ( nlev.ne.nl ) then
	  print*, ' ERROR - Number of levels wrong in co2_data file'
	  stop
      end if
c     Check that the sigma levels are the same
c     Note that the radiation data has the levels in the reverse order
      read(7,*) (sigin(k),k=nl,1,-1)
      do k=1,nl
	  if ( abs(sigma(k)-sigin(k)) .gt. sigtol ) then
	      print*, ' ERROR - sigma level wrong in co2_data file'
	      print*, k, sigma(k), sigin(k)
	      stop
          end if
      end do
      read(7,*) rrvco2
      write(*,*) ' CO2 mixing ratio is ', rrvco2*1e6,' ppmv'
      read(7,*) stemp 
      read(7,*) gtemp 
      read(7,*) cdt51
      read(7,*) co251
      read(7,*) c2d51
      read(7,*) cdt58
      read(7,*) co258
      read(7,*) c2d58
      read(7,*) cdtm51
      read(7,*) co2m51
      read(7,*) c2dm51
      read(7,*) cdtm58
      read(7,*) co2m58
      read(7,*) c2dm58
      read(7,*) cdt31
      read(7,*) co231
      read(7,*) c2d31
      read(7,*) cdt38
      read(7,*) co238
      read(7,*) c2d38
      read(7,*) cdt71
      read(7,*) co271
      read(7,*) c2d71
      read(7,*) cdt78
      read(7,*) co278
      read(7,*) c2d78
      read(7,*) co211
      read(7,*) co218
      close(7)

      return
      end
