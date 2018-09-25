c Removed unnecessary "include 'netcdf.inc'".
c SJP 2001/11/16
c
c $Log: hist_cld.f,v $
c Revision 1.14  2001/02/12 05:39:53  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.13  2000/11/14 03:11:39  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.12  1999/12/20 23:17:41  dix043
c Fix double row indexing of cloud level arrays.
c
c Revision 1.11  1999/05/20 06:23:55  rot032
c HBG changes to V5-2
c
c Revision 1.10  1998/12/10  00:55:50  ldr
c HBG changes to V5-1-21
c
c Revision 1.9  1997/12/17  23:22:54  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.8  1996/10/24  01:03:25  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.7  1996/06/13  02:08:46  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.6  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.3.2.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.5  1994/09/09  14:14:43  mrd
c Added possbility of saving two history files at different frequencies.
c
c Revision 1.4  94/09/09  09:58:20  mrd
c Removed some unused array declarations.
c 
c Revision 1.3  93/10/07  17:22:41  mrd
c Removed path from netcdf.inc include file.
c 
c Revision 1.2  93/10/07  17:18:23  mrd
c Added flags to control writing of individual variables to daily history.
c 
c Revision 1.1  93/08/18  15:13:23  mrd
c Initial revision
c
c     INPUT/OUTPUT
c     Input:   from common/hist in HIST.f
c                  cldgrid  - total cloud
c                  cllgrid,clmgrid,clhgrid  - low,mid,high level cloud
c                  ichgrid,icmgrid,ictgrid,icbgrid - cloud levels
c                  alsgrid  - albedo
c
c     Input:   from arguments
c                  ihist - history file no, either 1 or 2 
c 
      subroutine hist_cld(ihist)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer ihist

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'HIST.f'

C Local work arrays and variables

C Local data, functions etc

C Start code : ----------------------------------------------------------

C  This history data is stored on radsteps only

c     Cloud amounts
      if ( all_hflg(ihist). or. cld_hflg(ihist) ) then
        call histw3('cld',cldgrid,ihist)
      end if

      if ( all_hflg(ihist). or. cll_hflg(ihist) ) then
        call histw3('cll',cllgrid,ihist)
      end if

      if ( all_hflg(ihist). or. clm_hflg(ihist) ) then
        call histw3('clm',clmgrid,ihist)
      end if

      if ( all_hflg(ihist). or. clh_hflg(ihist) ) then
        call histw3('clh',clhgrid,ihist)
      end if

c     Cloud levels (integer data)
      if ( all_hflg(ihist). or. ich_hflg(ihist) ) then
        call histw4('ich',ichgrid,ihist)
      end if

      if ( all_hflg(ihist). or. icm_hflg(ihist) ) then
        call histw4('icm',icmgrid,ihist)
      end if

      if ( all_hflg(ihist). or. icb_hflg(ihist) ) then
        call histw4('icb',icbgrid,ihist)
      end if

      if ( all_hflg(ihist). or. ict_hflg(ihist) ) then
        call histw4('ict',ictgrid,ihist)
      end if

c     Albedo
      if ( all_hflg(ihist). or. als_hflg(ihist) ) then
        call histw3('als',alsgrid,ihist)
      end if

      return
      end
