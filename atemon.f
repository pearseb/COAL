c Fixing a bug in the handling of netCDF errors.
c SJP 2009/04/27
c
c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c Modified for five-character experiment names.
c SJP 2009/04/22
c
c Added IMPLICIT NONE statement, plus variable declarations for TDT, NDADD,
c NRUNDAY, IDAYPSAV and MSTEPSAV.
c SJP 2007/05/28
c
c Modified for five-digit year numbers.
c SJP 2004/09/29
c
c Declaration added for IERR.
c SJP 2001/12/10
c
c Updated from v2.4 to v3 of netCDF and "include 'netcdf.inc'" added.
c SJP 2001/11/21
c
c $Log: atemon.f,v $
c Revision 1.16  2000/08/16 02:59:24  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.15  1999/05/20 06:23:48  rot032
c HBG changes to V5-2
c
c Revision 1.14  1998/05/26  04:48:54  ldr
c Final changes (inc. ACH salinity stuff) for 10k run.
c
c Revision 1.13  1996/11/19  06:05:53  ldr
c Make call filewr the last thing that's done.
c
c Revision 1.12  1996/10/24  01:03:24  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.11  1996/06/13  02:08:43  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.10  1996/02/19  04:10:08  ldr
c Generalize for 24 levels.
c
c Revision 1.9  1994/09/09  14:15:24  mrd
c Added possbility of saving two history files at different frequencies.
c
c Revision 1.8  94/07/11  12:31:25  ldr
c Rearrange variables in common block /atdat to avoid misalignment warnings
c at 64 bit on SGI.
c 
c Revision 1.7  94/03/30  12:35:56  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.6  93/11/29  14:59:02  ldr
c Use icu_sflg and icv_sflg for treatment of ice U,V diagnostics, so they
c are consistent with other "Stats".
c 
c Revision 1.5  93/11/03  13:08:04  ldr
c Make locean a namelist flag and Replace silly flag zfl_sflg with savezflux.
c 
c Revision 1.4  93/10/05  13:07:56  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.3  93/08/17  11:58:20  ldr
c Add extra flag saveicuv to control dumping of SPO's U,V diagnostics.
c 
c Revision 1.2  93/08/10  16:14:00  ldr
c Merge of HBG changes to V4-3-28l with other changes since then.
c 
c Revision 1.1  93/08/10  15:03:41  ldr
c Initial revision
c 
c Revision 1.38  93/06/18  14:43:30  ldr
c
c     INPUT/OUTPUT:
c     Input:   from common/atdat in this subroutine
c                  endofmonth - counter
c
c              from common/fewflags in FEWFLAGS.f
c                  cdmap - if T, prtcd saves drag coefficient file
c                  filewrflag- if T, write restart file at end of run, 
c                               default F
c                  idyn      - if T and semice=T, leads=T, use dynamical 
c                              sea-ice model
c                  lcouple   - flag to run model in coupled mode, default F
c                  ltrace    - if T, calls Langrangian tracer routine, tracera
c                  savegbrf  - if T saves various global means
c                  statsflag - if T, model save monthly means of various fields
c                              covering the whole globe
c
c              from common/glmean in GLMEAN.f
c                  grunof    - global mean runoff
c
c              from common/hist_control in HIST.f
c                  histid   - id of netCDF file
c                  savehist - if T, saves detailed history of global fields in
c                             netCDF format
c
c              from common/stflags in STFLAGS.f
c                  icu_sflg - flag to save monthly mean zonal velocity of
c                             sea-ice
c                  icv_sflg - as above for meridional velocity of sea-ice
c
c     In/Out:  from common/files in FILES.f
c                  str - character string YYYYYMM.XXXXX where XXXXX is runtype
c
c              from common/timex in TIMEX.f
c                  nsteps - time step counter
c
c              from arguments
c                  exptyp - restart file header
c                  ijdyn  - ice dynamics time step factor
c
c 
      subroutine atemon(exptyp,ijdyn)

      implicit none

C Global parameters
      include 'PARAMS.f'

C Argument list
      character exptyp*50
      integer ijdyn

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f'
      include 'FILES.f'
      include 'GLMEAN.f'
      include 'HIST.f'
      include 'STFLAGS.f'
      include 'TIMEX.f'
      include 'netcdf.inc'

      real tdt
      integer ndadd,nrunday,idaypsav,mstepsav
      logical iday,endofmonth
      common/atdat/tdt,ndadd,nrunday,idaypsav,mstepsav,iday,endofmonth

C Local work arrays and variables
      integer ierr

C Local data, functions etc

C Start code : ----------------------------------------------------------

c     CALL TIMER('ATEMON  ',1)
C****
C**** END OF RUN - PRINT STATS (USUALLY PER MONTH)
C****

C*     IF END OF RUN ,UPDATE RESTART FILE ON DISK,DUMP DATA TO DISK,
C*     AND, IF DOING TRACER EXPERIMENT, OUTPUT POSITION OF PARTICLES

      if(.not.lcouple.or.lcouple.and.endofmonth)then
        if(ltrace) call traceout(1,nsteps)
        if(statsflag.and.cdmap)call prtcd
c       if(idyn.and.statsflag.and.(icu_sflg.or.icv_sflg))
c    &          call flatme(ijdyn)
        if(idyn.and.statsflag)call flatme(ijdyn)
        if(statsflag) call filest(str)
        write(6,6767)grunof
        if(savegbrf)write(69,6769)grunof
        if(filewrflag) call filewr(exptyp)
 6769   format(1x,e13.6)
 6767   format(1h ,'global mean runoff for month (mms)=',e13.6)
      endif

      if (savehist(1)) then
        ierr = nf_close(histid(1))
        if (ierr .ne. nf_noerr) stop "***  netCDF error in atemon"
      end if
      if (savehist(2)) then
        ierr = nf_close(histid(2))
        if (ierr .ne. nf_noerr) stop "***  netCDF error in atemon"
      end if

c     CALL TIMER('ATEMON  ',4)

      return
      end
