c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c Modified for five-character experiment names.
c SJP 2009/04/22
c
c Modified to enable five-digit year numbers.
c SJP 2004/09/22
c
c $Log: openfl.f,v $
c Revision 1.28  2001/02/22 05:34:38  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.27  2001/02/12 05:39:47  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.26  2000/11/14 03:11:36  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.25  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.24  1998/05/27  02:07:37  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.23  1996/06/13  02:07:14  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.22  1996/04/17  01:32:52  ldr
c Changes from TIE to allow initialization of monthly mean netcdf files
c from within model.
c
c Revision 1.21  1996/03/21  03:28:07  ldr
c A couple of housekeeping tid-ups from LDR.
c
c Revision 1.20  1995/05/02  07:25:58  ldr
c Merge HBG V4-5-43h2 stuff into latest version.
c
c Revision 1.19  94/09/09  14:14:47  mrd
c Added possbility of saving two history files at different frequencies.
c 
c Revision 1.18  94/08/08  17:21:54  ldr
c Strip off excessive RCS comments at top of file.
c 
c Revision 1.17.1.1  95/05/02  07:14:37  ldr
c Other changes from HBG that went into coupled run but not into V4-5-43h.
c 
c Revision 1.17  93/12/17  15:33:18  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.16  93/09/30  12:40:27  ldr
c Changes to add option for automatic naming of output restart files.
c 
c Revision 1.15  93/08/10  16:13:51  ldr
c Merge of HBG changes to V4-3-28l with other changes since then.
c 
c Revision 1.14  93/08/10  12:32:07  mrd
c Modified history to use netCDF
c 
c Revision 1.13  93/08/06  11:43:01  ldr
c Introduce new flag savefcor to control saving of flux correction data.
c 
c Revision 1.12.1.1  93/08/10  15:27:44  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c 
c Revision 1.12  93/07/23  15:02:44  mrd
c Added option of writing a daily history file.
c 
c Revision 1.11  93/06/18  15:23:41  ldr
c Removed unnecessary check for presence of SST file.
c 
c Revision 1.10  92/08/31  17:00:13  ldr
c /files/ in include statement.
c 
      subroutine openfl(exptyp)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      character*50 exptyp

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f'
      include 'FILES.f'
      include 'HIST.f'
      include 'NEST.f'
      include 'TIMEX.f'

C Local work arrays and variables
      character*4 chtst
      character*17 fname

C Local data, functions etc

C Start code : ----------------------------------------------------------

c SET UP FILE NAMES

      write(str,'(i5.5,i2.2,a1,a5)')iyear,month,'.',runtype


c Set up output restart file names if using automatic option

      if ( filewrflag .and. autoname ) then
        print*,'Setting output restart filenames automatically'
        orfilename='rest'//str
        if(orfilename.eq.irfilename)then
          print*,'WARNING: input restart file will be overwritten'
        endif
      endif

C****  OPEN THE HISTORY FILES FOR THIS RUN

      if(nestflag)then
       open(20,file='nest'//str,form='unformatted') !For nested model
      endif

      chtst='    ' ! dummy variable

      if(saveqflux)then
C**      FILE55 : IMPLIED OCEAN HEATING DATA (WRITTEN AT END OF
C**               MONTH IN timet.f)
        open(55,file='qflx'//str,form='unformatted') !Qfluxes

        write(55)exptyp
        write(55)kday,kdayp,chtst
      endif

      if(rainflag)then
        open(56,file='rain'//str,form='formatted')  !Daily rain

        write(56,9030)exptyp
        write(56,9040)kday,kdayp,chtst
 9030 format(1x,a50)
 9040 format(1x,2i6,2x,a4)
      endif

      if(tempflag)then
        open(57,file='tgmx'//str,status='unknown')  !Daily max surftemps
        open(58,file='tgmn'//str,status='unknown')  !Daily min surftemps
        open(64,file='tsmx'//str,status='unknown')  !Daily max screentem
        open(65,file='tsmn'//str,status='unknown')  !Daily min screentem
      endif

      if(saveglmean)then
        open(59,file='glmn'//str,status='unknown') !Global means
      endif

      if(sdiagflag)then
        open(60,file='spts'//str,status='unknown')  !Surface diagnostics
      endif                                         !for selected points

      if(uvtflag)then
       open(61,file='uday'//str,form='unformatted') !Daily zonal winds
       open(62,file='vday'//str,form='unformatted') !Daily merid winds
       open(63,file='tday'//str,form='unformatted') !Daily temps
      endif

      if(savegbrf)then
        open(69,file='gbrf'//str,form='formatted') !Monthly special data
      endif

      if(savefcor)then
        open(70,file='fcor'//str,form='unformatted') !Flux corrections
      endif

      if(savehist(1))then
	call openhist(1, iyear, month)
      endif
      if(savehist(2))then
	call openhist(2, iyear, month)
      endif

      if ( statsflag )then
c Call ncinit to specify which data to be compressed at T63 and to set up
c netCDF files.
        call ncinit
      endif   

      if (sltrace) then
        fname='trac'//str
        open(unit=89,file=fname,form='formatted')
c       if(mw.eq.22)then ! R21 points only
c         fname='trcg'//str
c         open(unit=87,file=fname,form='formatted')
c         fname='tras'//str
c         open(unit=86,file=fname,form='formatted')
c       endif
      end if

      return
      end
