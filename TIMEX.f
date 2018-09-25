c nstepsa counter for execution step counting is added to the COMMON
c block
c AJA 2009/03/19
c
c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c $Log: TIMEX.f,v $
c Revision 1.7  2001/02/12 05:39:54  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.6  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.5  1997/07/24 05:59:17  ldr
c Introduce logical variable called "start".
c
c Revision 1.4  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.3  1993/12/17  15:31:35  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.2  91/06/11  14:13:21  ldr
c Bugfix for coupled runs so that output files are not overwritten each month.
c 
c Revision 1.1  91/03/14  09:49:54  ldr
c Initial revision
c 
      real dt, ratlm
      integer mstep, minw, kday, kdayp, nstep, ndays, kdays,
     &       month, ldays, nsteps, iyear, ncepstrt, nrad, nstepsa
      integer*8 mins
      logical start, dcflag
      common /timex/ dt, ratlm, mins, mstep, minw, kday, kdayp,
     &               nstep, ndays, kdays, month, ldays, nsteps,
     &               iyear, ncepstrt, nrad, start, dcflag, nstepsa
