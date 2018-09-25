c $Log: prtcd.f,v $
c Revision 1.9  2001/02/22 05:34:39  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.8  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.7  1998/01/30  04:36:15  ldr
c Fix from SPO for Cd diagnostic.
c
c Revision 1.6  1997/07/24  05:42:50  ldr
c Tuning of ukconv version and T63 fixes from HBG
c
c Revision 1.5  1993/12/03  16:42:51  ldr
c write cdmap to output file like in filest (SPO).
c
c Revision 1.4  92/04/23  12:12:11  ldr
c Diagnostic routines generalized by HBG for R42.
c 
c Revision 1.3  91/05/15  11:09:02  ldr
c Some diagnostic code transferred here from csiro9.f
c 
c Revision 1.2  91/03/13  12:59:33  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:37:49  ldr
c Initial release V3-0
c 
      subroutine prtcd

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FILES.f'
      include 'TIMEX.f'

      real cdavg,pcdavg
      common/cdmapx/cdavg(ln2,lat),pcdavg(ln2,lat)

      real statsice
      common/masiv5/statsice(ln2,14,lat)

C Local work arrays and variables
      real ajcdavg(ln2,lat)
      character*3 varname

      integer lg
      integer mg

      real avg
      real pl

C Local data, functions etc

C Start code : ----------------------------------------------------------

      avg=1000.0/nsteps
      do 10 lg=1,lat
      do 10 mg=1,ln2
         pcdavg(mg,lg)=pcdavg(mg,lg)*avg
 10       cdavg(mg,lg)= cdavg(mg,lg)*avg
      do 25 lg=1,lat
      do 25 mg=1,ln2
        pl=statsice(mg,8,lg)
        ajcdavg(mg,lg)=(1.0-pl)*cdavg(mg,lg)+pl*pcdavg(mg,lg)
 25   continue

c write cdmap to ouput file. See ncinit.f to create nc file
c ncinit now set up to initialize for average cd stats if cdmap eq true
c     varname='cdm'
c     call colldy(varname,1.0,str,cdavg)
c     varname='pcd'
c     call colldy(varname,1.0,str,pcdavg)
      varname='acd'
      call colldy(varname,1.0,str,ajcdavg)
      return
      end
