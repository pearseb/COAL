c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c $Log: prdaily.f,v $
c Revision 1.9  2001/06/04 02:27:06  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.8  2001/02/22 05:34:45  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.7  2000/06/20 02:08:35  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.5  1999/05/20 06:23:57  rot032
c HBG changes to V5-2
c
c Revision 1.4  1993/07/05  11:33:01  ldr
c Corrected indenting in do 60 loop.
c
c Revision 1.3  93/07/02  12:58:01  ldr
c Changes to add diagnostics of mean, daily max and daily min tgg and tgf.
c 
c Revision 1.2  91/03/13  12:59:28  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:37:47  ldr
c Initial release V3-0
c 
      subroutine prdaily

c Routine to print out daily fields: Daily rain, max and min temps.
c Arrays are printed out in convenient form for plotting i.e. west to ea
c and south to north.
c Called from csiro9 at end of each day.   (LDR 7/90)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f'
      include 'TIMEX.f'

      real devap,drain,dvgg
      common/prtar/devap(ln2,lat),drain(ln2,lat),dvgg(ln2,lat)

      real brain
      common/rblock/brain(ln2,lat)

      real tgmax,tgmin,tsmax,tsmin,htsmax,htsmin
     & ,extsmax,extsmin,tggmax,tggmin,tgfmax,tgfmin
     & ,htggmax,htggmin,htgfmax,htgfmin,hpmc
      common/maxmin/tgmax(ln2,lat),tgmin(ln2,lat),
     &              tsmax(ln2,lat),tsmin(ln2,lat),
     &              htsmax(ln2,lat),htsmin(ln2,lat),
     &              extsmax(ln2,lat),extsmin(ln2,lat),
     &              tggmax(ln2,lat),tggmin(ln2,lat),
     &              tgfmax(ln2,lat),tgfmin(ln2,lat),
     &              htggmax(ln2,lat),htggmin(ln2,lat),
     &              htgfmax(ln2,lat),htgfmin(ln2,lat)
     &             ,hpmc(ln2,lat)

C Local work arrays and variables
C***      character*20 lfn

      integer lg
      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

c Output daily rain to its own file

      if(rainflag)then
        write(56,10)((brain(mg+lon,lg),mg=1,lon),lg=1,lat)
        write(56,10)((brain(mg,lg),mg=1,lon),lg=lat,1,-1)
 10   format(1p8e10.3)
      endif
C***        write(lfn,1616)iyear,month,ldays
C*** 1616 format('rain.y',i4.4,'m',i2.2,'d',i2.2)
C***        open(unit=51,file=lfn,form='formatted')
C***        write(51,6606)((brain(mg+lon,lg),mg=1,lon),lg=1,lat)
C***        write(51,6606)((brain(mg,lg),mg=1,lon),lg=lat,1,-1)
C*** 6606 format(11f6.1)
C***        close(unit=51)
C***        write(lfn,1617)iyear,month,ldays
C*** 1617 format('pmsl.y',i4.4,'m',i2.2,'d',i2.2)
C***        open(unit=51,file=lfn,form='formatted')
C***        write(51,6607)((dvgg(mg+lon,lg),mg=1,lon),lg=1,lat)
C***        write(51,6607)((dvgg(mg,lg),mg=1,lon),lg=lat,1,-1)
C*** 6607 format(11f7.1)
C***        close(unit=51)

c Output daily max and min temps to own file

      if(tempflag)then
        write(57,10)((tgmax(mg+lon,lg),mg=1,lon),lg=1,lat)
        write(57,10)((tgmax(mg,lg),mg=1,lon),lg=lat,1,-1)
        write(58,10)((tgmin(mg+lon,lg),mg=1,lon),lg=1,lat)
        write(58,10)((tgmin(mg,lg),mg=1,lon),lg=lat,1,-1)
        write(64,10)((tsmax(mg+lon,lg),mg=1,lon),lg=1,lat)
        write(64,10)((tsmax(mg,lg),mg=1,lon),lg=lat,1,-1)
        write(65,10)((tsmin(mg+lon,lg),mg=1,lon),lg=1,lat)
        write(65,10)((tsmin(mg,lg),mg=1,lon),lg=lat,1,-1)
      endif

c Add daily max and min tscreen values to monthly totals
c Also update monthly extreme tscreen values
      
      if(statsflag)then
        do 70 lg=1,lat
        do 70 mg=1,ln2
          htsmax(mg,lg)=htsmax(mg,lg)+tsmax(mg,lg)
          htsmin(mg,lg)=htsmin(mg,lg)+tsmin(mg,lg)
          extsmax(mg,lg)=max(extsmax(mg,lg),tsmax(mg,lg))
          extsmin(mg,lg)=min(extsmin(mg,lg),tsmin(mg,lg))
c Add daily max and min bare ground and veg. ground values to monthly totals
          htggmax(mg,lg)=htggmax(mg,lg)+tggmax(mg,lg)
          htggmin(mg,lg)=htggmin(mg,lg)+tggmin(mg,lg)
          htgfmax(mg,lg)=htgfmax(mg,lg)+tgfmax(mg,lg)
          htgfmin(mg,lg)=htgfmin(mg,lg)+tgfmin(mg,lg)
 70     continue 
      endif

c Zero brain and tgmax arrays. Set tgmin to a nice large value.

      do 90 lg=1,lat
      do 90 mg=1,ln2
          brain(mg,lg)=0.0
          tgmax(mg,lg)=0.0
          tgmin(mg,lg)=999.
          tsmax(mg,lg)=0.0
          tsmin(mg,lg)=999.
 90   continue

      return
      end
