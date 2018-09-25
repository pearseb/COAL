c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Transferred COMMON blocks to separate header files, as follows:
c /MDAY/  ->  MDAY.f
c SJP 2007/05/29
c
c $Log: timet.f,v $
c Revision 1.22  2001/02/22 05:34:38  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.21  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.20  1998/12/10  00:55:34  ldr
c HBG changes to V5-1-21
c
c Revision 1.19  1994/08/08  17:22:53  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.18  94/03/30  12:35:33  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.17  93/12/17  15:34:12  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.16  93/12/17  11:51:52  ldr
c Changes to V4-4-45l from HBG for coupled model
c 
c Revision 1.15  93/11/29  16:09:47  ldr
c Corrected length of common block hm.
c 
c Revision 1.14  93/08/17  11:58:46  ldr
c Corrections to new handling of multi-month runs.
c 
c Revision 1.13  93/08/10  16:13:58  ldr
c Merge of HBG changes to V4-3-28l with other changes since then.
c 
c Revision 1.12  93/07/29  12:46:40  ldr
c Corrected calculation of ratlm from previous revision.
c 
c Revision 1.11.1.1  93/08/10  15:28:11  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c 
c Revision 1.11  93/07/27  17:13:50  ldr
c Changes to ensure coherent approach to pure restarts, monthly
c interpolation factors..
c 
      subroutine timet

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f' 
      include 'LSMI.f'
      include 'TIMEX.f'
      include 'MDAY.f'

      real ochf,occur,g50dt
      integer khf
      common/curent/ochf(ln2,lat),occur(ln2,lat),khf(ln2,lat)
     &,g50dt(ln2,lat)

C Local work arrays and variables
c     integer ioc(32)

c     integer iocx
c     integer kk
      integer lg
c     integer loops
c     integer maxnum
c      integer mcur
      integer mg
c     integer mga
c     integer mgb

      real stepsm

C Local data, functions etc

C Start code : ----------------------------------------------------------

C**   TO UPDATE THE DAY COUNT AT THE COMPLETION OF EACH DAY.

      kdays=kdays+1
      ldays=ldays+1
      if(ldays.lt.mdays(month))go to 50
C****
C**** AVERAGE THE OCEAN HEATING OVER 1 MONTH
C****
c      mcur=month
C****
C***      ldays=0
C***      month=month+1
C***      if(month.gt.12)then
C***        month=1
C***        kdays=1
C***      endif

c**   SAVE QFLUXES.

      if(.not.lcouple)then
        if(saveqflux)then
          stepsm=(1440.0/mstep)*mdays(month)
          do 30 lg=1,lat
          do 30 mg=1,ln2
C**   Average the Q-Fluxes (occur) gathered over the month
C**   and the surface heating (ochf).
        ochf (mg,lg)=ochf (mg,lg)/stepsm
   30   occur(mg,lg)=occur(mg,lg)/stepsm
C**** PRINT THE IMPLIED OCEAN CURRENT FIELD ****
c         maxnum=32
c         loops=(lon-1)/maxnum+1
c         maxnum=lon/loops
c         do 664 kk=1,loops
c          mgb=kk*maxnum
c          mga=mgb-maxnum+1
c          if(kk.eq.loops)mgb=lon
c          write(6,659)mcur,mga,mgb
c659      format(1h1,'ocean currents for month ',i2,' mg=',i2,',',i2)
c          do 660 lg=1,lat
c           do 550 mg=mga,mgb
c           iocx=nint(occur(mg,lg))
c           if(imsl(mg,lg).eq.4)iocx=99999
c550        ioc(mg-mga+1)=iocx
c660       write(6,661)lg,ioc
c661      format(1h ,i2,1x,32i4)
c          do 662 lg=lat,1,-1
c           do 551 mg=mga,mgb
c           iocx=nint(occur(mg+lon,lg))
c           if(imsl(mg+lon,lg).eq.4)iocx=99999
c551        ioc(mg-mga+1)=iocx
c662       write(6,661)lg,ioc
c664      continue

C**** DUMP THE OCEAN CURRENT VALUES TO DISK FILE
          write(55)ochf
          write(55)occur
        endif
      endif
C****
 50   continue 
      ratlm=float(ldays)/mdays(month)
      write(6,6676)kdays,ratlm
 6676 format(1h ,' timet : kdays=',i5,' ratlm=',f7.4)
      return
      end
