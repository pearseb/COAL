c $Log: prsmap.f,v $
c Revision 1.6  2001/02/22 05:34:41  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.5  1999/05/20 06:23:52  rot032
c HBG changes to V5-2
c
c Revision 1.4  1993/12/17  15:33:26  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.3  92/04/23  12:11:25  ldr
c Diagnostic routines generalized by HBG for R42.
c 
c Revision 1.2  91/03/13  12:59:30  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:37:48  ldr
c Initial release V3-0
c 
      subroutine prsmap

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CHMAP.f'

C Local work arrays and variables
      character*1 chx(lon,lat2)

      integer kl
      integer lg
      integer loops
      integer lx
      integer ma
      integer maxnum
      integer mg
      integer mga
      integer mgb
      integer ns

C Local data, functions etc
      character*1 code(60,8)
      character*60 co60(8)
      equivalence (co60(1),code(1,1))
      data co60(8)
     &/'  ice and mlo temp difference from sst display.             '/
      data co60(7)
     &/'     *** code for this ice/mlo map ***                      '/
      data co60(6)
     &/'f=freezing (ice point),m=melting (ice point).               '/
      data co60(5)
     &/'rest are mlo pts: -= no diff from interp sst (i.sst) ;      '/
      data co60(4)
     &/'1,2,3,4,5,6,7,8,9= pts colder by tenths of deg from i.sst ; '/
      data co60(3)
     &/'t=pnt colder by more than ten tenths of deg from i.sst ;    '/
      data co60(2)
     &/'a,b,c,d,e,g,h,i,j= pts warmer by tenths of deg from i.sst ; '/
      data co60(1)
     &/'w=pnt warmer by more than ten tenths of deg from i.sst .    '/

C Start code : ----------------------------------------------------------

      do 10 lx=1,lat2
         ns=2-(lx-1)/lat
         lg=(ns-1)*lx+(lat2+1-lx)*(2-ns)
         do 10 mg=1,lon
         ma=mg+(ns-1)*lon
 10      chx(mg,lx)=i10d(ma,lg)

      do 15 lx=lat-3,lat+4
         do 15 mg=1,60
 15      chx(mg,lx)=code(mg,lx-lat+4)

c---- print out has max 128 numbers per line
      maxnum=128
      loops=(lon-1)/maxnum+1
      maxnum=lon/loops

      do 60 kl=1,loops
         mgb=kl*maxnum
         mga=mgb-maxnum+1
         if(kl.eq.loops)mgb=lon
         print *
         write (6,20) mga,mgb
 20      format ('* (torig-tg)*10 mlo map * : mg=',i3,' to ',i3)
c**** PRINT MLO TEMP DIFFERENCE*10 ARRAY
         do 30 lx=lat2,1,-1
            lg=lx
            if(lx.gt.lat)lg=lat2+1-lx
 30         write (6,40) lg,(chx(mg,lx),mg=mga,mgb)
 40      format (i2,1x,128a1)
 60   continue

      return
      end
