c $Log: prmlomap.f,v $
c Revision 1.7  2001/02/22 05:34:45  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.6  1999/05/20 06:23:57  rot032
c HBG changes to V5-2
c
c Revision 1.5  1994/07/11  12:31:38  ldr
c Tidy up RCS header to avoid annoying warning on SGI.
c
c Revision 1.4  93/12/23  15:31:26  ldr
c Changes to V4-4-54l from HBG for coupled model of V4-5
c 
c Revision 1.3  92/04/23  12:12:32  ldr
c Diagnostic routines generalized by HBG for R42.
c 
c Revision 1.2  91/05/20  15:20:47  ldr
c *** empty log message ***
c 
      subroutine prmlomap

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f'
      include 'QFLXDAT.f'

C Local work arrays and variables
      character*2 itch(lon)

      integer i
      integer ka
      integer lg
      integer loops
      integer lx
      integer maxnum
      integer mg
      integer mga
      integer mgb

C Local data, functions etc
      character*2 nums(42)
      character*42 numm,nump
      equivalence (numm,nums(1)), (nump,nums(22))
      data numm/'-K-J-I-H-G-F-E-D-C-B-A-9-8-7-6-5-4-3-2-1 0'/
      data nump/' 1 2 3 4 5 6 7 8 9101112131415161718192021'/

C Start code : ----------------------------------------------------------

c---- max 64 numbers per line
      maxnum=64
      loops=(lon-1)/maxnum+1
      do 70 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         mgb=min(lon,mgb)
         if(qflux)then
           write (6,10) mga,mgb
 10      format (1h1,'T(MLO)-T*(INTERPOLATED) : mg=',i3,' to ',i3)
         else
           write (6,11) mga,mgb
         end if
 11      format (1h1,'T(OCEAN MODEL)-T*(BMO) : mg=',i3,' to ',i3)
         do 30 lg=1,lat
            ka=0
            do 20 mg=mga,mgb
               i=itdiff(mg,lg)
               ka=ka+1
               if (i.lt.149) itch(ka)=nums(i+21)
               if (i.eq.150) itch(ka)='II'
               if (i.eq.9999) itch(ka)='**'
 20         continue
 30         write (6,40) lg,(itch(mg),mg=1,ka)
 40      format (1h ,i2,1x,64a2)
         do 60 lg=lat,1,-1
            ka=0
            do 50 mg=mga,mgb
               i=itdiff(mg+lon,lg)
               ka=ka+1
               if (i.lt.149) itch(ka)=nums(i+21)
               if (i.eq.150) itch(ka)='II'
               if (i.eq.9999) itch(ka)='**'
 50         continue
 60         write (6,40) lg,(itch(mg),mg=1,ka)
 70   continue
      return
      end
