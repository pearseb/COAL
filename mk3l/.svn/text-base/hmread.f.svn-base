c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Note that this subroutine is now only called if QFLUX is .TRUE., and that
c much of it is therfore redundant.
c SJP 2003/050/29
c
c $Log: hmread.f,v $
c Revision 1.11  2001/02/22 05:34:38  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.10  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.9  1998/05/27  02:07:36  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.8  1997/12/03  05:35:42  ldr
c Put IGW's qflux code with variable MLO depth into main version of model.
c
c Revision 1.7  1993/12/17  11:51:37  ldr
c Changes to V4-4-45l from HBG for coupled model
c
c Revision 1.6  93/11/29  11:38:28  ldr
c Changes to V4-4-32l from HBG for coupled model
c 
c Revision 1.5  93/10/05  13:06:21  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.4  93/09/26  18:38:43  ldr
c Corrected setting up of hmo array.
c 
c Revision 1.3  93/07/20  10:01:52  ldr
c Use free format read so that it works for R42 too.
c 
c Revision 1.2  93/07/14  17:51:25  ldr
c Tidied up printed message from this routine.
c 
c Revision 1.1  93/07/12  14:13:25  ldr
c Initial revision
c 
c 
      subroutine hmread(mnth)

c subroutine to read in mixed layer depths for two months
c ready for interpolation.

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer mnth

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f'
      include 'TIMEX.f'

      real hmo,hmo2
      common/hm/hmo(ln2,lat),hmo2(lon,lat2,2)

C Local work arrays and variables
      real arlon(lon),arlat(lat2)
      character*80 h2,fm

      integer i
      integer ierr
      integer j
      integer lg
      integer m
      integer ma
      integer mon
      integer nlat
      integer nlon
      integer ns

C Local data, functions etc

C Start code : ----------------------------------------------------------

c Code for IGW's qflux version with variable MLO depth...
c Only set up for R21 at present

      if ( qflux ) then

        if(lw.ne.22) stop 'error in hmread'

        fm ='mlod'
        open(unit=7,file=fm,form='formatted',status='old',iostat=ierr)
        call errcheck(ierr,'mlo data file    ','hmread    ')
        read(7,*) nlat 
        read(7,81) arlat
        read(7,*) nlon 
        read(7,81) arlon
        read(7,'(80a)') h2
        print *,h2
        read(7,81) ((hmo2(i,j,1),i=1,nlon),j=1,nlat)
 81     format(6e12.5)

        close(7)

        do 681 j=1,lat2
          ns=2-(j-1)/lat
          lg=(ns-1)*j+(lat2+1-j)*(2-ns)
          do 681 i=1,lon
          ma=i+(ns-1)*lon
 681      hmo(ma,lg)=hmo2(i,j,1)

      elseif(.not.qflux)then


c.... read in two  consecutive months (last month may be cyclic)
      do 500 mon=1,2
      m=mnth+mon-1
      if(m.gt.12)m=1
      write(6,*)'In hmread: Reading MLO depths for month ',m
      if(lw.eq.22)then
        write(fm,87)m
c87     format('ohm',i2.2)
 87     format('hhm',i2.2)
      elseif(lw.eq.43)then
        write(fm,88)m
 88     format('ohm',i2.2,'.R42')
      elseif(lw.eq.64)then
        write(fm,89)m
 89     format('ohmT63',i2.2)
      endif

      open(unit=7,file=fm,form='formatted',status='old',iostat=ierr)
      call errcheck(ierr,'mlo data file    ','hmread    ')

      read(7,50)h2
   50 format(80a)
      read(7,51)
   51 format(1x)
      read(7,769)(arlon(i),i=1,lon)
      do 767 j=1,lat2
  767 read(7,*)arlat(j),(hmo2(i,j,mon),i=1,lon)
  769 format(6x,8f8.2)

      close(7)
  500 continue

      endif

      return

      entry hmread1

c.... interpolate the 2 months of hmo data to the correct day

      do 780 j=1,lat2
      ns=2-(j-1)/lat
      lg=(ns-1)*j+(lat2+1-j)*(2-ns)
      do 781 i=1,lon
      ma=i+(ns-1)*lon
  781 hmo(ma,lg)=(1.0-ratlm)*hmo2(i,j,1)+ratlm*hmo2(i,j,2)
  780 continue 

      return
      end
