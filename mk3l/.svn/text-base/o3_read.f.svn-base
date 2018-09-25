c $Log: o3_read.f,v $
c Revision 1.6  2001/02/22 05:34:42  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.5  1998/05/27 02:07:38  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.4  1992/12/09  14:44:02  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c
c Revision 1.3  92/08/31  16:59:16  ldr
c Use variable filename for o3 read.
c 
c Revision 1.2  92/05/11  15:13:37  ldr
c Put Include PARAMS.f in source file rather than in RDPARM.f, 
c to avoid nested includes.
c 
c Revision 1.1  92/04/15  12:28:32  mrd
c Initial revision
c 
      subroutine o3_read(sigma)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'RDPARM.f'
      real sigtol
      parameter(sigtol=1.e-3)

C Argument list
      real sigma(nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FILES.f'

      real dduo3n,ddo3n2 ,ddo3n3,ddo3n4
      common /o3dat/ dduo3n(37,l),ddo3n2(37,l),ddo3n3(37,l),ddo3n4(37,l)

C Local work arrays and variables
      real sigin(nl)

      integer ierr
      integer k
      integer nlev

C Local data, functions etc

C Start code : ----------------------------------------------------------
c
c     Reads the ozone data from the o3_datafile (filename set in namelist)
c

      open(unit=7,file=o3_datafile,status='old',iostat=ierr)
      call errcheck(ierr,'o3_datafile      ','o3_read   ')

      read(7,*) nlev
      if ( nlev.ne.nl ) then
	  print*, ' ERROR - Number of levels wrong in o3_data file'
	  stop
      end if
c     Check that the sigma levels are the same
c     Note that the radiation data has the levels in the reverse order
      read(7,*) (sigin(k),k=nl,1,-1)
      do k=1,nl
	  if ( abs(sigma(k)-sigin(k)) .gt. sigtol ) then
	      print*, ' ERROR - sigma level wrong in o3_data file'
	      print*, k, sigma(k), sigin(k)
	      stop
          end if
      end do

c  Note that the data is written as MAM, SON, DJF, JJA. The arrays in
c  o3dat are in the order DJF, MAM, JJA, SON
      read(7,1000) ddo3n2
      read(7,1000) ddo3n4
      read(7,1000) dduo3n
      read(7,1000) ddo3n3
      close(unit=7)
 1000 format(9f8.5)

      return
      end
