c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c $Log: rads_l.f,v $
c Revision 1.1  2001/02/22 05:36:47  rot032
c Initial revision
c
      subroutine rads_zer(lg)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer lg

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'QCLOUD1.f'  ! cfrad,qcrad,qlrad,qfrad

C Local work arrays and variables
      integer k
      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

c---- Zero the arrays for summing cloud variables for use 
c----  by the radiation code radfs
      do k=1,nl
        do mg=1,ln2
          cfrad(mg,k,lg)=0.
          qcrad(mg,k,lg)=0.
          qlrad(mg,k,lg)=0.
          qfrad(mg,k,lg)=0.
        enddo
      enddo

      return
      end
C---------------------------------------------------------------------
      subroutine rads_acc(lg,cfsav,qccon,qlgsav,qfgsav)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer lg
      real cfsav(ln2,nl)
      real qccon(ln2,nl)
      real qlgsav(ln2,nl)
      real qfgsav(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'QCLOUD1.f'  ! cfrad,qcrad,qlrad,qfrad

C Local work arrays and variables
      integer k
      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

c---- Accumulate quantities for passing into radiation
      do k=1,nl-1
        do mg=1,ln2
          cfrad(mg,k,lg)=cfrad(mg,k,lg)+cfsav(mg,k)
          qcrad(mg,k,lg)=qcrad(mg,k,lg)+qccon(mg,k)
          qlrad(mg,k,lg)=qlrad(mg,k,lg)+qlgsav(mg,k)
          qfrad(mg,k,lg)=qfrad(mg,k,lg)+qfgsav(mg,k)
        enddo
      enddo

      return
      end
C---------------------------------------------------------------------
      subroutine rads_avg(lg,cfsav,qccon,qlgsav,qfgsav)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer lg
      real cfsav(ln2,nl)
      real qccon(ln2,nl)
      real qlgsav(ln2,nl)
      real qfgsav(ln2,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'QCLOUD1.f'  ! cfrad,qcrad,qlrad,qfrad
      include 'TIMEX.f'    ! mins,nrad

C Local work arrays and variables
      integer k
      integer mg
      integer nradcv

C Local data, functions etc

C Start code : ----------------------------------------------------------

c---- Average quantities for passing into radiation
      if(mod(mins,1440_8).eq.0_8)then
        nradcv=1
      else
        nradcv=nrad
      endif

      do k=1,nl-1
        do mg=1,ln2
          cfsav(mg,k)=cfrad(mg,k,lg)/nradcv
          if(cfsav(mg,k).lt.0.001)cfsav(mg,k)=0.
          qccon(mg,k)=qcrad(mg,k,lg)/nradcv
          qlgsav(mg,k)=qlrad(mg,k,lg)/nradcv
          qfgsav(mg,k)=qfrad(mg,k,lg)/nradcv
        enddo
      enddo

      return
      end
