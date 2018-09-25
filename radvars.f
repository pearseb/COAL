c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: radvars.f,v $
c Revision 1.1  2000/11/14 04:20:09  rot032
c Initial revision
c
      subroutine radvars(csqr,u,v,qtg,dqgdt,ipass,pttg,exptyp)

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /RVARSG/ )
!$OMP THREADPRIVATE ( /WORK1/ )

C Global parameters
      include 'PARAMS.f'
      real t00ec,t1ec,t0ec
      parameter (t00ec=288.0,t1ec=0.6652*t00ec,t0ec=t00ec-t1ec)

C Argument list
      real csqr
      real u(ln2,nl)
      real v(ln2,nl)
      real qtg(ln2,nl)
      real dqgdt(ln2,nl)
      integer ipass
      real pttg(ln2,nl) ! Only defined when ipass=2
      character*50 exptyp 

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'
      include 'RVARSG.f'
      include 'WORK1.f'

C Global data blocks
      include 'FEWFLAGS.f'
      include 'WORKA.f'

C Local work arrays and variables
      real xttg(ln2,nl) ! Dummy for ttg or pttg depending upon ipass

      integer k
      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

      do 40 k=1,nl
         do 40 mg=1,ln2
         u(mg,k)=un(mg,k)*csqr/muf(mg,k)
 40      v(mg,k)=vn(mg,k)*csqr/muf(mg,k)

c Convert from new (Chen) temperature variable ten to physical temperature ttg
c if chenflag is set and we are not on the first timestep of a restart from
c an old restart file.

CHEN
      If(chenflag.and.exptyp(42:45).eq.'CHEN')Then
      do 42 k=1,nl
         do 42 mg=1,ln2
 42      xttg(mg,k)=ten(mg,k)/muf(mg,k)+t0ec+t1ec*pdp00k(mg,k)
CHEN
      Else
      do 44 k=1,nl
         do 44 mg=1,ln2
 44      xttg(mg,k)=ten(mg,k)/muf(mg,k)+tmean(k)
      End If

      if(ipass.eq.1)then
         do k=1,nl
         do mg=1,ln2
            ttg(mg,k)=xttg(mg,k) ! common rvarsg
         enddo
         enddo
      else
         do k=1,nl
         do mg=1,ln2
            pttg(mg,k)=xttg(mg,k) ! argument list
         enddo
         enddo
      endif
       
      do 50 k=1,nl
         do 50 mg=1,ln2
         dqgdt(mg,k)=rmnt(mg,k)
 50      qtg(mg,k)=rmn(mg,k)

      return

      entry radvars2(qtg)

c Transfer back from ttg() into ten(), and from qtg() into rmn()

      If(chenflag)Then
      do 690 k=1,nl
      do 690 mg=1,ln2
  690   ten(mg,k)=(ttg(mg,k)-t0ec-t1ec*pdp00k(mg,k))*muf(mg,k)
      Else
      do 710 k=1,nl
      do 710 mg=1,ln2
  710   ten(mg,k)=(ttg(mg,k)-tmean(k))*muf(mg,k)
      End If

      do 620 k=1,nl
      do 620 mg=1,ln2
 620    rmn(mg,k)=qtg(mg,k)

      return
      end
