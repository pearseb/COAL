c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: radhyb.f,v $
c Revision 1.1  2000/11/14 04:20:09  rot032
c Initial revision
c
      subroutine radhyb

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )
!$OMP THREADPRIVATE ( /WORK1/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      real p00
      parameter (p00=1000.0)

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'
      include 'WORK1.f'

C Global data blocks
      include 'CNSTA.f'
      include 'FEWFLAGS.f'
      include 'HYBARR.f'

C Local work arrays and variables
      real pnk(ln2)

      integer k
      integer mg

C Local data, functions etc

C Start code : ----------------------------------------------------------

c---- hybrid stuff
c---- compute prh = pressure at half levels
      do mg=1,ln2
        prh(mg,1)=pn(mg)
        prh(mg,nlp)=0.0
      enddo
      do k=2,nl
      do mg=1,ln2
        prh(mg,k)=anh(k)+(bnh(k)*pn(mg))
      enddo
      enddo
c---- compute muf = dp/d(coord) for hybrid coords (=P* if sigma)
c---- compute prf = pressure at full levels
c---- compute dprf = pressure thickness of each level
      do k=1,nl
      do mg=1,ln2
        muf(mg,k)=dadnf(k)+(dbdnf(k)*pn(mg))
        prf(mg,k)=anf(k)+(bnf(k)*pn(mg))
        dprf(mg,k)=prh(mg,k)-prh(mg,k+1)
      enddo
      enddo
c---- compute pdpsk = (prf/P*)**cappa
c---- compute pdp00k = (prf/P00)**cappa
      if(hybrid)then
        do k=1,nl
        do mg=1,ln2
          pdpsk(mg,k)=(prf(mg,k)/pn(mg))**cappa
          pdp00k(mg,k)=(prf(mg,k)/p00)**cappa
        enddo
        enddo
      else
          do 41 mg=1,ln2
 41       pnk(mg)=(pn(mg)/p00)**cappa
        do k=1,nl
        do mg=1,ln2
          pdpsk(mg,k)=sigk(k)
          pdp00k(mg,k)=sigk(k)*pnk(mg)
        enddo
        enddo
      endif
c---- end of hybrid stuff

      return
      end
