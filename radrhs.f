c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: radrhs.f,v $
c Revision 1.1  2000/11/14 04:20:09  rot032
c Initial revision
c
      subroutine radrhs(pg,eg,fg,ttg,qtg,tscrn,rhscrn)

      implicit none

!$OMP THREADPRIVATE ( /HYBRPR/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      real pg(ln2)
      real eg(ln2)
      real fg(ln2)
      real ttg(ln2,nl)
      real qtg(ln2,nl)
      real tscrn(ln2)
      real rhscrn(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'HYBRPR.f'

C Global data blocks

C Local work arrays and variables
      integer mg

      real qs
      real qscrn
      real theta1

C Local data, functions etc
      include 'ESTABL.f'

C Start code : ----------------------------------------------------------

C To calculate sreen level RH

      do mg=1,ln2
        if(abs(fg(mg)).gt.1.e-8)then
          theta1=ttg(mg,1)/pdpsk(mg,1)
          qscrn=qtg(mg,1)+(tscrn(mg)-theta1)*eg(mg)/(hlcp*fg(mg))
        else
          qscrn=qtg(mg,1)
        endif
        qs=qsat(100.0*pg(mg),tscrn(mg))
        rhscrn(mg)=min(max(0., 100.*qscrn/qs), 100.)
      enddo

      return
      end
