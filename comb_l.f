c $Log: comb_l.f,v $
c Revision 1.1  2001/02/22 05:36:47  rot032
c Initial revision
c
      subroutine comb_1l(ompl,pl,arr,parr)
      
      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      real ompl(ln2)
      real pl(ln2)
      real arr(ln2)
      real parr(ln2)

C Local work arrays and variables
      integer mg

C Start code : ----------------------------------------------------------

C Combine the non-leads and leads components (1-level)

      do mg=1,ln2
         arr(mg)=ompl(mg)*arr(mg)+pl(mg)*parr(mg)
      enddo

      return
      end
C---------------------------------------------------------------------
      subroutine comb_nl(ompl,pl,arr,parr)
      
      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      real ompl(ln2)
      real pl(ln2)
      real arr(ln2,nl)
      real parr(ln2,nl)

C Local work arrays and variables
      integer k
      integer mg

C Start code : ----------------------------------------------------------

C Combine the non-leads and leads components (nl levels)

      do k=1,nl
      do mg=1,ln2
         arr(mg,k)=ompl(mg)*arr(mg,k)+pl(mg)*parr(mg,k)
      enddo
      enddo

      return
      end
