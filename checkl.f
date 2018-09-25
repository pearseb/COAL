c Added IMPLICIT NONE statement.
c SJP 2007/05/28
c
c $log$
      subroutine checkl(logarr,nel,flag)
c
c To check if any element of a logical array is true
c
      implicit none

      integer nel,mg
      logical logarr(nel),flag

      do mg=1,nel
       if(logarr(mg))go to 10
      enddo
      flag=.false.
      return

   10 flag=.true.
      return
      end
