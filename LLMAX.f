c $Log: LLMAX.f,v $
c Revision 1.3  1995/10/17 05:08:27  ldr
c Added a descriptive comment.
c
c Revision 1.2  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.1  1993/12/06  16:54:20  ldr
c Initial revision
c
c Statement function giving the max LL for rhomboidal/triangular truncation

      integer mm, llmax
      llmax(mm)=lw       !Rhomboidal
c     llmax(mm)=lw1-mm   !Triangular
c Note: On the vector machines, use rhomboidal setting even for T63.
