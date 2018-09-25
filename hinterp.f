c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c $Log: hinterp.f,v $
c Revision 1.10  2001/10/12 02:13:46  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.9  2001/06/04 02:27:05  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.8  1997/12/17 23:23:06  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.7  1996/10/24  01:02:50  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.6  1996/06/13  02:06:40  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.5  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.4  93/12/23  15:26:41  ldr
c Fix up Bermejo-Staniforth for cloud water/tracers.
c 
c Revision 1.3  93/12/17  15:32:42  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.2  93/09/15  17:09:44  ldr
c Changes to get inital version of cloud water scheme going.
c 
c Revision 1.1  93/02/03  11:17:28  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:  from arguments
c                 idel - Integer gridpoint no. of x departure point. 
c                 jdel - Integer gridpoint no. of y departure point.
c                 reuse - if T, the factors fx1, fx2, fx3, fx4 can be reused
c                         in subsequent calls
c                 xg - Fractional displacement of x departure pt. from idel.
c                 yg - Fractional displacement of y departure pt. from idel.
c
c     In/Out: from arguments
c                 s - New field after horizontal advection
c 
c 
      subroutine hinterp(nwig,reuse,xg,yg,idel,jdel,
     &     s)

!$OMP THREADPRIVATE ( /SAVEFACS/ )

c Interpolate to find value of s at departure point
c This has y interp via grid number. Bi-cubic method.

c INPUTS:
c nwig - Controls scheme used to dealt with wiggles (0 for none, 1 for BS, 2 for LDR)
c idel - Integer gridpoint no. of x departure point.
c xg   - Fractional displacement of x departure pt. from idel.
c jdel - Integer gridpoint no. of y departure point.
c yg   - Fractional displacement of y departure pt. from idel.

c OUTPUT:
c s     - New field after horizontal advection

      include 'PARAMS.f'
      parameter(il=lon,jl=lat2,il2=il+2,jl2=jl+2)
      parameter(il4=il+4,jl4=jl+4)
      dimension s(il4,jl4),temp(il4,jl4)
      real r(4,il4,jl4)
      dimension xg(il4,jl4),yg(il4,jl4),idel(il4,jl4),jdel(il4,jl4)
      common / savefacs / fx1(il4,jl4),fx2(il4,jl4),fx3(il4,jl4),
     &                    fx4(il4,jl4)
      logical reuse

      call extend(s)

c Lagrangian interpolation: see Abramowitz and Stegun (eqn. 25.2.13)
c Use the expanded form for the x interpolation so that the factors
c fx1,fx2,fx3,fx4 can be reused in subsequent calls (i.e. if reuse=.true.).
c This could also be done for the y interpolation, but it is perhaps not worth
c the memory, as it does only 25% as much work.

      if ( .not. reuse ) then
        do 3 j=3,jl2
        do 3 i=3,il2
        fx1(i,j) = -xg(i,j) * (xg(i,j)-1) * (xg(i,j)-2) / 6
        fx2(i,j) = (xg(i,j)**2-1) * (xg(i,j)-2) / 2
        fx3(i,j) = -xg(i,j) * (xg(i,j)+1) * (xg(i,j)-2) / 2
 3      fx4(i,j) = xg(i,j) * (xg(i,j)**2-1) / 6
      endif

      do 4 j=3,jl2
      do 4 i=3,il2
      do 4 nn=1,4
      c1=s(idel(i,j)-1,jdel(i,j)+nn-2)
      c2=s(idel(i,j)  ,jdel(i,j)+nn-2)
      c3=s(idel(i,j)+1,jdel(i,j)+nn-2)
      c4=s(idel(i,j)+2,jdel(i,j)+nn-2)
 4    r(nn,i,j) = fx1(i,j)*c1 + fx2(i,j)*c2 + fx3(i,j)*c3 + fx4(i,j)*c4

c     r       ={(1-x     )*{(2-x     )*[(1+x     )*c2-x     *c1/3]
c      -x     *(1+x     )*c4/3}
c       +x    *(1+x     )*(2-x     )*c3}/2
      do 7 j=3,jl2
      do 7 i=3,il2
 7    temp(i,j)=((1-yg(i,j))*((2-yg(i,j))*((1+yg(i,j))*r(2,i,j)
     & -yg(i,j)*r(1,i,j)/3)-yg(i,j)*(1+yg(i,j))*r(4,i,j)/3)
     & +yg(i,j)*(1+yg(i,j))*(2-yg(i,j))*r(3,i,j))/2


c Bermejo-Staniforth positive definite scheme (not used for water vapour): 
c Force s to lie between min(c1,c2,c3,c4) and max(c1,c2,c3,c4)

      if(nwig.eq.1)then
        do 8 j=3,jl2
        do 8 i=3,il2
        c1=s(idel(i,j),jdel(i,j))
        c2=s(idel(i,j)+1,jdel(i,j))
        c3=s(idel(i,j),jdel(i,j)+1)
        c4=s(idel(i,j)+1,jdel(i,j)+1)
 8      temp(i,j)=min(max(temp(i,j),min(c1,c2,c3,c4)),max(c1,c2,c3,c4))

c Alternative linear scheme... Not used by default.

      elseif(nwig.eq.2)then
        do 805 j=3,jl2
        do 805 i=3,il2
        c1=s(idel(i,j),jdel(i,j))
        c2=s(idel(i,j)+1,jdel(i,j))
        c3=s(idel(i,j),jdel(i,j)+1)
        c4=s(idel(i,j)+1,jdel(i,j)+1)
          temp(i,j) = c1 + (c2-c1)*xg(i,j) + (c3-c1)*yg(i,j)
     &                + (c4-c3-c2+c1)*xg(i,j)*yg(i,j)
 805    continue  

      endif

      do 81 j=3,jl2
      do 81 i=3,il2
 81   s(i,j)= temp(i,j)

      return

      end
