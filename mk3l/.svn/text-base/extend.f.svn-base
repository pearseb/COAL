c $Log: extend.f,v $
c Revision 1.1  1992/08/06 16:36:16  ldr
c Initial revision
c
      subroutine extend(x)

c Extend arrays beyond edge of grid for Semi-Lagrangian advection

      include 'PARAMS.f'
      parameter(il=lon,jl=lat2,ilp=il+1,jlp=jl+1,il2=il+2,jl2=jl+2)
      parameter(il3=il+3,jl3=jl+3,il4=il+4,jl4=jl+4)
      parameter(ilh=il/2)
      dimension x(il4,jl4)
c     impose boundary symmetries first north-south
      do 28 i=3,ilh+2
      x(i,2)=x(i+ilh,3)
      x(i+ilh,2)=x(i,3)
      x(i,1)=x(i+ilh,4)
      x(i+ilh,1)=x(i,4)
      x(i,jl3)=x(i+ilh,jl2)
      x(i+ilh,jl3)=x(i,jl2)
      x(i,jl4)=x(i+ilh,jlp)
28    x(i+ilh,jl4)=x(i,jlp)
c     then east-west
      do 29 j=1,jl4
      x(2,j)=x(il2,j)
      x(il3,j)=x(3,j)
      x(1,j)=x(ilp,j)
29    x(il4,j)=x(4,j)
      return
      end
