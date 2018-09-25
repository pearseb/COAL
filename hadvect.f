c Removed unnecessary "include 'MACHINE.f'"
c SJP 2001/11/22
c
c $Log: hadvect.f,v $
c Revision 1.18  1999/05/20 06:23:54  rot032
c HBG changes to V5-2
c
c Revision 1.17  1998/01/30  05:05:37  ldr
c Further parallelization stuff for NEC from MRD.
c
c Revision 1.16  1996/10/24  01:02:49  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.15  1996/06/13  02:06:39  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.14  1996/03/21  03:18:45  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.13  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.11.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.12  1994/08/08  17:21:25  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.11  93/12/17  15:32:41  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.10  93/10/07  12:09:49  ldr
c Move machine parameter statement into new include file MACHINE.f.
c 
c Revision 1.9  93/08/19  18:20:26  ldr
c Moved position of lock on initialization code.
c 
c Revision 1.8  93/08/19  15:08:12  ldr
c Minor cosmetic changes.
c 
c Revision 1.7  93/07/20  10:08:45  ldr
c Use block data for initializing data in common, to keep the VP happy.
c 
c Revision 1.6  93/02/03  11:16:29  ldr
c Split hinterp routine off into separate file.
c 
c Revision 1.5  93/01/26  16:26:43  ldr
c Changes to implement Semi-Lagrangian transport of tracers.
c 
c Revision 1.4  92/12/09  14:43:29  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.3  92/11/10  16:55:23  ldr
c Changes for SGI and some tidying up.
c 
      subroutine hadvect(u,v,dt,
     &                   xg,yg,idel,jdel)

c Horizontal component of JMcG's semi-Lagrangian advection scheme.
c Adapted for CSIRO9 model by LDR.
c Note: Indexing of 2D arrays is different in here - Grid runs from 3
c to il+2 and from 3 to jl+2, rather than from 1 to il and 1 to jl.

c INPUTS:
c u     - Zonal wind field in m/s (Grid runs from 3 to il2, 3 to jl2)
c v     - Meridional wind field in m/s ( " )
c s     - Field to be advected (Moisture for CSIRO9)
c dt    - Timestep for scheme (Passed in as tdt for CSIRO9 3 time level scheme)

c OUTPUTS:
c idel - Integer gridpoint no. of x departure point.
c xg   - Fractional displacement of x departure pt. from idel.
c jdel - Integer gridpoint no. of y departure point.
c yg   - Fractional displacement of y departure pt. from idel.

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      parameter(il=lon,jl=lat2,jlp=jl+1,il2=il+2,jl2=jl+2)
      parameter(jl3=jl+3,il4=il+4,jl4=jl+4)

C Argument list
      real u(il,jl)
      real v(il,jl)
      real dt
      real xg(il4,jl4)
      real yg(il4,jl4)
      integer idel(il4,jl4)
      integer jdel(il4,jl4)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      common/savetrig/rlat(jl4),rlong(il4),cs(jl4),sn(jl4),dlong
     & ,snlo(il4),cslo(il4),ucon(jl4),vcona(jl4),vconb(jl4),vconc(jl4)

C Local work arrays and variables
      real derx(il4,jl4),dery(il4,jl4),derz(il4,jl4)
     &  ,b1(il4,jl4),b2(il4,jl4),b3(il4,jl4)
      dimension x(il4,jl4),y(il4,jl4),z(il4,jl4),alat(il4,jl4)
      dimension along(il4,jl4)

C Local data, functions etc
      logical diag
      data diag,id,jd /.false.,20,10 /
      data nit/3/  !Order of scheme used to calculate departure points

C Start code : ----------------------------------------------------------

c Initialization: Extend latitude and longitude arrays by two at each end
c   rlat etc now set up in initax.f

c Calculate departure point

      do 2 j=3,jl2
      do 2 i=3,il2
      x(i,j)=erad*cs(j)*cslo(i)
      y(i,j)=erad*cs(j)*snlo(i)
      z(i,j)=erad*sn(j)
      derx(i,j)=( u(i-2,j-2)*snlo(i)+v(i-2,j-2)*sn(j)*cslo(i))*dt
      dery(i,j)=(-u(i-2,j-2)*cslo(i)+v(i-2,j-2)*sn(j)*snlo(i))*dt
   2  derz(i,j)=-v(i-2,j-2)*cs(j)*dt

      do 25 j=3,jl2
      do 25 i=3,il2
      x(i,j)=x(i,j)+derx(i,j)
      y(i,j)=y(i,j)+dery(i,j)
   25 z(i,j)=z(i,j)+derz(i,j)

c     if(diag)then
c       print*,'25: id, jd, x, y, z ',id,jd,x(id,jd),y(id,jd),z(id,jd)
c     endif

      do 32 itn=2,nit

c     impose boundary symmetries
      call extend(derx)
      call extend(dery)
      call extend(derz)
c     following is easier vectorized
      do 30 j=3,jl2
      do 30 i=3,il2
      b1(i,j)=-(u(i-2,j-2)*(derx(i+1,j)-derx(i-1,j))*ucon(j)
     & +v(i-2,j-2)*(derx(i,j-1)*vcona(j)+derx(i,j)*vconb(j)
     & +derx(i,j+1)*vconc(j) ))/itn
      b2(i,j)=-(u(i-2,j-2)*(dery(i+1,j)-dery(i-1,j))*ucon(j)
     & +v(i-2,j-2)*(dery(i,j-1)*vcona(j)+dery(i,j)*vconb(j)
     & +dery(i,j+1)*vconc(j) ))/itn
      b3(i,j)=-(u(i-2,j-2)*(derz(i+1,j)-derz(i-1,j))*ucon(j)
     & +v(i-2,j-2)*(derz(i,j-1)*vcona(j)+derz(i,j)*vconb(j)
     & +derz(i,j+1)*vconc(j) ))/itn
   30 continue

      do 31 j=3,jl2
      do 31 i=3,il2
      derx(i,j)=b1(i,j)
      dery(i,j)=b2(i,j)
      derz(i,j)=b3(i,j)
      x(i,j)=x(i,j)+derx(i,j)
      y(i,j)=y(i,j)+dery(i,j)
   31 z(i,j)=z(i,j)+derz(i,j)

   32 continue

c     bring back to earth's surface
      do 34 j=3,jl2
      do 34 i=3,il2
      rdiv=sqrt(x(i,j)**2+y(i,j)**2+z(i,j)**2)
      alat(i,j)=asin(z(i,j)/rdiv)
   34 along(i,j)=atan2(y(i,j),x(i,j))

c     if(diag)then
c       print*,'34: id, jd, x, y, z ',id,jd,x(id,jd),y(id,jd),z(id,jd)
c     endif

      do 38 j=3,jl2
      do 38 i=3,il2
c     xg first    ! (x departure point)
      xg(i,j)=(along(i,j)-rlong(3))/dlong +3.
      idel(i,j)=min(il2,max(int(xg(i,j)),2))
   38 xg(i,j)=xg(i,j)-idel(i,j)

c     there exists choice of methods for yg (y departure point)
        do 58 j=3,jl2
        do 58 i=3,il2
c       redefines y shortly, to be j grid number;
c       can then use original jlm bicubic interp (intsb)
c       define approximate grid number (adding 2 for boundary rows)
        yy=2. + .5*(jl+1 + (2*jl+1)*alat(i,j)/pi)
        jdel(i,j)=int(yy)
        yg(i,j)=(alat(i,j)-rlat(jdel(i,j)))/
     &               (rlat(jdel(i,j)+1)-rlat(jdel(i,j)))
   58   continue

c     if(diag)then
c       print*,'id, jd, xg, yg ',id,jd,xg(id,jd),yg(id,jd)
c     endif

      return
      end
