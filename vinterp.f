c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Removed DEFAULT(NONE) clauses from the OpenMP directives, in order to fix
c problems on both the SC and blizzard.
c SJP 2004/03/16
c
c Modified to make use of DYNAMIC scheduling, which leads to a considerable
c improvement in performance
c SJP 2003/03/29
c
c Modified so that CBRT function is not used for machine type 'ALPH',
c removing the need to link the model against the SCIPORT library.
c SJP 2003/03/22
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c MACHINE added to list of SHARED variables in parallel section of subroutine
c VINTERP3. The code will not compile on a CRAY otherwise.
c SJP 2001/12/05
c
c Changes to add machine type 'ALPH'.
c SJP 2001/11/22
c
c $Log: vinterp.f,v $
c Revision 1.17  2000/11/14 06:47:45  rot032
c Changes from LDR for aerosol model
c
c Revision 1.16  1999/05/20 06:23:57  rot032
c HBG changes to V5-2
c
c Revision 1.15  1998/12/10  00:56:01  ldr
c HBG changes to V5-1-21
c
c Revision 1.14  1998/01/30  05:05:37  ldr
c Further parallelization stuff for NEC from MRD.
c
c Revision 1.13  1997/03/06  23:31:48  ldr
c Add j to private directive (for T90).
c
c Revision 1.12  1996/10/24  00:31:00  ldr
c Corrections to makefile and vinterp.f for seca.
c
c Revision 1.11  1996/08/08  02:33:40  ldr
c Update qcloud to 24P: No anvils, usual vadv and tightened-up conservation.
c
c Revision 1.10  1996/02/19  23:23:07  ldr
c Do vertical SLT properly for 24 level version.
c
c Revision 1.9  1996/02/19  04:10:07  ldr
c Generalize for 24 levels.
c
c Revision 1.8  1995/07/11  02:52:21  ldr
c Split excessively long doall directive into two lines.
c
c Revision 1.7  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.6  94/03/22  15:37:46  ldr
c Do the B-Staniforth after multiplying through by sig^3.
c 
c Revision 1.5  93/12/17  15:34:27  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.4  93/09/15  17:09:49  ldr
c Changes to get inital version of cloud water scheme going.
c 
c Revision 1.3  93/07/12  11:14:26  ldr
c Use general Smagorinsky formula for sigma levels.
c 
c Revision 1.2  93/03/12  12:02:08  ldr
c Use general N-level formula for sigma values.
c 
c Revision 1.1  93/02/03  11:16:45  ldr
c Initial revision
c 
      subroutine vinterp(st,nfield,
     &                   qq)

c Interpolate to find qq at departure point.
c Part of vertical Semi-Lagrangian advection scheme.

c INPUTS:
c st       - array of departure points
c moisture - Set to 1 if doing water vapour, 2 for tracers, 3 for cloud water.
c qq       - field to be advected

c OUTPUTS:
c qq       - new value of advected field

C Global parameters
      include 'PARAMS.f'
      parameter(il=lon,jl=lat2,kl=nl)

C Argument list
      real st(il,jl,kl)
      integer nfield
      real qq(-1:il+2,-1:jl+2,kl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'

C Local work arrays and variables
      real zg(il,kl)
      integer kdel(il,kl),kkdel(il,kl)
      real temp(il,kl),sigdep(il,kl)

C Local data, functions etc
      data ntest/1/

C Start code : ----------------------------------------------------------

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (c1, c2, c3, c4, i, j, k, kdel, kkdel, sigdep, temp, zg)
!$OMP& SHARED  (nfield, ntest, qq, sig, st)

      do 999 j=1,jl
      do 38 k=1,kl
      do 38 i=1,il
c  Convert from k space to sigma
c  zg is fractional st displacement from point kkdel

      sigdep(i,k) = ((kl-st(i,j,k))**2*(kl+2*st(i,j,k))
     &            + (kl+1-st(i,j,k))**2*(kl-2+2*st(i,j,k)))/(2*kl**3)

      zg(i,k)=min( max(st(i,j,k),1.) ,float(kl))
      kdel(i,k)=min(kl-2,max(ifix(zg(i,k)),2))
      kkdel(i,k)=min(kl-1,max(ifix(zg(i,k)),1))  !used by B&S
38    zg(i,k)=zg(i,k)-kdel(i,k)   ! usual code
c38    zg(i,k)=zg(i,k)-kkdel(i,k)    ! this one for W&R tests

c     Convert qq to qq/(sigma^3) (moisture only)

      if(nfield.eq.1)then
        do 39 k=1,kl
        do 39 i=1,il
 39     qq(i,j,k)=qq(i,j,k)/sig(k)**3
      endif

      do 4 k=1,kl
      do 4 i=1,il
      c1=qq(i,j,kdel(i,k)-1)
      c2=qq(i,j,kdel(i,k))
      c3=qq(i,j,kdel(i,k)+1)
      c4=qq(i,j,kdel(i,k)+2)
 4    temp(i,k)=((1.-zg(i,k))*((2.-zg(i,k))*((1.+zg(i,k))*c2
     &  -zg(i,k)*c1/3.)-zg(i,k)*(1.+zg(i,k))*c4/3.)
     &     +zg(i,k)*(1.+zg(i,k))*(2.-zg(i,k))*c3)*.5

      if(nfield.eq.3)then !Set qc=0 if advecting upwards at k=1
        do 45 i=1,il
c          if(st(i,j,1).lt.1)temp(i,1)=(st(i,j,1)-.5)*qq(i,j,1)
          if(st(i,j,1).lt.1)temp(i,1)=0.
 45     continue
      endif

c Bermejo-Staniforth positive definite scheme: 
c Force qq to lie between c2 and c3. No need to do this if departure point st
c lies below k=1 or above k=nl, as in this case the field value is taken from
c the level itself, or handled specially in the case of qcloud below k=1.

      if(nfield.eq.1.and.ntest.eq.1)then  ! convert back to full qg first
        do k=1,kl
         do i=1,il
          temp(i,k)= temp(i,k)*sigdep(i,k)**3
          qq(i,j,k)= qq(i,j,k)*sig(k)**3
         enddo
        enddo
      endif

      do 5 k=1,kl
      do 5 i=1,il
      if(st(i,j,k).gt.1.and.st(i,j,k).lt.nl)then
        c2=qq(i,j,kkdel(i,k))
        c3=qq(i,j,kkdel(i,k)+1)
        temp(i,k)= min ( max (temp(i,k),min(c2,c3)) , max(c2,c3) )
      endif
 5    continue

      if(nfield.eq.1.and.ntest.eq.0)then  ! original code used ntest=0
        do k=1,kl
         do i=1,il
          qq(i,j,k)= temp(i,k)*sigdep(i,k)**3
         enddo
        enddo
      else
        do 52 k=1,kl
        do 52 i=1,il
            qq(i,j,k)=temp(i,k)
 52     continue
      endif

 999  continue 

!$OMP END PARALLEL DO

      return
      end

c-----------------------------------------------------------------------

      subroutine vinterp3(st,kdel,nfield,qq)

c Akima spline method
c
c Interpolate to find qq at departure point.
c Part of vertical Semi-Lagrangian advection scheme.

c INPUTS:
c st       - array of departure points
c kdel     - array of departure points in k space
c nfield   - Set to 1 if doing water vapour, 2 for tracers, 3 for cloud water.
c qq       - field to be advected

c OUTPUTS:
c qq       - new value of advected field

C Global parameters
      include 'PARAMS.f'
      parameter(il=lon,jl=lat2,kl=nl)

C Argument list
      real st(il,jl,kl)
      integer kdel(il,jl,kl)
      integer nfield
      real qq(-1:il+2,-1:jl+2,kl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'

C Local work arrays and variables
      real temp(il,kl),rhs(il,kl)
      real der(il,0:kl+1),x(il,-1:kl+1)

C Local data, functions etc
      logical B_S
      parameter (B_S=.true.) ! for Bermejo-Staniforth

C Start code : ----------------------------------------------------------

C Akima spline method for moisture variables only

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (aa, alfa, bb, beta, c2, c3, der, i, j, k, rhs, temp, x)
!$OMP& SHARED  (kdel, nfield, qq, st)

      do 999 j=1,jl

c     Convert qq to qq**(1/3) (moisture only)

      if(nfield.eq.1)then

CSJP  Former machine dependence at this point
c           following preserves sign of qq after cube root
          do k=1,kl
          do i=1,il
            qq(i,j,k)=sign( abs(qq(i,j,k))**(1./3.) , qq(i,j,k) )
          enddo
          enddo

      endif

      do k=1,kl-1
       do i=1,il
        x(i,k)=qq(i,j,k+1)-qq(i,j,k)  ! 1-sided derivative here
       enddo
      enddo

      do i=1,il
c     put top & bottom values into rhs(1) and rhs(kl)
        rhs(i,1)=qq(i,j,1)
        rhs(i,kl)=0.0
        x(i,-1)=0.      ! 1-sided here
        x(i,0)=0.       ! 1-sided here
        x(i,kl+1)=0.    ! 1-sided here
        x(i,kl)=rhs(i,kl)-qq(i,j,kl)   ! 1-sided here
      enddo

      do k=1,kl
       do i=1,il
        alfa=1.e-25+abs(x(i,k+1)-x(i,k))
        beta=1.e-25+abs(x(i,k-1)-x(i,k-2))
        der(i,k)=(alfa*x(i,k-1)+beta*x(i,k))/(alfa+beta)
       enddo
      enddo

      do i=1,il
        der(i,0)=0.
        der(i,kl+1)=0.
        x(i,0)=rhs(i,1)
        x(i,kl+1)=rhs(i,kl)
      enddo

c     to help vectorize move qq array into x array
      do k=1,kl
       do i=1,il
        x(i,k)=qq(i,j,k)
       enddo
      enddo

!     now interpolate at the departure points
      do k=1,kl
       do i=1,il
        c2=x(i,kdel(i,j,k))
        c3=x(i,kdel(i,j,k)+1)
        bb=st(i,j,k)
        aa=1.-bb
        temp(i,k)=aa*aa*(aa*c2 +bb*(3.*c2+der(i,kdel(i,j,k))))
     &              +bb*bb*(bb*c3 +aa*(3.*c3-der(i,kdel(i,j,k)+1)))
        if(B_S)then    ! B_S=.true. for Bermejo-Staniforth
          temp(i,k)=min( max(min(c2,c3),temp(i,k)) , max(c2,c3) )
        endif
       enddo
      enddo

      if(nfield.eq.1)then  ! convert back to full qq first
        do k=1,kl
          do i=1,il
            qq(i,j,k)=temp(i,k)**3
           enddo
        enddo
      else
        do k=1,kl
          do i=1,il
            qq(i,j,k)=temp(i,k)
          enddo
        enddo
      endif

 999  continue 

!$OMP END PARALLEL DO

      return
      end
