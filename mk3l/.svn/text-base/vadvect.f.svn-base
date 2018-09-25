c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Removed DEFAULT(NONE) clauses from the OpenMP directives, in order to fix
c problems on both the SC and blizzard.
c SJP 2004/03/16
c
c Modified to make use of DYNAMIC scheduling, which leads to a considerable
c improvement in performance
c SJP 2003/03/29
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: vadvect.f,v $
c Revision 1.22  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.21  2001/06/04 02:26:57  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.20  1998/12/10 00:55:45  ldr
c HBG changes to V5-1-21
c
c Revision 1.19  1998/01/30  05:05:37  ldr
c Further parallelization stuff for NEC from MRD.
c
c Revision 1.18  1997/03/06  23:31:48  ldr
c Add j to private directive (for T90).
c
c Revision 1.17  1996/02/19  23:23:07  ldr
c Do vertical SLT properly for 24 level version.
c
c Revision 1.16  1996/02/19  04:10:05  ldr
c Generalize for 24 levels.
c
c Revision 1.15  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.13.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.14  1994/08/08  17:23:22  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.13  93/12/17  15:34:26  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.12  93/08/19  15:11:10  ldr
c Minor cosmetic changes.
c 
c Revision 1.11  93/07/12  10:39:53  ldr
c  Use general Smagorinsky formula for sigma levels.
c 
c Revision 1.10  93/03/22  16:21:28  ldr
c Pass sdot into vadvect as argument.
c 
c Revision 1.9  93/03/15  15:46:05  ldr
c HBG's changes to combine uvpgd and sltjmcg common blocks.
c 
c Revision 1.8  93/02/03  11:15:11  ldr
c Split vinterp and hinterp routines off into separate files.
c 
c Revision 1.7  93/02/03  11:13:28  ldr
c Add zg to list of local variables for do 999 on SGI.
c 
c Revision 1.6  93/01/26  17:25:52  ldr
c Merge of tracer changes to V4-1 with other changes since V4-1.
c 
c Revision 1.5  92/12/09  14:44:47  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
      subroutine vadvect(tdt,sdot,           !inputs
     &                   st)                 !output

c Semi-Lagrangian vertical advection of moisture: find departure point st.
c Gridpoints are numbered 1 to jl and 1 to il in this routine.

c INPUTS:
c sdot   - sdot on full levels
c tdt    - timestep (2*dt for csiro9 leapfrog scheme)
c sig    - sigma levels
c dshdk  - d(sigh)/dk at bottom of layer k

c OUTPUT:
c st     - array of departure points

C Global parameters
      include 'PARAMS.f'
      parameter(il=lon,jl=lat2,kl=nl)

C Argument list
      real tdt
      real sdot(il,kl,jl)
      real st(il,jl,kl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'

C Local work arrays and variables
      real bb(il,kl),sd(il,0:kl+1),sdd(il,0:kl+1)
      real sdh(il,kl+1)

C Local data, functions etc
C     logical diag
C     data diag, id, jd /.false.,25,27 /

C Start code : ----------------------------------------------------------

C     Determine departure point st for all i,j,k:

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (bb, i, j, k, sd, sdd, sddd, sdh)
!$OMP& SHARED  (dshdk, sdot, st, tdt)

      do 34 j=1,jl
c     convert to units of grid-steps/(2*timestep)
c---- use sdot (full level data) to get sdh (half level values)
      do 200 i=1,il
  200 sdh(i,2)=2.0*sdot(i,1,j)
      do 210 k=3,kl
      do 210 i=1,il
  210 sdh(i,k)=2.0*sdot(i,k-1,j)-sdh(i,k-1)

c Divide through by the general dsigh/dk formula (dshdk set up in inital)
      do 21 k=2,kl
      do 21 i=1,il
c 21   sdh(i,k)=sdh(i,k)*tdt/(sig(k)-sig(k-1))
c 21   sdh(i,k)=-sdh(i,k)*tdt*(kl**3)/(6.*(k-1)*(kl+1-k)) !general formula
 21   sdh(i,k)=sdh(i,k)*tdt/dshdk(k)
c     &        /(-1.5*(alf/kl - (3.*alf-2)*(kl-2*k+2.)**2/kl**3))
      do 20 i=1,il
      sd(i,0)=0.
      sdd(i,0)=0.
      sd(i,kl+1)=0.
      sdd(i,kl+1)=0.
      sdh(i,1)=0.
 20   sdh(i,kl+1)=0.

      do 22 k=1,kl
      do 22 i=1,il
c     interpolate sdh to sd with cubic polynomials
c***  assumes equal dsig at present
      sd(i,k)=(sdh(i,k)+sdh(i,k+1))*.5
22    bb(i,k)=sdh(i,k+1)-sdh(i,k)
      do 27 k=2,kl-1
      do 27 i=1,il
27    sd(i,k)=sd(i,k)-(bb(i,k+1)-bb(i,k-1))/16.

c     build in 1/2 & 1/6 into sdd & sddd
c     OK for sd, sdd to go off ends of array as sdh(1)=0., sdh(kl)=0.
      do 28 k=1,kl
      do 28 i=1,il
c28    sdd(i,k)=(sdh(i,k)*(sd(i,k)-sd(i,k-1))  !Old formula
c     .       +sdh(i,k+1)*(sd(i,k+1)-sd(i,k)))/4.
28    sdd(i,k)=(sdh(i,k+1)**2-sdh(i,k)**2)/4.  !New improved (29/9/92)

      do 29 k=1,kl
      do 29 i=1,il
      sddd=(sdh(i,k)*(sdd(i,k)-sdd(i,k-1))
     &      +sdh(i,k+1)*(sdd(i,k+1)-sdd(i,k)))/12.
29    st(i,j,k)=k-sd(i,k)+sdd(i,k)-sddd

C***      if(diag.and.i.eq.id.and.j.eq.jd)then
C***        print *,'sig ',sig
C***        print *,'vadv30  i,j= ',i,j
C***        print *,'sdoth ',(sdh(i,k),k=1,kl+1)
C***        print *,'sd ', sd
C***        print *,'sdd ',sdd
C***        print *,'dep. level ',st
C***      endif

 34   continue 

!$OMP END PARALLEL DO

      return
      end

c----------------------------------------------------------------

      subroutine vadvect3(tdt,sdot,           !inputs
     &                    st,kdel)            !output

c Akima spline method.
c
c Semi-Lagrangian vertical advection of moisture: find departure point st.
c Gridpoints are numbered 1 to jl and 1 to il in this routine.

c INPUTS:
c sdot   - sdot on full levels
c tdt    - timestep (2*dt for csiro9 leapfrog scheme)
c sig    - sigma levels
c dshdk  - d(sigh)/dk at bottom of layer k

c OUTPUT:
c st     - array of departure points

C Global parameters
      include 'PARAMS.f'
      parameter(il=lon,jl=lat2,kl=nl)

C Argument list
      real tdt
      real sdot(il,kl,jl)
      real st(il,jl,kl)
      integer kdel(il,jl,kl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'

C Local work arrays and variables
      real bb(il,kl),sd(il,0:kl+1),sdd(il,0:kl+1)
      real sdh(il,kl+1)

C Local data, functions etc
C     logical diag
C     data diag, id, jd /.false.,25,27 /

C Start code : ----------------------------------------------------------

C     Determine departure point st for all i,j,k:

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (bb, i, j, k, sd, sdd, sddd, sdh)
!$OMP& SHARED  (dshdk, kdel, sdot, st, tdt)

      do 34 j=1,jl
c     convert to units of grid-steps/(2*timestep)
c---- use sdot (full level data) to get sdh (half level values)
      do 200 i=1,il
  200 sdh(i,2)=2.0*sdot(i,1,j)
      do 210 k=3,kl
      do 210 i=1,il
  210 sdh(i,k)=2.0*sdot(i,k-1,j)-sdh(i,k-1)

c Divide through by the general dsigh/dk formula (dshdk set up in inital)
      do 21 k=2,kl
      do 21 i=1,il
c 21   sdh(i,k)=sdh(i,k)*tdt/(sig(k)-sig(k-1))
c 21   sdh(i,k)=-sdh(i,k)*tdt*(kl**3)/(6.*(k-1)*(kl+1-k)) !general formula
 21   sdh(i,k)=sdh(i,k)*tdt/dshdk(k)
c     &        /(-1.5*(alf/kl - (3.*alf-2)*(kl-2*k+2.)**2/kl**3))
      do 20 i=1,il
      sd(i,0)=0.
      sdd(i,0)=0.
      sd(i,kl+1)=0.
      sdd(i,kl+1)=0.
      sdh(i,1)=0.
 20   sdh(i,kl+1)=0.

      do 22 k=1,kl
      do 22 i=1,il
c     interpolate sdh to sd with cubic polynomials
c***  assumes equal dsig at present
      sd(i,k)=(sdh(i,k)+sdh(i,k+1))*.5
22    bb(i,k)=sdh(i,k+1)-sdh(i,k)
      do 27 k=2,kl-1
      do 27 i=1,il
27    sd(i,k)=sd(i,k)-(bb(i,k+1)-bb(i,k-1))/16.

c     build in 1/2 & 1/6 into sdd & sddd
c     OK for sd, sdd to go off ends of array as sdh(1)=0., sdh(kl)=0.
      do 28 k=1,kl
      do 28 i=1,il
c28    sdd(i,k)=(sdh(i,k)*(sd(i,k)-sd(i,k-1))  !Old formula
c     .       +sdh(i,k+1)*(sd(i,k+1)-sd(i,k)))/4.
28    sdd(i,k)=(sdh(i,k+1)**2-sdh(i,k)**2)/4.  !New improved (29/9/92)

      do 29 k=1,kl
      do 29 i=1,il
      sddd=(sdh(i,k)*(sdd(i,k)-sdd(i,k-1))
     &      +sdh(i,k+1)*(sdd(i,k+1)-sdd(i,k)))/12.
29    st(i,j,k)=k-sd(i,k)+sdd(i,k)-sddd

C***      if(diag.and.i.eq.id.and.j.eq.jd)then
C***        print *,'sig ',sig
C***        print *,'vadv30  i,j= ',i,j
C***        print *,'sdoth ',(sdh(i,k),k=1,kl+1)
C***        print *,'sd ', sd
C***        print *,'sdd ',sdd
C***        print *,'dep. level ',st
C***      endif

C     here transform array st and create kdel
      do k=1,kl
      do i=1,il
c              In vinterp3 :
c  Upper boundary condition is qq(kl+1)=-qq(kl) etc.
C   i.e. qq(kl+0.5)~0.0  - so allow st() to go above level k=kl :
         st(i,j,k)=min( max(st(i,j,k),.5) ,(real(kl)+.5))
         kdel(i,j,k)=min(kl,max(int(st(i,j,k)),0))
         st(i,j,k)=st(i,j,k)-kdel(i,j,k)
      enddo
      enddo

 34   continue 

!$OMP END PARALLEL DO

      return
      end

c----------------------------------------------------------------

      subroutine vadvtvd (tdt,sdot,           !inputs
     &                    tr)                 !In & out

c This is a TVD scheme for vertical tracer advection, using the
c "MC" flux limiter of van Leer (1977). See chapter 5 of Durran (1999): 
c Numerical Methods for Wave Equations in Geophysical Fluid Dynamics.
 
c INPUTS:
c sdot   - sdot on full levels
c tdt    - timestep (2*dt for csiro9 leapfrog scheme)
c sig    - sigma levels (common block CNSTA)
c dshdk  - d(sigh)/dk at bottom of layer k (common block CNSTA)
c tr     - tracer array 

c OUTPUT:
c tr     - Modified tracer array

C Global parameters
      include 'PARAMS.f'
      parameter(il=lon,jl=lat2,kl=nl)

C Argument list
      real tdt
      real sdot(il,kl,jl)
      real tr(-1:il+2,-1:jl+2,kl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'

C Local work arrays and variables
      real sdh(il,kl+1)
      real delt(il,jl,0:kl)
      real fluxh(il,0:kl)

C Local data, functions etc
C     logical diag
C     data diag, id, jd /.false.,25,27 /

C Start code : ----------------------------------------------------------

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (fluxh, fluxhi, fluxlo, hdsdot, i, j, k, kp, kx, phi,
!$OMP&          rat, sdh)
!$OMP& SHARED  (delt, dshdk, sdot, tdt, tr)

      do 34 j=1,jl
c     convert to units of grid-steps/(2*timestep)
c---- use sdot (full level data) to get sdh (half level values)
      do 200 i=1,il
  200 sdh(i,2)=2.0*sdot(i,1,j)
      do 210 k=3,kl
      do 210 i=1,il
  210 sdh(i,k)=2.0*sdot(i,k-1,j)-sdh(i,k-1)

c Divide through by the general dsigh/dk formula (dshdk set up in inital)
c This makes sdh +ve upwards.
c According to JMcG, it is also OK to do sdh(i,k)=sdh(i,k)*tdt/(sig(k)-sig(k-1))

      do 21 k=2,kl
      do 21 i=1,il
 21   sdh(i,k)=sdh(i,k)*tdt/dshdk(k)
      do 20 i=1,il
      sdh(i,1)=0.
 20   sdh(i,kl+1)=0.

c Following code is adapted from JMcG's vadvtvd routine.

      do k=1,kl-1
       do i=1,il
         delt(i,j,k)=tr(i,j,k+1)-tr(i,j,k)
       enddo
      enddo
      do i=1,il
        delt(i,j,0)=min(delt(i,j,1),tr(i,j,1)) ! for non-negative tt
c        delt(i,j,0)=0.
        delt(i,j,kl)=0.                        ! safer
        fluxh(i,0)=0.
        fluxh(i,kl)=0.
      enddo

      do k=1,kl-1
        do i=1,il
          kp=sign(1.,sdh(i,k+1))
          kx=k+(1-kp)/2         !  k for sdh +ve,  k+1 for sdh -ve
          rat=delt(i,j,k-kp)/(delt(i,j,k)+sign(1.e-30,delt(i,j,k)))
          phi=max(0.,min(2.*rat,.5+.5*rat,2.))
          fluxlo=tr(i,j,kx)
          fluxhi=.5*(tr(i,j,k)+tr(i,j,k+1))
          fluxh(i,k)=sdh(i,k+1)*(fluxlo+phi*(fluxhi-fluxlo))
        enddo
      enddo
      do k=1,kl
       do i=1,il
         hdsdot=.5*(sdh(i,k+1)-sdh(i,k))
         tr(i,j,k)=(tr(i,j,k)+fluxh(i,k-1)-fluxh(i,k)
     &               +hdsdot*tr(i,j,k))/(1.-hdsdot)
       enddo
      enddo

 34   continue 

!$OMP END PARALLEL DO

      return
      end
