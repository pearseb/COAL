# 1 "matset.F"
c (1) Fixing syntax errors identified by the g95 Fortran compiler.
c (2) Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Modified to use DGETRF and DGETRI, rather than SGETRF and SGETRI, in order
c to enable the use of 32-bit integers and 64-bit real values.
c SJP 2003/03/22
c
c Changes to add machine type 'ALPH', which uses LAPACK routines SGETRF and
c SGETRI to compute inverses of matrices. These routines are probably also
c available on other machine types.
c SJP 2001/11/22
c
c $Log: matset.f,v $
c Revision 1.22  2001/10/12 02:06:57  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.21  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.20  1999/05/20 06:23:54  rot032
c HBG changes to V5-2
c
c Revision 1.19  1997/10/20  06:23:43  ldr
c Oops - NEC change wasn't actually checked in to V5-1, so here it is!
c
c Revision 1.18  1997/10/06  06:33:45  ldr
c Final corrections for V5-1.
c
c Revision 1.17  1996/10/24  01:03:02  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.17  1996/10/24  01:03:02  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.16  1996/06/13  02:07:08  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.15  1996/03/21  03:18:55  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.14  1996/02/19  04:09:54  ldr
c Generalize for 24 levels.
c
c Revision 1.13  1994/08/08  17:21:42  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.12  93/12/17  15:33:06  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.11  93/10/15  14:17:07  ldr
c Changes to create special 'FUJI' version of model with fast FFTs and
c Legendre transform.
c
c Revision 1.10  93/10/07  12:09:51  ldr
c Move machine parameter statement into new include file MACHINE.f.
c
c Revision 1.9  93/10/05  13:06:39  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c
c Revision 1.8  93/08/19  15:08:36  ldr
c Minor cosmetic changes.
c
c Revision 1.7  93/08/02  15:56:35  ldr
c Speedups for SGI.
c
c Revision 1.6  93/07/12  10:38:00  ldr
c Increase dimension of iagji to allow array bound checks.
c
c Revision 1.5  93/06/16  14:18:05  ldr
c HBG changes for Chen temperature variable.
c
c     INPUT/OUTPUT
c     Input:   from common/cnsta in CNSTA.f
c                  sig  - sigma values
c                  dsk - sigma thicknesses (1 to nl)
c
c              from common/const in CONST.f
c                  cappa - specific gas const dry air/spec heat dry air at
c                          constant pressure
c
c              from common/fewflags in FEWFLAGS.f
c                  chenflag - flag for CHEN version of model
c
c              from common/si in SI.f
c                  am - geopotential height matrix
c
c              from common/timex in TIMEX.f
c                  mins - current model time in mins
c
c              from common/worka in WORKA.f
c                  tmean - preset isothermal temperature
c
c              from arguments
c                  tdt - 2 X timestep (dt)
c
c
c
      subroutine matset(tdt)

      implicit none
C Global parameters
      include 'PARAMS.f'
      integer lmx
      parameter (lmx=lw+mw-2)
      include 'PHYSPARAMS.f'
c Used by DGETRF and DGETRI
      integer lwork
      parameter (lwork=nl)

C Argument list
      real tdt

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'
      include 'CNSTE.f'
      include 'FEWFLAGS.f'
      include 'SI.f'
      include 'TIMEX.f'
      include 'WORKA.f'

      real uk2di(nl,nl),uk2i(nl,nl),jm(nl,nl),gmd(nl,nl),ag(nl,nl)
     & ,amd(nl,nl),ada(nl,nl),iagji(0:lmx,nl,nl)
      common/matzer/uk2di,uk2i,jm,gmd,ag,amd,ada,iagji

C Local work arrays and variables
      real k2m(nl,nl),k2md(nl,nl),iagj(nl,nl)
      REAL WRK1(Nl,Nl),WRK2(Nl,Nl),WRK3(Nl,Nl)

# 141


C Local data, functions etc
      real akd
      real d
      real fllp
      real g1
      real g2
      real g1g2
      real temp
      real t00ec

      integer ierr
      integer j
      integer k
      integer l

C Start code : ----------------------------------------------------------

CHEN
      t00ec=288.0
c---- t00ec replaces tmean(k)
CHEN

C**** ROUTINE TO SET SEMI-IMPLICIT MATRICES
      g2=tdt*0.5
      g1=g2/eradsq
      g1g2=g1*g2
C**** NOTE THE MATRICES NOW HAVE THE COVENTIONAL SUBSCIPT-ARRAY LABELING
C****          !A11 A12 A13 ..! !B1!=!C1!
C****          !A21 A22 A23 ..! !B2! !C2!
C****          !A31 A32 A33 ..! !B3! !C3!
C****
      do 810 k=1,nl
      do 810 j=1,nl
      k2m(j,k)=0.0
      k2md(j,k)=0.0
      uk2i(j,k)=um(j,k)+k2m(j,k)
  810 uk2di(j,k)=um(j,k)+k2md(j,k)
      If(chenflag)Then
CHEN
      do 812 k=1,nl
      do 812 j=1,nl
  812 jm(j,k)=rdry*t00ec*dsk(k)*g1g2
CHEN
      Else
      do 814 k=1,nl
      do 814 j=1,nl
  814 jm(j,k)=rdry*tmean(j)*dsk(k)*g1g2
      End If
C**   IN MATINV, WRK1,WRK2,WRK3 ARE WORK ARRAYS

# 198

        call matinv(uk2i,nl,nl,wrk1,0,d,ierr,wrk2,wrk3)
        call matinv(uk2di,nl,nl,wrk1,0,d,ierr,wrk2,wrk3)


C**   FORM MATRIX AMD=AM*UK2DI
      call mtxmtx(am,uk2di,amd,1.0)
C**   FORM MATRIX ADA=AM+AMD
      do 820 k=1,nl
      do 820 j=1,nl
  820 ada(j,k)=am(j,k)+amd(j,k)
C**** CREATE GMD(,) MATRIX FOR CONSTANT TMEAN
      do 830 k=1,nl
      do 830 j=1,nl
  830 gmd(j,k)=0.0
      do 840 j=1,nl
      temp=cappa*tmean(j)/sig(j)
CHEN
      if(chenflag)temp=cappa*t00ec/sig(j)
CHEN
      do 840 k=j,nl
  840 gmd(j,k)=temp*dsk(k)
      do 850 k=1,nl
  850 gmd(k,k)=0.5*gmd(k,k)
C**   FROM THE REAL MATRICES AMD,GMD COMPUTE THE REAL
C**   MATRIX AG=(AMD*GMD).G1G2
      call mtxmtx(amd,gmd,ag,g1g2)
C**** COMPUTE THE INVERSE MATRIX FOR SOLVING THE SEMI-IMPLICIT
C**** DIVERGENCE EQUATION (L DEPENDENT)
C     AKD=0.0
      akd=5.0e-04
      if(mins.eq.0)then
       write(6,*)'WARNING: Using extra divergence diss. for 1st 3 days'
       write(0,*)'WARNING: Using extra divergence diss. for 1st 3 days'
      endif
      if(mins.ge.1440_8)akd=5.0e-05
      if(mins.ge.2880_8)akd=5.0e-06
      if(mins.ge.4320_8)akd=0.0

C---- SEE NEXT (FOR STARTUP ONLY) (+ C---- ABOVE)
      if((mod(mins,1440_8).eq.0_8).or.(ncepstrt.ge.-1))then

      do 300 l=1,lmx
      fllp=l*(l+1.0)
      do 165 k=1,nl
      do 165 j=1,nl
      iagj(j,k)=um(j,k)*(1.0+tdt*akd)+k2m(j,k)
     & +fllp*(ag(j,k)+jm(j,k))
  165 continue

# 251

        call matinv(iagj,nl,nl,wrk1,0,d,ierr,wrk2,wrk3)


      do 166 k=1,nl
      do 166 j=1,nl
  166 iagji(l,j,k)=iagj(j,k)
  300 continue
      do 167 k=1,nl
      do 167 j=1,nl
 167  iagji(0,j,k)=0.

      endif

      return
      end
