c $Log: gaussv.f,v $
c Revision 1.5  1996/10/24 01:02:47  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.4  1993/08/19  15:08:10  ldr
c Minor cosmetic changes.
c
c Revision 1.3  93/07/27  15:09:35  ldr
c Make scaling of SGI FFTs consistent with Cray FFTs. Now impvor option
c should work on SGI.
c 
c Revision 1.2  92/04/22  10:12:34  ldr
c Generalize for 64/32 bit machines.
c 
c Revision 1.1  91/02/22  16:37:27  ldr
c Initial release V3-0
c
c     INPUT/OUTPUT:
c     Output:  from arguments
c                  fs   - cosine of colatitudes
c                  rads - angle of Gaussian colatitude in radians
c                  sias - cosine latitude, double precision
c                  wts  - cosine of weights, double precision
c
c     In/Out:  from arguments
c                  nzero - no. Gauss latitudes per hemisphere
c
      subroutine gaussv(wts,fs,sias,rads,nzero)
C
C     THIS ROUTINE CALCULATES THE COSINE OF THE COLATITUDES "FS" AND
C     THE WEIGHTS "WTS" , FOR THE GAUSSIAN QUADRATURE WITH "NZERO" POINT
C     BETWEEN POLE AND EQUATOR
C     THIS VERSION WORKS IN DOUBLE PRECISION, BUT COMMUNICATES
C     OUTSIDE IN SINGLE
C
      implicit double precision (a-h,o-z)
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      parameter (isize=56)

C**** NOTE : TOTAL NUMBER OF GAUSS LATS BETWEEN NP AND SP IS
C****  GIVEN BY L > (5*M+1)/2 WHERE M=21 FOR AN R21 MODEL
C****  FOR M=21 , L>53 . FOR HISTORICAL REASONS (BOURKE MODEL)
C****  L WAS SET AT 56. THIS GIVES THE NUMBER OF GAUSS LATITUDES
C****  PER HEMISPHERE 'LAT'=L/2=28 (NZERO).
C****  PARAMETER 'ISIZE' MUST BE >= 'LAT' (NZERO).

      dimension f(isize), wt(isize), sia(isize), rad(isize)
      real fs(*), wts(*), sias(*), rads(*)
C
C----------------------------------------------------------------------
C
C
      if(nzero.gt.isize)then
        print *,' increase isize in gaussg '
        stop
      end if
      xlim=1.e-12
      ir = nzero+nzero
      fi=ir
      fi1=fi+1.0
      piov2 = pi * 0.5
      fn=piov2/float(nzero)
      do 10 i=1,nzero
   10 wt(i)=i-0.5
      do 20 i=1,nzero
   20 f(i)=sin(wt(i)*fn+piov2)
      dn = fi/sqrt(4.0  *fi*fi-1.0  )
      dn1=fi1/sqrt(4.0  *fi1*fi1-1.0  )
      a = dn1*fi
      b = dn*fi1
      irp = ir + 1
      irm = ir -1
      do 2 i=1,nzero
    5 call ordleg(g,f(i),ir)
      call ordleg(gm,f(i),irm)
      call ordleg(gp,f(i),irp)
      gt = (f(i)*f(i)-1.0  ) / (a*gp-b*gm)
      ftemp = f(i) - g*gt
      gtemp = f(i) - ftemp
      f(i) = ftemp
      if(abs(gtemp).gt.xlim) go to 5
    2 continue
      do 6 i=1,nzero
      a=2.0  *(1.0  -f(i)*f(i))
      call ordleg(b,f(i),irm)
      b = b*b*fi*fi
      wt(i)=a*(fi-0.5)/b
      rad(i) =acos(f(i))
    6 sia(i) =sin(rad(i))

c     now convert to single precision
 
      do 30 i=1,nzero
      fs(i) =      f(i)
      wts(i) =      wt(i)
      sias(i) =      sia(i)
      rads(i) =      rad(i)
   30 continue

      return
      end
