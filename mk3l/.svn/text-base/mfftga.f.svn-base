c Changes to ensure that integer arguments passed to FFTW routines have the
c correct precision - the modified routine must be compiled with default
c integer precision of 32 bits.
c SJP 2003/03/21
c
c New subroutine that carries out complex-to-real Fourier Transforms for
c machine type 'ALPH', using the super-fast FFTW FFT library.
c SJP 2001/11/22
c
c     INPUT/OUTPUT
c     Input:   from common/fftw in this subroutine
c                  plan_backward - used by FFTW FFTs
c
c     Output:  from arguments
c                  fg - grid point values of fourier coefficients
c
c     In/Out:  from arguments
c                  fm - fourier coefficients
c                  nex - total no. of vectors
c
      subroutine mfftga(fg, fm, nex)

C Global parameters
      include 'PARAMS.f'
      parameter (mmax=mw-1,nlev=(5*nl+1)*2)
c Parameters used by FFTW FFTs
      integer istride, idist, ostride, odist
      parameter (istride=1, idist=lon, ostride=1, odist=lon)

C Argument list
      real fg(lon,nlev)
      complex fm(0:mmax,nlev)
      integer nex

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      integer*8 plan_forward, plan_backward
      common /fftw/ plan_forward, plan_backward

C Local work arrays and variables
      integer howmany
      real in(lon, nlev), out(lon, nlev)

C Local data, functions etc

C Start code : ----------------------------------------------------------

C**** MODIFIED FOR TRANSFORMING A MAXIMUM OF (4*nl+2)*2 VARIABLES AT ONE HIT.
C**** (SPECTRAL "M" SPACE TO GRID)
C**** THIS IS SPECIFICALLY FOR THE Mk3 CSIRO MODEL.

C STORE COMPLEX ARRAY OF FOURIER COEFFS FM IN REAL ARRAY IN FOR INPUT TO
C RFFTW_F77

      do ih = 1, nex
        do mm = 1, lon
          in(mm, ih) = 0.0
        enddo
      enddo
      do 20 ih = 1, nex
         in(1, ih) = real(fm(0, ih))
         do m = 1, mmax
            in(m+1, ih) = real(fm(m, ih))
         enddo
         do m = 1, mmax
            in(lon+1-m, ih) = aimag(fm(m, ih))
         enddo
 20   continue

C CALL RFFTW_F77 THEN STORE OUTPUT ARRAY OUT IN FG

      howmany = nex
      call rfftw_f77(plan_backward, howmany, in, istride, idist, out,
     &               ostride, odist) 

      do 30 ih = 1, nex
        do 25 j = 1, lon
          fg(j, ih) = out(j, ih)
 25     continue
 30   continue

      return
      end
