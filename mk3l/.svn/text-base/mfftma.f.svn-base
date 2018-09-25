c Changes to ensure that integer arguments passed to FFTW routines have the
c correct precision - the modified routine must be compiled with default
c integer precision of 32 bits.
c SJP 2003/03/21
c
c New subroutine that carries out real-to-complex Fourier Transforms for
c machine type 'ALPH', using the super-fast FFTW FFT library.
c SJP 2001/11/22
c
c     INPUT/OUTPUT
c     Input:   from arguments
c                  fg - grid point values of fourier coefficients
c
c     Output:  from arguments
c                  fmr, fmi - fourier coefficients, real and imaginary
c
c     In/Out:  from arguments
c                  nex - total no of vectors
c 
      subroutine mfftma(fg, fmr, fmi, nex, ipd)

C Global parameters
      include 'PARAMS.f'
      parameter (nlev=6*nl*2)  ! max no. of stacked ffts
      parameter (mmax=mw-1)
c Number of variables for physics grid-spectral ffts
      parameter (ngtos_p=1)
c Total number of physics grid-spectral ffts per lat
      parameter (ntgtos_p=nl*2)
c Number of variables for dynamics grid-spectral ffts
      parameter (ngtos_d=6)
c Total number of dynamics grid-spectral ffts per lat
      parameter (ntgtos_d=6*nl*2)
c Parameters used by FFTW FFTs
      integer istride, idist, ostride, odist
      parameter (istride=1, idist=lon, ostride=1, odist=lon)

C Argument list
      real fg(lon,nlev)
      real fmr(0:mmax,nlev)
      real fmi(0:mmax,nlev)
      integer nex       ! number of stacked ffts
      integer ipd       ! =1 (physics), =2 (dynamics)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      integer*8 plan_forward, plan_backward
      common /fftw/ plan_forward, plan_backward
      common/fftgtos/igtos_p(ntgtos_p),igtos_po(ngtos_p)
     &              ,igtos_d(ntgtos_d),igtos_do(ngtos_d)

C Local work arrays and variables
      integer howmany
      real in(lon, nlev), out(lon, nlev)

C Local data, functions etc

C Start code : ----------------------------------------------------------

C**** MODIFIED FOR TRANSFORMING A MAXIMUM OF (6*nl)*2 VARIABLES AT ONE HIT.
C**** (GRID TO "M" SPACE - SPECTRAL)
C**** THIS IS SPECIFICALLY FOR THE MK3 CSIRO MODEL.
C**** DATA PROCESSED AS PAIRS OF LATS (NH, THEN SH)

C STORE DATA POINTS FG IN ARRAY IN, SCALING BY A FACTOR OF LON AS THE FFTW
C FFTs ARE UN-NORMALISED

      do 20 ih = 1, nex
        do 10 j = 1, lon
          in(j, ih) = fg(j, ih) / real(lon)
 10     continue
 20   continue

C CALL RFFTW_F77 THEN CONVERT OUTPUT ARRAY OUT TO REAL ARRAYS FMR AND FMI

      howmany = nex
      call rfftw_f77(plan_forward, howmany, in, istride, idist, out,
     &               ostride, odist) 

      if (ipd .eq. 1) then

c  Physics loop transform
        do ih = 1, nex
          fmr(0, igtos_p(ih)) = out(1, ih)
          fmi(0, igtos_p(ih)) = 0.0
          do m = 1, mmax
            fmr(m, igtos_p(ih)) = out(m+1, ih)
            fmi(m, igtos_p(ih)) = out(lon+1-m, ih)
          enddo
        enddo

      else      

c  Dynamics loop transform
        do ih = 1, nex
          fmr(0, igtos_d(ih)) = out(1, ih)
          fmi(0, igtos_d(ih)) = 0.0
          do m = 1, mmax
            fmr(m, igtos_d(ih)) = out(m+1, ih)
            fmi(m, igtos_d(ih)) = out(lon+1-m, ih)
          enddo
        enddo

      endif

      return
      end
