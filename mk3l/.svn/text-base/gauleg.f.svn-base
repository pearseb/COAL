c Removing the declaration of the COMMON block /GRADNS/ to a header file.
c SJP 2009/04/25
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Changes to ensure that integer arguments passed to FFTW routines have the
c correct precision - the modified routine must be compiled with default
c integer precision of 32 bits.
c SJP 2003/03/21
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c Changes to add machine type 'ALPH', which makes use of the super-fast FFTW
c FFT library. These FFTs are installed on the APAC National Facility and may
c be freely downloaded from http://www.fftw.org/ for installation on other
c systems. They are considerably faster than other FFTs - hence the name FFTW,
c which stands for "Fastest Fourier Transform in the West"!
c SJP 2001/11/22
c
c $Log: gauleg.f,v $
c Revision 1.22  1999/05/20 06:23:53  rot032
c HBG changes to V5-2
c
c Revision 1.21  1997/12/17  23:22:52  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.20  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.19  1996/10/24  01:27:40  ldr
c Merge of TIE changes with LDR/MRD changes since V5-0-12.
c
c Revision 1.18  1996/08/12  01:51:22  mrd
c Generalise for triangular truncations other than T63
c
c Revision 1.17.1.1  1996/10/24  01:02:47  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.17  1996/06/13  02:06:34  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.16  1996/04/04  02:16:16  mrd
c Move code to get a cleaner split between spectral initialisation in gauleg
c and other initialisation in inital. gauleg can now be called before inital.
c
c Revision 1.15  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.12.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.14  1995/05/04  04:06:11  ldr
c Put plmf2 into common block with padding to avoid -ve indefinite error
c on Cray with -ei option. Need to add to -Xlocaldata in makefile.
c
c Revision 1.13  1994/08/08  17:21:24  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.12  94/05/31  14:19:35  ldr
c Remove spurious call fftfax on Fujitsu.
c 
c Revision 1.11  93/12/17  15:32:39  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.10  93/10/15  14:17:05  ldr
c Changes to create special 'FUJI' version of model with fast FFTs and
c Legendre transform.
c 
c Revision 1.9  93/10/07  12:09:48  ldr
c Move machine parameter statement into new include file MACHINE.f.
c 
c Revision 1.8  93/10/05  13:06:18  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.7  93/08/19  15:08:09  ldr
c Minor cosmetic changes.
c 
c Revision 1.6  93/07/06  16:30:29  ldr
c      Setting plmn,cplm to zero initially (overshoot on plmn in ptogcray)
c 
c Revision 1.5  92/12/09  14:43:28  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c
c     INPUT/OUTPUT:
c     Output:  from common/gausl in GAUSL.f
c                  wocs - gauss weights/cos(lat)**2
c
c              from common/glegnd in this subroutine
c                  cplmg - global array of gradients of plmg(see below)
c
c              from common/legnd in this subroutine
c                  cplm - gradients of plm, where plm holds the Legendre
c                         polys for (lw,mw)
c                  pad  - pads array on the end of /legnd/ with zeroes
c                  plm  - Legendre polynomials for (lw,mw)
c                  plmx - as plm with extra element
c
c              from common/legnd in this subroutine
c                  padf - pads array with zeroes
c
c              from common/ramp in this subroutine
c                  rampm - ramp in m for all k levels
c
c     In/Out:  from common/cnste in CNSTE.f
c                  epsi - ((l**2-m**2)/((4*l**2)-1))**0.5
c                  flm  - l part of spectral form
c                  flm2 - l(l+1) part of spectral form
c                  flm2i - reciprocal of flm2 (see below)
c
c              from common/fft in this subroutine
c                  ifax, trigs - used for Cray FFTs only
c
c              from common/fftw in this subroutine
c                  plan_forward, plan_backward - used by FFTW FFTs only
c
c              from common/gausl in GAUSL.f
c                  acsq - 1.0/(sia(lg)*sia(lg))
c                  coa  - sin latitude      sia - cosine latitude
c                  w    - normal guass weights
c
c              from common/glegnd in this subroutine
c                  plmg - global array of Legendre polynomials
c
c              from common/gradns in GRADNS.f
c                  glats - physical value of latitudesc     
c 
      subroutine gauleg

!$OMP THREADPRIVATE ( /LEGND/ )
!$OMP THREADPRIVATE ( /LEGNDF/ )

C Global parameters
      include 'PARAMS.f'
c Number of variables for physics grid-spectral ffts
      parameter (ngtos_p=1)
c Total number of physics grid-spectral ffts per lat
      parameter (ntgtos_p=nl*2)
c Number of variables for dynamics grid-spectral ffts
      parameter (ngtos_d=6)
c Total number of dynamics grid-spectral ffts per lat
      parameter (ntgtos_d=6*nl*2)
c Parameters used by FFTW FFTs
      integer fftw_real_to_complex, fftw_complex_to_real
      parameter (fftw_real_to_complex=-1, fftw_complex_to_real=1)
      integer fftw_measure, fftw_threadsafe
      parameter (fftw_measure=1, fftw_threadsafe=128)

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      common/legnd/plmx(lw1,mw),plm(lw,mw),cplm(lw,mw),pad(2,mw)
      common/legndf/plmf2(lw,mw),padf(2,mw)

C Global data blocks
      include 'CNSTE.f'
      include 'GAUSL.f'
      include 'GRADNS.f'
      integer*8 plan_forward, plan_backward
      common /fftw/ plan_forward, plan_backward
      common /fft/ ifax(13),trigs(lon+2,2)  !Cray FFTs only
      common/fftgtos/igtos_p(ntgtos_p),igtos_po(ngtos_p)
     &             ,igtos_d(ntgtos_d),igtos_do(ngtos_d)
      parameter (lwh=(lw+1)/2)
      common/glegnd/plmg(lw1,mw,lat),cplmg(lw,mw,lat)
     &           ,plmgo(lat,lwh,mw),plmge(lat,lwh,mw)
     &           ,cplmgo(lat,lwh,mw),cplmge(lat,lwh,mw)
     &           ,rampm(mw,nl)

C Local work arrays and variables
      real delta(lat)
      double precision dplm(lw1,mw)

C Local data, functions etc
c Number of grid-spectral ffts in physics loop
      integer ngtos_pa(ngtos_p)
      data ngtos_pa/nl/
c Number of grid-spectral ffts in dynamics loop
      integer ngtos_da(ngtos_d)
      data ngtos_da/nl,nl,nl,nl,nl,nl/

C Start code : ----------------------------------------------------------

C**** CREATE GAUSSIAN WEIGHTS ETC & SOME FFT PARAMS 

C**** COSINE(GAUSSIAN CO-LATITUDES) : COA()
C**** GAUSS WEIGHTS : W()
      call gaussv(w,coa,sia,delta,lat)
      do lg=1,lat
         glats(lg)=-acos(sia(lg))
         glats(lat2+1-lg)=-glats(lg)
         acsq(lg)=1.0/(sia(lg)*sia(lg))
         wocs(lg)=w(lg)*acsq(lg)
      end do

C**** PARAMETERS FOR FFT ROUTINES
      write (*, *)
      write (*, *) "Initialising FFTW FFTs..."
      write (*, *)
        call rfftw_f77_create_plan(plan_backward, lon,
     &    fftw_complex_to_real, fftw_measure+fftw_threadsafe)
        call rfftw_f77_create_plan(plan_forward, lon,
     &    fftw_real_to_complex, fftw_measure+fftw_threadsafe)      

c------------------------------------------------------------------
c   Generalized code for stacked or sequential treatment of
c    grid to spectral fft transforms
c    in (phys.f & physseca.f, and dynm.f & dynmseca.f).
c    The number levels of data per group is specified by
c    data statements above (ngtos_pa & ngtos_da). These have to match
c    the stacked variables in the common blocks in dynm.f & phys.f
c------------------------------------------------------------------

c   Physics grid to spectral FFTs transforming 1 variables
c    (nl) for NH,SH at same time : nl*2
c
c   The input data may be stacked (ten etc) in a common block (work1)
c    to process as many "levels" of data at one time (phys.f).
c   The return data (tfnr,tfni etc) is also stacked in a common block
c    (workns). However, to facilitate vectorizing, a change of
c    position of the last two arguments in (tfnr,tfni) relative to (ten)
c    etc is required. The change from [ten(lon,2,nl) to tfnr(mw,nl,2) etc]
c    is achieved through the precomputed pointer "igtos_p".
c   The offset per variable (igtos_po) is also retained for
c    sequential (physseca.f) ffts.

      ih=0
      ix=0
      do nf=1,ngtos_p
        ka=ngtos_pa(nf)
        do k=1,ka
          ih=ih+2
          igtos_p(ih-1)=ix+k
          igtos_p(ih  )=ix+ka+k
        enddo
        igtos_po(nf)=ix
        ix=ix+ka*2
      enddo

c   Dynamics  grid to spectral FFTs transforming 6 variables 
c    (nl,nl,nl,nl,nl,nl) for NH,SH at same time : (6*nl)*2
c
c   The input data is stacked (aon,bon etc) in a common block (giant1)
c    to process as many "levels" of data at one time (dynm.f).
c   The return data (afnr,afni etc) is also stacked in a common block
c    (worknsd). However, to facilitate vectorizing, a change of
c    position of the last two arguments in (afnr,afni) relative to (aon)
c    etc is required. The change from [aon(lon,2,nl) to afnr(mw,nl,2) etc]
c    is achieved through the precomputed pointer "igtos_d".
c   The offset per variable (igtos_do) is also retained for
c   sequential (dynmseca.f) ffts.

      ih=0
      ix=0
      do nf=1,ngtos_d
        ka=ngtos_da(nf)
        do k=1,ka
          ih=ih+2
          igtos_d(ih-1)=ix+k
          igtos_d(ih  )=ix+ka+k
        enddo
        igtos_do(nf)=ix
        ix=ix+ka*2
      enddo


c------------------------------------------------------------------
c   Other spectral stuff
c------------------------------------------------------------------

      do mm=1,mw
         m=mm-1
         sqm=m*m
         epsi(1,mm)=0.0
         do ll=2,lw1
            l=m+ll-1
            sql=l*l
            epsi(ll,mm) = sqrt((sql-sqm) / abs(4*sql-1))
         end do
      end do

      fm = -1.0
      do mm=1,mw
         fm = fm+1.0
         fl = fm-1.0
         do ll=1,lw
            fl = fl+1.0
            flm(ll,mm) = fl
            flm2(ll,mm) = fl*(fl+1.0)
         end do
      end do
c     create 1/flm2 (take care of (1,1) which is 1/0 & not needed)
      flm2(1,1) = 1.0
      do mm=1,mw
         do ll=1,lw
            flm2i(ll,mm) = 1.0/flm2(ll,mm)
         end do
      end do
      flm2(1,1)=0.0

C**** 10: THE FOLLOWING COMPUTES P(L,M) & COS(LAT)D(PLM)/D(LAT)
C****  FOR EACH NH G-LAT
      do 10 lg=1,lat
      call lgndre(dplm,coa(lg),sia(lg),delta(lg),lw-1,lw1)
      if(trunc.eq.'T')then
c.... set top part of rhomboid to zero for triangular truncation
         do mm=2,mw
            do ll=lw1+2-mm,lw1
               dplm(ll,mm)=0.0
            end do
         end do
      endif
      do 22 mm=1,mw 
      do 22 ll=1,lw1
   22 plmg(ll,mm,lg)=dplm(ll,mm)
      do 16 mm=1,mw 
      cplmg(1,mm,lg)=-flm(1,mm)*epsi(2,mm)*plmg(2,mm,lg)
      do 16 ll=2,lw
      cplmg(ll,mm,lg)=(flm(ll,mm)+1.0)*epsi(ll,mm)*plmg(ll-1,mm,lg)
     & -flm(ll,mm)*epsi(ll+1,mm)*plmg(ll+1,mm,lg) 
   16 continue
      if(trunc.eq.'T')then
c.... set top part of rhomboid to zero for triangular truncation
         do mm=2,mw
            do ll=lw+2-mm,lw
               cplmg(ll,mm,lg)=0.0
            end do
         end do
      endif
   10 continue

      do 32 k=1,nl
      do 32 mm=1,mw
   32 rampm(mm,k)=mm-1.0

      do 400 mm=1,mw
      do 410 ll=1,lw
      plm(ll,mm)=0.0
      cplm(ll,mm)=0.0
  410 plmx(ll,mm)=0.0
  400 plmx(lw1,mm)=0.0

c Preset odd and even plm,cplm arrays
      do mm=1,mw
         if ( trunc .eq. 'R' ) then
            ll_lim = lw
         else
            ll_lim = lw1-mm
         end if

         do ll=1,(ll_lim+1)/2
            do lg=1,lat
               plmgo(lg,ll,mm)=plmg(ll*2-1,mm,lg)
               cplmgo(lg,ll,mm)=cplmg(ll*2-1,mm,lg)
            enddo
         enddo
         if(ll_lim.gt.1)then
         do ll=1,ll_lim/2
            do lg=1,lat
               plmge(lg,ll,mm)=plmg(ll*2,mm,lg)
               cplmge(lg,ll,mm)=cplmg(ll*2,mm,lg)
            enddo
         enddo
         endif
      enddo

c Pad array on the end of /legnd/ with zeroes.

      do mm=1,mw
        do ll=1,2
          pad(ll,mm)=0.
          padf(ll,mm)=0.
        enddo
      enddo

      return
      end 
