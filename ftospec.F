c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Fixed so that SGEMM is still called for machine type 'CRAY'.
c SJP 2003/06/20
c
c Modified to use DGEMM, rather than SGEMM, in order to enable the use of
c 32-bit integers and 64-bit real values.
c SJP 2003/03/22
c
c Changes to add machine type 'ALPH'.
c SJP 2001/11/22
c
c $Log: ftospec.f,v $
c Revision 1.5  1999/05/20 06:23:51  rot032
c HBG changes to V5-2
c
c Revision 1.4  1998/12/10  00:55:37  ldr
c HBG changes to V5-1-21
c
c Revision 1.3  1997/12/19  02:03:14  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.2  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.1  1996/10/24  01:33:41  ldr
c Initial revision
c

      subroutine ftospec(llmax,mm,plmo,plme,fro,foi,fer,fei
     &                   ,spectr,specti)

c     (Written by Tracey Elliott 9/96)
c
c     This subroutine carries out the Legendre transform to re-create 
c     spectral variables.  It should be used where there is only one
c     matrix multiply for each variable.
c
c     The subroutine takes the two input vectors for each variable,eg fro/foi,
c     merges them into a matrix, fm,  and carries out the matrix multiply.  
c     This is done once for the odd values of ll, and again for the even 
c     values.
c
c     The equations are of the form: C=A*B where
c            C=spectral variable, real or imaginary(spectr or specti)
c            A=array of Legendre polynomials (plmo or plme)
c            B=array of Fourier values, real AND imaginary(fm)
c
c     nstack - currently 2*nl
c     Input: plmo,plme - matrix of Legendre polynomials
c            fro - vector of Fourier values (real,odd)  )to be
c            foi - vector of Fourier values (imag,odd)  )multiplied by
c            fer - vector of Fourier values (real,even) )the matrix
c            fei - vector of Fourier values (imag,even) )poly
c
c     Output:spectr - spectral variable, real part
c            specti - spectral variable, imaginary part
c
c(hbg) Reordered matrix elements for vectorization purposes.
c(hbg) Replaced calls to matrix multiply "sgemm" with calls to "mmtx"
c(hbg) which will generate a call to an inhouse matrix
c(hbg) multiply routine (on NEC).

C Global parameters
      include 'PARAMS.f'
      parameter (nstack=2*nl)

C Argument list
      integer llmax
      integer mm
      real plmo(lat,*)
      real plme(lat,*)
c-- Dimension following at mw+1 to prevent Bank Conflicts at T63
      real fro(mw+1,nl,lat)
      real foi(mw+1,nl,lat)
      real fer(mw+1,nl,lat)
      real fei(mw+1,nl,lat)
      real spectr(lw,mw,nl)
      real specti(lw,mw,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      real prodo(nstack,(llmax+1)/2),prode(nstack,llmax/2)
      real fm(nstack,lat)

C Local data, functions etc

C Start code : ----------------------------------------------------------

c
c The input Legendre matrix in odd parts : plmo(lat,*)
c

c
c     Merge the two odd input vectors, real and imaginary, into a matrix
c

      do lg=1,lat
         do k=1,nl
            fm(k,lg)=fro(mm,k,lg)
            fm(k+nl,lg)=foi(mm,k,lg)
         enddo
      enddo

c
c     Carry out the matrix multiply, for the odd matrices
c           

#ifdef ALPH
        call dgemm('n','n',nstack,(llmax+1)/2,lat,1.0,fm,nstack
     &             ,plmo,lat,0.,prodo,nstack)
#else
        call mmtx(fm,plmo,prodo,nstack,lat,(llmax+1)/2)
#endif

c
c     Split the resultant matrix into a real output vector and an 
c     imaginary output vector 
c
      do k=1,nl
         do ll=1,(llmax+1)/2
            spectr(ll*2-1,mm,k)=prodo(k,ll)
            specti(ll*2-1,mm,k)=prodo(k+nl,ll)
         enddo
      enddo

c     With triangular truncation llmax is variable. If llmax=1 there
c     are no even components so just return.
      if ( llmax .eq. 1 ) return

c
c The input Legendre matrix in even parts : plme(lat,*)
c

c
c     Merge the two even input vectors,  real and imaginary, 
c     into a matrix
c

      do lg=1,lat
         do k=1,nl
            fm(k,lg)=fer(mm,k,lg)
            fm(k+nl,lg)=fei(mm,k,lg)
         enddo
      enddo

c
c     Carry out the matrix multiply, for the even matrices
c           

#ifdef ALPH
        call dgemm('n','n',nstack,llmax/2,lat,1.0,fm,nstack
     &             ,plme,lat,0.,prode,nstack)
#else
        call mmtx(fm,plme,prode,nstack,lat,llmax/2)
#endif

c
c     Split the resultant matrix into a real output vector and an 
c     imaginary output vector 
c
      do k=1,nl
         do ll=1,llmax/2
            spectr(ll*2,mm,k)=prode(k,ll)
            specti(ll*2,mm,k)=prode(k+nl,ll)
         enddo
      enddo

      return
      end
