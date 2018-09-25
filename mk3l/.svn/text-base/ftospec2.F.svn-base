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
c $Log: ftospec2.f,v $
c Revision 1.5  1999/05/20 06:23:56  rot032
c HBG changes to V5-2
c
c Revision 1.4  1998/12/10  00:55:53  ldr
c HBG changes to V5-1-21
c
c Revision 1.3  1997/12/19  02:03:18  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.2  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.1  1996/10/24  01:33:41  ldr
c Initial revision
c

      subroutine ftospec2(llmax,mm,ivi,ivr,plmo,plme,cplmo,cplme
     &           ,fro,cfro,foi,cfoi,fer,cfer,fei,cfei,spectr,specti)

c     (Written by Tracey Elliott 9/96)
c
c     This subroutine carries out the Legendre transform to re-create 
c     spectral variables.  It allows for the case where the re-creation of
c     the spectral variables requires 2 matrix multiplies for each variable.
c
c     The subroutine takes the two sets of input vectors for each variable, 
c     eg fro/foi and cfro/cfoi, merges them into a 2 matrices, fm and cfm,   
c     and carries out the matrix multiply.  This is done twice per variable, 
c     once for odd values of ll, and again for the even values
c
c     The equations are all of the form: C=D*E + iv* A*B 
c             C=spectral variable, real or imaginary(spectr or specti)
c             A=array of Legendre polynomials (plmo or plme)
c             B=array of Fourier values, real AND imaginary(fm)
c             D=array of gradients of Legendre polynomials (cplmo or cplme)
c             E=array of Fourier values, real AND imaginary(cfm)
c             iv=ivi or ivr.  This takes the value 1 when the eqn is of the
c             form C=D*E + A*B or -1 when the eqn is of the form C=D*E - A*B
c
c     nstack - currently 2*nl
c     Input:  ivi - integer value for the multiplication of the imaginary
c             part of fm, either 1 or -1.
c             ivr - same as above for the real part of fm
c             plmo,plme - matrix of Legendre polynomials
c             cplmo,cplme - matrix of gradients of Legendre polynomials
c             fro - vector of Fourier values (real,odd)  )to be
c             foi - vector of Fourier values (imag,odd)  )multiplied by
c             fer - vector of Fourier values (real,even) )the matrix
c             fei - vector of Fourier values (imag,even) )poly
c             cfro,cfoi,cfer,cfei - as above but multiplied by cpoly
c            
c     Output: spectr - spectral variable, real part
c             specti - spectral variable, imaginary part
c
c(hbg) Reordered matrix elements for vectorization purposes.
c(hbg) Replaced calls to matrix multiply "sgemm" with calls to "mmtx"
c(hbg) which will generate a call to an inhouse matrix
c(hbg) multiply routine (on NEC).
c!@#$

C Global parameters
      include 'PARAMS.f'
      parameter (nstack=2*nl)

C Argument list
      integer llmax
      integer mm
      integer ivi
      integer ivr
      real plmo(lat,*)
      real plme(lat,*)
      real cplmo(lat,*)
      real cplme(lat,*)
c-- Dimension following at mw+1 to prevent Bank Conflicts at T63
      real fro(mw+1,nl,lat)
      real cfro(mw+1,nl,lat)
      real foi(mw+1,nl,lat)
      real cfoi(mw+1,nl,lat)
      real fer(mw+1,nl,lat)
      real cfer(mw+1,nl,lat)
      real fei(mw+1,nl,lat)
      real cfei(mw+1,nl,lat)
      real spectr(lw,mw,nl)
      real specti(lw,mw,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      real prodo(nstack,(llmax+1)/2),prode(nstack,llmax/2)
      real cprodo(nstack,(llmax+1)/2),cprode(nstack,llmax/2)
      real fm(nstack,lat),cfm(nstack,lat)

C Local data, functions etc

C Start code : ----------------------------------------------------------

c
c     Merge the two odd input vectors, real and imaginary, into a matrix
c

      do lg=1,lat
         do k=1,nl
            fm(k,lg)=ivr*fro(mm,k,lg)
            fm(k+nl,lg)=ivi*foi(mm,k,lg)
            cfm(k,lg)=cfro(mm,k,lg)
            cfm(k+nl,lg)=cfoi(mm,k,lg)
         enddo
      enddo
c
c     Carry out the matrix multiplies , for the odd matrices with
c the input Legendre matrix in odd parts : plmo(lat,*) & cplmo(lat,*)
c
            
#ifdef ALPH
        call dgemm('n','n',nstack,(llmax+1)/2,lat,1.0,fm,nstack
     &             ,plmo,lat,0.,prodo,nstack)
        call dgemm('n','n',nstack,(llmax+1)/2,lat,1.0,cfm,nstack
     &             ,cplmo,lat,0.,cprodo,nstack)
#else
        call mmtx( fm, plmo, prodo,nstack,lat,(llmax+1)/2)
        call mmtx(cfm,cplmo,cprodo,nstack,lat,(llmax+1)/2)
#endif

c
c     Add the resultant matrices and then split into a real output 
c     vector and an imaginary output vector.
c

      do k=1,nl
         do ll=1,(llmax+1)/2
            spectr(ll*2-1,mm,k)=prodo(k,ll)+cprodo(k,ll)
            specti(ll*2-1,mm,k)=prodo(k+nl,ll)+cprodo(k+nl,ll)
         enddo
      enddo

c     With triangular truncation llmax is variable. If llmax=1 there
c     are no even components so just return.
      if ( llmax .eq. 1 ) return

c
c     Merge the two even input vectors,  real and imaginary, 
c     into a matrix
c

      do lg=1,lat
         do k=1,nl
            fm(k,lg)=ivr*fer(mm,k,lg)
            fm(k+nl,lg)=ivi*fei(mm,k,lg)
            cfm(k,lg)=cfer(mm,k,lg)
            cfm(k+nl,lg)=cfei(mm,k,lg)
         enddo
      enddo
c
c     Carry out the matrix multiply, for the even matrices with
c the input Legendre matrix in even parts : plme(lat,*) & cplme(lat,*)
c

#ifdef ALPH
        call dgemm('n','n',nstack,llmax/2,lat,1.0,fm,nstack
     &             ,plme,lat,0.,prode,nstack)
        call dgemm('n','n',nstack,llmax/2,lat,1.0,cfm,nstack
     &             ,cplme,lat,0.,cprode,nstack)
#else
        call mmtx( fm, plme, prode,nstack,lat,llmax/2)
        call mmtx(cfm,cplme,cprode,nstack,lat,llmax/2)
#endif

c
c     Add the resultant matrices and then split into a real output 
c     vector and an imaginary output vector 
c
      do k=1,nl
         do ll=1,llmax/2
            spectr(ll*2,mm,k)=prode(k,ll)+cprode(k,ll)
            specti(ll*2,mm,k)=prode(k+nl,ll)+cprode(k+nl,ll)
         enddo
      enddo

      return
      end
CXXXXXX
      subroutine ftospec3(llmax,mm,plmo,plme,cplmo,cplme
     &           ,ffor,gfor,hfoi,gfoi,ffer,gfer,hfei,gfei,itr,iti
     &           ,paor,pbor,paoi,pboi,paer,pber,paei,pbei,ipr,ipi
     &           ,xbor,xaor,xboi,xaoi,xber,xaer,xbei,xaei,ixr,ixi)

c     (Written by Tracey Elliott 9/96, & HBG 98/99)
c
c     This subroutine carries out the Legendre transform to re-create 
c     spectral variables.  It allows for the case where the re-creation of
c     the spectral variables requires 2 matrix multiplies for each variable.
c
c     The subroutine takes the two sets of input vectors for each variable, 
c     eg ffor/hfoi and gfor/gfoi, merges them into a 2 matrices, fm and cfm,   
c     and carries out the matrix multiply.  This is done twice per variable, 
c     once for odd values of ll, and again for the even values
c
c     The equations are all of the form: C=D*E + iv* A*B 
c             C=spectral variable, real or imaginary(itr or iti)
c             A=array of Legendre polynomials (plmo or plme)
c             B=array of Fourier values, real AND imaginary(fm)
c             D=array of gradients of Legendre polynomials (cplmo or cplme)
c             E=array of Fourier values, real AND imaginary(cfm)
c             iv=ivi or ivr.  This takes the value 1 when the eqn is of the
c             form C=D*E + A*B or -1 when the eqn is of the form C=D*E - A*B
c
c     nstack - currently 2*nl*3
c     Input:  plmo,plme - matrix of Legendre polynomials
c             cplmo,cplme - matrix of gradients of Legendre polynomials
c             ffor - vector of Fourier values (real,odd)  )to be
c             hfoi - vector of Fourier values (imag,odd)  )multiplied by
c             ffer - vector of Fourier values (real,even) )the matrix
c             hfei - vector of Fourier values (imag,even) )poly
c             gfor,gfoi,gfer,gfei - as above but multiplied by cpoly
c            
c     Output: itr - spectral temperature variable, real part
c             iti - spectral temperature variable, imaginary part
c
c     Similarly for other two variables (ipr,ipi) and (ixr,ixi).
c
c(hbg) Reordered matrix elements for vectorization purposes.
c(hbg) Replaced calls to matrix multiply "sgemm" with calls to "mmtx"
c(hbg) which will generate a call to an inhouse matrix
c(hbg) multiply routine (on NEC).
c(hbg) Reworked to do 3 calls of ftospec2 at same time (ftospec3)

C Global parameters
      include 'PARAMS.f'
      parameter (nstack=2*nl*3)

C Argument list
      integer llmax
      integer mm
      real plmo(lat,*)
      real plme(lat,*)
      real cplmo(lat,*)
      real cplme(lat,*)
c-- Dimension following at mw+1 to prevent Bank Conflicts at T63
      real ffor(mw+1,nl,lat)
      real gfor(mw+1,nl,lat)
      real hfoi(mw+1,nl,lat)
      real gfoi(mw+1,nl,lat)
      real ffer(mw+1,nl,lat)
      real gfer(mw+1,nl,lat)
      real hfei(mw+1,nl,lat)
      real gfei(mw+1,nl,lat)
      real itr(lw,mw,nl)
      real iti(lw,mw,nl)
      real paor(mw+1,nl,lat)
      real pbor(mw+1,nl,lat)
      real paoi(mw+1,nl,lat)
      real pboi(mw+1,nl,lat)
      real paer(mw+1,nl,lat)
      real pber(mw+1,nl,lat)
      real paei(mw+1,nl,lat)
      real pbei(mw+1,nl,lat)
      real ipr(lw,mw,nl)
      real ipi(lw,mw,nl)
      real xbor(mw+1,nl,lat)
      real xaor(mw+1,nl,lat)
      real xboi(mw+1,nl,lat)
      real xaoi(mw+1,nl,lat)
      real xber(mw+1,nl,lat)
      real xaer(mw+1,nl,lat)
      real xbei(mw+1,nl,lat)
      real xaei(mw+1,nl,lat)
      real ixr(lw,mw,nl)
      real ixi(lw,mw,nl)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      real prodo(nstack,(llmax+1)/2),prode(nstack,llmax/2)
      real cprodo(nstack,(llmax+1)/2),cprode(nstack,llmax/2)
      real fm(nstack,lat),cfm(nstack,lat)

C Local data, functions etc
c     integer ivi(3)
c     data ivi/1,-1,1/
c     integer ivr(3)
c     data ivr/1,1,-1/

C Start code : ----------------------------------------------------------

c
c Merge the pairs of odd input vectors, real and imaginary, into a matrix
c

      do lg=1,lat
         do k=1,nl
            fm (k,lg)= ffor(mm,k,lg)
            cfm(k,lg)= gfor(mm,k,lg)
            m=k+nl
            fm (m,lg)= hfoi(mm,k,lg)
            cfm(m,lg)= gfoi(mm,k,lg)
            m=k+2*nl
            fm (m,lg)= paor(mm,k,lg)
            cfm(m,lg)= pbor(mm,k,lg)
            m=k+3*nl
            fm (m,lg)=-paoi(mm,k,lg)
            cfm(m,lg)= pboi(mm,k,lg)
            m=k+4*nl
            fm (m,lg)=-xbor(mm,k,lg)
            cfm(m,lg)= xaor(mm,k,lg)
            m=k+5*nl
            fm (m,lg)= xboi(mm,k,lg)
            cfm(m,lg)= xaoi(mm,k,lg)
         enddo
      enddo
c
c     Carry out the matrix multiplies , for the odd matrices with
c the input Legendre matrix in odd parts : plmo(lat,*) & cplmo(lat,*)
c
            
#ifdef ALPH
        call dgemm('n','n',nstack,(llmax+1)/2,lat,1.0,fm,nstack
     &             ,plmo,lat,0.,prodo,nstack)
        call dgemm('n','n',nstack,(llmax+1)/2,lat,1.0,cfm,nstack
     &             ,cplmo,lat,0.,cprodo,nstack)
#else
        call mmtx( fm, plmo, prodo,nstack,lat,(llmax+1)/2)
        call mmtx(cfm,cplmo,cprodo,nstack,lat,(llmax+1)/2)
#endif

c
c     Add the resultant matrices and then split into a real output 
c     vector and an imaginary output vector.
c

      do k=1,nl
         do ll=1,(llmax+1)/2
            itr(ll*2-1,mm,k)=prodo(k,ll)+cprodo(k,ll)
            m=k+nl
            iti(ll*2-1,mm,k)=prodo(m,ll)+cprodo(m,ll)
            m=k+2*nl
            ipr(ll*2-1,mm,k)=prodo(m,ll)+cprodo(m,ll)
            m=k+3*nl
            ipi(ll*2-1,mm,k)=prodo(m,ll)+cprodo(m,ll)
            m=k+4*nl
            ixr(ll*2-1,mm,k)=prodo(m,ll)+cprodo(m,ll)
            m=k+5*nl
            ixi(ll*2-1,mm,k)=prodo(m,ll)+cprodo(m,ll)
         enddo
      enddo

c     With triangular truncation llmax is variable. If llmax=1 there
c     are no even components so just return.
      if ( llmax .eq. 1 ) return

c
c Merge the pairs of even input vectors, real and imaginary, into a matrix
c

      do lg=1,lat
         do k=1,nl
            fm (k,lg)= ffer(mm,k,lg)
            cfm(k,lg)= gfer(mm,k,lg)
            m=k+nl
            fm (m,lg)= hfei(mm,k,lg)
            cfm(m,lg)= gfei(mm,k,lg)
            m=k+2*nl
            fm (m,lg)= paer(mm,k,lg)
            cfm(m,lg)= pber(mm,k,lg)
            m=k+3*nl
            fm (m,lg)=-paei(mm,k,lg)
            cfm(m,lg)= pbei(mm,k,lg)
            m=k+4*nl
            fm (m,lg)=-xber(mm,k,lg)
            cfm(m,lg)= xaer(mm,k,lg)
            m=k+5*nl
            fm (m,lg)= xbei(mm,k,lg)
            cfm(m,lg)= xaei(mm,k,lg)
         enddo
      enddo
c
c     Carry out the matrix multiply, for the even matrices with
c the input Legendre matrix in even parts : plme(lat,*) & cplme(lat,*)
c

#ifdef ALPH
        call dgemm('n','n',nstack,llmax/2,lat,1.0,fm,nstack
     &             ,plme,lat,0.,prode,nstack)
        call dgemm('n','n',nstack,llmax/2,lat,1.0,cfm,nstack
     &             ,cplme,lat,0.,cprode,nstack)
#else
        call mmtx( fm, plme, prode,nstack,lat,llmax/2)
        call mmtx(cfm,cplme,cprode,nstack,lat,llmax/2)
#endif

c
c     Add the resultant matrices and then split into a real output 
c     vector and an imaginary output vector 
c
      do k=1,nl
         do ll=1,llmax/2
            itr(ll*2,mm,k)=prode(k,ll)+cprode(k,ll)
            m=k+nl
            iti(ll*2,mm,k)=prode(m,ll)+cprode(m,ll)
            m=k+2*nl
            ipr(ll*2,mm,k)=prode(m,ll)+cprode(m,ll)
            m=k+3*nl
            ipi(ll*2,mm,k)=prode(m,ll)+cprode(m,ll)
            m=k+4*nl
            ixr(ll*2,mm,k)=prode(m,ll)+cprode(m,ll)
            m=k+5*nl
            ixi(ll*2,mm,k)=prode(m,ll)+cprode(m,ll)
         enddo
      enddo

      return
      end
