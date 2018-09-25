c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c $Log: dynmvo.f,v $
c Revision 1.12  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.11  1999/06/16 06:21:54  rot032
c HBG changes to V5-3
c
c Revision 1.10  1999/05/20 06:23:53  rot032
c HBG changes to V5-2
c
c Revision 1.9  1998/12/10  00:55:42  ldr
c HBG changes to V5-1-21
c
c Revision 1.8  1997/12/17  23:22:52  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.7  1996/10/24  01:02:41  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.6  1996/06/13  02:06:26  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.5  1996/03/21  03:18:38  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.4  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.3  1994/03/30  12:34:16  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c
c Revision 1.2  93/12/17  15:32:17  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.1  93/07/14  15:20:01  ldr
c Initial revision
c
c     INPUT/OUTPUT:
c     Input:   from common/gausl in GAUSL.f
c                  acsq - 1.0/(sia(lg)*sia(lg))
c                  coa  - sin latitude
c
c              from common/ubarvo in this subroutine
c                  ubx  - mean weighted zonal wind
c                  dubx - gradient of above
c
c              from common/workf in this subroutine
c                  vof - fourier series vorticity at current time
c                  vofm - vof at time, t-1
c
c              from arguments
c                  lg  - latitude index
c                  tdt - 2 X timestep (dt)
c
c     Output:  from common/worknsvo in this subroutine
c                  real and imaginary Fourier coefficients for implicit
c                  vorticity version, n=northern hemisphere, s=southern hem.
c
c                  apni, apnr, bpni, bpnr
c
c     In/Out:  from common/worknsd in this subroutine
c                  afni, afnr, bfni, bfnr) - vorticity and divergence
c
c 
      subroutine dynmvo(lg,tdt)

!$OMP THREADPRIVATE ( /UBARVO/ )
!$OMP THREADPRIVATE ( /WORKF/ )
!$OMP THREADPRIVATE ( /WORKNSD/ )

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
       integer lg
       real tdt

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      common/ubarvo/ubx(nl,2),dubx(nl,2)
      complex dvf,vof,tef,cpf,plf,vofm
      common/workf/dvf(mw,2,nl),vof(mw,2,nl),tef(mw,2,nl)
     & ,cpf(mw,2),plf(mw,2),vofm(mw,2,nl)
      common/worknsd/
     &  afnr(mw,nl,2),bfnr(mw,nl,2),efnr(mw,nl,2),ffnr(mw,nl,2)
     & ,gfnr(mw,nl,2),hfnr(mw,nl,2)
     & ,afni(mw,nl,2),bfni(mw,nl,2),efni(mw,nl,2),ffni(mw,nl,2)
     & ,gfni(mw,nl,2),hfni(mw,nl,2)
     & ,apnr(mw,nl,2),bpnr(mw,nl,2)
     & ,apni(mw,nl,2),bpni(mw,nl,2)

C Global data blocks
      include 'CNSTA.f'
      include 'GAUSL.f'
      common/dystb/tpa(nl),tpb(nl),tpc(nl),wpa(nl),wpb(nl),wpc(nl)
     &,qpa(nl),qpb(nl),qpc(nl),plev(nl),fred(nl)

C Local work arrays and variables
      complex bm(mw,nl)
      complex amh(mw,nl),amhp(mw,nl),bmhp(mw,nl)
      complex cx1
      real alpha(nl),cfda(nl)

C Local data, functions etc

C Start code : ----------------------------------------------------------
c
c---- To implement implicit vorticity by changing
c---- the Fourier components for the vorticity tendencies
c---- This requires complex arithmetic.
c
      do 60 ins=1,2

      csqlg=acsq(lg)
c.... ubx=ubar.cos(lat)/erad
      sinphi=coa(lg)*(3.0-2.0*ins)
      do 15 k=1,nl
c optional height reduction factor fred(k) :
        ubx(k,ins)=ubx(k,ins)*fred(k)
        dubx(k,ins)=dubx(k,ins)*fred(k)
      alpha(k)=ubx(k,ins)*csqlg
   15 cfda(k)=tdt*csqlg*(dubx(k,ins)+2.0*sinphi*ubx(k,ins))

      do 10 k=1,nl
      do 10 mm=1,mw
      tma=tdt*(mm-1.0)*alpha(k)
      brm=1.0/(1.0+tma**2)
   10 bm(mm,k)=cmplx(brm,-tma*brm)

      do 20 k=1,nl
      do 20 mm=1,mw
c.... Form Am and Bm (complex)
        amh(mm,k)=cmplx(afnr(mm,k,ins),afni(mm,k,ins))
c.... Form Bm"
        bmhp(mm,k)=bm(mm,k)*
     &            cmplx(bfnr(mm,k,ins),bfni(mm,k,ins))
   20 continue

c.... Form Am"
      do 40 k=1,nl
      tubx=2.0*ubx(k,ins)
      do 40 mm=1,mw
        cx1=amh(mm,k)+tubx*(vofm(mm,ins,k)-vof(mm,ins,k))
     &               +cfda(k)*bmhp(mm,k)
        amhp(mm,k)=bm(mm,k)*cx1
   40 continue
   
c.... Put the components into arrays for use in the vorticity
c.... tendency (spectral) evaluation only.
c.... (apnr,apni,bpnr,bpni must not be used for divergence)
      do 50 k=1,nl
      do 50 mm=1,mw
        apnr(mm,k,ins)=real(amhp(mm,k))
        apni(mm,k,ins)=aimag(amhp(mm,k))
        bpnr(mm,k,ins)=real(bmhp(mm,k))
   50   bpni(mm,k,ins)=aimag(bmhp(mm,k))

   60 continue
      
      return
      end
