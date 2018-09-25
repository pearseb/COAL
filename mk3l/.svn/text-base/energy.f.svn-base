c Add calls to REAL to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/15
c
c Fixing syntax errors identified by the g95 Fortran compiler.
c SJP 2009/04/14
c
c $Log: energy.f,v $
c Revision 1.15  1999/06/16 06:21:53  rot032
c HBG changes to V5-3
c
c Revision 1.14  1996/10/24 01:02:42  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.13  1996/03/21  03:18:39  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.12  1994/08/08  17:21:11  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.11  93/12/17  15:32:22  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.10  93/10/05  13:06:01  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.9  93/08/10  16:13:47  ldr
c Merge of HBG changes to V4-3-28l with other changes since then.
c 
c Revision 1.8  93/08/03  11:25:45  ldr
c Replace intrinsic functions with generic forms where possible.
c 
c Revision 1.7.1.1  93/08/10  15:27:08  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c 
c Revision 1.7  93/07/06  16:17:45  ldr
c      Added small extension array to common block elmuv
c      (for overshoot in ptogcray). Added code at end of uvharm.f to set zero
c 
c Revision 1.6  92/12/10  09:55:26  ldr
c Minor fixes.
c 
c Revision 1.5  92/12/09  15:57:53  ldr
c More ogcm common block into include files.
c
c     INPUT/OUTPUT:
c     Input:   from common/cnsta in CNSTA.f
c                  dsk - sigma thicknesses (1 to nl) 
c
c              from common/const in CONST.f
c                  sq2 - sqrt(2.0)    tomg - 2 X omega
c
c              from common/fldri in FLDRI.f
c                  psii - spectral stream function, imaginary part
c                  psir - spectral stream function, real part
c                  xhii - velocity potental, imaginary part
c                  xhir - velocity potental, real part
c
c              from common/si in SI.f
c                  am - geopotential height matrix
c
c              from common/timex in TIMEX.f
c                  mins  - current model time in mins
c                  mstep - timestep in minutes 
c
c              from common/worka in WORKA.f
c                  emean - preset global mean energy const/model level
c                  pbar  - global mean surface pressure
c                  phisb - global mean geopotential height
c                  phnb  - approx global mean model level geopotential height
c                  tmean - preset isothermal temperature
c
c              from arguments
c                  iener   - if true, generate detailed energy diagnost
c                  ind     - if 0 = first pass
c
c     In/Out:  from common/fldri in FLDRI.f
c                  pi  - spectral pressure, imaginary part
c                  pr  - spectral pressure, real part
c                  tei - spectral temp, imaginary part
c                  ter - spectral temp, real part
c
c              from common/glmean in GLMEAN.f
c                  totkei - total kinetic energy
c                  totkzi - total kinetic energy, zonal average
c
c              from common/mountn in this subroutine
c                  phisi - spectral surface geopotential height, imaginary part
c                  phsir - spectral surface geopotential height, real part
c
c 
      subroutine energy(iener,ind)

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      logical iener
      integer ind

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'
      include 'FLDRI.f'
      include 'GLMEAN.f'   !Contains common block glmean used for global 
      include 'SI.f'
      include 'TIMEX.f'
      include 'WORKA.f'
      include 'WORKRID.f'
      common/mountn/phisr(lw,mw),phisi(lw,mw)

C Local work arrays and variables
      real totk(nl),totkz(nl),totp(nl)
      real tote(nl),totm(nl),psitot(nl),xhitot(nl),totke(nl),tottem(nl)
     & ,totsf(nl),totvp(nl)
      real sumat(nl)
      complex tet,phig,pgam

C Local data, functions etc

C Start code : ----------------------------------------------------------

      r2=eradsq
      ww=tomg*0.5
      poo=prr(1,1)
      pmb=poo/sq2
      pmb2=pmb*pmb
      psoo=log(pmb)
      cons=2.0/(poo*sqrt(3.0))
      cons2=ww*prr(3,1)/sqrt(15.0)
      cv=cp-rdry
C**** P () HOLDS P* SPECTRALLY IN MBS (*100 FOR MKS UNITS)
      vpress=0.0
      do 1 k=1,nl
      totk(k)=0.0
      totp(k)=0.0
      totkz(k)=0.0
      totke(k)=0.0
      totsf(k)=0.0
      totvp(k)=0.0
      tottem(k)=0.0
      psitot(k)=0.0
      xhitot(k)=0.0
    1 continue
      psphis=0.0
C**
C**
      do 20 mm=1,mw
      mm2=mm-2
      do 20 ll=1,lw
      l=ll+mm2
      fllp=l*(l+1)
      fllp2=fllp*fllp
      pgam=cmplx(prr(ll,mm),-pri(ll,mm))
      pressk=real(pgam*cmplx(prr(ll,mm),pri(ll,mm))*2.0, 8)
      if(mm.gt.1)go to 10
      pressk=pressk*0.5
      if(l.eq.0)go to 11
   10 vpress=vpress+pressk
   11 continue
      phig=pgam*cmplx(phisr(ll,mm),phisi(ll,mm))
      phig=phig+conjg(phig)
      if(mm.eq.1)phig=phig*0.5
      psphis=real(psphis+phig, 8)
      do 20 k=1,nl
      psisq=psir(ll,mm,k)*psir(ll,mm,k)+psii(ll,mm,k)*psii(ll,mm,k)
      psisq=psisq+psisq
      xhisq=xhir(ll,mm,k)*xhir(ll,mm,k)+xhii(ll,mm,k)*xhii(ll,mm,k)
      xhisq=xhisq+xhisq
      tet=cmplx(ter(ll,mm,k),tei(ll,mm,k))
     & +tmean(k)*cmplx(prr(ll,mm),pri(ll,mm))
      temsq=real(tet*conjg(tet), 8)
      temsq=temsq+temsq
      if(mm.gt.1)go to 15
      psisq=psisq*0.5
      xhisq=xhisq*0.5
      temsq=temsq*0.5
   15 if(l.eq.0)go to 20
      psitot(k)=psitot(k)+psisq*fllp2
      xhitot(k)=xhitot(k)+xhisq*fllp2
      tottem(k)=tottem(k)+temsq
      totsf(k)=totsf(k)+psisq
      totvp(k)=totvp(k)+xhisq
      tk=(psisq+xhisq)*fllp
      if(mm.gt.1)go to 18
      totkz(k)=totkz(k)+tk
      go to 20
   18 totke(k)=totke(k)+tk
   20 continue
      vpress=sqrt(vpress*0.5e+04)
      psphis=psphis*0.5/pmb
C**
C**
      if(ind.eq.0)then
c.... elr not defined on first pass (ind=0)
        do 24 k=1,nl
   24   elr(1,1,k)=0.0
      end if

      do 55 k=1,nl
      sumat(k)=0.0
      do 55 j=1,nl
   55 sumat(k)=sumat(k)+am(k,j)*ter(1,1,j)

      do 25 k=1,nl
C**   TOTM=ANG MMTM/UNIT MASS
      totm(k)=(2.0*ww/3.0-cons*(cons2+psir(1,2,k)))*r2
      totk(k)=emean(k)*r2+phisb+elr(1,1,k)*r2/poo-psphis
      totp(k)=cv*(ter(1,1,k)/poo+tmean(k))+phnb(k)+sumat(k)/poo+psphis
      totkz(k)=totkz(k)*r2*0.25/pmb2
      totke(k)=totke(k)*r2*0.25/pmb2
      psitot(k)=sqrt(psitot(k)*0.5)/pmb
      xhitot(k)=sqrt(xhitot(k)*0.5)/pmb
      totsf(k)=sqrt(totsf(k)*0.5)*r2 *1.0e-06/pmb
      totvp(k)=sqrt(totvp(k)*0.5)*r2 *1.0e-06/pmb
      tottem(k)=sqrt(tottem(k)*0.5)/pmb
   25 tote(k)=totk(k)+totp(k)
   
      totmi=0.0
      totpi=0.0
      totki=0.0
      totkei=0.0
      totkzi=0.0
      do 251 k=1,nl
      delsig=dsk(k)
      totmi=totmi+delsig*totm(k)
      totpi=totpi+delsig*totp(k)
      totki=totki+delsig*totk(k)
      totkei=totkei+delsig*totke(k)
  251 totkzi=totkzi+delsig*totkz(k)
      totei=totki+totpi
C**
C**
      if(mod(mins+int(mstep, 8),1440_8).eq.0_8)then
        if(iener)then
          write(6,30)psoo,pbar,vpress,totki,totpi,totei,totmi,totkzi,
     &         totkei,(k,totk(k),k,totp(k),k,tote(k),k,totm(k),
     &         k,psitot(k),k,xhitot(k),k=1,nl)
          write(6,31)(k,totkz(k),k,totke(k),k,tottem(k),k,totsf(k),k,
     &         totvp(k),k=1,nl)
        endif
      endif

 30     format(1h0,'log surface pressure=',e14.7,3x,
     &'surface pressure=',e14.7,3x,'press variance=',e14.7,/,
     &' ','     k=',e14.7,'     p=',e14.7,'     k+p=',e14.7,'     mom=',
     &e14.7,' kzon=',e14.7,' kedd=',e14.7,/,
     & (' ',' k(',i1,')=',e14.7,'  p(',i1,')=',e14.7,'  k+p(',i1,')=',
     &e14.7,'  mom(',i1,')=',e14.7,' v(',i1,')=',e14.7,
     &    ' d(',i1,')=',e14.7))

   31 format(1h ,'kz(',i1,')=',e14.7,' ke(',i1,')=',e14.7,
     & ' t(',i1,')=',e14.7,' sf(',i1,')=',e14.7,' vp(',i1,')=',e14.7)
      return
      end
