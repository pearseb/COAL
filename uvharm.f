c Removed DEFAULT(NONE) clauses from the OpenMP directives, in order to fix
c problems on both the SC and blizzard.
c SJP 2004/03/16
c
c Manually parallelised.
c SJP 2003/04/24
c
c $Log: uvharm.f,v $
c Revision 1.8  1993/12/17 15:34:25  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.7  93/12/06  16:55:37  ldr
c Changes to V4-4-45l to get T63 running on SGI (Seca).
c 
c Revision 1.6  93/10/05  13:07:49  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.5  93/07/06  16:21:40  ldr
c      Added small extension array to common block elmuv
c      (for overshoot in ptogcray). Added code at end of uvharm.f to set zero
c 
c Revision 1.4  92/12/09  14:44:47  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.3  92/04/22  12:00:38  ldr
c Put /fldri in include file.
c 
c Revision 1.2  91/03/13  13:01:20  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:38:23  ldr
c Initial release V3-0
c 
C$MP_SCHEDTYPE=INTERLEAVE
      subroutine uvharm
      include 'PARAMS.f'
      include 'CNSTE.f'
      include 'FLDRI.f'
      include 'DIFM.f'
      common/elmuv/ulmr(lw1,mw,nl),ulmi(lw1,mw,nl)
     & ,vlmr(lw1,mw,nl),vlmi(lw1,mw,nl)
     & ,elm4x(4)
      include 'WORKRID.f'
      include 'LLMAX.f'

C**** TO CALCULATE U(L,M),V(L,M) and the frictional dissipation terms
C**** Transfer values from dpsir,dpsii,dhxir,dxhii into ipr,ipi,ixr,ixi
C**** before creating frur,frui,frvr,frvi (since same arrays are used
C**** for both to save space)
C**** The arrys ipr,ipi,ixr,ixi will be reset to zero by ZEROGI later

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (k, ll, mm)
!$OMP& SHARED  (dpsii, dpsir, dxhii, dxhir, ipi, ipr, ixi, ixr)

      do 50 k=1,nl
      do 50 mm=1,mw
      do 50 ll=1,llmax(mm)
      ipr(ll,mm,k)=dpsir(ll,mm,k)
      ipi(ll,mm,k)=dpsii(ll,mm,k)
      ixr(ll,mm,k)=dxhir(ll,mm,k)
   50 ixi(ll,mm,k)=dxhii(ll,mm,k)

!$OMP END PARALLEL DO

C**** create values to lw+1 note.

!$OMP  PARALLEL DO SCHEDULE(DYNAMIC)
!$OMP& PRIVATE (fm, k, ll, mm, x1, x2)
!$OMP& SHARED  (epsi, flm, frui, frur, frvi, frvr, ipi, ipr, ixi, ixr,
!$OMP&          psii, psir, ulmi, ulmr, vlmi, vlmr, xhii, xhir)

      do 100 k=1,nl
      do 100 mm=1,mw
      fm=mm-1.0
      x2=(flm(1,mm)+2.0)*epsi(2,mm)
      ulmr(1,mm,k)=-x2*psir(2,mm,k)-fm*xhii(1,mm,k)
      ulmi(1,mm,k)=-x2*psii(2,mm,k)+fm*xhir(1,mm,k)
      vlmr(1,mm,k)= x2*xhir(2,mm,k)-fm*psii(1,mm,k)
      vlmi(1,mm,k)= x2*xhii(2,mm,k)+fm*psir(1,mm,k)
      frur(1,mm,k)=-x2*ipr(2,mm,k)-fm*ixi(1,mm,k)
      frui(1,mm,k)=-x2*ipi(2,mm,k)+fm*ixr(1,mm,k)
      frvr(1,mm,k)= x2*ixr(2,mm,k)-fm*ipi(1,mm,k)
      frvi(1,mm,k)= x2*ixi(2,mm,k)+fm*ipr(1,mm,k)
      do 105 ll=2,lw-1
      x1=(flm(ll,mm)-1.0)*epsi(ll,mm)
      x2=(flm(ll,mm)+2.0)*epsi(ll+1,mm) 
      ulmr(ll,mm,k)=x1*psir(ll-1,mm,k)-x2*psir(ll+1,mm,k)
     & -fm*xhii(ll,mm,k)
      ulmi(ll,mm,k)=x1*psii(ll-1,mm,k)-x2*psii(ll+1,mm,k)
     & +fm*xhir(ll,mm,k)
      vlmr(ll,mm,k)=x2*xhir(ll+1,mm,k)-x1*xhir(ll-1,mm,k)
     & -fm*psii(ll,mm,k)
      vlmi(ll,mm,k)=x2*xhii(ll+1,mm,k)-x1*xhii(ll-1,mm,k)
     & +fm*psir(ll,mm,k)
      frur(ll,mm,k)=x1*ipr(ll-1,mm,k)-x2*ipr(ll+1,mm,k) 
     & -fm*ixi(ll,mm,k)
      frui(ll,mm,k)=x1*ipi(ll-1,mm,k)-x2*ipi(ll+1,mm,k) 
     & +fm*ixr(ll,mm,k)
      frvr(ll,mm,k)=x2*ixr(ll+1,mm,k)-x1*ixr(ll-1,mm,k) 
     & -fm*ipi(ll,mm,k)
      frvi(ll,mm,k)=x2*ixi(ll+1,mm,k)-x1*ixi(ll-1,mm,k) 
     & +fm*ipr(ll,mm,k)
  105 continue
      x1=(flm(lw,mm)-1.0)*epsi(lw,mm)
      ulmr(lw,mm,k)= x1*psir(lw-1,mm,k)-fm*xhii(lw,mm,k)
      ulmi(lw,mm,k)= x1*psii(lw-1,mm,k)+fm*xhir(lw,mm,k)
      vlmr(lw,mm,k)=-x1*xhir(lw-1,mm,k)-fm*psii(lw,mm,k)
      vlmi(lw,mm,k)=-x1*xhii(lw-1,mm,k)+fm*psir(lw,mm,k)
      frur(lw,mm,k)= x1*ipr(lw-1,mm,k)-fm*ixi(lw,mm,k)
      frui(lw,mm,k)= x1*ipi(lw-1,mm,k)+fm*ixr(lw,mm,k)
      frvr(lw,mm,k)=-x1*ixr(lw-1,mm,k)-fm*ipi(lw,mm,k)
      frvi(lw,mm,k)=-x1*ixi(lw-1,mm,k)+fm*ipr(lw,mm,k)
c     x1=(flm(lw1,mm)-1.0)*epsi(lw1,mm) 
      x1=flm(lw,mm)*epsi(lw1,mm)
      ulmr(lw1,mm,k)=x1*psir(lw,mm,k)
      ulmi(lw1,mm,k)=x1*psii(lw,mm,k)
      vlmr(lw1,mm,k)=-x1*xhir(lw,mm,k)
      vlmi(lw1,mm,k)=-x1*xhii(lw,mm,k)
      frur(lw1,mm,k)=x1*ipr(lw,mm,k)
      frui(lw1,mm,k)=x1*ipi(lw,mm,k)
      frvr(lw1,mm,k)=-x1*ixr(lw,mm,k) 
      frvi(lw1,mm,k)=-x1*ixi(lw,mm,k) 
  100 continue

!$OMP END PARALLEL DO

      do 400 k=1,4
  400 elm4x(k)=0.0
      return
      end 
