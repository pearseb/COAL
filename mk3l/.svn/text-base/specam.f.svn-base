c Add calls to REAL to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/15
c
c $Log: specam.f,v $
c Revision 1.8  1999/06/16 06:21:53  rot032
c HBG changes to V5-3
c
c Revision 1.7  1996/10/24 01:03:16  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.6  1996/03/21  03:19:05  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.5  1993/12/17  15:33:47  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.4  93/10/05  13:07:29  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.3  92/12/09  14:44:32  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.2  91/03/13  13:00:37  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:38:09  ldr
c Initial release V3-0
c
c     INPUT/OUTPUT
c     Input:   from common/cnsta in CNSTA.f
c                  sig - sigma values
c
c              from common/const in CONST.f
c                  sq2  - sqrt(2.0)    
c                  tomg - 2 X omega
c
c              from common/worka in WORKA.f
c                  tmean - preset isothermal temperature
c
c
      subroutine specam

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'
      include 'FLDRI.f'
      include 'WORKA.f'

C Local work arrays and variables
      real ampre(mw),ampsf(mw),ampvor(mw),ampens(mw),ampk(mw)
     & ,ampmke(mw),ampte(mw)
      complex stf,surfp,tem,vor,div

C Local data, functions etc

C Start code : ----------------------------------------------------------

      r2=eradsq
      ww=tomg*0.5
      pmb=prr(1,1)/sq2
      pmb2=pmb*pmb
      do 3333 k=1,nl
      do 600 mm=1,mw
      ampre(mm)=0.0
      ampsf(mm)=0.0
      ampk(mm)=0.0
      ampmke(mm)=0.0
      ampens(mm)=0.0
      ampte(mm)=0.0
  600 ampvor(mm)=0.0
      do 440 mm=1,mw
      m=mm-1
      fllpi=0.0
      do 440 ll=1,lw
      l=m+ll-1
      fllp=l*(l+1)
      if(l.gt.0)fllpi=1.0/fllp
      stf=cmplx(psir(ll,mm,k),psii(ll,mm,k))
      vor=-stf*fllp
      enstsq=real(vor*conjg(vor), 8)
      surfp=cmplx(prr(ll,mm),pri(ll,mm))
      div=-cmplx(xhir(ll,mm,k),xhii(ll,mm,k))*fllp
      divsq=real(div*conjg(div), 8)
      pcsq=enstsq+divsq
      if(m.gt.0)enstsq=enstsq+enstsq
      ampens(mm)=ampens(mm)+enstsq
      amp=real(stf*conjg(stf), 8)
      if(m.gt.0)amp=amp+amp
      ampsf(mm)=ampsf(mm)+amp
      amp=real(vor*conjg(vor), 8)
      if(m.gt.0)amp=amp+amp
      ampvor(mm)=ampvor(mm)+amp
      amp=real(surfp*conjg(surfp), 8)
      if(m.gt.0)amp=amp+amp
      ampre(mm)=ampre(mm)+amp
      amp=divsq
      if(m.gt.0)amp=amp+amp
      ampk(mm)=ampk(mm)+amp
      if(m.gt.0)pcsq=pcsq+pcsq
      ampmke(mm)=ampmke(mm)+pcsq*fllpi
      tem=cmplx(ter(ll,mm,k),tei(ll,mm,k))+tmean(k)*surfp
      amp=real(tem*conjg(tem), 8)
      if(m.gt.0)amp=amp+amp
      ampte(mm)=ampte(mm)+amp
  440 continue
C**** REMOVE P* WEIGHTING APPROXIMATELY
      do 650 mm=1,mw
      amp       =ampens(mm)/(ww*ww) *0.5e02
      ampens(mm)=amp/pmb2
      amp      =sqrt(ampsf(mm)) *r2 *1.0e-06
      ampsf(mm)=amp/pmb
      ampre(mm)=sqrt(ampre(mm))
      amp       =sqrt(ampvor(mm))/ww
      ampvor(mm)=amp/pmb
      amp     =sqrt(ampk(mm))/ww *10.0
      ampk(mm)=amp/pmb
      amp      =sqrt(ampte(mm))
      ampte(mm)=amp/pmb
      amp       =ampmke(mm)*r2*0.25
  650 ampmke(mm)=amp/pmb2
      if(k.ne.1)go to 660
      write(6,500)
  500 format(1h0,'spectral decomposition -summed l amplitudes for each',
     &' m ')
      sigsf=1.0
      write(6,501)sigsf
  501 format(1h0,'surface pressure amplitudes (mbar)',14x,',sigma=',
     & f6.3)
      write(6,502)( ampre(mm),mm ,mm=1,mw)
  502 format(1h , 8(f8.3 ,'(mm=',i2,')',1x),/,
     &1x, 8(f8.3,'(mm=',i2,')',1x),/,1x,6(f8.3,'(mm=',i2,')',1x))
  660 write(6,503)sig(k)
  503 format(1h ,'stream function amplitudes (km**2/sec)',10x,',sigma=',
     & f6.3)
      write(6,502)( ampsf(mm),mm ,mm=1,mw)
      print 504,sig(k)
  504 format(1h ,'vorticity amplitudes (relative to omega)',8x,',sigma='
     &,f6.3)
      write(6,502)( ampvor(mm),mm ,mm=1,mw)
      write(6,505)sig(k)
  505 format(1h ,'divergence amplitudes (relative to omega/10)    ,sigma
     &=',f6.3)
      write(6,502)( ampk(mm),mm ,mm=1,mw)
      write(6,506)sig(k)
  506 format(1h ,'m.k.e. amplitudes in m.k.s',22x,',sigma=',f6.3)
      write(6,502)( ampmke(mm),mm ,mm=1,mw)
      write(6,507)sig(k)
  507 format(1h ,'enstrophy amplitudes relative to omega# 2x1.0e02,sigma
     &=',f6.3)
      write(6,502)( ampens(mm),mm ,mm=1,mw)
      write(6,508)sig(k)
  508 format(1h ,'temperature amplitudes (deg)',20x,',sigma=',f6.3)
      write(6,502)( ampte(mm),mm ,mm=1,mw)
 3333 continue
      return
      end
