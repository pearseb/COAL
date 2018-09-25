c $Log: diffn.f,v $
c Revision 1.22  2001/02/12 05:39:52  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.21  2000/06/20 02:08:33  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.20  1999/05/20 06:23:53  rot032
c HBG changes to V5-2
c
c Revision 1.19  1998/12/10  00:55:40  ldr
c HBG changes to V5-1-21
c
c Revision 1.18  1996/10/24  01:02:36  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.17  1996/03/21  03:18:33  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.16  1995/08/31  04:30:41  ldr
c Tidy-ups from HBG affecting character arrays and generality of NL
c
c Revision 1.15  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.14  1994/08/08  17:21:00  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.13  94/03/22  15:36:08  ldr
c Added extra diffusion at top in 18 level version.
c 
c Revision 1.12  93/12/06  16:55:24  ldr
c Changes to V4-4-45l to get T63 running on SGI (Seca).
c 
c Revision 1.11  93/10/05  13:05:48  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.10  93/08/06  11:48:00  ldr
c Remove extra diffusion on upper level divergence for R42.
c 
c Revision 1.9  93/07/06  16:27:49  ldr
c Only use split P* term for temp diffusion if not using Chen temp variable
c 
c Revision 1.8  93/06/16  14:10:30  ldr
c Changes for Chen temp variable merged with LDR changes to V4-3.
c 
c Revision 1.7  93/06/15  12:28:28  ldr
c Increase diffusion coefficient for divergence only at top two levels at R42.
c
c Revision 1.6.1.1  93/06/15  16:30:27  ldr
c Hal's changes to V4-3 for new Chen temperature variable.
c 
c Revision 1.6  92/12/10  09:56:15  ldr
c Replace np's with nlp's for compatibility with ogcm.
c
c     INPUT/OUTPUT:
c     Input:   from common/cnsta in CNSTA.f
c                  algf - log(sig(k))
c                  sig  - sigma values    sigk - sigma**cappa
c
c              from common/cnste in CNSTE.f
c                  flm  - l part of spectral form
c                  flm2 - l(l+1) part of spectral form
c                  toph - top half of diffusion indicator matrix
c
c              from common/fewflags in FEWFLAGS.f
c                  chenflag - flag for CHEN version of model
c
c              from common/fldri in FLDRI.f
c                  pi - spectral pressure, imaginary part
c                  pr - spectral pressure, real part
c
c              from common/timex in TIMEX.f
c                  mins - current model time in mins
c
c              from common/worka in WORKA.f
c                  tbar  - global mean temp at each model level
c                  tmean - preset isothermal temperature
c
c              from arguments
c                  tdt - 2 X timestep (dt)
c
c     In/Out:  from common/fldri in FLDRI.f
c                  psii - spectral stream function, imaginary part
c                  psir - spectral stream function, real part
c                  tei  - spectral temp, imaginary part
c                  ter  - spectral temp, real part
c                  xhii - velocity potental, imaginary part
c                  xhir - velocity potental, real part
c 
C$MP_SCHEDTYPE=INTERLEAVE

      subroutine diffn(tdt)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      real tdt

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'
      include 'CNSTE.f'
      include 'DIFM.f'
      include 'FEWFLAGS.f'
      include 'FLDRI.f'
      include 'HYBARR.f'
      include 'RKLF.f'
      include 'TIMEX.f'
      include 'WORKA.f'

C Local work arrays and variables
      real dftr(lw,mw),dfps(lw,mw),dfxh(lw,mw,nl)
      real dftri(lw,mw)
      real bij(nlp),dstds(nl),akhd(nl)

      integer k
      integer ll
c     integer mm ! Defined in LLMAX.f

      real akh
      real akhdr
      real bndsig
      real dfxhi
      real rx
      real t00ec
      real t1ec
      real t0ec

C Local data, functions etc
      include 'LLMAX.f'

C Start code : ----------------------------------------------------------

C**** TO ADD DIFFUSION
C**** FOR WAVES IN UPPER PART OF RHOMBOID (SELECT IN INITAL)
C**** ADD DIFFUSION ( AT T+1 ) VIA A SPLIT TIME SCHEME

c Horizontal diffusion coefficient/rad**2
      akh=2.5e+05/eradsq

C**** akhd is multiplier for divergence diffusion coefficient
      akhdr=1.0
      if(mins.le.4320)akhdr=0.0 ! When using divergence dissipation
      do k=1,nl
        akhd(k)=akhdr
      enddo
      if(nl.ge.18)then !Increase diffusion on xhi at top 2 levels 
        akhd(nl-1)=1.25*akhdr
        akhd(nl)  =2.50*akhdr
      endif
      do 1010 k=2,nl
 1010 bij(k)=(tbar(k-1)-tbar(k))*0.5/(algf(k-1)-algf(k))
      bij(1)=bij(2)
      bij(nlp)=bij(nl)

      If(chenflag)Then
CHEN
      t00ec=288.0
      t1ec=0.6652*t00ec
      t0ec=t00ec-t1ec
C**** DSTDS = T + d(T)/d(ln(sigma)) USING GLOBAL MEAN T
c****  - T0 -T1(1+K).sig**K
c**** This is same as T" + d(T")/d(ln(sigma)) where T" is
c****  the (global mean) Chen type temp variable.
      do 1011 k=1,nl
      bndsig=(bnf(k)/sig(k))
      dstds(k)=dbdnf(k)*(tbar(k)-t0ec)
     & -t1ec*(dbdnf(k)+cappa*bndsig)*sigk(k)
     & +bndsig*(bij(k)+bij(k+1))
 1011 continue

      Else

C**** DSTDS = T + d(T)/d(ln(sigma)) USING GLOBAL MEAN T
      do 1012 k=1,nl
      dstds(k)=(tbar(k)-tmean(k))*dbdnf(k)
     & +(bnf(k)/sig(k))*(bij(k)+bij(k+1))
 1012 continue

      End If

C**** SET UP THE DIFFUSION MULTIPLIER ARRAYS : SAME FOR EACH
C**** VERTICAL LEVEL
C**** UPPER PART OF SPECTRUM FORM, OR JSF'S DIFFUSION.

      if(jsfdiff)then
        do 13 mm=1,mw
        do 13 ll=1,llmax(mm)
        dftr(ll,mm)=akh*rklfred(1,ll,mm)
   13   dfps(ll,mm)=akh*rklfred(2,ll,mm)
      else
        do 14 mm=1,mw
        do 14 ll=1,llmax(mm)
        dftr(ll,mm)=akh*flm2(ll,mm)*toph(ll,mm)
   14   dfps(ll,mm)=akh*(flm(ll,mm)-1.0)*(flm(ll,mm)+2.0)
      endif

      if(nl.ge.18)then

      do 151 k=1,nl-2
      do 151 mm=1,mw
      do 151 ll=1,llmax(mm)
  151 dfxh(ll,mm,k)=akhd(k)*dfps(ll,mm)*toph(ll,mm)
c---- Top 2 levels have all wavenumber diffusion for Psi & Xhi
c---- Psi & Xhi both may have increased Kh values (akhd factor)
      do 152 k=nl-1,nl
      do 152 mm=1,mw
      do 152 ll=1,llmax(mm)
  152 dfxh(ll,mm,k)=akhd(k)*dfps(ll,mm)

      else

      do 15 k=1,nl
      do 15 mm=1,mw
      do 15 ll=1,llmax(mm)
   15 dfxh(ll,mm,k)=akhd(k)*dfps(ll,mm)*toph(ll,mm)

      endif

      If(.not.chenflag)Then
C**** ADD THE P TO SIGMA CONVERSION TERM FOR TEMP DIFFUSION
C****  INVOLVING P* (SPLIT) if not using Chen temp variable
      do 155 k=1,nl
      do 155 mm=1,mw
      do 155 ll=1,llmax(mm)
      rx=tdt*dstds(k)*dftr(ll,mm)
      ter(ll,mm,k)=ter(ll,mm,k)+rx*prr(ll,mm)
  155 tei(ll,mm,k)=tei(ll,mm,k)+rx*pri(ll,mm)
      End If

      do 16 mm=1,mw
      do 16 ll=1,llmax(mm)
   16 dftri(ll,mm)=1.0/(1.0+tdt*dftr(ll,mm))

C**** COMPUTE THE SPLIT FORWARD IMPLICIT DIFFUSION
      do 156 k=1,nl
      do 156 mm=1,mw
      do 156 ll=1,llmax(mm)
      ter(ll,mm,k) = ter(ll,mm,k)*dftri(ll,mm)
      tei(ll,mm,k) = tei(ll,mm,k)*dftri(ll,mm)
      dfxhi=1.0/(1.0+tdt*dfxh(ll,mm,k))
      psir(ll,mm,k)=psir(ll,mm,k)*dfxhi
      psii(ll,mm,k)=psii(ll,mm,k)*dfxhi
      xhir(ll,mm,k)=xhir(ll,mm,k)*dfxhi
      xhii(ll,mm,k)=xhii(ll,mm,k)*dfxhi
c**** retain spectral diffusion of psi (vorticity) and
c****  xhi (divergence) for calculating energy addition
c****  to temp equation
      dpsir(ll,mm,k)=dfxh(ll,mm,k)*psir(ll,mm,k)
      dpsii(ll,mm,k)=dfxh(ll,mm,k)*psii(ll,mm,k)
      dxhir(ll,mm,k)=dfxh(ll,mm,k)*xhir(ll,mm,k)
  156 dxhii(ll,mm,k)=dfxh(ll,mm,k)*xhii(ll,mm,k)
      return
      end
