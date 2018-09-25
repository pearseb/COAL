c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c $Log: lwr88.f,v $
c Revision 1.10  1998/12/10 00:56:00  ldr
c HBG changes to V5-1-21
c
c Revision 1.9  1997/12/17  23:23:05  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.8  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.7  93/12/17  15:33:03  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.6  93/08/19  15:08:33  ldr
c Minor cosmetic changes.
c 
c Revision 1.5  93/06/23  14:30:33  mrd
c Made padding of vtemp common block resolution independent
c 
c Revision 1.4  92/05/11  15:13:34  ldr
c Put Include PARAMS.f in source file rather than in RDPARM.f, 
c to avoid nested includes.
c 
c Revision 1.3  92/04/16  16:41:48  ldr
c 
c Reinstated common vtemp (needed for SGI.)
c 
c Revision 1.2  92/04/15  12:27:29  mrd
c Restructured radiation code include files and data input
c 
c Revision 1.1  91/02/22  16:37:37  ldr
c Initial release V3-0
c 
      subroutine lwr88

!$OMP THREADPRIVATE ( /KDACOM/ )
!$OMP THREADPRIVATE ( /RADISW/ )
!$OMP THREADPRIVATE ( /TFCOM/ )
!$OMP THREADPRIVATE ( /VTEMP/ )

c     subroutine lwr88 computes temperature-corrected co2 transmission
c   functions and also computes the pressure grid and layer optical 
c   paths.
c          inputs:                (common blocks) 
c      press,temp,rh2o,qo3             radisw 
c      co251,co258,cdt51,cdt58         co2bd3 
c      c2d51,c2d58,co2m51,co2m58       co2bd3 
c      cdtm51,cdtm58,c2dm51,c2dm58     co2bd3 
c      stemp,gtemp                     co2bd3 
c      co231,co238,cdt31,cdt38         co2bd2 
c      c2d31,c2d38                     co2bd2 
c      co271,co278,cdt71,cdt78         co2bd4 
c      c2d71,c2d78                     co2bd4 
c      betinw                          bdwide 
c          outputs: 
c      co21,co2nbl,co2sp1,co2sp2       tfcom
c      var1,var2,var3,var4,cntval      kdacom 
c      qh2o,p,delp,delp2,t             kdacom 
c          called by: 
c      radmn or input routine of model
c          calls: 
c      fst88
c 

C Global parameters
      include 'HCON.f'
      include 'PARAMS.f'
      include 'RDPARM.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'KDACOM.f'
      include 'RADISW.f'
      include 'TFCOM.f'
      common / vtemp / co2r(imax,lp1,lp1),dift(imax,lp1,lp1)
      common / vtemp / co2r1(imax,lp1),dco2d1(imax,lp1)
      common / vtemp / d2cd21(imax,lp1),d2cd22(imax,lp1)
      common / vtemp / co2r2(imax,lp1),dco2d2(imax,lp1)
      common / vtemp / co2mr(imax,l),co2md(imax,l),co2m2d(imax,l)
      common / vtemp / tdav(imax,lp1),tstdav(imax,lp1),
     & vv(imax,l),vsum3(imax,lp1),vsum1(imax),vsum2(imax)
      common / vtemp / a1(imax),a2(imax)
      common /vtemp/ dummy(imax*4)
      dimension diftd(imax,lp1,lp1) 
      dimension dco2dt(imax,lp1,lp1),d2cdt2(imax,lp1,lp1) 
      dimension texpsl(imax,lp1),tlsqu(imax,lp1)
      dimension vsum4(imax,l) 
      dimension dift1d(imax,lp1m) 
      equivalence (diftd,co2r)
c     equivalence (dco2dt,d2cdt2,co2r)
      equivalence (vsum3,tlsqu,texpsl)
      equivalence (vv,vsum4)
      equivalence (dift1d,dift) 

C Global data blocks
      include 'CO2DTA.f'
      include 'RNDDTA.f'

C Local work arrays and variables

C Local data, functions etc

C Start code : ----------------------------------------------------------

c****compute flux pressures (p) and differences (delp2,delp)
c****compute flux level temperatures (t) and continuum temperature
c    corrections (texpsl) 
c    changed to use supplied half level pressure and temperatures 
c    m dix 17/8/90
      do 103 k=1,lp1
      do 103 i=1,imax 
      p(i,k)=press2(i,k)
      t(i,k)=temp2(i,k)
103   continue
      do 107 k=1,l
      do 107 i=1,imax 
      delp2(i,k)=p(i,k+1)-p(i,k)
      delp(i,k)=one/delp2(i,k)
107   continue
c****compute argument for cont.temp.coeff.
c    (this is 1800.(1./temp-1./296.)) 
      do 125 k=1,lp1
      do 125 i=1,imax 
      texpsl(i,k)=h18e3/temp(i,k)-h6p08108
c...then take exponential 
      texpsl(i,k)=exp(texpsl(i,k))
125   continue
c***compute optical paths for h2o and o3, using the diffusivity 
c   approximation for the angular integration (1.66). obtain the
c   unweighted values(var1,var3) and the weighted values(var2,var4).
c   the quantities h3m4(.0003) and h3m3(.003) appearing in the var2 and 
c   var4 expressions are the approximate voigt corrections for h2o and
c   o3,respectively.
c 
      do 131 k=1,l
      do 131 i=1,imax 
      qh2o(i,k)=rh2o(i,k)*diffctr 
c---vv is the layer-mean pressure (in atm),which is not the same as 
c   the level pressure (press)
      vv(i,k)=haf*(p(i,k+1)+p(i,k))*p0inv 
      var1(i,k)=delp2(i,k)*qh2o(i,k)*ginv
      var3(i,k)=delp2(i,k)*qo3(i,k)*diffctr*ginv
      var2(i,k)=var1(i,k)*(vv(i,k)+h3m4)
      var4(i,k)=var3(i,k)*(vv(i,k)+h3m3)
c  compute optical path for the h2o continuum, using roberts coeffs.
c  (betinw),and temp. correction (texpsl). the diffusivity factor 
c  (which cancels out in this expression) is assumed to be 1.66. the
c  use of the diffusivity factor has been shown to be a significant 
c  source of error in the continuum calcs.,but the time penalty of
c  an angular integration is severe.
c 
      cntval(i,k)=texpsl(i,k)*rh2o(i,k)*var2(i,k)*betinw/
     &             (rh2o(i,k)+rath2omw)
131   continue
c***compute weighted temperature (tdav) and pressure (tstdav) integrals 
c   for use in obtaining temp. difference bet. sounding and std.
c   temp. sounding (dift) 
      do 161 i=1,imax 
      tstdav(i,1)=zero
      tdav(i,1)=zero
161   continue
      do 162 k=1,lp1
      do 162 i=1,imax 
      vsum3(i,k)=temp(i,k)-stemp(k) 
162   continue
      do 163 k=1,l
      do 165 i=1,imax 
      vsum2(i)=gtemp(k)*delp2(i,k)
      vsum1(i)=vsum2(i)*vsum3(i,k)
      tstdav(i,k+1)=tstdav(i,k)+vsum2(i)
      tdav(i,k+1)=tdav(i,k)+vsum1(i)
165   continue
163   continue
c***compute dift
      do 204 k=1,l
      do 204 kp=k+1,lp1
      do 204 i=1,imax
      dift(i,kp,k)=(tdav(i,kp)-tdav(i,k))/
     &              (tstdav(i,kp)-tstdav(i,k))
204   continue
      do 205 i=1,imax
      dift(i,1,1)=0.
205   continue
      do 206 k=1,l
      do 206 i=1,imax 
      dift(i,k+1,k+1)=haf*(vsum3(i,k+1)+vsum3(i,k))
206   continue
      do 207 k=2,lp1
      do 207 kp=1,k-1
      do 207 i=1,imax
      dift(i,kp,k)=dift(i,k,kp)
207   continue
c 
c****evaluate coefficients for co2 pressure interpolation (a1,a2) 
      do 171 i=1,imax 
      a1(i)=(press(i,lp1)-p0xzp8)/p0xzp2
      a2(i)=(p0-press(i,lp1))/p0xzp2
171   continue
c***perform co2 pressure interpolation on all inputted transmission 
c   functions and temp. derivatives 
c---successively computing co2r,dco2dt and d2cdt2 is done to save 
c   storage (at a slight loss in computation time)
      do 184 k=1,lp1
      do 184 i=1,imax 
        co2r1(i,k)=a1(i)*co231(k)+a2(i)*co238(k)
        d2cd21(i,k)=h1m3*(a1(i)*c2d31(k)+a2(i)*c2d38(k))
        dco2d1(i,k)=h1m2*(a1(i)*cdt31(k)+a2(i)*cdt38(k))
        co2r2(i,k)=a1(i)*co271(k)+a2(i)*co278(k)
        d2cd22(i,k)=h1m3*(a1(i)*c2d71(k)+a2(i)*c2d78(k))
        dco2d2(i,k)=h1m2*(a1(i)*cdt71(k)+a2(i)*cdt78(k))
184   continue
      do 190 k=1,l
      do 190 i=1,imax 
        co2mr(i,k)=a1(i)*co2m51(k)+a2(i)*co2m58(k)
        co2md(i,k)=h1m2*(a1(i)*cdtm51(k)+a2(i)*cdtm58(k))
        co2m2d(i,k)=h1m3*(a1(i)*c2dm51(k)+a2(i)*c2dm58(k))
190   continue
c***compute co2 temperature interpolations for all bands,using dift 
      do 240 k=1,lp1
      do 240 kp=1,lp1 
      do 240 i=1,imax 
      co2r(i,kp,k)=a1(i)*co251(kp,k)+a2(i)*co258(kp,k)
      dco2dt(i,kp,k)=h1m2*(a1(i)*cdt51(kp,k)+a2(i)*cdt58(kp,k))
      d2cdt2(i,kp,k)=h1m3*(a1(i)*c2d51(kp,k)+a2(i)*c2d58(kp,k))
      co21(i,kp,k)=co2r(i,kp,k)+dift(i,kp,k)*(dco2dt(i,kp,k)+
     &             haf*dift(i,kp,k)*d2cdt2(i,kp,k))
240   continue
c***compute transmission fctns used in spa88
c---(in the 250 loop,dift really should be (i,1,k), but dift is 
c    invariant with respect to k,kp,and so (i,1,k)=(i,k,1)) 
      do 250 k=1,lp1
      do 250 i=1,imax 
      co2sp1(i,k)=co2r1(i,k)+dift(i,k,1)*(dco2d1(i,k)+haf*dift(i,k,1)*
     & d2cd21(i,k)) 
      co2sp2(i,k)=co2r2(i,k)+dift(i,k,1)*(dco2d2(i,k)+haf*dift(i,k,1)*
     & d2cd22(i,k)) 
250   continue
c--- we aren't doing nbl tfs on the 100 cm-1 bands .
      do 260 k=1,l
      do 260 i=1,imax 
      co2nbl(i,k)=co2mr(i,k)+dift(i,k,k+1)*(co2md(i,k)+haf* 
     & dift(i,k,k+1)*co2m2d(i,k)) 
260   continue
c***compute temp. coefficient based on t(k) (see ref.2) 
      do 264 k=1,lp1
      do 264 i=1,imax
      if (t(i,k).le.h25e2) then
         tlsqu(i,k)=b0+(t(i,k)-h25e2)*
     &                      (b1+(t(i,k)-h25e2)*
     &                   (b2+b3*(t(i,k)-h25e2))) 
      else 
         tlsqu(i,k)=b0 
      endif
264   continue
c***apply to all co2 tfs
      do 280 k=1,lp1
      do 282 kp=1,lp1 
      do 282 i=1,imax 
      co21(i,kp,k)=co21(i,kp,k)*(one-tlsqu(i,kp))+tlsqu(i,kp) 
282   continue
280   continue
      do 284 k=1,lp1
      do 286 i=1,imax 
      co2sp1(i,k)=co2sp1(i,k)*(one-tlsqu(i,1))+tlsqu(i,1) 
      co2sp2(i,k)=co2sp2(i,k)*(one-tlsqu(i,1))+tlsqu(i,1) 
286   continue
284   continue
      do 288 k=1,l
      do 290 i=1,imax 
      co2nbl(i,k)=co2nbl(i,k)*(one-tlsqu(i,k))+tlsqu(i,k) 
290   continue
288   continue
      call fst88
      return
      end 
