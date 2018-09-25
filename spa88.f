c Loops re-rolled so that code passes array-bounds checking.
c SJP 2003/03/08
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c $Log: spa88.f,v $
c Revision 1.15  1998/12/10 00:56:03  ldr
c HBG changes to V5-1-21
c
c Revision 1.14  1997/12/17  23:23:09  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.13  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.11.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.12  1994/08/08  17:22:33  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.11  94/05/27  16:35:33  ldr
c Make some variables double precision for Seca to avoid underflows.
c 
c Revision 1.10  93/12/17  15:33:45  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.9  93/08/19  15:10:47  ldr
c Minor cosmetic changes.
c 
c Revision 1.8  93/07/27  14:58:50  ldr
c Merge contents of common blocks clrtemp and lwoutclr into lwout.
c 
c Revision 1.7  93/07/22  10:17:25  mrd
c Add diagnostics for clear sky net surface radiation and downward radiation
c at surface.
c 
c Revision 1.6  93/06/23  14:30:35  mrd
c Made padding of vtemp common block resolution independent
c 
c Revision 1.5  92/12/09  14:44:30  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.4  92/05/11  15:13:48  ldr
c Put Include PARAMS.f in source file rather than in RDPARM.f, 
c to avoid nested includes.
c 
c     subroutine spa88 computes exact cts heating rates and fluxes and
c  corresponding cts emissivity quantities for h2o,co2 and o3.
c          inputs:                (common blocks) 
c       acomb,bcomb,apcm,bpcm                  bdcomb 
c       atpcm,btpcm,betacm                     bdcomb 
c       betinw                                 bdwide 
c       temp,press                             radisw 
c       var1,var2,p,delp,delp2                 kdacom 
c       totvo2,to3spc,co2sp1,co2sp2,co2sp      tfcom
c       cldfac                                 cldcom 
c       sko2d                                  tabcom 
c       sorc,csour,osour                       srccom 
c           outputs:  
c       excts,exctsn,ctso3                     tfcom
c        gxcts,fctsg                           rdflux 
c           called by:  
c       fst88 
c            calls: 
 
      subroutine spa88
 
!$OMP THREADPRIVATE ( /CLDCOM/ )
!$OMP THREADPRIVATE ( /KDACOM/ )
!$OMP THREADPRIVATE ( /LWOUT/ )
!$OMP THREADPRIVATE ( /RADISW/ )
!$OMP THREADPRIVATE ( /RDFLUX/ )
!$OMP THREADPRIVATE ( /CLRFLX/ )
!$OMP THREADPRIVATE ( /SRCCOM/ )
!$OMP THREADPRIVATE ( /TFCOM/ )
!$OMP THREADPRIVATE ( /VTEMP/ )

C Global parameters
      include 'HCON.f'
      include 'PARAMS.f'
      include 'RDPARM.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      include 'CLDCOM.f'
      include 'KDACOM.f'
      include 'LWOUT.f'
      include 'RADISW.f'
      include 'RDFLUX.f' ! includes CLRFLX
      include 'SRCCOM.f'
      include 'TFCOM.f'
      common / vtemp / phitmp(imax,l),psitmp(imax,l),tt(imax,l),
     &          fac1(imax,l),fac2(imax,l),
     &          ctmp(imax,lp1),x(imax,l),y(imax,l),
     &          topm(imax,l),topphi(imax,l),
     &          ctmp3(imax,lp1),ctmp2(imax,lp1)
      common /vtemp/ dummy(imax*(2*l*l+5*l+16))
      double precision tt,fac1
      dimension f(imax,l),ag(imax,l)
      equivalence (f,ag,phitmp) 
      dimension ff(imax,l),agg(imax,l)
      equivalence (ff,agg,psitmp)

C Global data blocks
      include 'RNDDTA.f'

C Local work arrays and variables

C Local data, functions etc
      integer i, j

C Start code : ----------------------------------------------------------

c---compute temperature quantities for use in program
      do 101 k=1,l
      do 101 i=1,imax 
      x(i,k)=temp(i,k)-h25e2
      y(i,k)=x(i,k)*x(i,k)
101   continue
c---initialize ctmp(i,1),ctmp2(i,1),ctmp3(i,1) to unity; these are 
c   transmission fctns at the top.
      do 345 i=1,imax 
      ctmp(i,1)=one 
      ctmp2(i,1)=1.
      ctmp3(i,1)=1.
345   continue
c***begin loop on frequency bands (1)***
c 
c---calculation for band 1 (combined band 1)  
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 301 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(1)*x(i,1)+bpcm(1)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(1)*x(i,1)+btpcm(1)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 301 j = 1, l
        do 301 i = 1, imax
          f(i,j)=h44194m2*(apcm(1)*x(i,j)+bpcm(1)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(1)*x(i,j)+btpcm(1)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
301   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 315 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
315   continue
      do 319 k=2,l
      do 317 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
317   continue
319   continue
c---tt is the cloud-free cts transmission function
CSJP      do 321 i=1,imax*l
CSJP      fac1(i,1)=acomb(1)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(1)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*fac1(i,1)/sqrt(1.+fac2(i,1)))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 321 j = 1, l
        do 321 i = 1, imax
          fac1(i,j)=acomb(1)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(1)*topphi(i,j))
          tt(i,j)=exp(hm1ez*fac1(i,j)/sqrt(1.+fac2(i,j)))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
321   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 353 i=1,imax*l
CSJP      exctsn(i,1,1)=radcon*delp(i,1)*sorc(i,1,1)*(ctmp(i,2)-ctmp(i,1)) 
CSJP      excts(i,1)=exctsn(i,1,1) 
CSJP      exctsclr(i,1)=radcon*delp(i,1)*sorc(i,1,1)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 353 j = 1, l
        do 353 i = 1, imax
          exctsn(i,j,1)=radcon*delp(i,j)*sorc(i,j,1)*(ctmp(i,j+1)-
     &                                                ctmp(i,j))
          excts(i,j)=exctsn(i,j,1)
          exctsclr(i,j)=radcon*delp(i,j)*sorc(i,j,1)*(ctmp2(i,j+1)-
     &                                                ctmp2(i,j))
353   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 361 i=1,imax 
      fctsg(i,1)=tt(i,l)*sorc(i,l,1)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,1)-sorc(i,l,1))
      gxctsclr(i)=fctsg(i,1)
      fctsg(i,1)=cldfac(i,lp1,1)*fctsg(i,1)
      gxcts(i)=fctsg(i,1)
361   continue
c 
c 
c-----calculation for band 2 (combined band 2)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 401 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(2)*x(i,1)+bpcm(2)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(2)*x(i,1)+btpcm(2)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 401 j = 1, l
        do 401 i = 1, imax
          f(i,j)=h44194m2*(apcm(2)*x(i,j)+bpcm(2)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(2)*x(i,j)+btpcm(2)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
401   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 415 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
415   continue
      do 419 k=2,l
      do 417 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
417   continue
419   continue
c---tt is the cloud-free cts transmission function
CSJP      do 421 i=1,imax*l
CSJP      fac1(i,1)=acomb(2)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(2)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*fac1(i,1)/sqrt(1.+fac2(i,1)))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 421 j = 1, l
        do 421 i = 1, imax
          fac1(i,j)=acomb(2)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(2)*topphi(i,j))
          tt(i,j)=exp(hm1ez*fac1(i,j)/sqrt(1.+fac2(i,j)))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
421   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 453 i=1,imax*l
CSJP      exctsn(i,1,2)=radcon*delp(i,1)*sorc(i,1,2)*(ctmp(i,2)-ctmp(i,1)) 
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,2) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,2)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 453 j = 1, l
        do 453 i = 1, imax
          exctsn(i,j,2)=radcon*delp(i,j)*sorc(i,j,2)*(ctmp(i,j+1)-
     &                                                ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,2)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,2)*(ctmp2(i,j+1)-ctmp2(i,j))
453   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 461 i=1,imax 
      fctsg(i,2)=tt(i,l)*sorc(i,l,2)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,2)-sorc(i,l,2))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,2)
      fctsg(i,2)=cldfac(i,lp1,1)*fctsg(i,2)
      gxcts(i)=gxcts(i)+fctsg(i,2)
461   continue
c 
c-----calculation for band 3 (combined band 3)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 501 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(3)*x(i,1)+bpcm(3)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(3)*x(i,1)+btpcm(3)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 501 j = 1, l
        do 501 i = 1, imax
          f(i,j)=h44194m2*(apcm(3)*x(i,j)+bpcm(3)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(3)*x(i,j)+btpcm(3)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
501   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 515 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
515   continue
      do 519 k=2,l
      do 517 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
517   continue
519   continue
c---tt is the cloud-free cts transmission function
CSJP      do 521 i=1,imax*l
CSJP      fac1(i,1)=acomb(3)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(3)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*fac1(i,1)/sqrt(1.+fac2(i,1)))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 521 j = 1, l
        do 521 i = 1, imax
          fac1(i,j)=acomb(3)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(3)*topphi(i,j))
          tt(i,j)=exp(hm1ez*fac1(i,j)/sqrt(1.+fac2(i,j)))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
521   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 553 i=1,imax*l
CSJP      exctsn(i,1,3)=radcon*delp(i,1)*sorc(i,1,3)*(ctmp(i,2)-ctmp(i,1)) 
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,3) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,3)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 553 j = 1, l
        do 553 i = 1, imax
          exctsn(i,j,3)=radcon*delp(i,j)*sorc(i,j,3)*(ctmp(i,j+1)-
     &                                                ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,3)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,3)*(ctmp2(i,j+1)-ctmp2(i,j))
553   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 561 i=1,imax 
      fctsg(i,3)=tt(i,l)*sorc(i,l,3)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,3)-sorc(i,l,3))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,3)
      fctsg(i,3)=cldfac(i,lp1,1)*fctsg(i,3)
      gxcts(i)=gxcts(i)+fctsg(i,3)
561   continue
c 
c-----calculation for band 4 (combined band 4)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 601 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(4)*x(i,1)+bpcm(4)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(4)*x(i,1)+btpcm(4)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 601 j = 1, l
        do 601 i = 1, imax
          f(i,j)=h44194m2*(apcm(4)*x(i,j)+bpcm(4)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(4)*x(i,j)+btpcm(4)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
601   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 615 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
615   continue
      do 619 k=2,l
      do 617 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
617   continue
619   continue
c---tt is the cloud-free cts transmission function
CSJP      do 621 i=1,imax*l
CSJP      fac1(i,1)=acomb(4)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(4)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*fac1(i,1)/sqrt(1.+fac2(i,1)))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 621 j = 1, l
        do 621 i = 1, imax
          fac1(i,j)=acomb(4)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(4)*topphi(i,j))
          tt(i,j)=exp(hm1ez*fac1(i,j)/sqrt(1.+fac2(i,j)))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
621   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 653 i=1,imax*l
CSJP      exctsn(i,1,4)=radcon*delp(i,1)*sorc(i,1,4)*(ctmp(i,2)-ctmp(i,1)) 
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,4) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,4)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 653 j = 1, l
        do 653 i = 1, imax
          exctsn(i,j,4)=radcon*delp(i,j)*sorc(i,j,4)*(ctmp(i,j+1)-
     &                                                ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,4)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,4)*(ctmp2(i,j+1)-ctmp2(i,j))
653   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 661 i=1,imax 
      fctsg(i,4)=tt(i,l)*sorc(i,l,4)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,4)-sorc(i,l,4))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,4)
      fctsg(i,4)=cldfac(i,lp1,1)*fctsg(i,4)
      gxcts(i)=gxcts(i)+fctsg(i,4)
661   continue
c 
c-----calculation for band 5 (combined band 5)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 701 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(5)*x(i,1)+bpcm(5)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(5)*x(i,1)+btpcm(5)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 701 j = 1, l
        do 701 i = 1, imax
          f(i,j)=h44194m2*(apcm(5)*x(i,j)+bpcm(5)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(5)*x(i,j)+btpcm(5)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
701   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 715 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
715   continue
      do 719 k=2,l
      do 717 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
717   continue
719   continue
c---tt is the cloud-free cts transmission function
CSJP      do 721 i=1,imax*l
CSJP      fac1(i,1)=acomb(5)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(5)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
CSJP     &           betacm(5)*totvo2(i,2)*sko2d))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 721 j = 1, l
        do 721 i = 1, imax
          fac1(i,j)=acomb(5)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(5)*topphi(i,j))
          tt(i,j)=exp(hm1ez*(fac1(i,j)/sqrt(one+fac2(i,j))+
     &            betacm(5)*totvo2(i,j+1)*sko2d))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
721   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 753 i=1,imax*l
CSJP      exctsn(i,1,5)=radcon*delp(i,1)*sorc(i,1,5)*(ctmp(i,2)-ctmp(i,1)) 
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,5) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,5)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 753 j = 1, l
        do 753 i = 1, imax
          exctsn(i,j,5)=radcon*delp(i,j)*sorc(i,j,5)*(ctmp(i,j+1)-
     &                                                ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,5)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,5)*(ctmp2(i,j+1)-ctmp2(i,j))
753   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 761 i=1,imax 
      fctsg(i,5)=tt(i,l)*sorc(i,l,5)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,5)-sorc(i,l,5))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,5)
      fctsg(i,5)=cldfac(i,lp1,1)*fctsg(i,5)
      gxcts(i)=gxcts(i)+fctsg(i,5)
761   continue
c 
c-----calculation for band 6 (combined band 6)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 801 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(6)*x(i,1)+bpcm(6)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(6)*x(i,1)+btpcm(6)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 801 j = 1, l
        do 801 i = 1, imax
          f(i,j)=h44194m2*(apcm(6)*x(i,j)+bpcm(6)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(6)*x(i,j)+btpcm(6)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
801   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 815 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
815   continue
      do 819 k=2,l
      do 817 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
817   continue
819   continue
c---tt is the cloud-free cts transmission function
CSJP      do 821 i=1,imax*l
CSJP      fac1(i,1)=acomb(6)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(6)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
CSJP     &           betacm(6)*totvo2(i,2)*sko2d))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 821 j = 1, l
        do 821 i = 1, imax
          fac1(i,j)=acomb(6)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(6)*topphi(i,j))
          tt(i,j)=exp(hm1ez*(fac1(i,j)/sqrt(one+fac2(i,j))+
     &            betacm(6)*totvo2(i,j+1)*sko2d))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
821   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 853 i=1,imax*l
CSJP      exctsn(i,1,6)=radcon*delp(i,1)*sorc(i,1,6)*(ctmp(i,2)-ctmp(i,1)) 
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,6) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,6)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 853 j = 1, l
        do 853 i = 1, imax
          exctsn(i,j,6)=radcon*delp(i,j)*sorc(i,j,6)*(ctmp(i,j+1)-
     &                                                ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,6)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,6)*(ctmp2(i,j+1)-ctmp2(i,j))
853   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 861 i=1,imax 
      fctsg(i,6)=tt(i,l)*sorc(i,l,6)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,6)-sorc(i,l,6))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,6)
      fctsg(i,6)=cldfac(i,lp1,1)*fctsg(i,6)
      gxcts(i)=gxcts(i)+fctsg(i,6)
861   continue
c 
c-----calculation for band 7 (combined band 7)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 901 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(7)*x(i,1)+bpcm(7)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(7)*x(i,1)+btpcm(7)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 901 j = 1, l
        do 901 i = 1, imax
          f(i,j)=h44194m2*(apcm(7)*x(i,j)+bpcm(7)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(7)*x(i,j)+btpcm(7)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
901   continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 915 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
915   continue
      do 919 k=2,l
      do 917 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
917   continue
919   continue
c---tt is the cloud-free cts transmission function
CSJP      do 921 i=1,imax*l
CSJP      fac1(i,1)=acomb(7)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(7)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
CSJP     &           betacm(7)*totvo2(i,2)*sko2d))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 921 j = 1, l
        do 921 i = 1, imax
          fac1(i,j)=acomb(7)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(7)*topphi(i,j))
          tt(i,j)=exp(hm1ez*(fac1(i,j)/sqrt(one+fac2(i,j))+
     &            betacm(7)*totvo2(i,j+1)*sko2d))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
921   continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 953 i=1,imax*l
CSJP      exctsn(i,1,7)=radcon*delp(i,1)*sorc(i,1,7)*(ctmp(i,2)-ctmp(i,1)) 
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,7) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,7)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 953 j = 1, l
        do 953 i = 1, imax
          exctsn(i,j,7)=radcon*delp(i,j)*sorc(i,j,7)*(ctmp(i,j+1)-
     &                                                ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,7)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,7)*(ctmp2(i,j+1)-ctmp2(i,j))
953   continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 961 i=1,imax 
      fctsg(i,7)=tt(i,l)*sorc(i,l,7)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,7)-sorc(i,l,7))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,7)
      fctsg(i,7)=cldfac(i,lp1,1)*fctsg(i,7)
      gxcts(i)=gxcts(i)+fctsg(i,7)
961   continue
c 
c-----calculation for band 8 (combined band 8)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 1001 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(8)*x(i,1)+bpcm(8)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(8)*x(i,1)+btpcm(8)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 1001 j = 1, l
        do 1001 i = 1, imax
          f(i,j)=h44194m2*(apcm(8)*x(i,j)+bpcm(8)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(8)*x(i,j)+btpcm(8)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
1001  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1015 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1015  continue
      do 1019 k=2,l
      do 1017 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1017  continue
1019  continue
c---tt is the cloud-free cts transmission function
CSJP      do 1021 i=1,imax*l
CSJP      fac1(i,1)=acomb(8)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(8)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
CSJP     &           betacm(8)*totvo2(i,2)*sko2d))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 1021 j = 1, l
        do 1021 i = 1, imax
          fac1(i,j)=acomb(8)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(8)*topphi(i,j))
          tt(i,j)=exp(hm1ez*(fac1(i,j)/sqrt(one+fac2(i,j))+
     &            betacm(8)*totvo2(i,j+1)*sko2d))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
1021  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 1053 i=1,imax*l
CSJP      exctsn(i,1,8)=radcon*delp(i,1)*sorc(i,1,8)*(ctmp(i,2)-ctmp(i,1)) 
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,8) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,8)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 1053 j = 1, l
        do 1053 i = 1, imax
          exctsn(i,j,8)=radcon*delp(i,j)*sorc(i,j,8)*(ctmp(i,j+1)
     &                                                -ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,8)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,8)*(ctmp2(i,j+1)-ctmp2(i,j))
1053  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1061 i=1,imax 
      fctsg(i,8)=tt(i,l)*sorc(i,l,8)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,8)-sorc(i,l,8))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,8)
      fctsg(i,8)=cldfac(i,lp1,1)*fctsg(i,8)
      gxcts(i)=gxcts(i)+fctsg(i,8)
1061  continue
c 
c-----calculation for band 9 ( 560-670 cm-1; includes co2)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 1101 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(9)*x(i,1)+bpcm(9)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(9)*x(i,1)+btpcm(9)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 1101 j = 1, l
        do 1101 i = 1, imax
          f(i,j)=h44194m2*(apcm(9)*x(i,j)+bpcm(9)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(9)*x(i,j)+btpcm(9)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
1101  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1115 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1115  continue
      do 1119 k=2,l
      do 1117 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1117  continue
1119  continue
c---tt is the cloud-free cts transmission function
CSJP      do 1121 i=1,imax*l
CSJP      fac1(i,1)=acomb(9)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(9)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
CSJP     &           betacm(9)*totvo2(i,2)*sko2d))*co2sp1(i,2)
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 1121 j = 1, l
        do 1121 i = 1, imax
          fac1(i,j)=acomb(9)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(9)*topphi(i,j))
          tt(i,j)=exp(hm1ez*(fac1(i,j)/sqrt(one+fac2(i,j))+
     &            betacm(9)*totvo2(i,j+1)*sko2d))*co2sp1(i,j+1)
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
1121  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 1153 i=1,imax*l
CSJP      exctsn(i,1,9)=radcon*delp(i,1)*sorc(i,1,9)*(ctmp(i,2)-ctmp(i,1)) 
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,9) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,9)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 1153 j = 1, l
        do 1153 i = 1, imax
          exctsn(i,j,9)=radcon*delp(i,j)*sorc(i,j,9)*(ctmp(i,j+1)-
     &                                                ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,9)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,9)*(ctmp2(i,j+1)-ctmp2(i,j))
1153  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1161 i=1,imax 
      fctsg(i,9)=tt(i,l)*sorc(i,l,9)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,9)-sorc(i,l,9))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,9)
      fctsg(i,9)=cldfac(i,lp1,1)*fctsg(i,9)
      gxcts(i)=gxcts(i)+fctsg(i,9)
1161  continue
c 
c-----calculation for band 10 (670-800 cm-1; includes co2)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 1201 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(10)*x(i,1)+bpcm(10)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(10)*x(i,1)+btpcm(10)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 1201 j = 1, l
        do 1201 i = 1, imax
          f(i,j)=h44194m2*(apcm(10)*x(i,j)+bpcm(10)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(10)*x(i,j)+btpcm(10)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
1201  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1215 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1215  continue
      do 1219 k=2,l
      do 1217 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1217  continue
1219  continue
c---tt is the cloud-free cts transmission function
CSJP      do 1221 i=1,imax*l
CSJP      fac1(i,1)=acomb(10)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(10)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
CSJP     &           betacm(10)*totvo2(i,2)*sko2d))*co2sp2(i,2)
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 1221 j = 1, l
        do 1221 i = 1, imax
          fac1(i,j)=acomb(10)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(10)*topphi(i,j))
          tt(i,j)=exp(hm1ez*(fac1(i,j)/sqrt(one+fac2(i,j))+
     &            betacm(10)*totvo2(i,j+1)*sko2d))*co2sp2(i,j+1)
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
1221  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 1253 i=1,imax*l
CSJP      exctsn(i,1,10)=radcon*delp(i,1)*sorc(i,1,10)*
CSJP     &               (ctmp(i,2)-ctmp(i,1))
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,10) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,10)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 1253 j = 1, l
        do 1253 i = 1, imax
          exctsn(i,j,10)=radcon*delp(i,j)*sorc(i,j,10)*
     &                   (ctmp(i,j+1)-ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,10)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,10)*(ctmp2(i,j+1)-ctmp2(i,j))
1253  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1261 i=1,imax 
      fctsg(i,10)=tt(i,l)*sorc(i,l,10)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,10)-sorc(i,l,10))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,10)
      fctsg(i,10)=cldfac(i,lp1,1)*fctsg(i,10)
      gxcts(i)=gxcts(i)+fctsg(i,10)
1261  continue
c 
c-----calculation for band 11 (800-900 cm-1)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 1301 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(11)*x(i,1)+bpcm(11)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(11)*x(i,1)+btpcm(11)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 1301 j = 1, l
        do 1301 i = 1, imax
          f(i,j)=h44194m2*(apcm(11)*x(i,j)+bpcm(11)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(11)*x(i,j)+btpcm(11)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
1301  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1315 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1315  continue
      do 1319 k=2,l
      do 1317 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1317  continue
1319  continue
c---tt is the cloud-free cts transmission function
CSJP      do 1321 i=1,imax*l
CSJP      fac1(i,1)=acomb(11)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(11)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
CSJP     &           betacm(11)*totvo2(i,2)*sko2d))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 1321 j = 1, l
        do 1321 i = 1, imax
          fac1(i,j)=acomb(11)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(11)*topphi(i,j))
          tt(i,j)=exp(hm1ez*(fac1(i,j)/sqrt(one+fac2(i,j))+
     &            betacm(11)*totvo2(i,j+1)*sko2d))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
1321  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 1353 i=1,imax*l
CSJP      exctsn(i,1,11)=radcon*delp(i,1)*sorc(i,1,11)*
CSJP     &               (ctmp(i,2)-ctmp(i,1))
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,11) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,11)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 1353 j = 1, l
        do 1353 i = 1, imax
          exctsn(i,j,11)=radcon*delp(i,j)*sorc(i,j,11)*
     &                   (ctmp(i,j+1)-ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,11)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,11)*(ctmp2(i,j+1)-ctmp2(i,j))
1353  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1361 i=1,imax 
      fctsg(i,11)=tt(i,l)*sorc(i,l,11)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,11)-sorc(i,l,11))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,11)
      fctsg(i,11)=cldfac(i,lp1,1)*fctsg(i,11)
      gxcts(i)=gxcts(i)+fctsg(i,11)
1361  continue
c 
c-----calculation for band 12 (900-990 cm-1)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 1401 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(12)*x(i,1)+bpcm(12)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(12)*x(i,1)+btpcm(12)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 1401 j = 1, l
        do 1401 i = 1, imax
          f(i,j)=h44194m2*(apcm(12)*x(i,j)+bpcm(12)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(12)*x(i,j)+btpcm(12)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
1401  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1415 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1415  continue
      do 1419 k=2,l
      do 1417 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1417  continue
1419  continue
c---tt is the cloud-free cts transmission function
CSJP      do 1421 i=1,imax*l
CSJP      fac1(i,1)=acomb(12)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(12)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
CSJP     &           betacm(12)*totvo2(i,2)*sko2d))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 1421 j = 1, l
        do 1421 i = 1, imax
          fac1(i,j)=acomb(12)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(12)*topphi(i,j))
          tt(i,j)=exp(hm1ez*(fac1(i,j)/sqrt(one+fac2(i,j))+
     &            betacm(12)*totvo2(i,j+1)*sko2d))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
1421  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 1453 i=1,imax*l
CSJP      exctsn(i,1,12)=radcon*delp(i,1)*sorc(i,1,12)*
CSJP     &               (ctmp(i,2)-ctmp(i,1))
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,12) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,12)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 1453 j = 1, l
        do 1453 i = 1, imax
          exctsn(i,j,12)=radcon*delp(i,j)*sorc(i,j,12)*
     &                   (ctmp(i,j+1)-ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,12)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,12)*(ctmp2(i,j+1)-ctmp2(i,j))
1453  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1461 i=1,imax 
      fctsg(i,12)=tt(i,l)*sorc(i,l,12)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,12)-sorc(i,l,12))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,12)
      fctsg(i,12)=cldfac(i,lp1,1)*fctsg(i,12)
      gxcts(i)=gxcts(i)+fctsg(i,12)
1461  continue
c 
c-----calculation for band 13 (990-1070 cm-1; includes o3))
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 1501 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(13)*x(i,1)+bpcm(13)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(13)*x(i,1)+btpcm(13)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 1501 j = 1, l
        do 1501 i = 1, imax
          f(i,j)=h44194m2*(apcm(13)*x(i,j)+bpcm(13)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(13)*x(i,j)+btpcm(13)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
1501  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1515 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1515  continue
      do 1519 k=2,l
      do 1517 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1517  continue
1519  continue
c---tt is the cloud-free cts transmission function
CSJP      do 1521 i=1,imax*l
CSJP      fac1(i,1)=acomb(13)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(13)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
CSJP     &           betacm(13)*totvo2(i,2)*sko2d +to3spc(i,1)))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 1521 j = 1, l
        do 1521 i = 1, imax
          fac1(i,j)=acomb(13)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(13)*topphi(i,j))
          tt(i,j)=exp(hm1ez*(fac1(i,j)/sqrt(one+fac2(i,j))+
     &            betacm(13)*totvo2(i,j+1)*sko2d +to3spc(i,j)))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
1521  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 1553 i=1,imax*l
CSJP      exctsn(i,1,13)=radcon*delp(i,1)*sorc(i,1,13)*
CSJP     &               (ctmp(i,2)-ctmp(i,1))
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,13) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,13)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 1553 j = 1, l
        do 1553 i = 1, imax
          exctsn(i,j,13)=radcon*delp(i,j)*sorc(i,j,13)*
     &                   (ctmp(i,j+1)-ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,13)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,13)*(ctmp2(i,j+1)-ctmp2(i,j))
1553  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1561 i=1,imax 
      fctsg(i,13)=tt(i,l)*sorc(i,l,13)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,13)-sorc(i,l,13))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,13)
      fctsg(i,13)=cldfac(i,lp1,1)*fctsg(i,13)
      gxcts(i)=gxcts(i)+fctsg(i,13)
1561  continue
c 
c-----calculation for band 14 (1070-1200 cm-1)
c 
c
c---obtain temperature correction (capphi,cappsi),then multiply
c   by optical path (var1,var2) to compute temperature-corrected
c   optical path and mean pressure for a layer (phitmp,psitmp)
CSJP      do 1601 i=1,imax*l
CSJP      f(i,1)=h44194m2*(apcm(14)*x(i,1)+bpcm(14)*y(i,1)) 
CSJP      ff(i,1)=h44194m2*(atpcm(14)*x(i,1)+btpcm(14)*y(i,1))
CSJP      ag(i,1)=(h1p41819+f(i,1))*f(i,1)+one
CSJP      agg(i,1)=(h1p41819+ff(i,1))*ff(i,1)+one 
CSJP      phitmp(i,1)=var1(i,1)*(((( ag(i,1)*ag(i,1))**2)**2)**2)
CSJP      psitmp(i,1)=var2(i,1)*(((( agg(i,1)*agg(i,1))**2)**2)**2)
      do 1601 j = 1, l
        do 1601 i = 1, imax
          f(i,j)=h44194m2*(apcm(14)*x(i,j)+bpcm(14)*y(i,j))
          ff(i,j)=h44194m2*(atpcm(14)*x(i,j)+btpcm(14)*y(i,j))
          ag(i,j)=(h1p41819+f(i,j))*f(i,j)+one
          agg(i,j)=(h1p41819+ff(i,j))*ff(i,j)+one
          phitmp(i,j)=var1(i,j)*(((( ag(i,j)*ag(i,j))**2)**2)**2)
          psitmp(i,j)=var2(i,j)*(((( agg(i,j)*agg(i,j))**2)**2)**2)
1601  continue
c---obtain optical path,mean pressure from the top to the pressure
c   p(k) (topm,topphi)
      do 1615 i=1,imax 
      topm(i,1)=phitmp(i,1) 
      topphi(i,1)=psitmp(i,1) 
1615  continue
      do 1619 k=2,l
      do 1617 i=1,imax 
      topm(i,k)=topm(i,k-1)+phitmp(i,k) 
      topphi(i,k)=topphi(i,k-1)+psitmp(i,k) 
1617  continue
1619  continue
c---tt is the cloud-free cts transmission function
CSJP      do 1621 i=1,imax*l
CSJP      fac1(i,1)=acomb(14)*topm(i,1)
CSJP      fac2(i,1)=fac1(i,1)*topm(i,1)/(bcomb(14)*topphi(i,1))
CSJP      tt(i,1)=exp(hm1ez*(fac1(i,1)/sqrt(one+fac2(i,1))+
CSJP     &           betacm(14)*totvo2(i,2)*sko2d))
CSJP      ctmp(i,2)=tt(i,1)*cldfac(i,2,1) 
CSJP      ctmp2(i,2)=tt(i,1)
      do 1621 j = 1, l
        do 1621 i = 1, imax
          fac1(i,j)=acomb(14)*topm(i,j)
          fac2(i,j)=fac1(i,j)*topm(i,j)/(bcomb(14)*topphi(i,j))
          tt(i,j)=exp(hm1ez*(fac1(i,j)/sqrt(one+fac2(i,j))+
     &            betacm(14)*totvo2(i,j+1)*sko2d))
          ctmp(i,j+1)=tt(i,j)*cldfac(i,j+1,1)
          ctmp2(i,j+1)=tt(i,j)
1621  continue
c---excts is the cts cooling rate (exctsn is the rate for band nl)
CSJP      do 1653 i=1,imax*l
CSJP      exctsn(i,1,14)=radcon*delp(i,1)*sorc(i,1,14)*
CSJP     &               (ctmp(i,2)-ctmp(i,1))
CSJP      excts(i,1)=excts(i,1)+exctsn(i,1,14) 
CSJP      exctsclr(i,1)=exctsclr(i,1)+
CSJP     & radcon*delp(i,1)*sorc(i,1,14)*(ctmp2(i,2)-ctmp2(i,1)) 
      do 1653 j = 1, l
        do 1653 i = 1, imax
          exctsn(i,j,14)=radcon*delp(i,j)*sorc(i,j,14)*
     &                   (ctmp(i,j+1)-ctmp(i,j))
          excts(i,j)=excts(i,j)+exctsn(i,j,14)
          exctsclr(i,j)=exctsclr(i,j)+
     &      radcon*delp(i,j)*sorc(i,j,14)*(ctmp2(i,j+1)-ctmp2(i,j))
1653  continue
c---gxcts is the exact cts top flux (fctsg is the flux for band nl)
      do 1661 i=1,imax 
      fctsg(i,14)=tt(i,l)*sorc(i,l,14)+
     &   (haf*delp(i,l)*(tt(i,lm1)*(p(i,lp1)-press(i,l)) + 
     &   tt(i,l)*(p(i,lp1)+press(i,l)-two*p(i,l)))) *
     &   (sorc(i,lp1,14)-sorc(i,l,14))
      gxctsclr(i)=gxctsclr(i)+fctsg(i,14)
      fctsg(i,14)=cldfac(i,lp1,1)*fctsg(i,14)
      gxcts(i)=gxcts(i)+fctsg(i,14)
1661  continue
c 
c---this is the end of the exact cts computations; at this point
c   excts has its appropriate value. the exact cts flux at the 
c   ground is now evaluated by using the (now available) cts cooling 
c   rate and the exact cts top flux.
c*** compute approximate cts heating rates for 15um and 9.6 um bands
c     (ctso3) 
CSJP      do 1711 i=1,imax*l
CSJP      ctmp2(i,2)=co2sp(i,2)*cldfac(i,2,1)
CSJP      ctmp3(i,2)=to3(i,2,1)*cldfac(i,2,1)
      do 1711 j = 1, l
        do 1711 i = 1, imax
          ctmp2(i,j+1)=co2sp(i,j+1)*cldfac(i,j+1,1)
          ctmp3(i,j+1)=to3(i,j+1,1)*cldfac(i,j+1,1)
1711  continue
      do 1701 k=1,l
      do 1701 i=1,imax 
      gxcts(i)=gxcts(i)-excts(i,k)*delp2(i,k)*radcon1
      gxctsclr(i)=gxctsclr(i)-exctsclr(i,k)*delp2(i,k)*radcon1
c     if ( i.eq.1 ) print*, ' ctso3 ',
c    &     radcon*delp(i,k)*csour(i,k)*(ctmp2(i,k+1)-ctmp2(i,k)) ,
c    &     radcon*delp(i,k)*osour(i,k)*(ctmp3(i,k+1)-ctmp3(i,k)) 
      ctso3(i,k)=radcon*delp(i,k)*
     &     (csour(i,k)*(ctmp2(i,k+1)-ctmp2(i,k)) +
     &      osour(i,k)*(ctmp3(i,k+1)-ctmp3(i,k)))
      ctso3clr(i,k)=radcon*delp(i,k)*
     &     (csour(i,k)*(co2sp(i,k+1)-co2sp(i,k)) +
     &      osour(i,k)*(to3(i,k+1,1)-to3(i,k,1)))
1701  continue
      return
      end 
