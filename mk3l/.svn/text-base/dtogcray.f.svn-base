c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Loops re-rolled so that code passes array-bounds checking.
c SJP 2003/03/08
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/12
c
c Changes to add machine type 'ALPH', which makes use of the super-fast FFTW
c FFT library.
c SJP 2001/11/22
c
c $Log: dtogcray.f,v $
c Revision 1.17  1999/05/20 06:23:52  rot032
c HBG changes to V5-2
c
c Revision 1.16  1997/12/17  23:22:49  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.15  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.14  1996/10/24  01:02:38  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.13  1996/06/13  02:06:11  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.12  1996/03/21  03:18:34  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.11  1995/08/14  07:03:25  ldr
c Tidy up indices in declaration of rotated arrays.
c
c Revision 1.10  1995/07/05  06:14:32  ldr
c Changes from Dave Micklethwaite for Cray parallel execution tidied by LDR
c and merged into V4-7-2l.
c
c Revision 1.8.1.1  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.9  1995/05/04  04:06:11  ldr
c Put plmf2 into common block with padding to avoid -ve indefinite error
c on Cray with -ei option. Need to add to -Xlocaldata in makefile.
c
c Revision 1.8  1993/12/17  15:32:06  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.7  93/10/15  14:16:59  ldr
c Changes to create special 'FUJI' version of model with fast FFTs and
c Legendre transform.
c 
c Revision 1.6  93/10/05  13:05:53  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.5  93/07/14  14:51:04  ldr
c  ECMWF implicit treatment of vorticity eqn is now an option (HBG).
c 
c Revision 1.4  92/12/09  14:43:13  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.3  92/12/07  15:25:40  ldr
c HBG speedups.
c 
c Revision 1.2  92/08/06  16:29:37  ldr
c Turn off moisture stuff if nsemilag.ne.0.
c 
c Revision 1.1  92/04/16  16:44:48  ldr
c Initial revision
c 
c Revision 1.4  92/03/19  11:38:48  ldr
c Added 2 lines which were omitted from general (commented out) version.
c 
c Revision 1.3  91/03/13  12:57:36  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.2  91/03/01  17:16:21  hbg
c Miklethwaite/HBG super fast loop unrolled version
c 
c Revision 1.1  91/02/22  16:37:09  ldr
c Initial release V3-0
c 
*VOCL TOTAL,REPEAT(999999)
      subroutine dtogcray

!$OMP THREADPRIVATE ( /LEGND/ )
!$OMP THREADPRIVATE ( /LEGNDF/ )
!$OMP THREADPRIVATE ( /UBARVO/ )
!$OMP THREADPRIVATE ( /WORK1X/ )
!$OMP THREADPRIVATE ( /WORKF/ )

C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      common/legnd/plmx(lw1,mw),plm(lw,mw),cplm(lw,mw),pad(2,mw)
      common/legndf/plmf2(lw,mw),padf(2,mw)
      common/ubarvo/ubx(nl,2),dubx(nl,2)
      common/work1x/un(ln2,nl),vn(ln2,nl),pn(ln2),dvn(ln2,nl)
     & ,von(ln2,nl),ten(ln2,nl),cpn(ln2),pln(ln2),rmn(ln2,nl)
      complex dvf,vof,tef,cpf,plf,vofm
      common/workf/dvf(mw,2,nl),vof(mw,2,nl),tef(mw,2,nl)
     & ,cpf(mw,2),plf(mw,2),vofm(mw,2,nl)

C Global data blocks
      include 'CNSTE.f'
      include 'FEWFLAGS.f'
      include 'FLDRI.f'
      include 'FLDMRI.f'
      common/ubplm/utr(nl,lw),ur(lw,nl)

C Local work arrays and variables
      real psumr(lw,mw),cpsumr(lw,mw)
      real psumi(lw,mw),cpsumi(lw,mw)
      real tsumr(lw,mw),vosumr(lw,mw),dvsumr(lw,mw)
      real tsumi(lw,mw),vosumi(lw,mw),dvsumi(lw,mw)
      real vmsumr(lw,mw),vmsumi(lw,mw)

C Local data, functions etc

C Start code : ----------------------------------------------------------

c**** Now general for R21/R42 (LDR 4/92)      
c**** DTOG code specifically for R21 model - loop unrolled (HBG Dec90)
c**** recoding by Micklethwaite Jan,1991.
c**** more recoding by HBG Feb,1991.

C**** NOTE THE COMPLEX ARITHMETIC HAS BEEN SPLIT INTO REAL & IMAGINARY
C**** PARTS.THE SUMMATIONS ARE CARRIED OUT SEPARATELY OVER ODD AND EVEN
C**** L COMPONENTS.
C**** TRANSFORM THE SPECTRAL COMPONENTS TO FOURIER COMPONENTS AT A
C**** PARTICULAR LATITUDE.

C**** ADD THE ODD AND EVEN COMPONENTS SEPARATELY
C**   ADD ODD COMPONENTS IN COMPUTER NOTATION (EVEN SPECTRAL L VALUES)
C**   AND EVEN COMPONENTS (ODD SPECTRAL L VALUES)

C     CALL TIMER('DTOG ADD',3)

      do 80 mm=1,mw
        do 80 ll=1,lw-2
CSJP      mm=1
CSJP      do 80 ll=1,lw*mw
          psumr(ll,mm)=prr(ll,mm)*plm(ll,mm)+prr(ll+2,mm)*plm(ll+2,mm)
          cpsumr(ll,mm)=prr(ll,mm)*cplm(ll,mm)+prr(ll+2,mm)*
     &    cplm(ll+2,mm)
          psumi(ll,mm)=pri(ll,mm)*plm(ll,mm)+pri(ll+2,mm)*plm(ll+2,mm)
          cpsumi(ll,mm)=pri(ll,mm)*cplm(ll,mm)+pri(ll+2,mm)*
     &    cplm(ll+2,mm)
   80 continue

      if(lw.eq.22)then
        do mm = 1, mw
          por=psumr(1,mm)+psumr(5,mm)+psumr(9,mm)+psumr(13,mm)+
     &         psumr(17,mm)+prr(21,mm)*plm(21,mm)
          per=psumr(2,mm)+psumr(6,mm)+psumr(10,mm)+psumr(14,mm)+
     &         psumr(18,mm)+prr(22,mm)*plm(22,mm)
          cpor=cpsumr(1,mm)+cpsumr(5,mm)+cpsumr(9,mm)+cpsumr(13,mm)+
     &         cpsumr(17,mm)+prr(21,mm)*cplm(21,mm)
          cper=cpsumr(2,mm)+cpsumr(6,mm)+cpsumr(10,mm)+cpsumr(14,mm)+
     &         cpsumr(18,mm)+prr(22,mm)*cplm(22,mm)
          poi=psumi(1,mm)+psumi(5,mm)+psumi(9,mm)+psumi(13,mm)+
     &         psumi(17,mm)+pri(21,mm)*plm(21,mm)
          pei=psumi(2,mm)+psumi(6,mm)+psumi(10,mm)+psumi(14,mm)+
     &         psumi(18,mm)+pri(22,mm)*plm(22,mm)
          cpoi=cpsumi(1,mm)+cpsumi(5,mm)+cpsumi(9,mm)+cpsumi(13,mm)+
     &         cpsumi(17,mm)+pri(21,mm)*cplm(21,mm)
          cpei=cpsumi(2,mm)+cpsumi(6,mm)+cpsumi(10,mm)+cpsumi(14,mm)+
     &         cpsumi(18,mm)+pri(22,mm)*cplm(22,mm)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
        plf(mm,1)=cmplx(-(mm-1)*poi,(mm-1)*por)
     &           +cmplx(-(mm-1)*pei,(mm-1)*per)
        cpf(mm,1)=cmplx(cpor+cper,cpoi+cpei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD(REVERSE COS.D(PLM) TERM)
        plf(mm,2)=cmplx(-(mm-1)*poi,(mm-1)*por)
     &           -cmplx(-(mm-1)*pei,(mm-1)*per)
        cpf(mm,2)=cmplx(cper-cpor,cpei-cpoi)
        enddo

      elseif(lw.eq.43)then
        do  mm = 1, mw
          por=psumr(1,mm)+psumr(5,mm)+psumr(9,mm)+psumr(13,mm)+
     &         psumr(17,mm)+psumr(21,mm)+psumr(25,mm)+psumr(29,mm)+
     &         psumr(33,mm)+psumr(37,mm)+psumr(41,mm)
          poi=psumi(1,mm)+psumi(5,mm)+psumi(9,mm)+psumi(13,mm)+
     &         psumi(17,mm)+psumi(21,mm)+psumi(25,mm)+psumi(29,mm)+
     &         psumi(33,mm)+psumi(37,mm)+psumi(41,mm)
          per=psumr(2,mm)+psumr(6,mm)+psumr(10,mm)+psumr(14,mm)+
     &         psumr(18,mm)+psumr(22,mm)+psumr(26,mm)+psumr(30,mm)+
     &         psumr(34,mm)+psumr(38,mm)+prr(42,mm)*plm(42,mm)
          pei=psumi(2,mm)+psumi(6,mm)+psumi(10,mm)+psumi(14,mm)+
     &         psumi(18,mm)+psumi(22,mm)+psumi(26,mm)+psumi(30,mm)+
     &         psumi(34,mm)+psumi(38,mm)+pri(42,mm)*plm(42,mm)
          cpor=cpsumr(1,mm)+cpsumr(5,mm)+cpsumr(9,mm)+cpsumr(13,mm)+
     &         cpsumr(17,mm)+cpsumr(21,mm)+cpsumr(25,mm)+cpsumr(29,mm)+
     &         cpsumr(33,mm)+cpsumr(37,mm)+cpsumr(41,mm)
          cpoi=cpsumi(1,mm)+cpsumi(5,mm)+cpsumi(9,mm)+cpsumi(13,mm)+
     &         cpsumi(17,mm)+cpsumi(21,mm)+cpsumi(25,mm)+cpsumi(29,mm)+
     &         cpsumi(33,mm)+cpsumi(37,mm)+cpsumi(41,mm)
          cper=cpsumr(2,mm)+cpsumr(6,mm)+cpsumr(10,mm)+cpsumr(14,mm)+
     &         cpsumr(18,mm)+cpsumr(22,mm)+cpsumr(26,mm)+cpsumr(30,mm)+
     &         cpsumr(34,mm)+cpsumr(38,mm)+prr(42,mm)*cplm(42,mm)
          cpei=cpsumi(2,mm)+cpsumi(6,mm)+cpsumi(10,mm)+cpsumi(14,mm)+
     &         cpsumi(18,mm)+cpsumi(22,mm)+cpsumi(26,mm)+cpsumi(30,mm)+
     &         cpsumi(34,mm)+cpsumi(38,mm)+pri(42,mm)*cplm(42,mm)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
        plf(mm,1)=cmplx(-(mm-1)*poi,(mm-1)*por)
     &           +cmplx(-(mm-1)*pei,(mm-1)*per)
        cpf(mm,1)=cmplx(cpor+cper,cpoi+cpei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD(REVERSE COS.D(PLM) TERM)
        plf(mm,2)=cmplx(-(mm-1)*poi,(mm-1)*por)
     &           -cmplx(-(mm-1)*pei,(mm-1)*per)
        cpf(mm,2)=cmplx(cper-cpor,cpei-cpoi)
        enddo
      endif

      do 300 mm=1,mw
        do 300 ll=1,lw
          plmf2(ll,mm) = plm(ll,mm)*flm2(ll,mm)
 300  continue

      do 125 k=1,nl
        do 100 mm=1,mw
          do 100 ll=1,lw-2
CSJP      mm=1
CSJP      do 100 ll=1,lw*mw
            tsumr(ll,mm)=ter(ll,mm,k)*plm(ll,mm)+
     &			ter(ll+2,mm,k)*plm(ll+2,mm)
            vosumr(ll,mm)=psir(ll,mm,k)*plmf2(ll,mm)+
     &			psir(ll+2,mm,k)*plmf2(ll+2,mm)
            dvsumr(ll,mm)=xhir(ll,mm,k)*plmf2(ll,mm)+
     &			xhir(ll+2,mm,k)*plmf2(ll+2,mm)
          tsumi(ll,mm)=tei(ll,mm,k)*plm(ll,mm)+
     &			tei(ll+2,mm,k)*plm(ll+2,mm)
          vosumi(ll,mm)=psii(ll,mm,k)*plmf2(ll,mm)+
     &			psii(ll+2,mm,k)*plmf2(ll+2,mm)
          dvsumi(ll,mm)=xhii(ll,mm,k)*plmf2(ll,mm)+
     &			xhii(ll+2,mm,k)*plmf2(ll+2,mm)
  100 continue

      if(impvor)then
CSJP      mm=1
CSJP      do 102 ll=1,lw*mw
        do 102 mm = 1, mw
          do 102 ll = 1, lw-2
            vmsumr(ll,mm)=psimr(ll,mm,k)*plmf2(ll,mm)+
     &                    psimr(ll+2,mm,k)*plmf2(ll+2,mm)
            vmsumi(ll,mm)=psimi(ll,mm,k)*plmf2(ll,mm)+
     &                    psimi(ll+2,mm,k)*plmf2(ll+2,mm)
  102 continue
          uor=0.0
          ucor=0.0
          do 106 ll=1,lw,2
            uor=uor+ur(ll,k)*plm(ll,1)
  106       ucor=ucor+ur(ll,k)*cplm(ll,1)
          uer=0.0
          ucer=0.0
          do 108 ll=2,lw,2
            uer=uer+ur(ll,k)*plm(ll,1)
  108       ucer=ucer+ur(ll,k)*cplm(ll,1)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          ubx(k,1)=uor+uer
          dubx(k,1)=ucor+ucer
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD(REVERSE COS.D(PLM) TERM)
          ubx(k,2)=uor-uer
          dubx(k,2)=ucer-ucor
      end if

      if(lw.eq.22)then
        do mm = 1, mw
          tor=tsumr(1,mm)+tsumr(5,mm)+tsumr(9,mm)+tsumr(13,mm)+
     &         tsumr(17,mm)+ter(21,mm,k)*plm(21,mm)
          tre=tsumr(2,mm)+tsumr(6,mm)+tsumr(10,mm)+tsumr(14,mm)+
     &         tsumr(18,mm)+ter(22,mm,k)*plm(22,mm)
          voor=vosumr(1,mm)+vosumr(5,mm)+vosumr(9,mm)+vosumr(13,mm)+
     &         vosumr(17,mm)+psir(21,mm,k)*plmf2(21,mm)
          voer=vosumr(2,mm)+vosumr(6,mm)+vosumr(10,mm)+vosumr(14,mm)+
     &         vosumr(18,mm)+psir(22,mm,k)*plmf2(22,mm)
          dvor=dvsumr(1,mm)+dvsumr(5,mm)+dvsumr(9,mm)+dvsumr(13,mm)+
     &         dvsumr(17,mm)+xhir(21,mm,k)*plmf2(21,mm)
          dver=dvsumr(2,mm)+dvsumr(6,mm)+dvsumr(10,mm)+dvsumr(14,mm)+
     &         dvsumr(18,mm)+xhir(22,mm,k)*plmf2(22,mm)
          toi=tsumi(1,mm)+tsumi(5,mm)+tsumi(9,mm)+tsumi(13,mm)+
     &         tsumi(17,mm)+tei(21,mm,k)*plm(21,mm)
          tie=tsumi(2,mm)+tsumi(6,mm)+tsumi(10,mm)+tsumi(14,mm)+
     &         tsumi(18,mm)+tei(22,mm,k)*plm(22,mm)
          vooi=vosumi(1,mm)+vosumi(5,mm)+vosumi(9,mm)+vosumi(13,mm)+
     &         vosumi(17,mm)+psii(21,mm,k)*plmf2(21,mm)
          voei=vosumi(2,mm)+vosumi(6,mm)+vosumi(10,mm)+vosumi(14,mm)+
     &         vosumi(18,mm)+psii(22,mm,k)*plmf2(22,mm)
          dvoi=dvsumi(1,mm)+dvsumi(5,mm)+dvsumi(9,mm)+dvsumi(13,mm)+
     &         dvsumi(17,mm)+xhii(21,mm,k)*plmf2(21,mm)
          dvei=dvsumi(2,mm)+dvsumi(6,mm)+dvsumi(10,mm)+dvsumi(14,mm)+
     &         dvsumi(18,mm)+xhii(22,mm,k)*plmf2(22,mm)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          tef(mm,1,k)=cmplx(tor+tre,toi+tie)
          vof(mm,1,k)=cmplx(-voor,-vooi)+cmplx(-voer,-voei)
          dvf(mm,1,k)=cmplx(-dvor,-dvoi)+cmplx(-dver,-dvei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD(REVERSE COS.D(PLM) TERM)
          tef(mm,2,k)=cmplx(tor-tre,toi-tie)
          vof(mm,2,k)=cmplx(-voor,-vooi)-cmplx(-voer,-voei)
          dvf(mm,2,k)=cmplx(-dvor,-dvoi)-cmplx(-dver,-dvei)
        enddo

      if(impvor)then
        do mm = 1, mw
          voor=vmsumr(1,mm)+vmsumr(5,mm)+vmsumr(9,mm)+vmsumr(13,mm)+
     &         vmsumr(17,mm)+psimr(21,mm,k)*plmf2(21,mm)
          voer=vmsumr(2,mm)+vmsumr(6,mm)+vmsumr(10,mm)+vmsumr(14,mm)+
     &         vmsumr(18,mm)+psimr(22,mm,k)*plmf2(22,mm)
          vooi=vmsumi(1,mm)+vmsumi(5,mm)+vmsumi(9,mm)+vmsumi(13,mm)+
     &         vmsumi(17,mm)+psimi(21,mm,k)*plmf2(21,mm)
          voei=vmsumi(2,mm)+vmsumi(6,mm)+vmsumi(10,mm)+vmsumi(14,mm)+
     &         vmsumi(18,mm)+psimi(22,mm,k)*plmf2(22,mm)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          vofm(mm,1,k)=cmplx(-voor,-vooi)+cmplx(-voer,-voei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD(REVERSE COS.D(PLM) TERM)
          vofm(mm,2,k)=cmplx(-voor,-vooi)-cmplx(-voer,-voei)
        enddo
      endif

      elseif(lw.eq.43)then
        do  mm = 1, mw
          tor=tsumr(1,mm)+tsumr(5,mm)+tsumr(9,mm)+tsumr(13,mm)+
     &         tsumr(17,mm)+tsumr(21,mm)+tsumr(25,mm)+tsumr(29,mm)+
     &         tsumr(33,mm)+tsumr(37,mm)+tsumr(41,mm)
          toi=tsumi(1,mm)+tsumi(5,mm)+tsumi(9,mm)+tsumi(13,mm)+
     &         tsumi(17,mm)+tsumi(21,mm)+tsumi(25,mm)+tsumi(29,mm)+
     &         tsumi(33,mm)+tsumi(37,mm)+tsumi(41,mm)
          tre=tsumr(2,mm)+tsumr(6,mm)+tsumr(10,mm)+tsumr(14,mm)+
     &         tsumr(18,mm)+tsumr(22,mm)+tsumr(26,mm)+tsumr(30,mm)+
     &         tsumr(34,mm)+tsumr(38,mm)+ter(42,mm,k)*plm(42,mm)
          tie=tsumi(2,mm)+tsumi(6,mm)+tsumi(10,mm)+tsumi(14,mm)+
     &         tsumi(18,mm)+tsumi(22,mm)+tsumi(26,mm)+tsumi(30,mm)+
     &         tsumi(34,mm)+tsumi(38,mm)+tei(42,mm,k)*plm(42,mm)
          voor=vosumr(1,mm)+vosumr(5,mm)+vosumr(9,mm)+vosumr(13,mm)+
     &         vosumr(17,mm)+vosumr(21,mm)+vosumr(25,mm)+vosumr(29,mm)+
     &         vosumr(33,mm)+vosumr(37,mm)+vosumr(41,mm)
          vooi=vosumi(1,mm)+vosumi(5,mm)+vosumi(9,mm)+vosumi(13,mm)+
     &         vosumi(17,mm)+vosumi(21,mm)+vosumi(25,mm)+vosumi(29,mm)+
     &         vosumi(33,mm)+vosumi(37,mm)+vosumi(41,mm)
          voer=vosumr(2,mm)+vosumr(6,mm)+vosumr(10,mm)+vosumr(14,mm)+
     &         vosumr(18,mm)+vosumr(22,mm)+vosumr(26,mm)+vosumr(30,mm)+
     &         vosumr(34,mm)+vosumr(38,mm)+psir(42,mm,k)*plmf2(42,mm)
          voei=vosumi(2,mm)+vosumi(6,mm)+vosumi(10,mm)+vosumi(14,mm)+
     &         vosumi(18,mm)+vosumi(22,mm)+vosumi(26,mm)+vosumi(30,mm)+
     &         vosumi(34,mm)+vosumi(38,mm)+psii(42,mm,k)*plmf2(42,mm)
          dvor=dvsumr(1,mm)+dvsumr(5,mm)+dvsumr(9,mm)+dvsumr(13,mm)+
     &         dvsumr(17,mm)+dvsumr(21,mm)+dvsumr(25,mm)+dvsumr(29,mm)+
     &         dvsumr(33,mm)+dvsumr(37,mm)+dvsumr(41,mm)
          dvoi=dvsumi(1,mm)+dvsumi(5,mm)+dvsumi(9,mm)+dvsumi(13,mm)+
     &         dvsumi(17,mm)+dvsumi(21,mm)+dvsumi(25,mm)+dvsumi(29,mm)+
     &         dvsumi(33,mm)+dvsumi(37,mm)+dvsumi(41,mm)
          dver=dvsumr(2,mm)+dvsumr(6,mm)+dvsumr(10,mm)+dvsumr(14,mm)+
     &         dvsumr(18,mm)+dvsumr(22,mm)+dvsumr(26,mm)+dvsumr(30,mm)+
     &         dvsumr(34,mm)+dvsumr(38,mm)+xhir(42,mm,k)*plmf2(42,mm)
          dvei=dvsumi(2,mm)+dvsumi(6,mm)+dvsumi(10,mm)+dvsumi(14,mm)+
     &         dvsumi(18,mm)+dvsumi(22,mm)+dvsumi(26,mm)+dvsumi(30,mm)+
     &         dvsumi(34,mm)+dvsumi(38,mm)+xhii(42,mm,k)*plmf2(42,mm)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          tef(mm,1,k)=cmplx(tor+tre,toi+tie)
          vof(mm,1,k)=cmplx(-voor,-vooi)+cmplx(-voer,-voei)
          dvf(mm,1,k)=cmplx(-dvor,-dvoi)+cmplx(-dver,-dvei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD(REVERSE COS.D(PLM) TERM)
          tef(mm,2,k)=cmplx(tor-tre,toi-tie)
          vof(mm,2,k)=cmplx(-voor,-vooi)-cmplx(-voer,-voei)
          dvf(mm,2,k)=cmplx(-dvor,-dvoi)-cmplx(-dver,-dvei)
        enddo

      if(impvor)then
        do  mm = 1, mw
          voor=vmsumr(1,mm)+vmsumr(5,mm)+vmsumr(9,mm)+vmsumr(13,mm)+
     &         vmsumr(17,mm)+vmsumr(21,mm)+vmsumr(25,mm)+vmsumr(29,mm)+
     &         vmsumr(33,mm)+vmsumr(37,mm)+vmsumr(41,mm)
          vooi=vmsumi(1,mm)+vmsumi(5,mm)+vmsumi(9,mm)+vmsumi(13,mm)+
     &         vmsumi(17,mm)+vmsumi(21,mm)+vmsumi(25,mm)+vmsumi(29,mm)+
     &         vmsumi(33,mm)+vmsumi(37,mm)+vmsumi(41,mm)
          voer=vmsumr(2,mm)+vmsumr(6,mm)+vmsumr(10,mm)+vmsumr(14,mm)+
     &         vmsumr(18,mm)+vmsumr(22,mm)+vmsumr(26,mm)+vmsumr(30,mm)+
     &         vmsumr(34,mm)+vmsumr(38,mm)+psimr(42,mm,k)*plmf2(42,mm)
          voei=vmsumi(2,mm)+vmsumi(6,mm)+vmsumi(10,mm)+vmsumi(14,mm)+
     &         vmsumi(18,mm)+vmsumi(22,mm)+vmsumi(26,mm)+vmsumi(30,mm)+
     &         vmsumi(34,mm)+vmsumi(38,mm)+psimi(42,mm,k)*plmf2(42,mm)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          vofm(mm,1,k)=cmplx(-voor,-vooi)+cmplx(-voer,-voei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD(REVERSE COS.D(PLM) TERM)
          vofm(mm,2,k)=cmplx(-voor,-vooi)-cmplx(-voer,-voei)
        enddo
      endif

      endif

 125  continue


C     CALL TIMER('DTOG ADD',4)

C**   THE FOURIER ARRAYS DVF() ETC HAVE NOW BEEN FILLED WITH THE
C**   FOURIER COMPONENTS FOR A PARTICULAR LATITUDE.
C**   TRANSFORM THE FOURIER COMPONENTS ONTO THE LONGITUDE GRID
C**    ,THE RESULTS BEING PLACED IN DVN() ETC.

C     CALL TIMER('DTOG MFF',3)

      nex=(3*nl+2)*2

CSJP  Former machine dependence at this point
      call mfftga(dvn, dvf, nex)

C     CALL TIMER('DTOG MFF',4)
C**   NOTE CPN() HOLDS COS(LAT)D(P*)/D(LAT)
C**        PLN() HOLDS D(P*)/D(LONG)
      return
      end
C----
C---- The next part is for T63 form of dtogcray
C----
      subroutine dtogcr63

!$OMP THREADPRIVATE ( /LEGND/ )
!$OMP THREADPRIVATE ( /UBARVO/ )
!$OMP THREADPRIVATE ( /WORK1X/ )
!$OMP THREADPRIVATE ( /WORKF/ )

C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      common/legnd/plmx(lw1,mw),plm(lw,mw),cplm(lw,mw),pad(2,mw)
      common/ubarvo/ubx(nl,2),dubx(nl,2)
      common/work1x/un(ln2,nl),vn(ln2,nl),pn(ln2),dvn(ln2,nl)
     & ,von(ln2,nl),ten(ln2,nl),cpn(ln2),pln(ln2),rmn(ln2,nl)
      complex dvf,vof,tef,cpf,plf,vofm
      common/workf/dvf(mw,2,nl),vof(mw,2,nl),tef(mw,2,nl)
     & ,cpf(mw,2),plf(mw,2),vofm(mw,2,nl)

C Global data blocks
      include 'CNSTE.f'
      include 'FEWFLAGS.f'
      include 'FLDRI.f'
      include 'FLDMRI.f'
      common/ubplm/utr(nl,lw),ur(lw,nl)

C Local work arrays and variables
      real plmf2(lw,mw)
c --- Note following arrays at dimension mw+1 to prevent "bank conflicts"
      real psumr(mw+1,lw),cpsumr(mw+1,lw)
      real psumi(mw+1,lw),cpsumi(mw+1,lw)
      real tsumr(mw+1,lw),vosumr(mw+1,lw),dvsumr(mw+1,lw)
      real tsumi(mw+1,lw),vosumi(mw+1,lw),dvsumi(mw+1,lw)
      real vmsumr(mw+1,lw),vmsumi(mw+1,lw)

C Local data, functions etc

C Start code : ----------------------------------------------------------

c**** DTOG code specifically for R/T63 model - loop unrolled (HBG Dec90)
c**** recoding by Micklethwaite Jan,1991.
c**** more recoding by HBG Feb,1991.

c---- This code will only work at R63,T63 : If MW is not equal to
c---- 63 + 1 = 64 then this code will not work. Do not change next line.
      if(mw.ne.64)then
        print *,'Wrong resolution for fast DTOG R63,T63 code'
        stop
      end if

C**** NOTE THE COMPLEX ARITHMETIC HAS BEEN SPLIT INTO REAL & IMAGINARY
C**** PARTS.THE SUMMATIONS ARE CARRIED OUT SEPARATELY OVER ODD AND EVEN
C**** L COMPONENTS.
C**** TRANSFORM THE SPECTRAL COMPONENTS TO FOURIER COMPONENTS AT A
C**** PARTICULAR LATITUDE.

C**** ADD THE ODD AND EVEN COMPONENTS SEPARATELY
C**   ADD ODD COMPONENTS IN COMPUTER NOTATION (EVEN SPECTRAL L VALUES)
C**   AND EVEN COMPONENTS (ODD SPECTRAL L VALUES)

C     CALL TIMER('DTOG ADD',3)

C
C-- Compute gradients of surface pressure
C
      do 80 mm=1,mw
       do 80 ll=1,lw-2
          psumr(mm,ll)=prr(ll,mm)*plm(ll,mm)+prr(ll+2,mm)*plm(ll+2,mm)
          cpsumr(mm,ll)=prr(ll,mm)*cplm(ll,mm)+prr(ll+2,mm)*
     &    cplm(ll+2,mm)
          psumi(mm,ll)=pri(ll,mm)*plm(ll,mm)+pri(ll+2,mm)*plm(ll+2,mm)
          cpsumi(mm,ll)=pri(ll,mm)*cplm(ll,mm)+pri(ll+2,mm)*
     &    cplm(ll+2,mm)
   80 continue

c---- Do not change the following loop unrolled coding
      do 85 mm = 1, mw
      por=psumr(mm, 1)+psumr(mm, 5)+psumr(mm, 9)+psumr(mm,13)+
     &    psumr(mm,17)+psumr(mm,21)+psumr(mm,25)+psumr(mm,29)+
     &    psumr(mm,33)+psumr(mm,37)+psumr(mm,41)+psumr(mm,45)+
     &    psumr(mm,49)+psumr(mm,53)+psumr(mm,57)+psumr(mm,61)
      poi=psumi(mm, 1)+psumi(mm, 5)+psumi(mm, 9)+psumi(mm,13)+
     &    psumi(mm,17)+psumi(mm,21)+psumi(mm,25)+psumi(mm,29)+
     &    psumi(mm,33)+psumi(mm,37)+psumi(mm,41)+psumi(mm,45)+
     &    psumi(mm,49)+psumi(mm,53)+psumi(mm,57)+psumi(mm,61)
      per=psumr(mm, 2)+psumr(mm, 6)+psumr(mm,10)+psumr(mm,14)+
     &    psumr(mm,18)+psumr(mm,22)+psumr(mm,26)+psumr(mm,30)+
     &    psumr(mm,34)+psumr(mm,38)+psumr(mm,42)+psumr(mm,46)+
     &    psumr(mm,50)+psumr(mm,54)+psumr(mm,58)+psumr(mm,62)
      pei=psumi(mm, 2)+psumi(mm, 6)+psumi(mm,10)+psumi(mm,14)+
     &    psumi(mm,18)+psumi(mm,22)+psumi(mm,26)+psumi(mm,30)+
     &    psumi(mm,34)+psumi(mm,38)+psumi(mm,42)+psumi(mm,46)+
     &    psumi(mm,50)+psumi(mm,54)+psumi(mm,58)+psumi(mm,62)
      cpor=cpsumr(mm, 1)+cpsumr(mm, 5)+cpsumr(mm, 9)+cpsumr(mm,13)+
     &    cpsumr(mm,17)+cpsumr(mm,21)+cpsumr(mm,25)+cpsumr(mm,29)+
     &    cpsumr(mm,33)+cpsumr(mm,37)+cpsumr(mm,41)+cpsumr(mm,45)+
     &    cpsumr(mm,49)+cpsumr(mm,53)+cpsumr(mm,57)+cpsumr(mm,61)
      cpoi=cpsumi(mm, 1)+cpsumi(mm, 5)+cpsumi(mm, 9)+cpsumi(mm,13)+
     &    cpsumi(mm,17)+cpsumi(mm,21)+cpsumi(mm,25)+cpsumi(mm,29)+
     &    cpsumi(mm,33)+cpsumi(mm,37)+cpsumi(mm,41)+cpsumi(mm,45)+
     &    cpsumi(mm,49)+cpsumi(mm,53)+cpsumi(mm,57)+cpsumi(mm,61)
      cper=cpsumr(mm, 2)+cpsumr(mm, 6)+cpsumr(mm,10)+cpsumr(mm,14)+
     &    cpsumr(mm,18)+cpsumr(mm,22)+cpsumr(mm,26)+cpsumr(mm,30)+
     &    cpsumr(mm,34)+cpsumr(mm,38)+cpsumr(mm,42)+cpsumr(mm,46)+
     &    cpsumr(mm,50)+cpsumr(mm,54)+cpsumr(mm,58)+cpsumr(mm,62)
      cpei=cpsumi(mm, 2)+cpsumi(mm, 6)+cpsumi(mm,10)+cpsumi(mm,14)+
     &    cpsumi(mm,18)+cpsumi(mm,22)+cpsumi(mm,26)+cpsumi(mm,30)+
     &    cpsumi(mm,34)+cpsumi(mm,38)+cpsumi(mm,42)+cpsumi(mm,46)+
     &    cpsumi(mm,50)+cpsumi(mm,54)+cpsumi(mm,58)+cpsumi(mm,62)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
        plf(mm,1)=cmplx(-(mm-1)*poi,(mm-1)*por)
     &           +cmplx(-(mm-1)*pei,(mm-1)*per)
        cpf(mm,1)=cmplx(cpor+cper,cpoi+cpei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD(REVERSE COS.D(PLM) TERM)
        plf(mm,2)=cmplx(-(mm-1)*poi,(mm-1)*por)
     &           -cmplx(-(mm-1)*pei,(mm-1)*per)
        cpf(mm,2)=cmplx(cper-cpor,cpei-cpoi)
   85 continue


      do 300 mm=1,mw
        do 300 ll=1,lw
          plmf2(ll,mm) = plm(ll,mm)*flm2(ll,mm)
 300  continue

C
C-- Compute fields on model levels (T,Vort,Divg)
C
      do 125 k=1,nl

       do 100 mm=1,mw
         do 100 ll=1,lw-2
          tsumr(mm,ll)=ter(ll,mm,k)*plm(ll,mm)+
     &			ter(ll+2,mm,k)*plm(ll+2,mm)
          vosumr(mm,ll)=psir(ll,mm,k)*plmf2(ll,mm)+
     &			psir(ll+2,mm,k)*plmf2(ll+2,mm)
          dvsumr(mm,ll)=xhir(ll,mm,k)*plmf2(ll,mm)+
     &			xhir(ll+2,mm,k)*plmf2(ll+2,mm)
          tsumi(mm,ll)=tei(ll,mm,k)*plm(ll,mm)+
     &			tei(ll+2,mm,k)*plm(ll+2,mm)
          vosumi(mm,ll)=psii(ll,mm,k)*plmf2(ll,mm)+
     &			psii(ll+2,mm,k)*plmf2(ll+2,mm)
          dvsumi(mm,ll)=xhii(ll,mm,k)*plmf2(ll,mm)+
     &			xhii(ll+2,mm,k)*plmf2(ll+2,mm)
  100   continue

c---- Do not change the following loop unrolled coding
        do 105 mm = 1, mw
      tor=tsumr(mm, 1)+tsumr(mm, 5)+tsumr(mm, 9)+tsumr(mm,13)+
     &    tsumr(mm,17)+tsumr(mm,21)+tsumr(mm,25)+tsumr(mm,29)+
     &    tsumr(mm,33)+tsumr(mm,37)+tsumr(mm,41)+tsumr(mm,45)+
     &    tsumr(mm,49)+tsumr(mm,53)+tsumr(mm,57)+tsumr(mm,61)
      toi=tsumi(mm, 1)+tsumi(mm, 5)+tsumi(mm, 9)+tsumi(mm,13)+
     &    tsumi(mm,17)+tsumi(mm,21)+tsumi(mm,25)+tsumi(mm,29)+
     &    tsumi(mm,33)+tsumi(mm,37)+tsumi(mm,41)+tsumi(mm,45)+
     &    tsumi(mm,49)+tsumi(mm,53)+tsumi(mm,57)+tsumi(mm,61)
      tre=tsumr(mm, 2)+tsumr(mm, 6)+tsumr(mm,10)+tsumr(mm,14)+
     &    tsumr(mm,18)+tsumr(mm,22)+tsumr(mm,26)+tsumr(mm,30)+
     &    tsumr(mm,34)+tsumr(mm,38)+tsumr(mm,42)+tsumr(mm,46)+
     &    tsumr(mm,50)+tsumr(mm,54)+tsumr(mm,58)+tsumr(mm,62)
      tie=tsumi(mm, 2)+tsumi(mm, 6)+tsumi(mm,10)+tsumi(mm,14)+
     &    tsumi(mm,18)+tsumi(mm,22)+tsumi(mm,26)+tsumi(mm,30)+
     &    tsumi(mm,34)+tsumi(mm,38)+tsumi(mm,42)+tsumi(mm,46)+
     &    tsumi(mm,50)+tsumi(mm,54)+tsumi(mm,58)+tsumi(mm,62)
      voor=vosumr(mm, 1)+vosumr(mm, 5)+vosumr(mm, 9)+vosumr(mm,13)+
     &    vosumr(mm,17)+vosumr(mm,21)+vosumr(mm,25)+vosumr(mm,29)+
     &    vosumr(mm,33)+vosumr(mm,37)+vosumr(mm,41)+vosumr(mm,45)+
     &    vosumr(mm,49)+vosumr(mm,53)+vosumr(mm,57)+vosumr(mm,61)
      vooi=vosumi(mm, 1)+vosumi(mm, 5)+vosumi(mm, 9)+vosumi(mm,13)+
     &    vosumi(mm,17)+vosumi(mm,21)+vosumi(mm,25)+vosumi(mm,29)+
     &    vosumi(mm,33)+vosumi(mm,37)+vosumi(mm,41)+vosumi(mm,45)+
     &    vosumi(mm,49)+vosumi(mm,53)+vosumi(mm,57)+vosumi(mm,61)
      voer=vosumr(mm, 2)+vosumr(mm, 6)+vosumr(mm,10)+vosumr(mm,14)+
     &    vosumr(mm,18)+vosumr(mm,22)+vosumr(mm,26)+vosumr(mm,30)+
     &    vosumr(mm,34)+vosumr(mm,38)+vosumr(mm,42)+vosumr(mm,46)+
     &    vosumr(mm,50)+vosumr(mm,54)+vosumr(mm,58)+vosumr(mm,62)
      voei=vosumi(mm, 2)+vosumi(mm, 6)+vosumi(mm,10)+vosumi(mm,14)+
     &    vosumi(mm,18)+vosumi(mm,22)+vosumi(mm,26)+vosumi(mm,30)+
     &    vosumi(mm,34)+vosumi(mm,38)+vosumi(mm,42)+vosumi(mm,46)+
     &    vosumi(mm,50)+vosumi(mm,54)+vosumi(mm,58)+vosumi(mm,62)
      dvor=dvsumr(mm, 1)+dvsumr(mm, 5)+dvsumr(mm, 9)+dvsumr(mm,13)+
     &    dvsumr(mm,17)+dvsumr(mm,21)+dvsumr(mm,25)+dvsumr(mm,29)+
     &    dvsumr(mm,33)+dvsumr(mm,37)+dvsumr(mm,41)+dvsumr(mm,45)+
     &    dvsumr(mm,49)+dvsumr(mm,53)+dvsumr(mm,57)+dvsumr(mm,61)
      dvoi=dvsumi(mm, 1)+dvsumi(mm, 5)+dvsumi(mm, 9)+dvsumi(mm,13)+
     &    dvsumi(mm,17)+dvsumi(mm,21)+dvsumi(mm,25)+dvsumi(mm,29)+
     &    dvsumi(mm,33)+dvsumi(mm,37)+dvsumi(mm,41)+dvsumi(mm,45)+
     &    dvsumi(mm,49)+dvsumi(mm,53)+dvsumi(mm,57)+dvsumi(mm,61)
      dver=dvsumr(mm, 2)+dvsumr(mm, 6)+dvsumr(mm,10)+dvsumr(mm,14)+
     &    dvsumr(mm,18)+dvsumr(mm,22)+dvsumr(mm,26)+dvsumr(mm,30)+
     &    dvsumr(mm,34)+dvsumr(mm,38)+dvsumr(mm,42)+dvsumr(mm,46)+
     &    dvsumr(mm,50)+dvsumr(mm,54)+dvsumr(mm,58)+dvsumr(mm,62)
      dvei=dvsumi(mm, 2)+dvsumi(mm, 6)+dvsumi(mm,10)+dvsumi(mm,14)+
     &    dvsumi(mm,18)+dvsumi(mm,22)+dvsumi(mm,26)+dvsumi(mm,30)+
     &    dvsumi(mm,34)+dvsumi(mm,38)+dvsumi(mm,42)+dvsumi(mm,46)+
     &    dvsumi(mm,50)+dvsumi(mm,54)+dvsumi(mm,58)+dvsumi(mm,62)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          tef(mm,1,k)=cmplx(tor+tre,toi+tie)
          vof(mm,1,k)=cmplx(-voor,-vooi)+cmplx(-voer,-voei)
          dvf(mm,1,k)=cmplx(-dvor,-dvoi)+cmplx(-dver,-dvei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD
          tef(mm,2,k)=cmplx(tor-tre,toi-tie)
          vof(mm,2,k)=cmplx(-voor,-vooi)-cmplx(-voer,-voei)
          dvf(mm,2,k)=cmplx(-dvor,-dvoi)-cmplx(-dver,-dvei)
  105   continue

      if(impvor)then
C
C-- If implicit vorticity, then calculate Vort at T-1 etc
C
       do 102 mm=1,mw
         do 102 ll=1,lw-2
          vmsumr(mm,ll)=psimr(ll,mm,k)*plmf2(ll,mm)+
     &                  psimr(ll+2,mm,k)*plmf2(ll+2,mm)
          vmsumi(mm,ll)=psimi(ll,mm,k)*plmf2(ll,mm)+
     &                  psimi(ll+2,mm,k)*plmf2(ll+2,mm)
  102 continue

        do mm = 1, mw
          voor=vmsumr(mm,1)+vmsumr(mm,5)+vmsumr(mm,9)+vmsumr(mm,13)+
     &      vmsumr(mm,17)+vmsumr(mm,21)+vmsumr(mm,25)+vmsumr(mm,29)+
     &      vmsumr(mm,33)+vmsumr(mm,37)+vmsumr(mm,41)+vmsumr(mm,45)+
     &      vmsumr(mm,49)+vmsumr(mm,53)+vmsumr(mm,57)+vmsumr(mm,61)
          voer=vmsumr(mm,2)+vmsumr(mm,6)+vmsumr(mm,10)+vmsumr(mm,14)+
     &      vmsumr(mm,18)+vmsumr(mm,22)+vmsumr(mm,26)+vmsumr(mm,30)+
     &      vmsumr(mm,34)+vmsumr(mm,38)+vmsumr(mm,42)+vmsumr(mm,46)+
     &      vmsumr(mm,50)+vmsumr(mm,54)+vmsumr(mm,58)+vmsumr(mm,62)
          vooi=vmsumi(mm,1)+vmsumi(mm,5)+vmsumi(mm,9)+vmsumi(mm,13)+
     &      vmsumi(mm,17)+vmsumi(mm,21)+vmsumi(mm,25)+vmsumi(mm,29)+
     &      vmsumi(mm,33)+vmsumi(mm,37)+vmsumi(mm,41)+vmsumi(mm,45)+
     &      vmsumi(mm,49)+vmsumi(mm,53)+vmsumi(mm,57)+vmsumi(mm,61)
          voei=vmsumi(mm,2)+vmsumi(mm,6)+vmsumi(mm,10)+vmsumi(mm,14)+
     &      vmsumi(mm,18)+vmsumi(mm,22)+vmsumi(mm,26)+vmsumi(mm,30)+
     &      vmsumi(mm,34)+vmsumi(mm,38)+vmsumi(mm,42)+vmsumi(mm,46)+
     &      vmsumi(mm,50)+vmsumi(mm,54)+vmsumi(mm,58)+vmsumi(mm,62)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          vofm(mm,1,k)=cmplx(-voor,-vooi)+cmplx(-voer,-voei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD
          vofm(mm,2,k)=cmplx(-voor,-vooi)-cmplx(-voer,-voei)
        enddo

          uor=0.0
          ucor=0.0
          do 106 ll=1,lw,2
            uor=uor+ur(ll,k)*plm(ll,1)
  106       ucor=ucor+ur(ll,k)*cplm(ll,1)
          uer=0.0
          ucer=0.0
          do 108 ll=2,lw,2
            uer=uer+ur(ll,k)*plm(ll,1)
  108       ucer=ucer+ur(ll,k)*cplm(ll,1)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          ubx(k,1)=uor+uer
          dubx(k,1)=ucor+ucer
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD(REVERSE COS.D(PLM) TERM)
          ubx(k,2)=uor-uer
          dubx(k,2)=ucer-ucor

      endif

 125  continue


C     CALL TIMER('DTOG ADD',4)

C**   THE FOURIER ARRAYS DVF() ETC HAVE NOW BEEN FILLED WITH THE
C**   FOURIER COMPONENTS FOR A PARTICULAR LATITUDE.
C**   TRANSFORM THE FOURIER COMPONENTS ONTO THE LONGITUDE GRID
C**    ,THE RESULTS BEING PLACED IN DVN() ETC.

C     CALL TIMER('DTOG MFF',3)

      nex=(3*nl+2)*2

CSJP  Former machine dependence at this point
      call mfftga(dvn, dvf, nex)

C     CALL TIMER('DTOG MFF',4)
C**   NOTE CPN() HOLDS COS(LAT)D(P*)/D(LAT)
C**        PLN() HOLDS D(P*)/D(LONG)
      return
      end
