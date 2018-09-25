c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c Loops re-rolled so that code passes array-bounds checking.
c SJP 2003/03/04
c
c Parallel compiler directives replaced with equivalent OpenMP instructions.
c SJP 2001/12/13
c
c Changes to add machine type 'ALPH', which makes use of the super-fast FFTW
c FFT library.
c SJP 2001/11/22
c
c $Log: ptogcray.f,v $
c Revision 1.16  1999/05/20 06:23:52  rot032
c HBG changes to V5-2
c
c Revision 1.15  1998/12/10  00:55:39  ldr
c HBG changes to V5-1-21
c
c Revision 1.14  1997/12/17  23:22:49  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.13  1997/10/03  05:32:08  ldr
c Changes for NEC from MRD and HBG
c
c Revision 1.12  1996/10/24  01:27:40  ldr
c Merge of TIE changes with LDR/MRD changes since V5-0-12.
c
c Revision 1.11  1996/08/13  07:30:05  ldr
c Correct length of /workoe/ to avoid warning on Cray.
c
c Revision 1.10.1.1  1996/10/24  01:03:09  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.10  1996/06/13  02:07:40  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.9  1996/03/21  03:18:59  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.8  1995/08/14  07:03:25  ldr
c Tidy up indices in declaration of rotated arrays.
c
c Revision 1.7  1995/07/03  05:48:41  ldr
c Mickles' changes to V4-5-30l for Cray multiprocessing, tidied by LDR.
c
c Revision 1.6  93/12/17  15:33:32  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.5  93/10/15  14:17:12  ldr
c Changes to create special 'FUJI' version of model with fast FFTs and
c Legendre transform.
c 
c Revision 1.4  93/10/05  13:07:07  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.3  93/07/06  16:21:29  ldr
c      Added small extension array to common block elmuv
c      (for overshoot in ptogcray). Added code at end of uvharm.f to set zero
c 
c Revision 1.2  92/12/09  14:44:13  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.1  92/04/22  14:45:50  ldr
c Initial revision
c 
c Revision 1.4  92/03/19  11:46:45  ldr
c Corrected general (commented out) version to use plmx rather than plm
c where necessary.
c 
c Revision 1.3  91/03/13  12:59:45  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.2  91/03/01  17:25:29  hbg
c Miklethwaite/HBG super fast loop unrolled version. Please include
c V3-0-1h for all of dtog,dynm,gauleg,mfftm,phys,ptog,zerogi
c 
c Revision 1.1  91/02/22  16:37:53  ldr
c Initial release V3-0
c 
*VOCL TOTAL,REPEAT(999999)
      subroutine ptogcray

!$OMP THREADPRIVATE ( /LEGND/ )
!$OMP THREADPRIVATE ( /WORK1/ )
!$OMP THREADPRIVATE ( /WORKG/ )

C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      common/legnd/plmx(lw1,mw),plm(lw,mw),cplm(lw,mw),pad(2,mw)
      include 'WORK1.f'
      complex uf,vf,pf,tef,fuf,fvf
      common/workg/uf(mw,2,nl),vf(mw,2,nl),pf(mw,2),tef(mw,2,nl)
     & ,fuf(mw,2,nl),fvf(mw,2,nl) 

C Global data blocks
      include 'DIFM.f'
      include 'FLDRI.f'
      common/elmuv/ulmr(lw1,mw,nl),ulmi(lw1,mw,nl)
     & ,vlmr(lw1,mw,nl),vlmi(lw1,mw,nl)
     & ,elm4x(4)

C Local work arrays and variables
      real psumr(lw,mw),tsumr(lw,mw)
      real psumi(lw,mw),tsumi(lw,mw)
      real usumr(lw1,mw),vsumr(lw1,mw),fusmr(lw1,mw),fvsmr(lw1,mw)
      real usumi(lw1,mw),vsumi(lw1,mw),fusmi(lw1,mw),fvsmi(lw1,mw)

C Local data, functions etc

C Start code : ----------------------------------------------------------

c**** Now general for R21/R42 (LDR 4/92)
c**** PTOG code specifically for R21 model - loop unrolled (HBG Dec90)
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

C     CALL TIMER('PTOG ADD',3)

      do 81 mm=1,mw
        do 81 ll=1,lw-2
CSJP      mm=1
CSJP      do 81 ll=1,lw*mw
          psumr(ll,mm)=prr(ll,mm)*plm(ll,mm)+prr(ll+2,mm)*plm(ll+2,mm)
  81      psumi(ll,mm)=pri(ll,mm)*plm(ll,mm)+pri(ll+2,mm)*plm(ll+2,mm)
      if(lw.eq.22)then
        do mm = 1, mw
          por=psumr(1,mm)+psumr(5,mm)+psumr(9,mm)+psumr(13,mm)+
     &		psumr(17,mm)+prr(21,mm)*plm(21,mm)
          per=psumr(2,mm)+psumr(6,mm)+psumr(10,mm)+psumr(14,mm)+
     &         psumr(18,mm)+prr(22,mm)*plm(22,mm)
          poi=psumi(1,mm)+psumi(5,mm)+psumi(9,mm)+psumi(13,mm)+
     &         psumi(17,mm)+pri(21,mm)*plm(21,mm)
          pei=psumi(2,mm)+psumi(6,mm)+psumi(10,mm)+psumi(14,mm)+
     &		psumi(18,mm)+pri(22,mm)*plm(22,mm)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          pf(mm,1)=cmplx(por+per,poi+pei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD
          pf(mm,2)=cmplx(por-per,poi-pei)
        enddo
      elseif(lw.eq.43)then
        do mm = 1, mw
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
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          pf(mm,1)=cmplx(por+per,poi+pei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD
          pf(mm,2)=cmplx(por-per,poi-pei)
        enddo
      endif
        
      do 200 k=1,nl
        do 110 mm=1,mw
          do 110 ll=1,lw-2
CSJP      mm=1
CSJP      do 110 ll=1,lw*mw
	    tsumr(ll,mm)=ter(ll,mm,k)*plm(ll,mm)+
     &		ter(ll+2,mm,k)*plm(ll+2,mm)
            tsumi(ll,mm)=tei(ll,mm,k)*plm(ll,mm)+
     &		tei(ll+2,mm,k)*plm(ll+2,mm)
  110   continue

        if(lw.eq.22)then
          do mm = 1, mw
            tor=tsumr(1,mm)+tsumr(5,mm)+tsumr(9,mm)+tsumr(13,mm)+
     &           tsumr(17,mm)+ter(21,mm,k)*plm(21,mm)
            tre=tsumr(2,mm)+tsumr(6,mm)+tsumr(10,mm)+tsumr(14,mm)+
     &           tsumr(18,mm)+ter(22,mm,k)*plm(22,mm)
            toi=tsumi(1,mm)+tsumi(5,mm)+tsumi(9,mm)+tsumi(13,mm)+
     &           tsumi(17,mm)+tei(21,mm,k)*plm(21,mm)
            tie=tsumi(2,mm)+tsumi(6,mm)+tsumi(10,mm)+tsumi(14,mm)+
     &           tsumi(18,mm)+tei(22,mm,k)*plm(22,mm)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
            tef(mm,1,k)=cmplx(tor+tre,toi+tie)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD
            tef(mm,2,k)=cmplx(tor-tre,toi-tie)
          enddo
        elseif(lw.eq.43)then
          do mm = 1, mw
            tor=tsumr(1,mm)+tsumr(5,mm)+tsumr(9,mm)+tsumr(13,mm)+
     &           tsumr(17,mm)+tsumr(21,mm)+tsumr(25,mm)+tsumr(29,mm)+
     &           tsumr(33,mm)+tsumr(37,mm)+tsumr(41,mm)
            toi=tsumi(1,mm)+tsumi(5,mm)+tsumi(9,mm)+tsumi(13,mm)+
     &           tsumi(17,mm)+tsumi(21,mm)+tsumi(25,mm)+tsumi(29,mm)+
     &           tsumi(33,mm)+tsumi(37,mm)+tsumi(41,mm)
            tre=tsumr(2,mm)+tsumr(6,mm)+tsumr(10,mm)+tsumr(14,mm)+
     &           tsumr(18,mm)+tsumr(22,mm)+tsumr(26,mm)+tsumr(30,mm)+
     &           tsumr(34,mm)+tsumr(38,mm)+ter(42,mm,k)*plm(42,mm)
            tie=tsumi(2,mm)+tsumi(6,mm)+tsumi(10,mm)+tsumi(14,mm)+
     &           tsumi(18,mm)+tsumi(22,mm)+tsumi(26,mm)+tsumi(30,mm)+
     &           tsumi(34,mm)+tsumi(38,mm)+tei(42,mm,k)*plm(42,mm)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
            tef(mm,1,k)=cmplx(tor+tre,toi+tie)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD
            tef(mm,2,k)=cmplx(tor-tre,toi-tie)
          enddo
        endif

        do 135 mm = 1, mw
          do 135 ll=1,lw1-2
CSJP      mm=1
CSJP      do 135 ll=1,lw1*mw
            usumr(ll,mm)=ulmr(ll,mm,k)*plmx(ll,mm)+
     &		ulmr(ll+2,mm,k)*plmx(ll+2,mm)
            vsumr(ll,mm)=vlmr(ll,mm,k)*plmx(ll,mm)+
     &		vlmr(ll+2,mm,k)*plmx(ll+2,mm)
            fusmr(ll,mm)=frur(ll,mm,k)*plmx(ll,mm)+
     &		frur(ll+2,mm,k)*plmx(ll+2,mm)
            fvsmr(ll,mm)=frvr(ll,mm,k)*plmx(ll,mm)+
     &		frvr(ll+2,mm,k)*plmx(ll+2,mm)
            usumi(ll,mm)=ulmi(ll,mm,k)*plmx(ll,mm)+
     &		ulmi(ll+2,mm,k)*plmx(ll+2,mm)
            vsumi(ll,mm)=vlmi(ll,mm,k)*plmx(ll,mm)+
     &		vlmi(ll+2,mm,k)*plmx(ll+2,mm)
            fusmi(ll,mm)=frui(ll,mm,k)*plmx(ll,mm)+
     &		frui(ll+2,mm,k)*plmx(ll+2,mm)
            fvsmi(ll,mm)=frvi(ll,mm,k)*plmx(ll,mm)+
     &		frvi(ll+2,mm,k)*plmx(ll+2,mm)
  135   continue

        if(lw.eq.22)then
          do  mm = 1, mw
            uor=usumr(1,mm)+usumr(5,mm)+usumr(9,mm)+usumr(13,mm)+
     &           usumr(17,mm)+usumr(21,mm)
            uer=usumr(2,mm)+usumr(6,mm)+usumr(10,mm)+usumr(14,mm)+
     &           usumr(18,mm)+ulmr(22,mm,k)*plmx(22,mm)
            vor=vsumr(1,mm)+vsumr(5,mm)+vsumr(9,mm)+vsumr(13,mm)+
     &           vsumr(17,mm)+vsumr(21,mm)
            ver=vsumr(2,mm)+vsumr(6,mm)+vsumr(10,mm)+vsumr(14,mm)+
     &           vsumr(18,mm)+vlmr(22,mm,k)*plmx(22,mm)
            fruor=fusmr(1,mm)+fusmr(5,mm)+fusmr(9,mm)+fusmr(13,mm)+
     &           fusmr(17,mm)+fusmr(21,mm)
            fruer=fusmr(2,mm)+fusmr(6,mm)+fusmr(10,mm)+fusmr(14,mm)+
     &	 	 fusmr(18,mm)+frur(22,mm,k)*plmx(22,mm)
            frvor=fvsmr(1,mm)+fvsmr(5,mm)+fvsmr(9,mm)+fvsmr(13,mm)+
     &           fvsmr(17,mm)+fvsmr(21,mm)
            frver=fvsmr(2,mm)+fvsmr(6,mm)+fvsmr(10,mm)+fvsmr(14,mm)+
     &           fvsmr(18,mm)+frvr(22,mm,k)*plmx(22,mm)
            uoi=usumi(1,mm)+usumi(5,mm)+usumi(9,mm)+usumi(13,mm)+
     &           usumi(17,mm)+usumi(21,mm)
            uei=usumi(2,mm)+usumi(6,mm)+usumi(10,mm)+usumi(14,mm)+
     &           usumi(18,mm)+ulmi(22,mm,k)*plmx(22,mm)
            voi=vsumi(1,mm)+vsumi(5,mm)+vsumi(9,mm)+vsumi(13,mm)+
     &           vsumi(17,mm)+vsumi(21,mm)
            vei=vsumi(2,mm)+vsumi(6,mm)+vsumi(10,mm)+vsumi(14,mm)+
     &           vsumi(18,mm)+vlmi(22,mm,k)*plmx(22,mm)
            fruoi=fusmi(1,mm)+fusmi(5,mm)+fusmi(9,mm)+
     &           fusmi(13,mm)+fusmi(17,mm)+fusmi(21,mm)
            fruei=fusmi(2,mm)+fusmi(6,mm)+fusmi(10,mm)+
     &           fusmi(14,mm)+fusmi(18,mm)+frui(22,mm,k)*plmx(22,mm)
            frvoi=fvsmi(1,mm)+fvsmi(5,mm)+fvsmi(9,mm)+
     &           fvsmi(13,mm)+fvsmi(17,mm)+fvsmi(21,mm)
            frvei=fvsmi(2,mm)+fvsmi(6,mm)+fvsmi(10,mm)+
     &           fvsmi(14,mm)+fvsmi(18,mm)+frvi(22,mm,k)*plmx(22,mm)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
            uf(mm,1,k)=cmplx(uor+uer,uoi+uei)
            vf(mm,1,k)=cmplx(vor+ver,voi+vei)
            fuf(mm,1,k)=cmplx(fruor+fruer,fruoi+fruei)
            fvf(mm,1,k)=cmplx(frvor+frver,frvoi+frvei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD
            uf(mm,2,k)=cmplx(uor-uer,uoi-uei)
            vf(mm,2,k)=cmplx(vor-ver,voi-vei)
            fuf(mm,2,k)=cmplx(fruor-fruer,fruoi-fruei)
            fvf(mm,2,k)=cmplx(frvor-frver,frvoi-frvei)
          enddo
        elseif(lw.eq.43)then
          do mm = 1, mw
            uor=usumr(1,mm)+usumr(5,mm)+usumr(9,mm)+usumr(13,mm)+
     &           usumr(17,mm)+usumr(21,mm)+usumr(25,mm)+usumr(29,mm)+
     &           usumr(33,mm)+usumr(37,mm)+usumr(41,mm)
            uoi=usumi(1,mm)+usumi(5,mm)+usumi(9,mm)+usumi(13,mm)+
     &           usumi(17,mm)+usumi(21,mm)+usumi(25,mm)+usumi(29,mm)+
     &           usumi(33,mm)+usumi(37,mm)+usumi(41,mm)
            uer=usumr(2,mm)+usumr(6,mm)+usumr(10,mm)+usumr(14,mm)+
     &           usumr(18,mm)+usumr(22,mm)+usumr(26,mm)+usumr(30,mm)+
     &           usumr(34,mm)+usumr(38,mm)+usumr(42,mm)
            uei=usumi(2,mm)+usumi(6,mm)+usumi(10,mm)+usumi(14,mm)+
     &           usumi(18,mm)+usumi(22,mm)+usumi(26,mm)+usumi(30,mm)+
     &           usumi(34,mm)+usumi(38,mm)+usumi(42,mm)
            vor=vsumr(1,mm)+vsumr(5,mm)+vsumr(9,mm)+vsumr(13,mm)+
     &           vsumr(17,mm)+vsumr(21,mm)+vsumr(25,mm)+vsumr(29,mm)+
     &           vsumr(33,mm)+vsumr(37,mm)+vsumr(41,mm)
            voi=vsumi(1,mm)+vsumi(5,mm)+vsumi(9,mm)+vsumi(13,mm)+
     &           vsumi(17,mm)+vsumi(21,mm)+vsumi(25,mm)+vsumi(29,mm)+
     &           vsumi(33,mm)+vsumi(37,mm)+vsumi(41,mm)
            ver=vsumr(2,mm)+vsumr(6,mm)+vsumr(10,mm)+vsumr(14,mm)+
     &           vsumr(18,mm)+vsumr(22,mm)+vsumr(26,mm)+vsumr(30,mm)+
     &           vsumr(34,mm)+vsumr(38,mm)+vsumr(42,mm)
            vei=vsumi(2,mm)+vsumi(6,mm)+vsumi(10,mm)+vsumi(14,mm)+
     &           vsumi(18,mm)+vsumi(22,mm)+vsumi(26,mm)+vsumi(30,mm)+
     &           vsumi(34,mm)+vsumi(38,mm)+vsumi(42,mm)
            fruor=fusmr(1,mm)+fusmr(5,mm)+fusmr(9,mm)+fusmr(13,mm)+
     &           fusmr(17,mm)+fusmr(21,mm)+fusmr(25,mm)+fusmr(29,mm)+
     &           fusmr(33,mm)+fusmr(37,mm)+fusmr(41,mm)
            fruoi=fusmi(1,mm)+fusmi(5,mm)+fusmi(9,mm)+fusmi(13,mm)+
     &           fusmi(17,mm)+fusmi(21,mm)+fusmi(25,mm)+fusmi(29,mm)+
     &           fusmi(33,mm)+fusmi(37,mm)+fusmi(41,mm)
            fruer=fusmr(2,mm)+fusmr(6,mm)+fusmr(10,mm)+fusmr(14,mm)+
     &           fusmr(18,mm)+fusmr(22,mm)+fusmr(26,mm)+fusmr(30,mm)+
     &           fusmr(34,mm)+fusmr(38,mm)+fusmr(42,mm)
            fruei=fusmi(2,mm)+fusmi(6,mm)+fusmi(10,mm)+fusmi(14,mm)+
     &           fusmi(18,mm)+fusmi(22,mm)+fusmi(26,mm)+fusmi(30,mm)+
     &           fusmi(34,mm)+fusmi(38,mm)+fusmi(42,mm)
            frvor=fvsmr(1,mm)+fvsmr(5,mm)+fvsmr(9,mm)+fvsmr(13,mm)+
     &           fvsmr(17,mm)+fvsmr(21,mm)+fvsmr(25,mm)+fvsmr(29,mm)+
     &           fvsmr(33,mm)+fvsmr(37,mm)+fvsmr(41,mm)
            frvoi=fvsmi(1,mm)+fvsmi(5,mm)+fvsmi(9,mm)+fvsmi(13,mm)+
     &           fvsmi(17,mm)+fvsmi(21,mm)+fvsmi(25,mm)+fvsmi(29,mm)+
     &           fvsmi(33,mm)+fvsmi(37,mm)+fvsmi(41,mm)
            frver=fvsmr(2,mm)+fvsmr(6,mm)+fvsmr(10,mm)+fvsmr(14,mm)+
     &           fvsmr(18,mm)+fvsmr(22,mm)+fvsmr(26,mm)+fvsmr(30,mm)+
     &           fvsmr(34,mm)+fvsmr(38,mm)+fvsmr(42,mm)
            frvei=fvsmi(2,mm)+fvsmi(6,mm)+fvsmi(10,mm)+fvsmi(14,mm)+
     &           fvsmi(18,mm)+fvsmi(22,mm)+fvsmi(26,mm)+fvsmi(30,mm)+
     &           fvsmi(34,mm)+fvsmi(38,mm)+fvsmi(42,mm)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
            uf(mm,1,k)=cmplx(uor+uer,uoi+uei)
            vf(mm,1,k)=cmplx(vor+ver,voi+vei)
            fuf(mm,1,k)=cmplx(fruor+fruer,fruoi+fruei)
            fvf(mm,1,k)=cmplx(frvor+frver,frvoi+frvei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD
            uf(mm,2,k)=cmplx(uor-uer,uoi-uei)
            vf(mm,2,k)=cmplx(vor-ver,voi-vei)
            fuf(mm,2,k)=cmplx(fruor-fruer,fruoi-fruei)
            fvf(mm,2,k)=cmplx(frvor-frver,frvoi-frvei)
          enddo
        endif

  200 continue


C     CALL TIMER('PTOG ADD',4)

C**   THE FOURIER ARRAYS UF() ETC HAVE NOW BEEN FILLED WITH THE 
C**   FOURIER COMPONENTS FOR A PARTICULAR LATITUDE.
C**   TRANSFORM THE FOURIER COMPONENTS ONTO THE LONGITUDE GRID
C**    ,THE RESULTS BEING PLACED IN UN() ETC.

C     CALL TIMER('PTOG MFF',3)

      nex=(5*nl+1)*2

CSJP  Former machine dependence at this point
      call mfftga(un, uf, nex)

C     CALL TIMER('PTOG MFF',4)
C**     UN() HOLDS P*.U.COS(LAT)/RAD
C**     VN() HOLDS P*.V.COS(LAT)/RAD
      return
      end
C----
C---- The next part is for T63 form of ptogcray
C----
      subroutine ptogcr63

!$OMP THREADPRIVATE ( /LEGND/ )
!$OMP THREADPRIVATE ( /WORK1/ )
!$OMP THREADPRIVATE ( /WORKG/ )

C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)
      common/legnd/plmx(lw1,mw),plm(lw,mw),cplm(lw,mw),pad(2,mw)
      include 'WORK1.f'
      complex uf,vf,pf,tef,fuf,fvf
      common/workg/uf(mw,2,nl),vf(mw,2,nl),pf(mw,2),tef(mw,2,nl)
     & ,fuf(mw,2,nl),fvf(mw,2,nl) 

C Global data blocks
      include 'DIFM.f'
      include 'FLDRI.f'
      common/elmuv/ulmr(lw1,mw,nl),ulmi(lw1,mw,nl)
     & ,vlmr(lw1,mw,nl),vlmi(lw1,mw,nl)
     & ,elm4x(4)

C Local work arrays and variables
c --- Note following arrays at dimension mw+1 to prevent "bank conflicts"
      real psumr(mw+1,lw),tsumr(mw+1,lw)
      real psumi(mw+1,lw),tsumi(mw+1,lw)
      real usumr(mw+1,lw),vsumr(mw+1,lw),fusmr(mw+1,lw),fvsmr(mw+1,lw)
      real usumi(mw+1,lw),vsumi(mw+1,lw),fusmi(mw+1,lw),fvsmi(mw+1,lw)

C Local data, functions etc

C Start code : ----------------------------------------------------------

c**** PTOG code specifically for R/T63 model - loop unrolled (HBG Dec90)
c**** recoding by Micklethwaite Jan,1991.
c**** more recoding by HBG Feb,1991.
c---- This code will only work at R63,T63 : If MW is not equal to
c---- 63 + 1 = 64 then this code will not work. Do not change next line.
      if(mw.ne.64)then
        print *,'Wrong resolution for fast PTOG R63,T63 code'
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

C     CALL TIMER('PTOG ADD',3)

C
C-- Compute surface pressure
C
      do 81 mm=1,mw
       do 81 ll=1,lw-2
          psumr(mm,ll)=prr(ll,mm)*plm(ll,mm)+prr(ll+2,mm)*plm(ll+2,mm)
          psumi(mm,ll)=pri(ll,mm)*plm(ll,mm)+pri(ll+2,mm)*plm(ll+2,mm)
  81  continue

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
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
        pf(mm,1)=cmplx(por+per,poi+pei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD
        pf(mm,2)=cmplx(por-per,poi-pei)
   85   continue

C
C-- Compute fields on model levels (T,U,V,Fu,Fv)
C
      do 200 k=1,nl

       do 110 mm=1,mw
         do 110 ll=1,lw-2
	    tsumr(mm,ll)=ter(ll,mm,k)*plm(ll,mm)+
     &		ter(ll+2,mm,k)*plm(ll+2,mm)
            tsumi(mm,ll)=tei(ll,mm,k)*plm(ll,mm)+
     &		tei(ll+2,mm,k)*plm(ll+2,mm)
  110   continue

c---- Do not change the following loop unrolled coding
        do 120 mm = 1, mw
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
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          tef(mm,1,k)=cmplx(tor+tre,toi+tie)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD
          tef(mm,2,k)=cmplx(tor-tre,toi-tie)
  120   continue
   
       do 135 mm = 1, mw
         do 135 ll=1,lw1-2
            usumr(mm,ll)=ulmr(ll,mm,k)*plmx(ll,mm)+
     &		ulmr(ll+2,mm,k)*plmx(ll+2,mm)
            vsumr(mm,ll)=vlmr(ll,mm,k)*plmx(ll,mm)+
     &		vlmr(ll+2,mm,k)*plmx(ll+2,mm)
            fusmr(mm,ll)=frur(ll,mm,k)*plmx(ll,mm)+
     &		frur(ll+2,mm,k)*plmx(ll+2,mm)
            fvsmr(mm,ll)=frvr(ll,mm,k)*plmx(ll,mm)+
     &		frvr(ll+2,mm,k)*plmx(ll+2,mm)
            usumi(mm,ll)=ulmi(ll,mm,k)*plmx(ll,mm)+
     &		ulmi(ll+2,mm,k)*plmx(ll+2,mm)
            vsumi(mm,ll)=vlmi(ll,mm,k)*plmx(ll,mm)+
     &		vlmi(ll+2,mm,k)*plmx(ll+2,mm)
            fusmi(mm,ll)=frui(ll,mm,k)*plmx(ll,mm)+
     &		frui(ll+2,mm,k)*plmx(ll+2,mm)
            fvsmi(mm,ll)=frvi(ll,mm,k)*plmx(ll,mm)+
     &		frvi(ll+2,mm,k)*plmx(ll+2,mm)
  135   continue

c---- Do not change the following loop unrolled coding
        do 140 mm = 1, mw
      uor=usumr(mm, 1)+usumr(mm, 5)+usumr(mm, 9)+usumr(mm,13)+
     &    usumr(mm,17)+usumr(mm,21)+usumr(mm,25)+usumr(mm,29)+
     &    usumr(mm,33)+usumr(mm,37)+usumr(mm,41)+usumr(mm,45)+
     &    usumr(mm,49)+usumr(mm,53)+usumr(mm,57)+usumr(mm,61)
     &    +ulmr(65,mm,k)*plmx(65,mm)
      uoi=usumi(mm, 1)+usumi(mm, 5)+usumi(mm, 9)+usumi(mm,13)+
     &    usumi(mm,17)+usumi(mm,21)+usumi(mm,25)+usumi(mm,29)+
     &    usumi(mm,33)+usumi(mm,37)+usumi(mm,41)+usumi(mm,45)+
     &    usumi(mm,49)+usumi(mm,53)+usumi(mm,57)+usumi(mm,61)
     &    +ulmi(65,mm,k)*plmx(65,mm)
      uer=usumr(mm, 2)+usumr(mm, 6)+usumr(mm,10)+usumr(mm,14)+
     &    usumr(mm,18)+usumr(mm,22)+usumr(mm,26)+usumr(mm,30)+
     &    usumr(mm,34)+usumr(mm,38)+usumr(mm,42)+usumr(mm,46)+
     &    usumr(mm,50)+usumr(mm,54)+usumr(mm,58)+usumr(mm,62)
      uei=usumi(mm, 2)+usumi(mm, 6)+usumi(mm,10)+usumi(mm,14)+
     &    usumi(mm,18)+usumi(mm,22)+usumi(mm,26)+usumi(mm,30)+
     &    usumi(mm,34)+usumi(mm,38)+usumi(mm,42)+usumi(mm,46)+
     &    usumi(mm,50)+usumi(mm,54)+usumi(mm,58)+usumi(mm,62)
      vor=vsumr(mm, 1)+vsumr(mm, 5)+vsumr(mm, 9)+vsumr(mm,13)+
     &    vsumr(mm,17)+vsumr(mm,21)+vsumr(mm,25)+vsumr(mm,29)+
     &    vsumr(mm,33)+vsumr(mm,37)+vsumr(mm,41)+vsumr(mm,45)+
     &    vsumr(mm,49)+vsumr(mm,53)+vsumr(mm,57)+vsumr(mm,61)
     &    +vlmr(65,mm,k)*plmx(65,mm)
      voi=vsumi(mm, 1)+vsumi(mm, 5)+vsumi(mm, 9)+vsumi(mm,13)+
     &    vsumi(mm,17)+vsumi(mm,21)+vsumi(mm,25)+vsumi(mm,29)+
     &    vsumi(mm,33)+vsumi(mm,37)+vsumi(mm,41)+vsumi(mm,45)+
     &    vsumi(mm,49)+vsumi(mm,53)+vsumi(mm,57)+vsumi(mm,61)
     &    +vlmi(65,mm,k)*plmx(65,mm)
      ver=vsumr(mm, 2)+vsumr(mm, 6)+vsumr(mm,10)+vsumr(mm,14)+
     &    vsumr(mm,18)+vsumr(mm,22)+vsumr(mm,26)+vsumr(mm,30)+
     &    vsumr(mm,34)+vsumr(mm,38)+vsumr(mm,42)+vsumr(mm,46)+
     &    vsumr(mm,50)+vsumr(mm,54)+vsumr(mm,58)+vsumr(mm,62)
      vei=vsumi(mm, 2)+vsumi(mm, 6)+vsumi(mm,10)+vsumi(mm,14)+
     &    vsumi(mm,18)+vsumi(mm,22)+vsumi(mm,26)+vsumi(mm,30)+
     &    vsumi(mm,34)+vsumi(mm,38)+vsumi(mm,42)+vsumi(mm,46)+
     &    vsumi(mm,50)+vsumi(mm,54)+vsumi(mm,58)+vsumi(mm,62)
      fruor=fusmr(mm, 1)+fusmr(mm, 5)+fusmr(mm, 9)+fusmr(mm,13)+
     &    fusmr(mm,17)+fusmr(mm,21)+fusmr(mm,25)+fusmr(mm,29)+
     &    fusmr(mm,33)+fusmr(mm,37)+fusmr(mm,41)+fusmr(mm,45)+
     &    fusmr(mm,49)+fusmr(mm,53)+fusmr(mm,57)+fusmr(mm,61)
     &    +frur(65,mm,k)*plmx(65,mm)
      fruoi=fusmi(mm, 1)+fusmi(mm, 5)+fusmi(mm, 9)+fusmi(mm,13)+
     &    fusmi(mm,17)+fusmi(mm,21)+fusmi(mm,25)+fusmi(mm,29)+
     &    fusmi(mm,33)+fusmi(mm,37)+fusmi(mm,41)+fusmi(mm,45)+
     &    fusmi(mm,49)+fusmi(mm,53)+fusmi(mm,57)+fusmi(mm,61)
     &    +frui(65,mm,k)*plmx(65,mm)
      fruer=fusmr(mm, 2)+fusmr(mm, 6)+fusmr(mm,10)+fusmr(mm,14)+
     &    fusmr(mm,18)+fusmr(mm,22)+fusmr(mm,26)+fusmr(mm,30)+
     &    fusmr(mm,34)+fusmr(mm,38)+fusmr(mm,42)+fusmr(mm,46)+
     &    fusmr(mm,50)+fusmr(mm,54)+fusmr(mm,58)+fusmr(mm,62)
      fruei=fusmi(mm, 2)+fusmi(mm, 6)+fusmi(mm,10)+fusmi(mm,14)+
     &    fusmi(mm,18)+fusmi(mm,22)+fusmi(mm,26)+fusmi(mm,30)+
     &    fusmi(mm,34)+fusmi(mm,38)+fusmi(mm,42)+fusmi(mm,46)+
     &    fusmi(mm,50)+fusmi(mm,54)+fusmi(mm,58)+fusmi(mm,62)
      frvor=fvsmr(mm, 1)+fvsmr(mm, 5)+fvsmr(mm, 9)+fvsmr(mm,13)+
     &    fvsmr(mm,17)+fvsmr(mm,21)+fvsmr(mm,25)+fvsmr(mm,29)+
     &    fvsmr(mm,33)+fvsmr(mm,37)+fvsmr(mm,41)+fvsmr(mm,45)+
     &    fvsmr(mm,49)+fvsmr(mm,53)+fvsmr(mm,57)+fvsmr(mm,61)
     &    +frvr(65,mm,k)*plmx(65,mm)
      frvoi=fvsmi(mm, 1)+fvsmi(mm, 5)+fvsmi(mm, 9)+fvsmi(mm,13)+
     &    fvsmi(mm,17)+fvsmi(mm,21)+fvsmi(mm,25)+fvsmi(mm,29)+
     &    fvsmi(mm,33)+fvsmi(mm,37)+fvsmi(mm,41)+fvsmi(mm,45)+
     &    fvsmi(mm,49)+fvsmi(mm,53)+fvsmi(mm,57)+fvsmi(mm,61)
     &    +frvi(65,mm,k)*plmx(65,mm)
      frver=fvsmr(mm, 2)+fvsmr(mm, 6)+fvsmr(mm,10)+fvsmr(mm,14)+
     &    fvsmr(mm,18)+fvsmr(mm,22)+fvsmr(mm,26)+fvsmr(mm,30)+
     &    fvsmr(mm,34)+fvsmr(mm,38)+fvsmr(mm,42)+fvsmr(mm,46)+
     &    fvsmr(mm,50)+fvsmr(mm,54)+fvsmr(mm,58)+fvsmr(mm,62)
      frvei=fvsmi(mm, 2)+fvsmi(mm, 6)+fvsmi(mm,10)+fvsmi(mm,14)+
     &    fvsmi(mm,18)+fvsmi(mm,22)+fvsmi(mm,26)+fvsmi(mm,30)+
     &    fvsmi(mm,34)+fvsmi(mm,38)+fvsmi(mm,42)+fvsmi(mm,46)+
     &    fvsmi(mm,50)+fvsmi(mm,54)+fvsmi(mm,58)+fvsmi(mm,62)
c     NORTH LAT BAND - ADD THE EVEN AND ODD COMPONENTS
          uf(mm,1,k)=cmplx(uor+uer,uoi+uei)
          vf(mm,1,k)=cmplx(vor+ver,voi+vei)
          fuf(mm,1,k)=cmplx(fruor+fruer,fruoi+fruei)
          fvf(mm,1,k)=cmplx(frvor+frver,frvoi+frvei)
c     SOUTH LAT BAND - SUBTRACT EVEN FROM ODD
          uf(mm,2,k)=cmplx(uor-uer,uoi-uei)
          vf(mm,2,k)=cmplx(vor-ver,voi-vei)
          fuf(mm,2,k)=cmplx(fruor-fruer,fruoi-fruei)
          fvf(mm,2,k)=cmplx(frvor-frver,frvoi-frvei)
  140   continue
  200 continue


C     CALL TIMER('PTOG ADD',4)

C**   THE FOURIER ARRAYS UF() ETC HAVE NOW BEEN FILLED WITH THE 
C**   FOURIER COMPONENTS FOR A PARTICULAR LATITUDE.
C**   TRANSFORM THE FOURIER COMPONENTS ONTO THE LONGITUDE GRID
C**    ,THE RESULTS BEING PLACED IN UN() ETC.

C     CALL TIMER('PTOG MFF',3)

      nex=(5*nl+1)*2

CSJP  Former machine dependence at this point
      call mfftga(un, uf, nex)

C     CALL TIMER('PTOG MFF',4)
C**     UN() HOLDS P*.U.COS(LAT)/RAD
C**     VN() HOLDS P*.V.COS(LAT)/RAD
      return
      end 
