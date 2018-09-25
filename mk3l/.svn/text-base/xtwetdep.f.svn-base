c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c $Log: xtwetdep.f,v $
c Revision 1.11  2001/10/12 02:13:46  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.10  2001/06/04 02:27:04  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.9  2001/03/07 04:29:02  rot032
c Changes from LDR to bring sulfur cycle to run N63
c
c Revision 1.8  2001/02/28 04:40:57  rot032
c Bring sulfur cycle to run N58.
c
c Revision 1.7  2001/02/22 05:56:38  rot032
c Changes from LDR to sulfur cycle, to bring it to run N52
c
c Revision 1.6  2000/12/08 03:58:54  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.5  2000/11/14 06:47:45  rot032
c Changes from LDR for aerosol model
c
c Revision 1.4  2000/08/16 02:59:26  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.3  2000/06/20 07:39:54  rot032
c Remove setting of wet dep (SO4) to zero (sens. test only).
c
c Revision 1.2  2000/06/20 02:27:11  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.1  2000/02/02 00:51:29  rot032
c Initial revision
c
      SUBROUTINE XTWETDEP
     * (KTRAC,krow,
     *  PTMST, PCONS2, PDTIME,
     *  PDPP1,
     *  PMRATEP, PFPREC, PFEVAP,
     *  PCLCOVER, PRHOP1, PSOLUB, psolubc, pmlwc,
     &  pfsnow,pfsubl,pcfcover,pmiwc,pmaccr,pfmelt,pfstay,pqfsed,plambs,
     &  prscav,pfconv,pclcon,pcscav,      !Inputs
     &  PXTP10, PXTP1C, PDEP3D, conwd,pxtp1con,   !In & Out
     &  prevap)                          !Outputs

C
C   *XTWETDEP* CALCULATES THE WET DEPOSITION OF TRACE GASES OR AEROSOLS
C
C   JOHANN FEICHTER              UNI HAMBURG            08-91
C
C   PURPOSE
C  ---------
C   TO CALCULATE THE WET SCAVENGING OF GASES OR AEROSOLS IN CLOUDS
C
C   INTERFACE
C  -------------
C   THIS ROUTINE IS CALLED FROM *XTCHEM*
C   C
C  METHOD
C  -------
C
C   NO EXTERNALS
C---------------
C

C Global parameters

      include 'PARAMS.f'
      include 'PHYSPARAMS.f' !Needed for CPARAMS.f
      include 'CPARAMS.f' !Input rhos 
      include 'ECPARM.f'  !Input ntrac
      include 'FEWFLAGS.f'  

C Argument list

      real pcscav(klon,klev) !Convective tracer scavenging fraction (b/w 0 and 1)
      INTEGER KTRAC
      REAL PTMST
      REAL PCONS2
      REAL PDTIME
      REAL PXTP10(KLON,KLEV) !Tracer m.r. outside liquid-water cloud (clear air/ice cloud)
      REAL PXTP1C(KLON,KLEV) !Tracer m.r.  inside liquid-water cloud
      real pxtp1con(klon,klev) !Tracer m.r.  inside convective cloud
      REAL PDPP1(KLON,KLEV)
      REAL PMRATEP(KLON,KLEV)
      REAL PFPREC(KLON,KLEV)
      REAL PFEVAP(KLON,KLEV)
      REAL PDEP3D(KLON,KLEV)
      REAL PCLCOVER(KLON,KLEV)
      REAL PRHOP1(KLON,KLEV)
      REAL PSOLUB(KLON,KLEV)
      real psolubc(klon,klev)
      real pmlwc(klon,klev)
      real pfsnow(klon,klev)
      real pfconv(klon,klev)
      real pclcon(klon,klev)
      real pfsubl(klon,klev)
      real pcfcover(klon,klev)
      real pmiwc(klon,klev)
      real pmaccr(klon,klev)
      real pfmelt(klon,klev)
      real pfstay(klon,klev)
      real pqfsed(klon,klev)
      real plambs(klon,klev)
      real prscav(klon,klev)
      real prevap(klon,klev)
      real conwd(klon,ntrac)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
c      REAL ZDEP(KLON) !Only needed for old code
      REAL ZDEPS(KLON),ZDEPR(KLON),
     *     ZMTOF(KLON),   ZFTOM(KLON),
     *     ZCLEAR(KLON), ZCLR0(KLON)
      real zcollefr(ntrac), zcollefs(ntrac), Ecols(ntrac), Rcoeff(ntrac)
      real Evfac(ntrac)
      integer kbase(klon)

C Local data, functions etc

      pow75(x)=sqrt(x*sqrt(x))

C Start code : ----------------------------------------------------------

C    PHYSICAL CONSTANTS

      KTOP=4     !Top level for wet deposition (counting from top)
c      PQTMST=1./PTMST
      ZMIN=1.E-20

      ZCOLLEFR(:)=0.05          !Below-cloud scav eff. for rain
c      ZCOLLEFR(:)=0.1          !Below-cloud scav eff. for rain

      ZCOLLEFS(:)=0.01          !Below-cloud scav eff. for snow
      zcollefs(itracso2)=0.0
 
      Ecols(:)=0.              !In-cloud scav. eff. (snow)

      Rcoeff(:)=1.             !Retention coeff. on riming
      Rcoeff(itracso2)=0.62    !Smaller value for SO2 (Iribarne et al, 1990)
      
      Evfac(:)=0.5             !Relative reevaporation rate < 1 for aerosols
      Evfac(itracso2)=1.       !Relative reevaporation rate = 1 for SO2

      prevap(:,:)=0.

      if(debug)then
        if(krow.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, start xtwetdep, JT= ',ktrac
          write(25,9)'pdep3d ',(pdep3d(mg,k),k=nl,1,-1)
          write(25,9)'pxtp1con ',(pxtp1con(mg,k),k=nl,1,-1)
          write(25,9)'pxtp1c ',(pxtp1c(mg,k),k=nl,1,-1)
          write(25,9)'pxtp10 ',(pxtp10(mg,k),k=nl,1,-1)
        endif
      endif

C
      CALL RESETR(ZDEPS,KLON,0.)
      CALL RESETR(ZDEPR,KLON,0.)

c Search for convective cloud base

      kbase(:)=9999
      do jk=ktop,klev
        do jl=1,klon
          if(pclcon(jl,jk).gt.zmin)kbase(jl)=jk
        enddo
      enddo


C     BEGIN OF VERTICAL LOOP
        DO 160 JK=KTOP,KLEV
C
         DO 110 JL=1,KLON
          ZCLEAR(JL)=1.-PCLCOVER(JL,JK)-pcfcover(jl,jk)-pclcon(jl,jk)
          ZCLR0(JL)=1.-PCLCOVER(JL,JK)-pclcon(jl,jk) !Clear air or ice cloud (applies to pxtp10)
          ZMTOF(JL)=PDPP1(JL,JK)*PCONS2
          ZFTOM(JL)=1./ZMTOF(JL)
          PXTP1C(JL,JK)=AMAX1(0.,PXTP1C(JL,JK))
          PXTP10(JL,JK)=AMAX1(0.,PXTP10(JL,JK))
  110    CONTINUE

c In-cloud ice scavenging (including vertical redistribution when snow
c evaporates or falls into a layer). Include accretion of ql by snow.

         if(Ecols(ktrac).gt.zmin)then
           do jl=1,klon
             ziicscav=0.
             if(pmiwc(jl,jk).gt.zmin)then
              ziicscav=Ecols(ktrac)*pqfsed(jl,jk) !qfsed is the fractional sedimentation in dt
              xdep=pxtp10(jl,jk)*ziicscav*pcfcover(jl,jk)
              pdep3d(jl,jk)=pdep3d(jl,jk)+xdep
              pxtp10(jl,jk)=pxtp10(jl,jk)
     &              *(zclear(jl)/(1.-pclcover(jl,jk))
     &              +(1.-ziicscav)*pcfcover(jl,jk)/(1.-pclcover(jl,jk)))
              zdeps(jl)=zdeps(jl)+xdep*zmtof(jl)
             endif
           enddo
         endif

         do jl=1,klon
           zilcscav=0.
           if(pmlwc(jl,jk).gt.zmin)then
             zilcscav=Rcoeff(ktrac)
     &               *psolub(jl,jk)*(pmaccr(jl,jk)*ptmst/pmlwc(jl,jk))
             xdep=pxtp1c(jl,jk)*zilcscav*pclcover(jl,jk)
             pdep3d(jl,jk)=pdep3d(jl,jk)+xdep
             pxtp1c(jl,jk)=pxtp1c(jl,jk)*(1.-zilcscav)
             zdeps(jl)=zdeps(jl)+xdep*zmtof(jl)
           endif
         enddo


c Below-cloud scavenging by snow

        do jl=1,klon
          if(pfsnow(jl,jk).gt.zmin.and.pclcover(jl,jk).lt.1.-zmin)then
           plambda=min(plambs(jl,jk),8.e3) !Cut it off at about -30 deg. C
           zbcscav=zcollefs(ktrac)*plambda*pfsnow(jl,jk)*ptmst/(2*rhos)
           zbcscav=min(1.,zbcscav/(1.+0.5*zbcscav)) !Time-centred
           xbcscav=zbcscav*pxtp10(jl,jk)*zclr0(jl)
c     &          +(pxtp1c(jl,jk)*pclcover(jl,jk))
           pdep3d(jl,jk)=pdep3d(jl,jk)+xbcscav
c           pxtp1c(jl,jk)=pxtp1c(jl,jk)*(1.-zbcscav)
           pxtp10(jl,jk)=pxtp10(jl,jk)*(1.-zbcscav)
           zdeps(jl)=zdeps(jl)+xbcscav*zmtof(jl)
          end if
        enddo


c Redistribution by snow that evaporates or stays in layer

        do jl=1,klon
          if (pfsnow(jl,jk).gt.zmin) then
            zstay=(pfsubl(jl,jk)+pfstay(jl,jk))/pfsnow(jl,jk)
            zstay=min(1., zstay)
            xstay=zdeps(jl)*zstay*zftom(jl)
            zdeps(jl)=zdeps(jl)*(1.-zstay)
            zdeps(jl)=max(0.,zdeps(jl))
            pdep3d(jl,jk)=pdep3d(jl,jk)-xstay
            if(zclr0(jl).gt.zmin)then
              pxtp10(jl,jk)=pxtp10(jl,jk)+xstay/zclr0(jl)
            else
              pxtp1c(jl,jk)=pxtp1c(jl,jk)+xstay/pclcover(jl,jk)
            endif
          end if
        enddo


c Melting of snow... 

        do jl=1,klon
          if (pfmelt(jl,jk).gt.zmin) then
            zdepr(jl)=zdepr(jl)+zdeps(jl)
            zdeps(jl)=0.
          endif
        enddo

C
C  In-cloud scavenging by warm-rain processes (autoconversion and collection)

        do jl=1,klon
          if(pmratep(jl,jk).gt.zmin) then
            zicscav=psolub(jl,jk)*(pmratep(jl,jk)*ptmst/pmlwc(jl,jk))
            xicscav=pxtp1c(jl,jk)*zicscav*pclcover(jl,jk) !gridbox mean
            pxtp1c(jl,jk)=pxtp1c(jl,jk)*(1.-zicscav)
            pdep3d(jl,jk)=pdep3d(jl,jk)+xicscav
            zdepr(jl)=zdepr(jl)+xicscav*zmtof(jl)
          end if
        enddo


c Below-cloud scavenging by stratiform rain (conv done below)

        if(zcollefr(ktrac).gt.zmin)then
          do jl=1,klon
            if(pfprec(jl,jk).gt.zmin.and.zclr0(jl).gt.zmin)then
              zbcscav=zcollefr(ktrac)*prscav(jl,jk)
              zbcscav=min(1.,zbcscav/(1.+0.5*zbcscav)) !Time-centred
              xbcscav=zbcscav*pxtp10(jl,jk)*zclr0(jl)
              pdep3d(jl,jk)=pdep3d(jl,jk)+xbcscav
              pxtp10(jl,jk)=pxtp10(jl,jk)*(1.-zbcscav)
              zdepr(jl)=zdepr(jl)+xbcscav*zmtof(jl)
            end if
          enddo
        endif

C
C Reevaporation of rain

        do jl=1,klon
          if (pfprec(jl,jk).gt.zmin.and.zclear(jl).gt.zmin) then
            zevap=pfevap(jl,jk)/pfprec(jl,jk)
            zevap=min(1., zevap)
            if(zevap.lt.1.)zevap=Evfac(ktrac)*zevap
            xevap=zdepr(jl)*zevap*zftom(jl) !xevap is the grid-box-mean m.r. change
            zdepr(jl)=max(0.,zdepr(jl)*(1.-zevap))
            pdep3d(jl,jk)=pdep3d(jl,jk)-xevap
            prevap(jl,jk)=xevap
            pxtp10(jl,jk)=pxtp10(jl,jk)+xevap/zclr0(jl)
          end if
        enddo
C
C   END OF VERTICAL LOOP
  160 CONTINUE


c Now do the convective below-cloud bit...

      do jk=ktop,klev
        do jl=1,klon
          zmtof(jl)=pdpp1(jl,jk)*pcons2
          zftom(jl)=1./zmtof(jl)
          zclr0(jl)=1.-pclcover(jl,jk)-pclcon(jl,jk)
          
c Below-cloud scavenging by convective precipitation (assumed to be rain)

          if(pfconv(jl,jk-1).gt.zmin.and.zclr0(jl).gt.zmin)then
            fracc=0.3           !Convective rain fraction
            Frc=max(0.,pfconv(jl,jk-1)/fracc)
            zbcscav=zcollefr(ktrac)*fracc*0.24*ptmst*pow75(Frc)
            zbcscav=min(1.,zbcscav/(1.+0.5*zbcscav)) !Time-centred
            xbcscav=zbcscav*pxtp10(jl,jk)*zclr0(jl)
            conwd(jl,ktrac)=conwd(jl,ktrac)+xbcscav*zmtof(jl)
            pdep3d(jl,jk)=pdep3d(jl,jk)+xbcscav
            pxtp10(jl,jk)=pxtp10(jl,jk)*(1.-zbcscav)
          endif

c Below-cloud reevaporation of convective rain

          if(jk.gt.kbase(jl).and.
     &       pfconv(jl,jk-1).gt.zmin.and.zclr0(jl).gt.zmin)then
            pcevap=pfconv(jl,jk-1)-pfconv(jl,jk)
            zevap=pcevap/pfconv(jl,jk-1)
            zevap=max(0.,min(1.,zevap))
            if(zevap.lt.1.)zevap=Evfac(ktrac)*zevap
            xevap=conwd(jl,ktrac)*zevap*zftom(jl) !xevap is the grid-box-mean m.r. change
            conwd(jl,ktrac)=max(0.,conwd(jl,ktrac)*(1.-zevap))
            pdep3d(jl,jk)=pdep3d(jl,jk)-xevap
            prevap(jl,jk)=prevap(jl,jk)+xevap
            pxtp10(jl,jk)=pxtp10(jl,jk)+xevap/zclr0(jl)
          end if
        enddo
      enddo

      if(debug)then
        if(krow.eq.lgdebug)then
          ins=insdebug
          mg=mgdebug+(ins-1)*lon
          write(25,'(a,3i3)')'IPASS=1, in xtwetdep, JT= ',ktrac
          write(25,9)'pdep3d ',(pdep3d(mg,k),k=nl,1,-1)
          write(25,9)'prevap ',(prevap(mg,k),k=nl,1,-1)
          write(25,9)'pxtp1con ',(pxtp1con(mg,k),k=nl,1,-1)
          write(25,9)'pclcon ',(pclcon(mg,k),k=nl,1,-1)
          write(25,9)'pcscav ',(pcscav(mg,k),k=nl,1,-1)
          write(25,9)'psolubc ',(psolubc(mg,k),k=nl,1,-1)
          write(25,9)'pxtp1c ',(pxtp1c(mg,k),k=nl,1,-1)
          write(25,9)'psolub ',(psolub(mg,k),k=nl,1,-1)
          write(25,9)'pxtp10 ',(pxtp10(mg,k),k=nl,1,-1)
        endif
      endif

 9    format(a,30g10.3)

      RETURN
      END

C Start code : ----------------------------------------------------------

C    PHYSICAL CONSTANTS
C***      KTOP=4     !Top level for wet deposition (counting from top)
C***      PQTMST=1./PTMST
C***      ZCONS1=5.2
C***      ZCOLLEFR(:)=0.1
C***c      ZCOLLEFR(2)=0.05
C***      ZMIN=1.E-20
C***      ZEP=1.
C***      prevap(:,:)=0.
C***C
C***C
C***      CALL RESETR(ZDEP,KLON,0.)
C***      if(ktrac.eq.3)psolub(:,:)=1. !SO4 totally dissolved
C***C
C***C     BEGIN OF VERTICAL LOOP
C***        DO 160 JK=KTOP,KLEV
C***C
C***         DO 110 JL=1,KLON
C***          pfprec(jl,jk)=pfprec(jl,jk)+pfsnow(jl,jk) !Total precip flux
C***          pmratep(jl,jk)=pmratep(jl,jk)+pmaccr(jl,jk) !Total removal of ql
C***          ZCLEAR(JL)=1.-PCLCOVER(JL,JK)
C***          ZMTOF(JL)=PDPP1(JL,JK)*PCONS2
C***          ZFTOM(JL)=1./ZMTOF(JL)
C***          PXTP1C(JL,JK)=AMAX1(0.,PXTP1C(JL,JK))
C***          PXTP10(JL,JK)=AMAX1(0.,PXTP10(JL,JK))
C***  110    CONTINUE
C***C
C***C   2. IN-CLOUD SCAVENGING (GIORGI + CHAMEIDES)
C***         DO 120 JL=1,KLON
C***          IF(PMRATEP(JL,JK).GT.ZMIN) THEN
C***            ZTAU=0.5E-03/PMRATEP(JL,JK) !Using fixed LWC!!!
C***c            ZTAU=pmlwc(jl,jk)/PMRATEP(JL,JK) !Using actual LWC!!!
C***            ZF=0.8/(1.+0.8E-04*ZTAU)
C***            ZBETA=1.E-04+1./(0.8*ZTAU)
C***            ZICSCAV=ZF*(1.-EXP(-ZBETA*PTMST))
C***            ZICSCAV=AMAX1(0.,ZICSCAV)
C***            ZICSCAV=AMIN1(1.,ZICSCAV)
C***            ZEPS=ZEP*PSOLUB(JL,JK)
C***            ZICSCAV=ZEPS*ZICSCAV
C***
C***c            ZICSCAV=0.
C***C         NO SCAVENGING IN ICE-CLOUDS
C***C
C***            PDEP3D(JL,JK)=PDEP3D(JL,JK)+PXTP1C(JL,JK)*ZICSCAV
C***     *             *PCLCOVER(JL,JK)
C***            zdep(jl)=zdep(jl)+pdep3d(jl,jk)*zmtof(jl)
C***            PXTP1C(JL,JK)=PXTP1C(JL,JK)*(1.-ZICSCAV)
C***          END IF
C***  120    CONTINUE
C***C
C***C   3. REEVAPORATION
C***C***         DO 130 JL=1,KLON
C***C***          IF (PFPREC(JL,JK).GT.ZMIN.AND.ZCLEAR(JL).GT.ZMIN) THEN
C***C***           ZEVAP=PFEVAP(JL,JK)/PFPREC(JL,JK)
C***C***           ZEVAP=AMAX1(0.,ZEVAP)
C***C***           ZEVAP=AMIN1(1.,ZEVAP)
C***C***           xevap=ZDEP(JL)*ZEVAP*ZFTOM(JL) !xevap is the grid-box-mean m.r. change
C***C***c           PXTP10(JL,JK)=PXTP10(JL,JK)+xevap/zclear(jl)
C***C***           ZDEP(JL)=ZDEP(JL)*(1.-ZEVAP)
C***C***           PDEP3D(JL,JK)=PDEP3D(JL,JK)-xevap
C***C***          END IF
C***C***  130    CONTINUE
C***C
C***C   4. BELOW CLOUD SCAVENGING
C***         DO 140 JL=1,KLON
C***          IF(PFPREC(JL,JK).GT.ZMIN.AND.PCLCOVER(JL,JK).LE.ZMIN) THEN
C***c          ZZEFF=ZCOLLEFF*3.
C***c          ZBCSCAV=ZCONS1*ZZEFF*PFPREC(JL,JK)*ZFTOM(JL)*PRHOP1(JL,JK)
C***          zbcscav=1.04*zcollefr(ktrac)*pfprec(jl,jk)*ptmst/prhop1(jl,jk) !V=5 m/s approx.
C***          ZBCSCAV=AMIN1(1.,ZBCSCAV)
C***          ZBCSCAV=AMAX1(0.,ZBCSCAV)
C***c          ZBCSCAV=0.
C***          PDEP3D(JL,JK)=PDEP3D(JL,JK)+ZBCSCAV*PXTP10(JL,JK)
C***          PXTP10(JL,JK)=PXTP10(JL,JK)*(1.-ZBCSCAV)
C***          END IF
C***  140    CONTINUE
C***C
C***C
C***C   END OF VERTICAL LOOP
C***  160   CONTINUE
C***C
C***      RETURN
C***      END
