c Further modifications to the ocean model output routines to improve memory
c usage. The features that were causing stack overflows have now been resolved.
c SJP 2009/05/11
c
c Modifying the ocean model output routines for improved memory usage.
c SJP 2009/05/06
c
c Adding density as an ocean model statistic.
c SJP 2009/04/21
c
c Enabling user control over the ocean model statistics that are saved to file.
c SJP 2009/04/19
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Implementing the equation of state of McDougall et al (2003).
c SJP 2009/04/08
c
c Implementing a Robert time filter.
c SJP 2008/12/17
c
c Modified so that flux adjustments are only applied when FLUXADJ=T.
c SJP 2008/03/08
c
c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/17
c
c Added IMPLICIT NONE statement, plus variable declarations as necessary.
c SJP 2007/05/30
c
c Transferred COMMON block declarations to separate header files.
c SJP 2007/05/29
c
c Re-inserted the application of flux adjustments to the wind stresses in the
c coupled model.
c SJP 2004/11/18
c
c Removed URAT, which is essentially redundant.
c SJP 2004/02/14
c
c Removed the call to CLINIC1, which is redundant.
c SJP 2004/01/06
c
c Array bounds violations in loops 460, 524 and 520 fixed.
c SJP 2004/01/04
c
c Subroutine CLINIC1 split off into a separate file.
c SJP 2003/12/30
c
c (1) Some tidying up.
c (2) Adjacent loops with identical outer bounds fused, where possible.
c SJP 2003/12/20
c
c (1) Modified to make use of /orestart/.
c (2) Reomved the application of flux adjustments to the wind stresses in the
c     coupled model, as I do not intend to apply such adjustments.
c SJP 2003/09/04
c
c Various loops modified so that code passes array bounds checking.
c SJP 2003/05/01
c
c Parameter definitions moved to inlude file OPARAMS.f.
c SJP 2003/04/29
c
c $Log: clinic.f,v $
c Revision 1.8  1997/12/19 01:25:41  ldr
c Changes from ACH for 21 OGCM levels and optinal GM mixing scheme...
c
c Revision 1.7  1994/03/30  12:33:59  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c
c Revision 1.6  93/12/17  15:31:49  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.5  93/11/29  11:38:20  ldr
c Changes to V4-4-32l from HBG for coupled model
c 
c Revision 1.4  93/11/03  11:44:06  ldr
c Changes to V4-4-24l from HBG for coupled model
c 
c Revision 1.3  93/10/05  13:05:26  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.2  93/08/10  15:26:01  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c 
c Revision 1.1  93/02/05  16:21:58  ldr
c Initial revision
c 
      SUBROUTINE CLINIC(J)
C
C=======================================================================
C                                                                    ===
C  CLINIC COMPUTES, FOR ONE ROW, THE INTERNAL MODE COMPONENT OF      ===
C         THE U AND V VELOCITIES, AS WELL AS THE VORTICITY DRIVING   ===
C         FUNCTION FOR USE BY "RELAX" LATER IN DETERMINING THE       ===
C         EXTERNAL MODES, WHERE:                                     ===
C           J=THE ROW NUMBER                                         ===
C                                                                    ===
C=======================================================================
C
      implicit none
C
C---------------------------------------------------------------------
C  DEFINE GLOBAL DATA
C---------------------------------------------------------------------
C
      include 'OPARAMS.f'
      include 'OCEAN_NML.f'
      include 'ORESTART.f'
      include 'FULLWD.f'
      include 'SCALAR.f'
      include 'ONEDIM.f'
      include 'FIELDS.f'
      include 'WORKSP.f'
      include 'LEVD.f'
      include 'COEFFS.f'
      include 'SFC1.f'
      include 'FEWFLAGS.f'
      include 'FCOR.f'
      include 'TIME1.f'
      include 'OHIST.f'
      include 'bulkf_variables.f' ! pjb
C
C---------------------------------------------------------------------
C  DEFINE AND EQUIVALENCE LOCAL DATA
C---------------------------------------------------------------------
C
      integer i, j, k, kz, ll, is, ie, jj, isave, ieave, l, iredo,
     &        im, m, n, ism1, iea, ieb, ii
      real dpdx, dpdy, ueng, veng, fx, diag1, diag2, fxa, fxb, bbuj,
     &     ccuj, dduj, gguj, hhuj, detmr, engtmp
      DIMENSION DPDX(IMT,KM),DPDY(IMT,KM),UENG(IMT,KM),VENG(IMT,KM)
      EQUIVALENCE (TDIF(1,1,1),DPDX(1,1),UENG(1,1)),
     &            (TDIF(1,1,2),DPDY(1,1),VENG(1,1))
      real sden(imt), tden(imt), rden(imt)
C
C---------------------------------------------------------------------
C  BEGIN EXECUTABLE CODE
C---------------------------------------------------------------------
C
C=======================================================================
C  BEGIN INTRODUCTORY SECTION, PREPARING VARIOUS     ===================
C  ARRAYS FOR THE COMPUTATION OF THE INTERNAL MODES  ===================
C=======================================================================
C
C---------------------------------------------------------------------
C  FIND ADVECTIVE COEFFICIENT 'FUW' FOR WEST  FACE OF U,V BOX
C                           & 'FVN' FOR NORTH FACE OF U,V BOX
C---------------------------------------------------------------------
C
C  1ST, CALCULATE EXT. MODE U AT WEST  FACE OF U,V BOX
C                         & V AT NORTH FACE OF U,V BOX
C

      DO 100 I=2,IMTM1
        SFU(I)=-(P(I  ,J+1)-P(I,J  ))*DYUR(J)*MIN(HR(I-1,J  ),HR(I,J))
        SFV(I)= (P(I+1,J+1)-P(I,J+1))*DXUR(I)*MIN(HR(I  ,J+1),HR(I,J))
     &                                       *CSTR(J+1)
 100  CONTINUE

      SFU(1) = SFU(IMTM1)
      SFV(1) = SFV(IMTM1)
      SFU(IMT) = SFU(2)
      SFV(IMT) = SFV(2)
C
C  2ND, CALCULATE INT. MODE U AT WEST  FACE OF U,V BOX
C                         & V AT NORTH FACE OF U,V BOX
C
C  3RD, COMBINE EXT. AND INT. MODES AND ADD GRID WGT. FACTOR
C
C---------------------------------------------------------------------
C  SAVE INTERNAL MODE VELOCITIES
C---------------------------------------------------------------------
C

      FX=DYU2R(J)*CSR(J)*CST(J+1)

      DO K=1,KM
        DO I=2,IMT
          FUW(I,K)=(UCLIN(I,K)+UCLIN(I-1,K))*0.5
          FVN(I,K)=(VP   (I,K)+VCLIN(I  ,K))*0.5
        END DO
        FUW(1,K)=FUW(IMTM1,K)
        FVN(1,K)=FVN(IMTM1,K)
        do i = 1, imt
          FUW(I,K)=(FUW(I,K)+SFU(I))*CSR(J)
          FVN(I,K)=(FVN(I,K)+SFV(I))*FX
          USAV(I,K)=UCLIN(I,K)
          VSAV(I,K)=VCLIN(I,K)
          UCLIN(I,K)=UP(I,K)
          VCLIN(I,K)=VP(I,K)
        end do
      END DO

C
C---------------------------------------------------------------------
C  COMPUTE EXTERNAL MODE VELOCITIES FOR ROW J+1
C---------------------------------------------------------------------
C
C  1ST, COMPUTE FOR TAU-1 TIME LEVEL
C
C  2ND, COMPUTE FOR TAU TIME LEVEL
C

      DO 150 I=1,IMTM1
        DIAG1=PB(I+1,J+2)-PB(I  ,J+1)
        DIAG2=PB(I  ,J+2)-PB(I+1,J+1)
        SFUB(I)=-(DIAG1+DIAG2)*DYU2R(J+1)*HR(I,J+1)
        SFVB(I)= (DIAG1-DIAG2)*DXU2R(I  )*HR(I,J+1)*CSR(J+1)
        DIAG1=P (I+1,J+2)-P (I  ,J+1)
        DIAG2=P (I  ,J+2)-P (I+1,J+1)
        SFU (I)=-(DIAG1+DIAG2)*DYU2R(J+1)*HR(I,J+1)
        SFV (I)= (DIAG1-DIAG2)*DXU2R(I  )*HR(I,J+1)*CSR(J+1)
 150  CONTINUE

C
C  3RD, SET CYCLIC BOUNDARY CONDITIONS
C
      SFUB(IMT)=SFUB(2)
      SFVB(IMT)=SFVB(2)
      SFU (IMT)=SFU (2)
      SFV (IMT)=SFV (2)

c...  If using Robert time filter, save external mode
      if (robert_time_filter) then
        do i = 1, imt
          ssfubp(i) = sfub(i)
          ssfvbp(i) = sfvb(i)
        end do
      end if

C
C---------------------------------------------------------------------
C  ADD EXTERNAL MODE TO INTERNAL MODE FOR ROW J+1  (OCEAN PTS. ONLY)
C---------------------------------------------------------------------
C

      DO 170 K=1,KM
      DO 170 I=1,IMU
        IF(KMU(I,J+1).GE.KAR(K)) THEN
          UBP(I,K)=UBP(I,K)+SFUB(I)
          VBP(I,K)=VBP(I,K)+SFVB(I)
          UP (I,K)=UP (I,K)+SFU (I)
          VP (I,K)=VP (I,K)+SFV (I)
        ENDIF
 170  CONTINUE

C
C---------------------------------------------------------------------
C  ACCUMULATE KINETIC ENERGY FROM ROW J+1 EVERY NTSI TIMESTEPS
C---------------------------------------------------------------------
C
      IF(MOD(ITT,NTSI).EQ.0) THEN
        FX=0.25*CS(J+1)*DYU(J+1)

        do k = 1, km
          do i = 1, imt
            UENG(I,K)=(FX*(UP(I,K)*UP(I,K)+VP(I,K)*VP(I,K)))
     &                *C2DZQ(I,K)*DXUQ(I,K)
          end do
          do i = 2, imum1
            EKTOT=EKTOT+UENG(I,K)
          end do
        end do

      ENDIF
      CONTINUE
C
C
C---------------------------------------------------------------------
C  COMPUTE DENSITY OF ROW J+1
C---------------------------------------------------------------------
C
      if (m2003_eos) then

        do k = 1, km
          if (omask(j+1, k)) then
            do i = 1, imt
              sden(i) = 1000.0*tp(i, k, 2)+35.0
              tden(i) = tp(i, k, 1)
            end do
            call m2003(sden, tden, zdb(k), rden)
            do i = 1, imt
              rhon(i, k) = rden(i)
            end do
          else
            do i = 1, imt
              rhon(i, k) = 1.0
            end do
          end if
        end do

c.....  Save the density if required
        if (save_rho .and. (mxp .eq. 0)) then
          do k = 1, km
            do i = 2, imtm1
              ohist_rho(i, j+1, k) = ohist_rho(i, j+1, k) + rhon(i, k)
            end do
          end do
        end if

      else

        CALL STATE(TP(1,1,1),TP(1,1,2),RHON,TDIF(1,1,1),TDIF(1,1,2))

      end if
C
C  SET CYCLIC BOUNDARY CONDITIONS
C

      DO 232 K=1,KM
        RHON(IMT,K)=RHON(2,K)
 232  CONTINUE

C
C---------------------------------------------------------------------
C  COMPUTE VERTICAL VELOCITY IN U,V COLUMNS
C---------------------------------------------------------------------
C
C 1ST, SET VERTICAL VELOCITY AT THE SURFACE TO ZERO (RIGID-LID)
C      SET VERTICAL VELOCITY AT KMP1 TO ZERO
C
      FX=0.0

      DO 240 I=1,IMT
        W(I,1)=FX
        W(I,KMP1)=FX
 240  CONTINUE

C
C  2ND, COMPUTE CHANGE OF W BETWEEN LEVELS
C

      DO K=1,KMM1
        DO I=1,IMTM1
          W(I,K+1)=C2DZQ(I,K)*((FUW(I+1,K)-FUW (I,K))*DXU2RQ(I,K)
     &                          +FVN(I  ,K)-FVSU(I,K))
        END DO
        W(IMT, K+1) = W(2, K+1)
      END DO

C
C  3RD, INTEGRATE DOWNWARD FROM THE SURFACE
C

      DO 255 K=1,KMM1
      DO 255 I=1,IMT
        W(I,K+1)=W(I,K)+W(I,K+1)
 255  CONTINUE

C
C---------------------------------------------------------------------
C  COMPUTE HYDROSTATIC PRESSURE GRADIENT
C---------------------------------------------------------------------
C
C  1ST, COMPUTE IT AT THE FIRST LEVEL
C
      FXA=GRAV*DZZ(1)*CSR(J)
      FXB=GRAV*DZZ(1)*DYU2R(J)

      DO 260 I=1,IMTM1
        UDIF(I,1)=RHON(I+1,1)-RHOS(I  ,1)
        VDIF(I,1)=RHON(I  ,1)-RHOS(I+1,1)
        DPDX(I,1)=((UDIF(I,1)-VDIF(I,1))*FXA)*DXU2R(I)
        DPDY(I,1)= (UDIF(I,1)+VDIF(I,1))*FXB
 260  CONTINUE

      UDIF(IMT,1)=UDIF(2,1)
      VDIF(IMT,1)=VDIF(2,1)
      DPDX(IMT,1)=DPDX(2,1)
      DPDY(IMT,1)=DPDY(2,1)
C
C  2ND, COMPUTE THE CHANGE IN PRESSURE GRADIENT BETWEEN LEVELS
C
      FXA=GRAV*CSR(J)*0.5
      FXB=GRAV*DYU4R(J)

      do k = 2, km
        do i = 1, imt
          DPDX(I,K)=RHON(I,K-1)+RHON(I,K)
          DPDY(I,K)=RHOS(I,K-1)+RHOS(I,K)
        end do
        do i = 1, imtm1
          UDIF(I,K)=DPDX(I+1,K)-DPDY(I  ,K)
          VDIF(I,K)=DPDX(I  ,K)-DPDY(I+1,K)
          DPDX(I,K)=(FXA*(UDIF(I,K)-VDIF(I,K)))*DZZQ(I,K)*DXU2RQ(I,K)
          DPDY(I,K)=(FXB*(UDIF(I,K)+VDIF(I,K)))*DZZQ(I,K)
        end do
        UDIF(IMT,K)=UDIF(2,K)
        VDIF(IMT,K)=VDIF(2,K)
        DPDX(IMT,K)=DPDX(2,K)
        DPDY(IMT,K)=DPDY(2,K)
      end do

C
C  3RD, INTEGRATE DOWNWARD FROM THE FIRST LEVEL
C

      DO 275 K=1,KMM1
      DO 275 I=1,IMT
        DPDX(I,K+1)=DPDX(I,K)+DPDX(I,K+1)
        DPDY(I,K+1)=DPDY(I,K)+DPDY(I,K+1)
 275  CONTINUE

C
C---------------------------------------------------------------------
C  SET BOUNDARY CONDITIONS FOR THE COMPUTATION OF
C       VERTICAL DIFFUSION OF MOMENTUM
C---------------------------------------------------------------------

      If(.not.lcouple)Then

      DO 277 I = 1,IMT
         WSX(I,J) = W1*STRX(I,J,lmon1) + W2*STRX(I,J,lmon2)
         WSY(I,J) = W1*STRY(I,J,lmon1) + W2*STRY(I,J,lmon2)

         if (bulk_force) then
            
           call bulkf_formula(
     I          uwnd(i,j),vwnd(i,j),TB(i,1,1),spechum(i,j),airtem(i,j),
     I          swdown(i,j),fclt(i,j),fice(i,j),lwdown(i,j),
     I          rain(i,j),snow(i,j),roff(i,j),
     O          swflx(i,j),lwup(i,j),lwflx(i,j),senflx(i,j),latflx(i,j),
     O          hflx(i,j),fflx(i,j),WSX(I,J),WSY(I,J),
     I          i,j, ncar_bulk_method, mit_bulk_method)

           ! WSX and WSY overwritten by above call to bulkf_formula

           ! Convert from N/m2 to Dynes
           WSX(I,J) = WSX(I,J)*10.0
           WSY(I,J) = WSY(I,J)*10.0

         endif

  277  CONTINUE

      Else

c.... Note wt1=W1 and wt2=W2 (see step.f)

        if (fluxadj) then

          do i = 2, imtm1
            wsx(i, j) = strx(i, j, 1) - wt1 * txcor(i, j, lmon1)
     &                                - wt2 * txcor(i, j, lmon2)
            wsy(i, j) = stry(i, j, 1) - wt1 * tycor(i, j, lmon1)
     &                                - wt2 * tycor(i, j, lmon2)
          end do

        else

          do i = 2, imtm1
            wsx(i, j) = strx(i, j, 1)
            wsy(i, j) = stry(i, j, 1)
          end do

        end if

      WSX(1,J)=WSX(IMT-1,J)
      WSY(1,J)=WSY(IMT-1,J)
      WSX(IMT,J)=WSX(2,J)
      WSY(IMT,J)=WSY(2,J)

      End If

c...  Save the wind stresses if required
      if (save_smfzon .and. (mix .eq. 0)) then
        do i = 2, imtm1
          ohist_smfzon(i, j) = ohist_smfzon(i, j) + wsx(i, j)
        end do
      end if
      if (save_smfmer .and. (mix .eq. 0)) then
        do i = 2, imtm1
          ohist_smfmer(i, j) = ohist_smfmer(i, j) + wsy(i, j)
        end do
      end if

C
C  1ST, TRANSFER INTERIOR POINTS INTO DIFFUSION COMPUTATION ARRAYS
C

      DO 280 K=1,KM
      DO 280 I=1,IMT
        UDIF(I,K)=UB(I,K)
        VDIF(I,K)=VB(I,K)
 280  CONTINUE

C
C  2ND, SET K=0 ELEMENTS OF DIFF. COMP. ARRAYS TO REFLECT WIND STRESS
C
C  3RD, SET FIRST LAND LEVEL IN EACH COLUMN TO REFLECT BOTTOM CONDITION
C
C*************  INSERT RAYLEIGH FRICTION BOTTOM DRAG  *****************

      FX=DZZ(1)/FKPMF
      FXA = CDRAG/FKPMF

      DO 290 I=1,IMT
        UOVER(I)=UB(I,1)+WSX(I,J)*FX
        VOVER(I)=VB(I,1)+WSY(I,J)*FX
        KZ=KMU(I,J)
        if (kz .gt. 0) then
          if (kz .lt. km) then
            UDIF(I,KZ+1) = UB(I,KZ) - FXA*UB(I,KZ)*DZZ(KZ+1)
            VDIF(I,KZ+1) = VB(I,KZ) - FXA*VB(I,KZ)*DZZ(KZ+1)
          else
            UUNDER(I) = UB(I,KZ) - FXA*UB(I,KZ)*DZZ(KZ+1)
            VUNDER(I) = VB(I,KZ) - FXA*VB(I,KZ)*DZZ(KZ+1)
          end if
        end if
 290  CONTINUE

C
C=======================================================================
C  END INTRODUCTORY SECTION  ===========================================
C=======================================================================
C
C=======================================================================
C  BEGIN COMPUTATION OF THE INTERNAL MODES.                 ============
C  THE NEW VALUES "UA" AND "VA", WILL FIRST BE LOADED WITH  ============
C  THE TIME RATE OF CHANGE, AND THEN UPDATED.               ============
C=======================================================================
C
C---------------------------------------------------------------------
C  COMPUTE TOTAL ADVECTION OF MOMENTUM
C---------------------------------------------------------------------
C
C  1ST, COMPUTE FLUX THROUGH WEST FACE OF U,V BOX
C
C  2ND, COMPUTE ZONAL FLUX DIVERGENCE
C
C  3RD, ADD IN MERIDIONAL FLUX DIVERGENCE
C

      DO K=1,KM
        DO I=2,IMT
          TEMPA(I,K)=FUW(I,K)*(U(I-1,K)+U(I,K))
          TEMPB(I,K)=FUW(I,K)*(V(I-1,K)+V(I,K))
        END DO
        TEMPA(1,K)=TEMPA(IMTM1,K)
        TEMPB(1,K)=TEMPB(IMTM1,K)
        DO I=1,IMTM1
          UA(I,K)=(TEMPA(I,K)-TEMPA(I+1,K))*DXU2RQ(I,K)
          VA(I,K)=(TEMPB(I,K)-TEMPB(I+1,K))*DXU2RQ(I,K)
          UA(I,K)=UA(I,K)-FVN (I,K)*(UP(I,K)+U (I,K))
     &                   +FVSU(I,K)*(U (I,K)+UM(I,K))
          VA(I,K)=VA(I,K)-FVN (I,K)*(VP(I,K)+V (I,K))
     &                   +FVSU(I,K)*(V (I,K)+VM(I,K))
        END DO
        UA(IMT,K)=UA(2,K)
        VA(IMT,K)=VA(2,K)
      END DO

C
C  4TH, COMPUTE FLUX THROUGH TOP OF U,V BOX
C

      DO 340 K=2,KM
      DO 340 I=1,IMT
        TEMPA(I,K)=W(I,K)*(U(I,K-1)+U(I,K))
        TEMPB(I,K)=W(I,K)*(V(I,K-1)+V(I,K))
 340  CONTINUE

      DO I = 1, IMT
        TEMPA(I, 1) = 0.0
        TEMPB(I, 1) = 0.0
        TEMPA(I, KMP1) = 0.0
        TEMPB(I, KMP1) = 0.0
      END DO

C
C  5TH, ADD IN VERTICAL FLUX DIVERGENCE
C

      DO 343 K=1,KM
      DO 343 I=1,IMT
        UA(I,K)=UA(I,K)+(TEMPA(I,K+1)-TEMPA(I,K))*DZ2RQ(I,K)
        VA(I,K)=VA(I,K)+(TEMPB(I,K+1)-TEMPB(I,K))*DZ2RQ(I,K)
 343  CONTINUE

C
C---------------------------------------------------------------------
C  ADD IN HORIZONTAL DIFFUSION OF MOMENTUM (EVAL. AT TAU-1 TSTEP)
C---------------------------------------------------------------------
C
C  1ST, COMPUTE SEVERAL COEFFICIENTS DEPENDENT ONLY ON LATITUDE
C
      BBUJ=8.0*AMF*CSR(J)*CSR(J)
      CCUJ=AMF*CST(J+1)*DYTR(J+1)*DYUR(J)*CSR(J)
      DDUJ=AMF*CST(J  )*DYTR(J  )*DYUR(J)*CSR(J)
      GGUJ=AMF*(1.0-TNG(J)*TNG(J))/(RADIUS*RADIUS)
      HHUJ=2.0*AMF*SINE(J)/(RADIUS*CS(J)*CS(J))
C ***** REDUCE HORIZ VISCOSITY OVER ARCTIC, TO SATISFY STABILITY CRITERI
      BBUJ = BBUJ * AMFAC(J)
      CCUJ = CCUJ * AMFAC(J)
      DDUJ = DDUJ * AMFAC(J)
      GGUJ = GGUJ * AMFAC(J)
      HHUJ = HHUJ * AMFAC(J)
C
C  2ND, COMPUTE GRADIENTS AT WEST FACE OF U,V BOX
C
C  3RD, ADD IN FINAL CONTRIBUTION FROM HOR. DIFF. OF MOMENTUM
C

      DO K=1,KM
        DO I=2,IMT
          TEMPA(I,K)=DXT4RQ(I,K)*(UB(I,K)-UB(I-1,K))
          TEMPB(I,K)=DXT4RQ(I,K)*(VB(I,K)-VB(I-1,K))
        END DO
        TEMPA(1,K)=TEMPA(IMTM1,K)
        TEMPB(1,K)=TEMPB(IMTM1,K)
        DO I=2,IMTM1
          UA(I,K)=UA(I,K)+BBUJ*(DXU2RQ(I,K)*(TEMPA(I+1,K)-TEMPA(I,K)))
     &            +CCUJ*(UBP(I,K)-UB(I,K))+DDUJ*(UBM(I,K)-UB(I,K))
     &            +GGUJ*UB(I,K)-HHUJ*DXU2RQ(I,K)*(VB(I+1,K)-VB(I-1,K))
          VA(I,K)=VA(I,K)+BBUJ*(DXU2RQ(I,K)*(TEMPB(I+1,K)-TEMPB(I,K)))
     &            +CCUJ*(VBP(I,K)-VB(I,K))+DDUJ*(VBM(I,K)-VB(I,K))
     &            +GGUJ*VB(I,K)+HHUJ*DXU2RQ(I,K)*(UB(I+1,K)-UB(I-1,K))
        END DO
        UA(1,K)=UA(IMTM1,K)
        VA(1,K)=VA(IMTM1,K)
        UA(IMT,K)=UA(2,K)
        VA(IMT,K)=VA(2,K)
      END DO

C
C---------------------------------------------------------------------
C  ADD IN VERTICAL DIFFUSION OF MOMENTUM
C---------------------------------------------------------------------
C
C  1ST, COMPUTE GRADIENTS AT TOP OF U,V BOX
C

      DO 345 K=2,KM
      DO 345 I=1,IMT
        TEMPA(I,K)=UDIF(I,K-1)-UDIF(I,K)
        TEMPB(I,K)=VDIF(I,K-1)-VDIF(I,K)
 345  CONTINUE

      DO I = 1, IMT
        TEMPA(I,1) = UOVER(I)-UDIF(I,1)
        TEMPB(I,1) = VOVER(I)-VDIF(I,1)
        TEMPA(I,KMP1) = UDIF(I,KM)-UUNDER(I) 
        TEMPB(I,KMP1) = VDIF(I,KM)-VUNDER(I)
      END DO

C
C  2ND, ADD IN FINAL CONTRIBUTION FROM VERT. DIFF. OF MOMENTUM
C

      DO 348 K=1,KM
      DO 348 I=1,IMT
        UA(I,K)=UA(I,K)+EEMQ(I,K)*TEMPA(I,K)-FFMQ(I,K)*TEMPA(I,K+1)
        VA(I,K)=VA(I,K)+EEMQ(I,K)*TEMPB(I,K)-FFMQ(I,K)*TEMPB(I,K+1)
 348  CONTINUE

C
C---------------------------------------------------------------------
C  ADD IN CORIOLIS FORCE (EVAL. ON TAU   TSTEP FOR EXPLICIT TRTMNT;
C                         EVAL. ON TAU-1 TSTEP FOR IMPLICIT TRTMNT
C                         WITH REMAINDER OF TERM TO BE ADDED LATER)
C---------------------------------------------------------------------
C
C---------------------------------------------------------------------
C  ADD IN PRESSURE TERM AND MASK OUT LAND
C---------------------------------------------------------------------
C
      FX=2.0*OMEGA*SINE(J)
      IF(ACORF.EQ.0.) THEN

        DO 357 K=1,KM
        DO 357 I=1,IMT
          UA(I,K)=UA(I,K)+FX*V(I,K)
          VA(I,K)=VA(I,K)-FX*U(I,K)
          UA(I,K)=GM(I,K)*(UA(I,K)-DPDX(I,K))
          VA(I,K)=GM(I,K)*(VA(I,K)-DPDY(I,K))
 357    CONTINUE

      ELSE

        DO 359 K=1,KM
        DO 359 I=1,IMT
          UA(I,K)=UA(I,K)+FX*VB(I,K)
          VA(I,K)=VA(I,K)-FX*UB(I,K)
          UA(I,K)=GM(I,K)*(UA(I,K)-DPDX(I,K))
          VA(I,K)=GM(I,K)*(VA(I,K)-DPDY(I,K))
 359    CONTINUE

      ENDIF
C
C---------------------------------------------------------------------
C  FORM TIME CHANGE OF VERTICALLY AVERAGED FORCING
C---------------------------------------------------------------------
C
C  1ST, INTEGRATE TIME CHANGE VERTICALLY
C
      FX=0.0

      DO 370 I=1,IMT
        ZUN(I)=FX
        ZVN(I)=FX
 370  CONTINUE

      DO 380 K=1,KM
        FX=C2DTSF*DZ(K)
      DO 380 I=1,IMT
        ZUN(I)=ZUN(I)+UA(I,K)*FX
        ZVN(I)=ZVN(I)+VA(I,K)*FX
 380  CONTINUE

C
C  2ND, FORM AVERAGE BY DIVIDING BY DEPTH
C

      DO 390 I=1,IMT
        ZUN(I)=ZUN(I)*HR(I,J)
        ZVN(I)=ZVN(I)*HR(I,J)
 390  CONTINUE

C
C---------------------------------------------------------------------
C  DO ANALYSIS OF INTERNAL MODE FORCING ON ENERGY TIMESTEP
C  ALSO, FORM VERT AVE. FOR USE LATER IN EXT. MODE ANALYSIS
C---------------------------------------------------------------------
C
      FX=0.0

      DO 395 LL=1,8
      DO 395 I=1,IMT
        ZUNENG(I,LL)=FX
        ZVNENG(I,LL)=FX
 395  CONTINUE

C
C  1ST, COMPUTE KE CHANGE DUE TO PRESSURE TERM
C
      IF(NERGY.EQ.0.OR.MXP.EQ.1) GO TO 550
      FX=CS(J)*DYU(J)

      DO 400 K=1,KM
      DO 400 I=2,IMUM1
        UENG(I,K)=GM(I,K)*(-DPDX(I,K))
        VENG(I,K)=GM(I,K)*(-DPDY(I,K))
        ENGINT(6)=ENGINT(6)+(USAV(I,K)*UENG(I,K)
     &                      +VSAV(I,K)*VENG(I,K))*FX*DXU(I)*DZ(K)
        ZUNENG(I,6)=ZUNENG(I,6)+UENG(I,K)*DZ(K)*HR(I,J)
        ZVNENG(I,6)=ZVNENG(I,6)+VENG(I,K)*DZ(K)*HR(I,J)
 400  CONTINUE

C
C  2ND, COMPUTE KE CHANGE DUE TO ADVECTION OF MOMENTUM
C

      DO 430 K=1,KM
      DO 430 I=2,IMUM1
        UENG(I,K)=GM(I,K)*((-FUW (I+1,K)*(U (I+1,K)+U (I  ,K))
     &                      +FUW (I  ,K)*(U (I  ,K)+U (I-1,K)))*DXU2R(I)
     &                      -FVN (I  ,K)*(UP(I  ,K)+U (I  ,K))
     &                      +FVSU(I  ,K)*(U (I  ,K)+UM(I  ,K)))
        VENG(I,K)=GM(I,K)*((-FUW (I+1,K)*(V (I+1,K)+V (I  ,K))
     &                      +FUW (I  ,K)*(V (I  ,K)+V (I-1,K)))*DXU2R(I)
     &                      -FVN (I  ,K)*(VP(I  ,K)+V (I  ,K))
     &                      +FVSU(I  ,K)*(V (I  ,K)+VM(I  ,K)))
        ENGINT(2)=ENGINT(2)+(USAV(I,K)*UENG(I,K)
     &                      +VSAV(I,K)*VENG(I,K))*FX*DXU(I)*DZ(K)
        ZUNENG(I,2)=ZUNENG(I,2)+UENG(I,K)*DZ(K)*HR(I,J)
        ZVNENG(I,2)=ZVNENG(I,2)+VENG(I,K)*DZ(K)*HR(I,J)
 430  CONTINUE

      DO 460 K=1,KM
      DO 460 I=2,IMUM1
        if (k .eq. 1) then
            UENG(I,1)=GM(I,1)*(-(W(I,1)*(U(I,1))
     &                          -W(I,2)*(U(I,1)+U(I,2)))*DZ2RQ(I,1))
            VENG(I,1)=GM(I,1)*(-(W(I,1)*(V(I,1))
     &                          -W(I,2)*(V(I,1)+V(I,2)))*DZ2RQ(I,1))
        else if (k .eq. km) then
          UENG(I,KM)=GM(I,KM)*(-(W(I,KM)*(U(I,KM-1)+U(I,KM))
     &                        -W(I,KM+1)*U(I,KM))*DZ2RQ(I,KM))
          VENG(I,KM)=GM(I,KM)*(-(W(I,KM)*(V(I,KM-1)+V(I,KM))
     &                        -W(I,KM+1)*V(I,KM))*DZ2RQ(I,KM))
        else
          UENG(I,K)=GM(I,K)*(-(W(I,K  )*(U(I,K-1)+U(I,K  ))
     &                        -W(I,K+1)*(U(I,K  )+U(I,K+1)))*DZ2RQ(I,K))
          VENG(I,K)=GM(I,K)*(-(W(I,K  )*(V(I,K-1)+V(I,K  ))
     &                        -W(I,K+1)*(V(I,K  )+V(I,K+1)))*DZ2RQ(I,K))
        end if
        ENGINT(3)=ENGINT(3)+(USAV(I,K)*UENG(I,K)
     &                      +VSAV(I,K)*VENG(I,K))*FX*DXU(I)*DZ(K)
        ZUNENG(I,3)=ZUNENG(I,3)+UENG(I,K)*DZ(K)*HR(I,J)
        ZVNENG(I,3)=ZVNENG(I,3)+VENG(I,K)*DZ(K)*HR(I,J)
 460  CONTINUE

C
C  3RD, COMPUTE KE CHANGE DUE TO HOR. DIFFUSION OF MOMENTUM
C

      DO 490 K=1,KM
      DO 490 I=2,IMUM1
        UENG(I,K)=GM(I,K)*(
     &            +BBUJ*DXU2R(I)*(DXT4R(I+1)*(UB(I+1,K)-UB(I,K))
     &                           +DXT4R(I  )*(UB(I-1,K)-UB(I,K)))
     &            +CCUJ*(UBP(I,K)-UB(I,K))+DDUJ*(UBM(I,K)-UB(I,K))
     &            +GGUJ*UB(I,K)-HHUJ*DXU2R(I)*(VB(I+1,K)-VB(I-1,K)))
        VENG(I,K)=GM(I,K)*(
     &            +BBUJ*DXU2R(I)*(DXT4R(I+1)*(VB(I+1,K)-VB(I,K))
     &                           +DXT4R(I  )*(VB(I-1,K)-VB(I,K)))
     &            +CCUJ*(VBP(I,K)-VB(I,K))+DDUJ*(VBM(I,K)-VB(I,K))
     &            +GGUJ*VB(I,K)+HHUJ*DXU2R(I)*(UB(I+1,K)-UB(I-1,K)))
        ENGINT(4)=ENGINT(4)+(USAV(I,K)*UENG(I,K)
     &                      +VSAV(I,K)*VENG(I,K))*FX*DXU(I)*DZ(K)
        ZUNENG(I,4)=ZUNENG(I,4)+UENG(I,K)*DZ(K)*HR(I,J)
        ZVNENG(I,4)=ZVNENG(I,4)+VENG(I,K)*DZ(K)*HR(I,J)
 490  CONTINUE

C
C  4TH, COMPUTE KE CHANGE DUE TO WIND STRESS
C
C  5TH, COMPUTE KE CHANGE DUE TO BOTTOM DRAG
C

      DO 522 I=2,IMUM1
        UENG(I,1)=GM(I,1)*EEM(1)*(UOVER(I)-UDIF(I,1))
        VENG(I,1)=GM(I,1)*EEM(1)*(VOVER(I)-VDIF(I,1))
        ENGINT(7)=ENGINT(7)+(USAV(I,1)*UENG(I,1)
     &                      +VSAV(I,1)*VENG(I,1))*FX*DXU(I)*DZ(1)
        ZUNENG(I,7)=ZUNENG(I,7)+UENG(I,1)*DZ(1)*HR(I,J)
        ZVNENG(I,7)=ZVNENG(I,7)+VENG(I,1)*DZ(1)*HR(I,J)
        KZ=KMU(I,J)
        IF(KZ.EQ.0) GO TO 522
        if (kz .eq. km) then
          UENG(I,KM)=GM(I,KM)*FFM(KM)*(0.0-UDIF(I,KM))
          VENG(I,KM)=GM(I,KM)*FFM(KM)*(0.0-VDIF(I,KM))
        else
          UENG(I,KZ)=GM(I,KZ)*FFM(KZ)*(UDIF(I,KZ+1)-UDIF(I,KZ))
          VENG(I,KZ)=GM(I,KZ)*FFM(KZ)*(VDIF(I,KZ+1)-VDIF(I,KZ))
        end if
        ENGINT(8)=ENGINT(8)+(USAV(I,KZ)*UENG(I,KZ)
     &                      +VSAV(I,KZ)*VENG(I,KZ))*FX*DXU(I)*DZ(KZ)
        ZUNENG(I,8)=ZUNENG(I,8)+UENG(I,KZ)*DZ(KZ)*HR(I,J)
        ZVNENG(I,8)=ZVNENG(I,8)+VENG(I,KZ)*DZ(KZ)*HR(I,J)
 522  CONTINUE

C
C  6TH, COMPUTE KE CHANGE DUE TO VERT. DIFFUSION OF MOMENTUM
C

      DO 520 I=2,IMUM1
        KZ=KMU(I,J)
        IF(KZ.EQ.0) GO TO 520
      DO 521 K=1,KZ
        if (k .eq. 1) then
          UENG(I,1)=GM(I,1)*(0.0-FFM(1)*(UDIF(I,1)-UDIF(I,2)))
          VENG(I,1)=GM(I,1)*(0.0-FFM(1)*(VDIF(I,1)-VDIF(I,2)))
        else if (k .eq. km) then
          UENG(I,KM)=GM(I,KM)*EEM(KM)*(UDIF(I,KM-1)-UDIF(I,KM))
          VENG(I,KM)=GM(I,KM)*EEM(KM)*(VDIF(I,KM-1)-VDIF(I,KM))
        else
          UENG(I,K)=GM(I,K)*(EEM(K)*(UDIF(I,K-1)-UDIF(I,K  ))
     &                      -FFM(K)*(UDIF(I,K  )-UDIF(I,K+1)))
          VENG(I,K)=GM(I,K)*(EEM(K)*(VDIF(I,K-1)-VDIF(I,K  ))
     &                      -FFM(K)*(VDIF(I,K  )-VDIF(I,K+1)))
        end if
        ENGINT(5)=ENGINT(5)+(USAV(I,K)*UENG(I,K)
     &                      +VSAV(I,K)*VENG(I,K))*FX*DXU(I)*DZ(K)
        ZUNENG(I,5)=ZUNENG(I,5)+UENG(I,K)*DZ(K)*HR(I,J)
        ZVNENG(I,5)=ZVNENG(I,5)+VENG(I,K)*DZ(K)*HR(I,J)
 521  CONTINUE
 520  CONTINUE

 550  CONTINUE
C
C---------------------------------------------------------------------
C  COMPUTE NEW VELOCITIES (WITH INCORRECT VERTICAL MEANS)
C  ALSO, ADD IN REMAINDER OF CORIOLIS TERM IF TREATED IMPLICITLY
C---------------------------------------------------------------------
C
      IF(ACORF.EQ.0.) THEN

        DO 560 K=1,KM
        DO 560 I=1,IMT
          UA(I,K)=UB(I,K)+C2DTUV*UA(I,K)
          VA(I,K)=VB(I,K)+C2DTUV*VA(I,K)
 560    CONTINUE

      ELSE
        FX=C2DTUV*ACORF*2.0*OMEGA*SINE(J)
        DETMR=1.0/(1.0+FX*FX)

        DO 565 K=1,KM
        DO 565 I=1,IMT
          UDIF(I,K)=(UA(I,K)+FX*VA(I,K))*DETMR
          VDIF(I,K)=(VA(I,K)-FX*UA(I,K))*DETMR
          UA(I,K)=UB(I,K)+C2DTUV*UDIF(I,K)
          VA(I,K)=VB(I,K)+C2DTUV*VDIF(I,K)
 565    CONTINUE

      ENDIF
C
C---------------------------------------------------------------------
C  DETERMINE THE INCORRECT VERTICAL MEANS OF THE NEW VELOCITIES
C---------------------------------------------------------------------
C
      FX=0.0

      DO 575 I=1,IMT
        SFU(I)=FX
        SFV(I)=FX
 575  CONTINUE

      DO 580 K=1,KM
      DO 580 I=1,IMT
        SFU(I)=SFU(I)+UA(I,K)*DZ(K)
        SFV(I)=SFV(I)+VA(I,K)*DZ(K)
 580  CONTINUE

      DO 590 I=1,IMT
        SFU(I)=SFU(I)*HR(I,J)
        SFV(I)=SFV(I)*HR(I,J)
 590  CONTINUE

C
C---------------------------------------------------------------------
C  SUBTRACT INCORRECT VERTICAL MEAN TO GET INTERNAL MODE
C---------------------------------------------------------------------
C

      DO 600 K=1,KM
      DO 600 I=1,IMT
        UA(I,K)=UA(I,K)-SFU(I)
        VA(I,K)=VA(I,K)-SFV(I)
        UA(I,K)=GM(I,K)*UA(I,K)
        VA(I,K)=GM(I,K)*VA(I,K)
 600  CONTINUE

C
C---------------------------------------------------------------------
C  COMPUTE TOTAL CHANGE OF K.E. OF INTERNAL MODE ON ENERGY TIMESTEP
C---------------------------------------------------------------------
C
      IF(NERGY.EQ.1 .AND. MXP.NE.1) THEN

      DO 605 K=1,KM
        FX=CS(J)*DYU(J)*DZ(K)/C2DTUV
        DO 605 I=2,IMUM1
          ENGINT(1)=ENGINT(1)+(USAV(I,K)*(UA(I,K)-UB(I,K))
     &                        +VSAV(I,K)*(VA(I,K)-VB(I,K)))*FX*DXU(I)
 605    CONTINUE

      ENDIF
C
C=======================================================================
C  END COMPUTATION OF INTERNAL MODES  ==================================
C=======================================================================
C
C=======================================================================
C  BEGIN COMPUTATION OF VORTICITY FOR INPUT TO "RELAX"  ================
C=======================================================================
C
C---------------------------------------------------------------------
C  SET CYCLIC BOUNDARY CONDITIONS ON EXT. MODE FORCING FUNCTIONS
C---------------------------------------------------------------------
C
      ZUN(1)=ZUN(IMUM1)
      ZVN(1)=ZVN(IMUM1)
      IF(NERGY.EQ.0.OR.MXP.EQ.1) GO TO 613

      DO 612 LL=2,8
      ZUNENG(1,LL)=ZUNENG(IMUM1,LL)
      ZVNENG(1,LL)=ZVNENG(IMUM1,LL)
 612  CONTINUE

 613  CONTINUE
C
C---------------------------------------------------------------------
C  FORM CURL OF TIME CHANGE IN VERTICALLY AVERAGED EQUATIONS
C---------------------------------------------------------------------
C
C  ALL VORTICITY PTS. ARE COMPUTED SO THAT THOSE NEEDED FOR THE LINE
C  INTEGRAL OF HOLE RELAXATION (IMMEDIATELY ADJACENT TO ISLANDS) WILL
C  BE DEFINED.
C
      IS=2
      IE=IMTM1

      DO 620 I=IS,IE
        ZTD(I,J)=((ZUN(I)*DXU(I)+ZUN(I-1)*DXU(I-1))*CS(J  )
     &           -(ZUS(I)*DXU(I)+ZUS(I-1)*DXU(I-1))*CS(J-1))
        ZTD(I,J)=(((ZVN(I)-ZVN(I-1))*DYU(J  )
     &            +(ZVS(I)-ZVS(I-1))*DYU(J-1)
     &            -ZTD(I,J))*DXT2R(I)*DYTR(J))*CSTR(J)
 620  CONTINUE

      CONTINUE
C
C---------------------------------------------------------------------
C  DO ANALYSIS OF EXTERNAL MODE FORCING ON ENERGY TIMESTEP
C---------------------------------------------------------------------
C
      IF(NERGY.EQ.1. AND .MXP.NE.1) THEN

        DO 630 LL=2,8
        DO 630 I=2,IMTM1
          ENGEXT(LL)=ENGEXT(LL)
     &     -P(I,J)*(((ZVNENG(I,LL)-ZVNENG(I-1,LL))*DYU(J)
     &              +(ZVSENG(I,LL)-ZVSENG(I-1,LL))*DYU(J-1))
     &            *DXT2R(I)*DYTR(J)
     &         -((ZUNENG(I,LL)*DXU(I)+ZUNENG(I-1,LL)*DXU(I-1))*CS(J)
     &          -(ZUSENG(I,LL)*DXU(I)+ZUSENG(I-1,LL)*DXU(I-1))*CS(J-1))
     &            *DYT2R(J)*DXTR(I))*DXT(I)*DYT(J)
 630    CONTINUE

        FX=CST(J)*DYT(J)/C2DTSF
        ENGTMP=0.

        DO 635 I=2,IMTM1
          ENGTMP   =ENGTMP   -P(I,J)*ZTD(I,J)*FX*DXT(I)
 635    CONTINUE

        ENGEXT(1)=ENGEXT(1)+ENGTMP
      ENDIF
C
C=======================================================================
C  END COMPUTATION OF VORTICITY  =======================================
C=======================================================================
C
C---------------------------------------------------------------------
C  FOURIER FILTER U AND V AT HIGH LATITUDES
C---------------------------------------------------------------------
C
      IF((J.GT.JFU1.AND.J.LT.JFU2).OR.J.LT.JFRST)GO TO 840
      JJ=J-JFRST+1
      IF (J.GE.JFU2) JJ=JJ-JSKPU+1
      FX=-1.0
      IF (PHI(J).GT.0.) FX=1.0
      ISAVE=0
      IEAVE=0

      DO 740 L=1,LSEGF
        DO 730 K=1,KM
          IF(ISUF(JJ,L,K).EQ.0) GO TO 730
          IS=ISUF(JJ,L,K)
          IE=IEUF(JJ,L,K)
          IREDO=1
          IF(IS.NE.ISAVE .OR. IE.NE.IEAVE) THEN
            IREDO=0
            IM=IE-IS+1
            ISAVE=IS
            IEAVE=IE
            IF(IM.NE.IMTM2) THEN
              M=2
              N=NINT(IM*CS(J)*CSR(JFU0))
            ELSE
              M=3
              N=NINT(IM*CS(J)*CSR(JFU0)*0.5)
            ENDIF
          ENDIF
          ISM1=IS-1
          IEA=IE
          IF(IE.GE.IMU) IEA=IMUM1
          DO 700 I=IS,IEA
            UDIF(I-ISM1,K)=-FX*UA(I,K)*SPSIN(I)-VA(I,K)*SPCOS(I)
            VDIF(I-ISM1,K)= FX*UA(I,K)*SPCOS(I)-VA(I,K)*SPSIN(I)
  700     CONTINUE
          IF(IE.GE.IMU) THEN
            IEB=IE-IMUM2
            II=IMUM1-IS
            DO 702 I=2,IEB
              UDIF(I+II,K)=-FX*UA(I,K)*SPSIN(I)-VA(I,K)*SPCOS(I)
              VDIF(I+II,K)= FX*UA(I,K)*SPCOS(I)-VA(I,K)*SPSIN(I)
 702        CONTINUE
          ENDIF
C
          CALL FILTER(UDIF(1,K),IM,M,N,IREDO)
C
          CALL FILTER(VDIF(1,K),IM,M,N,1)
C
          DO 720 I=IS,IEA
            UA(I,K)=FX*(-UDIF(I-ISM1,K)*SPSIN(I)
     &                  +VDIF(I-ISM1,K)*SPCOS(I))
            VA(I,K)=-UDIF(I-ISM1,K)*SPCOS(I)-VDIF(I-ISM1,K)*SPSIN(I)
  720     CONTINUE
          IF(IE.GE.IMT) THEN
            DO 722 I=2,IEB
              UA(I,K)=FX*(-UDIF(I+II,K)*SPSIN(I)
     &                    +VDIF(I+II,K)*SPCOS(I))
              VA(I,K)=-UDIF(I+II,K)*SPCOS(I)-VDIF(I+II,K)*SPSIN(I)
 722        CONTINUE
          ENDIF
  730   CONTINUE
  740 CONTINUE

      DO 750 I=1,IMT
        UOVER(I)=0.0
        VOVER(I)=0.0
  750 CONTINUE

      DO 760 K=1,KM
      DO 760 I=1,IMT
        UOVER(I)=UOVER(I)+UA(I,K)*DZ(K)
        VOVER(I)=VOVER(I)+VA(I,K)*DZ(K)
  760 CONTINUE

      DO 770 I=1,IMT
        UOVER(I)=UOVER(I)*HR(I,J)
        VOVER(I)=VOVER(I)*HR(I,J)
  770 CONTINUE

      DO 780 K=1,KM
      DO 780 I=1,IMT
        UA(I,K)=UA(I,K)-UOVER(I)
        VA(I,K)=VA(I,K)-VOVER(I)
        UA(I,K)=UA(I,K)*GM(I,K)
        VA(I,K)=VA(I,K)*GM(I,K)
  780 CONTINUE

C
C---------------------------------------------------------------------
C   FOURIER FILTER ZTD AT HIGH LATITUDES
C---------------------------------------------------------------------
C
      IF(J.EQ.JFU2) GO TO 840

      DO 830 L=1,LSEGF
        IS=ISZF(JJ,L)
        IF(IS.EQ.0) GO TO 840
        IE=IEZF(JJ,L)
        DO 800 II=IS,IE
          I=MOD(II-2,IMTM2)+2
          UDIF(II+1-IS,1)=ZTD(I,J)
  800   CONTINUE
        IM=IE-IS+1
        IF(IM.NE.IMTM2) THEN
           M=1
           N=NINT(IM*CST(J)*CSR(JFU0))
        ELSE
           M=3
           N=NINT(IM*CST(J)*CSR(JFU0)*0.5)
        ENDIF
C
        CALL FILTER(UDIF(1,1),IM,M,N,0)
C
        DO 820 II=IS,IE
          I=MOD(II-2,IMTM2)+2
          ZTD(I,J)=UDIF(II+1-IS,1)
  820   CONTINUE
  830 CONTINUE

  840 CONTINUE
C
C---------------------------------------------------------------------
C  TRANSFER QUANTITIES COMPUTED TO THE NORTH OF THE PRESENT ROW
C  TO BE DEFINED TO THE SOUTH IN THE COMPUTATION OF THE NEXT ROW
C---------------------------------------------------------------------
C
      FX=CS(J)*DYU(J)*CSR(J+1)*DYUR(J+1)

      DO 644 K=1,KM
      DO 644 I=1,IMT
        FVSU(I,K)=FVN(I,K)*FX
 644  CONTINUE

      DO 650 I=1,IMT
        ZUS(I)=ZUN(I)
        ZVS(I)=ZVN(I)
 650  CONTINUE

      IF(NERGY.EQ.1. AND .MXP.NE.1) THEN

        DO 660 LL=2,8
        DO 660 I=1,IMT
          ZUSENG(I,LL)=ZUNENG(I,LL)
          ZVSENG(I,LL)=ZVNENG(I,LL)
 660    CONTINUE

      ENDIF
C
C---------------------------------------------------------------------
C  SET CYCLIC BOUNDARY CONDITIONS ON NEWLY COMPUTED INTERNAL MODE
C---------------------------------------------------------------------
C

      DO 662 K=1,KM
        UA(1,K)=UA(IMUM1,K)
        VA(1,K)=VA(IMUM1,K)
        UA(IMU,K)=UA(2,K)
        VA(IMU,K)=VA(2,K)
 662  CONTINUE

C
      RETURN
      END
