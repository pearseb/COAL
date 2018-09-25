# 1 "relax.F"
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Implementing a Robert time filter.
c SJP 2008/12/17
c
c Modified so that a warning message is displayed whenever the maximum
c permitted number of relaxation scans is reached.
c SJP 2008/03/09
c
c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Modified to remove the fix to ISMASK in the Indonesian region when using the
c new ocean model grid, as this is no longer required.
c SJP 2007/10/10

c Major tidy-up of the ocean model source code.
c SJP 2007/06/16
c
c Further modification to enable the new ocean model grid.
c SJP 2007/06/11
c
c Modified to use OCEAN_LOW preprocessor macro, rather than OCEAN_DBL.
c SJP 2007/06/04
c
c Modified to enable the new ocean model grid.
c SJP 2007/06/01
c
c Added IMPLICIT NONE statement, plus variable declarations as necessary.
c SJP 2007/05/30
c
c Transferred COMMON block declarations to separate header files.
c SJP 2007/05/29
c
c Replaced calls to OGET/OPUT with data access via COMMON block /orestart/
c SJP 2003/09/02
c
c Loop 230 restructured to avoid array bounds violations.
c SJP 2003/05/02
c
c Calls to OFIND removed, as they did nothing. Parameter definitions moved
c to include file OPARAMS.f.
c SJP 2003/04/29
c
c $Log: relax.f,v $
c Revision 1.5  1997/12/19 01:25:42  ldr
c Changes from ACH for 21 OGCM levels and optinal GM mixing scheme...
c
c Revision 1.4  1994/03/30  12:34:56  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c
c Revision 1.3  93/12/17  15:33:41  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.2  93/10/05  13:07:27  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c
c Revision 1.1  93/02/05  16:22:25  ldr
c Initial revision
c
      SUBROUTINE RELAX
C
C=======================================================================
C                                                                    ===
C  RELAX TAKES AS INPUT THE VORTICITY DRIVING FUNCTION COMPUTED IN   ===
C        "CLINIC" (ZTD) AND, USING SEQUENTIAL OVER-RELAXATION,       ===
C        SOLVES THE LAPLACIAN EQUATION FOR THE EXTERNAL MODE OF      ===
C        VELOCITY IN TERMS OF A MASS TRANSPORT STREAM FUNCTION (P).  ===
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
C
C---------------------------------------------------------------------
C  DEMENSION AND EQUIVALENCE LOCAL DATA
C---------------------------------------------------------------------
C
      integer ii, jj, ismask, i, j, luptdb, luptd, l, is, ie, isle,
     &        js, je
      real ptdb, h, cof, cofis, cfn, cfs, cfe, cfw, ptd, res, cpf,
     &     fx, test1, test2, fxa, fxb, fxc, fxd, crtp, resmax, resis
      DIMENSION PTDB(IMT,JMT),H(IMT,JMT),ISMASK(IMT,JMT),COF(IMT,JMT),
     & COFIS(NIEVEN),CFN(IMT,JMT),CFS(IMT,JMT),CFE(IMT,JMT),
     & CFW(IMT,JMT),PTD(IMT,JMT),RES(IMT,JMT),CPF(IMT)
      EQUIVALENCE (RES,PTDB),(PTD,H)
C
C---------------------------------------------------------------------
C  BEGIN EXECUTABLE CODE
C---------------------------------------------------------------------
C
C=======================================================================
C  BEGIN INTRODUCTORY SECTION TO PREPARE FOR THE RELAXATION  ===========
C=======================================================================
C
C---------------------------------------------------------------------
C  INITIATE READIN OF RELAXATION SOLUTION OF 2 TIMESTEPS BACK
C  FOR THE PURPOSE OF COMPUTING AN INITIAL GUESS FOR THE PRESENT PASS
C  (INPUT UNIT ALTERNATES ON TSTEPS BETWEEN 5 & 6)
C---------------------------------------------------------------------
C
      LUPTDB=6-MOD(ITT,2)
      LUPTD =5+MOD(ITT,2)
C
C---------------------------------------------------------------------
C  INITIALIZE THE WORK AREA
C---------------------------------------------------------------------
C
      FX=0.0
      DO 90 J=1,JMT
      DO 90 I=1,IMT
        CFN(I,J)=FX
        CFS(I,J)=FX
        CFE(I,J)=FX
        CFW(I,J)=FX
        COF(I,J)=FX
        ISMASK(I,J)=0
 90   CONTINUE
C
C---------------------------------------------------------------------
C  COMPLETE READIN OF RELAXATION SOLUTION OF 2 TIMESTEPS BACK AND
C  INITIATE READIN OF RELAXATION SOLUTION OF PREVIOUS TIMESTEP
C  (INPUT UNIT ALTERNATES ON TIMESTEPS BETWEEN 6 & 5)
C---------------------------------------------------------------------
C
      do jj = 1, jmt
        do ii = 1, imt
          ptdb(ii, jj) = ptd2(ii, jj, luptdb-4)
        end do
      end do
C
C---------------------------------------------------------------------
C  FORM ISLAND MASK  (ISMASK=0 OVER INTERIOR OCEAN POINTS
C                     ISMASK=1 OVER PERIMETER OCEAN POINTS
C                     ISMASK=2 OVER LAND POINTS)
C---------------------------------------------------------------------
C
      DO  95 I=2,IMU
      DO  95 J=2,JMTM1
        TEST1=HR(I,J)*HR(I-1,J)*HR(I,J-1)*HR(I-1,J-1)
        TEST2=HR(I,J)+HR(I-1,J)+HR(I,J-1)+HR(I-1,J-1)
        IF(TEST1.EQ.0.0) ISMASK(I,J)=1
        IF(TEST2.EQ.0.0) ISMASK(I,J)=2
  95  CONTINUE

# 167


c...  Set IMSASK equal to zero over gridpoints (43:44, 55:56), so that these
c...  points are not included when calculating the line integral around
c...  Australia/New Guinea. This is necessary not because of a bug in the
c...  model as such, but because the method that it uses to define islands is
c...  excessively crude.
c...
c...  Note that any change to the positions of the coastlines on the new ocean
c...  model grid may require that any or all of the following changes be made
c...  to this section of code:
c...
c...  - that the following loop be modified
c...  - that the following loop be removed entirely
c...  - that additional loops be added, to mask out other regions
c...
c...  There is considerable scope to improve the model's treatment of islands,
c...  and to make it more user-friendly!
c...
c...  SJP 2007/06/11

CSJP      do j = 55, 56
CSJP        do i = 43, 44
CSJP          ismask(i, j) = 0
CSJP        end do
CSJP      end do



C---------------------------------------------------------------------
C  CALCULATE DEPTH FIELD FROM ITS RECIPROCAL
C  (..NOTE.. H IS EQUIVALENCED WITH PTD)
C---------------------------------------------------------------------
C
C  COMPUTE RECIPROCALS WITH AN EPSILON ADDED TO AVOID ZERO DIVIDE
C
      ! computes water column depths from HR(I,J), which is the inverse
      ! of depth in cm. It sets all land points where HR(I,J) equals
      ! zero to 1E+20
      FXA=1.0
      FXB=1.E-20
      DO 110 J=2,JMTM1
      DO 110 I=1,IMT
        H(I,J)=FXA/(HR(I,J)+FXB)
 110  CONTINUE
      
C
C  THEN RESET H OVER LAND TO ZERO
C
      FX=0.0
      DO 112 J=1,JMT
      DO 112 I=1,IMT
        IF(HR(I,J).EQ.FX)H(I,J)=FX
 112  CONTINUE
C
C  THEN SET CYCLIC BOUNDARY CONDITIONS ON H
C
      DO 117 J=1,JMT
        H(1,J)=H(IMUM1,J)
 117  CONTINUE
C
C---------------------------------------------------------------------
C  GENERATE ARRAYS OF COEFFICIENTS FOR RELAXATION
C---------------------------------------------------------------------
C
C  1ST, COMPUTE COEFFICIENTS OF THE LAPLACIAN STAR
C  (HOLD NON-INTERIOR POINTS TO ZERO USING START AND END INDICES)
C
      DO 130 J=3,JMTM1
      DO 130 L=1,LSEG
        IS=ISZ(J,L)
        IF(IS.EQ.0) GO TO 130
        IE=IEZ(J,L)
        FXA=2.0*CSTR(J)*CSTR(J)
        FXB=2.0*CS(J  )*CSTR(J)*DYTR(J)*DYUR(J  )
        FXC=2.0*CS(J-1)*CSTR(J)*DYTR(J)*DYUR(J-1)
        FX=1.0
        DO 120 I=IS,IE
          CFN(I,J)=FXB/(H(I-1,J)+H(I,J))
          CFS(I,J)=FXC/(H(I-1,J-1)+H(I,J-1))
          CFE(I,J)=FXA*DXUR(I)*DXTR(I)/(H(I,J)+H(I,J-1))
          CFW(I,J)=FXA*DXUR(I-1)*DXTR(I)/(H(I-1,J)+H(I-1,J-1))
          CPF(I)=FX/(CFN(I,J)+CFS(I,J)+CFE(I,J)+CFW(I,J))
 120    CONTINUE
C
C  2ND, AUGMENT COEFFICIENTS FOR IMPLICIT TREATMENT OF CORIOLIS TERM
C
        IF(ACORF.NE.0.0) THEN
          FX=-C2DTSF*ACORF*CSTR(J)*DYTR(J)*OMEGA
          DO 125 I=IS,IE
            CFN(I,J)=CFN(I,J)+(HR(I,J  )-HR(I-1,J  ))*SINE(J  )
     &                        *FX*DXTR(I)
            CFS(I,J)=CFS(I,J)-(HR(I,J-1)-HR(I-1,J-1))*SINE(J-1)
     &                        *FX*DXTR(I)
            CFE(I,J)=CFE(I,J)-(HR(I  ,J)*SINE(J)-HR(I  ,J-1)*SINE(J-1))
     &                        *FX*DXTR(I)
            CFW(I,J)=CFW(I,J)+(HR(I-1,J)*SINE(J)-HR(I-1,J-1)*SINE(J-1))
     &                        *FX*DXTR(I)
 125      CONTINUE
        ENDIF
C
C  3RD, NORMALIZE COEFFICIENTS AND FORCING TERM FOR EFFICIENCY
C
        DO 128 I=IS,IE
          CFN(I,J)=CFN(I,J)*CPF(I)
          CFS(I,J)=CFS(I,J)*CPF(I)
          CFE(I,J)=CFE(I,J)*CPF(I)
          CFW(I,J)=CFW(I,J)*CPF(I)
          ZTD(I,J)=ZTD(I,J)*CPF(I)
 128    CONTINUE
 130  CONTINUE
C
C  4TH, COMPUTE COEFFICIENTS ON ISLAND PERIMETER POINTS
C
      DO 180 ISLE=1,NISLE
        COFIS(ISLE)=0.0
        IS=ISIS(ISLE)
        IE=IEIS(ISLE)
        JS=JSIS(ISLE)
        JE=JEIS(ISLE)
        DO 175 J=JS,JE
          FXA=2.0*CSTR(J)*CSTR(J)
          FXB=2.0*CS(J  )*DYUR(J  )*DYTR(J)*CSTR(J)
          FXC=2.0*CS(J-1)*DYUR(J-1)*DYTR(J)*CSTR(J)
          FXD=-C2DTSF*ACORF*CSTR(J)*DYTR(J)*OMEGA
        DO 175 I=IS,IE
          IF(ISMASK(I,J).NE.1) GO TO 175
          IF(HR(I-1,J  ).NE.0.0 .OR. HR(I  ,J  ).NE.0.0)
     &      CFN(I,J)=FXB/(PTD(I-1,J)+PTD(I,J))
     &         +FXD*DXTR(I)*(HR(I,J)-HR(I-1,J))*SINE(J)
          IF(HR(I-1,J-1).NE.0.0 .OR. HR(I  ,J-1).NE.0.0)
     &      CFS(I,J)=FXC/(PTD(I-1,J-1)+PTD(I,J-1))
     &         -FXD*DXTR(I)*(HR(I,J-1)-HR(I-1,J-1))*SINE(J-1)
          IF(HR(I  ,J  ).NE.0.0 .OR. HR(I  ,J-1).NE.0.0)
     &      CFE(I,J)=FXA*DXTR(I)*DXUR(I)/(PTD(I,J)+PTD(I,J-1))
     &         -FXD*DXTR(I)*(HR(I,J)*SINE(J)-HR(I,J-1)*SINE(J-1))
          IF(HR(I-1,J  ).NE.0.0 .OR. HR(I-1,J-1).NE.0.0)
     &      CFW(I,J)=FXA*DXTR(I)*DXUR(I-1)/(PTD(I-1,J)+PTD(I-1,J-1))
     &         +FXD*DXTR(I)*(HR(I-1,J)*SINE(J)-HR(I-1,J-1)*SINE(J-1))
          COF(I,J)=1.0/(CFN(I,J)+CFS(I,J)+CFE(I,J)+CFW(I,J))
          CFN(I,J)=CFN(I,J)*COF(I,J)
          CFS(I,J)=CFS(I,J)*COF(I,J)
          CFE(I,J)=CFE(I,J)*COF(I,J)
          CFW(I,J)=CFW(I,J)*COF(I,J)
          ZTD(I,J)=ZTD(I,J)*COF(I,J)
          COF(I,J)=CST(J)*DXT(I)*DYT(J)/COF(I,J)
          COFIS(ISLE)=COFIS(ISLE)+COF(I,J)
 175    CONTINUE
        COFIS(ISLE)=1./COFIS(ISLE)
 180  CONTINUE
C
C  FINALLY, MULTIPLY THE COEFFICIENTS BY THE OVERRELAXATION FACTOR
C
      DO 182 J=1,JMT
      DO 182 I=1,IMT
        CFN(I,J)=CFN(I,J)*SORF
        CFS(I,J)=CFS(I,J)*SORF
        CFE(I,J)=CFE(I,J)*SORF
        CFW(I,J)=CFW(I,J)*SORF
 182  CONTINUE
C
C---------------------------------------------------------------------
C  COMPUTE A FIRST GUESS FOR THE RELAXATION BY EXTRAPOLATING THE TWO
C  PREVIOUS SOLUTIONS FORWARD IN TIME.
C---------------------------------------------------------------------
C
C  1ST, COMPLETE READIN OF RELAXATION SOLUTION OF PREVIOUS TIMESTEP
C
      do jj = 1, jmt
        do ii = 1, imt
          ptd(ii, jj) = ptd2(ii, jj, luptd-4)
        end do
      end do
C
C  2ND, PERFORM TIME EXTRAPOLATION (ACCOUNTING FOR MIXING TIMESTEP)
C
      FXA=1.0
      FXB=2.0
      IF(MIX.EQ.1.OR.MXP.EQ.1) FXA=0.5
      DO 135 J=1,JMT
      DO 135 I=1,IMT
        PTD(I,J)=FXA*(FXB*PTD(I,J)-PTDB(I,J))
 135  CONTINUE
C
C  COMPUTE CRITERION FOR CONVERGENCE OF RELAXATION AND SET RESIDUALS 0
C
      CRTP=CRITF*FXA*SORF
      FX=0.0
      DO 140 J=1,JMT
      DO 140 I=1,IMT
        RES(I,J)=FX
 140  CONTINUE
C
C=======================================================================
C  END INTRODUCTORY SECTION  ===========================================
C=======================================================================
C
C=======================================================================
C  BEGIN SECTION TO DO THE RELAXATION  =================================
C=======================================================================
C
      MSCAN=0
 300  MSCAN=MSCAN+1
C
C---------------------------------------------------------------------
C  COMPUTE ENTIRE FIELD OF RESIDUALS AS IN SIMULTANEOUS RELAXATION
C---------------------------------------------------------------------
C
      DO J=3,JSCAN
        DO I=2,IMTM1
          RES(I,J)=CFN(I,J)*PTD(I,J+1)
     &            +CFS(I,J)*PTD(I,J-1)
     &            +CFE(I,J)*PTD(I+1,J)
     &            +CFW(I,J)*PTD(I-1,J)
     &            -SORF*(PTD(I,J)+ZTD(I,J))
        END DO
        RES(1,J)=RES(IMTM1,J)
        RES(IMT,J)=RES(2,J)
      END DO
C
C---------------------------------------------------------------------
C  RESET RESIDUALS TO ZERO OVER LAND
C---------------------------------------------------------------------
C
      FX=0.0
      DO 224 J=3,JSCAN
      DO 220 I=1,IMT
        CPF(I)=FX
 220  CONTINUE
      DO 221 L=1,LSEG
        IS=ISZ(J,L)
        IF(IS.EQ.0) GO TO 222
        IE=IEZ(J,L)
      DO 221 I=IS,IE
        CPF(I)=RES(I,J)
 221  CONTINUE
 222  CONTINUE
      DO 223 I=1,IMT
        RES(I,J)=CPF(I)
 223  CONTINUE
 224  CONTINUE
C
C---------------------------------------------------------------------
C  SET CYCLIC BOUNDARY CONDITIONS ON THE RESIDUALS
C---------------------------------------------------------------------
C
      DO 226 J=3,JSCAN
        RES(1,J)=RES(IMUM1,J)
        RES(IMT,J)=RES(2,J)
 226  CONTINUE
C
C---------------------------------------------------------------------
C  PERFORM CORRECTION ON SOUTHERN POINT TO YIELD SEQUENTIAL RELAXATION
C---------------------------------------------------------------------
C
      DO 247 J=3,JSCAN
        DO 254 L=1,LSEG
          IS=ISZ(J,L)
          IF(IS.EQ.0) GO TO 255
          IE=IEZ(J,L)
        DO 254 I=IS,IE
          RES(I,J)=RES(I,J)+CFS(I,J)*RES(I,J-1)
 254    CONTINUE
 255    CONTINUE
C
C---------------------------------------------------------------------
C  PERFORM CORRECTION ON WESTERN POINT TO YIELD SEQUENTIAL RELAXATION
C---------------------------------------------------------------------
C
        DO 245 L=1,LSEG
          IS=ISZ(J,L)
          IF(IS.EQ.0) GO TO 247
          IE=IEZ(J,L)
          DO 237 I=IS,IE
            RES(I,J)=RES(I,J)+CFW(I,J)*RES(I-1,J)
 237      CONTINUE
 245    CONTINUE
 247  CONTINUE
C
C---------------------------------------------------------------------
C  MAKE A CORRECTION TO PTD BASED ON THE RESIDUALS
C---------------------------------------------------------------------
C
      DO 235 J=3,JSCAN
      DO 235 I=1,IMT
        PTD(I,J)=PTD(I,J)+RES(I,J)
 235  CONTINUE
C
C---------------------------------------------------------------------
C  FIND THE MAXIMUM ABSOLUTE RESIDUAL TO DETERMINE CONVERGENCE
C---------------------------------------------------------------------
C
      RESMAX=0.0
      DO 253 J=3,JSCAN
      DO 253 I=2,IMTM1
        RESMAX=MAX(ABS(RES(I,J)),RESMAX)
 253  CONTINUE
C
C---------------------------------------------------------------------
C  DO HOLE RELAXATION FOR EACH ISLAND
C---------------------------------------------------------------------
C
      FX=0.0
      DO 260 ISLE=1,NISLE
        IS=ISIS(ISLE)
        IE=IEIS(ISLE)
        JS=JSIS(ISLE)
        JE=JEIS(ISLE)
        RESIS=FX
        DO 262 J=JS,JE
        DO 262 I=IS,IE
          IF(ISMASK(I,J).EQ.1) THEN
            RESIS=RESIS+(CFN(I,J)*PTD(I  ,J+1)
     &                  +CFS(I,J)*PTD(I  ,J-1)
     &                  +CFE(I,J)*PTD(I+1,J  )
     &                  +CFW(I,J)*PTD(I-1,J  )
     &                  -SORF*(PTD(I,J)+ZTD(I,J)))*COF(I,J)
          ENDIF
262     CONTINUE
C
C  NORMALIZE THE ISLAND RESIDUAL AND UPDATE THE MAXIMUM
C  ABSOLUTE RESIDUAL OF THE RELAXATION IF NECESSARY
C
        RESIS=RESIS*COFIS(ISLE)
        RESMAX=MAX(ABS(RESIS),RESMAX)
C
C  MAKE A CORRECTION TO PTD OVER THE ISLAND AND ITS PERIMETER POINTS
C
        DO 250 J=JS,JE
        DO 250 I=IS,IE
          IF(ISMASK(I,J).GE.1) THEN
            PTD(I,J)=PTD(I,J)+RESIS
          ENDIF
 250    CONTINUE
 260  CONTINUE
C
C---------------------------------------------------------------------
C  SET CYCLIC BOUNDARY CONDITION
C---------------------------------------------------------------------
C
      DO 272 J=1,JMT
        PTD(1,J)=PTD(IMUM1,J)
        PTD(IMU,J)=PTD(2,J)
 272  CONTINUE
C
C---------------------------------------------------------------------
C  TEST MAXIMUM RESIDUAL FOR CONVERGENCE OF THE RELAXATION.
C  IF NOT CONVERGED, PROCEED WITH ANOTHER SCAN.
C  (..NOTE.. IF THE NUMBER OF SCANS REACHES MXSCAN, LEAVE THE LOOP)
C---------------------------------------------------------------------
C
      IF(RESMAX.GE.CRTP .AND. MSCAN.LT.MXSCAN) GO TO 300

      if (mscan .eq. mxscan) then
        write (*, *)
        write (*, *)
     &    "***  WARNING: Maximum number of relaxation scans reached"
        write (*, *) "***"
        write (*, *) "***  Number of scans performed = ", mscan
        write (*, *) "***"
        write (*, *) "***  Consider increasing the value of MXSCAN"
        write (*, *)
      end if

C
C=======================================================================
C  END OF THE RELAXATION  ==============================================
C=======================================================================
C
C---------------------------------------------------------------------
C  UPDATE THE STREAM FUNCTION BASED UPON THE RELAXATION SOLUTION
C---------------------------------------------------------------------
C
      if (robert_time_filter) then

        do j = 1, jmt
          do i = 1, imt
            ptdb(i, j) = pb(i, j) + ptd(i, j)
            pb(i, j) = pnu2m * p(i, j) + pnu * (pb(i, j) + ptdb(i, j))
            p(i, j) = ptdb(i,j)
          end do
        end do

      else

      IF(MXP.EQ.0) THEN
      DO 340 J=1,JMT
      DO 340 I=1,IMT
        PTDB(I,J)=PB(I,J)+PTD(I,J)
        PB(I,J)=P(I,J)
        P(I,J)=PTDB(I,J)
 340  CONTINUE
      ELSE
      DO 342 J=1,JMT
      DO 342 I=1,IMT
        P(I,J)=PB(I,J)+PTD(I,J)
 342  CONTINUE
      ENDIF

      end if
C
C---------------------------------------------------------------------
C  SAVE PTD TO COMPUTE 1ST GUESS FOR RELAXATION NEXT TIMESTEP
C  (..NOTE.. ON 1ST PASS OF EULER BACKWARD TIMESTEP, BYPASS THIS
C            SAVE, SINCE IT WILL BE DONE ON THE 2ND PASS)
C  (..NOTE.. ON A MIXING TIMESTEP, ALTER PTD TO BE CONSISTENT WITH
C            NORMAL, LEAP-FROG STEPPING)
C---------------------------------------------------------------------
C
      IF(MIX.EQ.1 .AND. EB) RETURN
      IF(MXP.NE.0. OR. MIX.NE.0) THEN
        FX=2.0
        DO 350 J=1,JMT
        DO 350 I=1,IMT
          PTD(I,J)=FX*PTD(I,J)
 350    CONTINUE
      ENDIF
      do jj = 1, jmt
        do ii = 1, imt
          ptd2(ii, jj, luptdb-4) = ptd(ii, jj)
        end do
      end do
      CONTINUE
      RETURN
      END
