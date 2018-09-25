c Add calls to NINT to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/15
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/16
c
c Added IMPLICIT NONE statement, plus variable declarations as necessary.
c SJP 2007/05/30
c
c Remove unnecessary variable declarations.
c SJP 2004/01/06
c
c Parameter definitions moved to include file OPARAMS.f.
c SJP 2003/04/29
c
c $Log: filter.f,v $
c Revision 1.6  1997/12/19 01:25:35  ldr
c Changes from ACH for 21 OGCM levels and optinal GM mixing scheme...
c
c Revision 1.5  1997/03/06  23:33:59  ldr
c Corrections from HBG to the recent tidy-ups to ocean routines.
c
c Revision 1.4  1996/10/24  01:02:45  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.3  1994/03/30  12:34:22  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c
c Revision 1.2  93/12/17  15:32:32  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.1  93/02/05  16:22:05  ldr
c Initial revision
c 
      SUBROUTINE FILTER(S,IM,MM,N,ISS)
C
C=======================================================================
C                                                                    ===
C  FILTER FOURIER ANALYSES THE ARRAYS OF VARIOUS                     ===
C         PHYSICAL QUANTITIES, THEN TRUNCATES THE SERIES AND         ===
C         RESYNTHESIZES THE FILTERED QUANTITIES WHERE:               ===
C             S  =THE STRING TO BE FILTERED                          ===
C             IM =THE LENGTH OF S                                    ===
C             MM =1 (COSINE SERIES, DERIV AT BNDRY PTS=0)            ===
C                =2 (  SINE SERIES,          BNDRY PTS=0)            ===
C                =3 (FULL SERIES, CYCLIC)                            ===
C             N  =NUMBER OF WAVES TO KEEP                            ===
C             ISS=0 (CANT USE FOURIER COEFS FROM PREVIOUS CALL)      ===
C             ISS>0 (CAN  USE FOURIER COEFS FROM PREVIOUS CALL)      ===
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
C
C  FILTER T TO YIELD EQUIV DX AT J=JFT0 FROM J=JFRST TO JFT1 AND
C       J=JFT2 TO JMTM1
C  FILTER U TO YIELD EQUIV DX AT J=JFU0 FROM J=JFRST TO JFU1 AND
C       J=JFU2 TO JMTM2
C
C---------------------------------------------------------------------
C  DEFINE LOCAL DATA AND DIMENSION ARGUMENT ARRAYS
C---------------------------------------------------------------------
C
      integer imtx2, ni, imtd2, lqmsum, lhsum, imtx4, imtx8, imtimt,
     &        imp1x2, icbase, idbase, ind, indx, im, mm, n, iss,
     &        imsave, i, ibase, jbase, imm1, imqc, nmax, nmaxp1,
     &        lcy, lh, lhm1, lqm, lcyp1, imx4, imx8, maxind, ncyc,
     &        maxndx, npwr, np, j, ioff1, ioff2, joff, ioff
      real temp, denmsv, cosnpi, circle, cof, cosine, ftarr, denom,
     &     s, sprime, pi, fimr, c1, c2, fnorm, ssum, fim, stemp,
     &     fxa, fact1, fact2, genadj, fxb, ssm
      PARAMETER (IMTX2=IMT*2,NI=IMT)
      PARAMETER (IMTD2=IMT/2,LQMSUM=IMTD2*(IMT-IMTD2),LHSUM=IMT*IMTP1/2)
      PARAMETER (IMTX4=IMT*4,IMTX8=IMT*8,IMTIMT=IMT*IMT)
      PARAMETER (IMP1X2=IMTP1*2)
C
C     COSSAV MUST REMAIN FULL PRECISION IF MOST OF FILTER IS MADE HALF-P
      REAL COSSAV
C
      DIMENSION ICBASE(IMTP1),IDBASE(IMTP1),IND(IMTX8),TEMP(IMTX4)
      DIMENSION COSSAV(LQMSUM),DENMSV(LHSUM),COSNPI(IMT)
      DIMENSION CIRCLE(4)
      DIMENSION INDX(IMTX8),COF(IMTX8)
      DIMENSION COSINE(IMTX8),FTARR(IMTIMT)
      DIMENSION DENOM(IMTX4)
      DIMENSION S(IMT),SPRIME(IMT)
C
C
      LOGICAL INITDN
      DATA INITDN/.FALSE./
C
      DATA PI/3.141592653589793/, CIRCLE/0.,-1.,0.,1./
C
C********************************  INSERT CRAY FIXES  ******************
      SAVE DENMSV,COSSAV,COSNPI,IDBASE,ICBASE,IND,INITDN,FTARR
C---------------------------------------------------------------------
C  BEGIN EXECUTABLE CODE
C---------------------------------------------------------------------
C
      IF (IM.LT.1 .OR. MM.LT.1 .OR. MM.GT.3 .OR. N.LT.0 .OR. ISS.LT.0)
     & GOTO 4000
C
      IF (INITDN) GOTO 90
C
C  THIS SECTION SETS UP TABLES FOR FILTER; IT MUST BE CALLED ONCE PER
C    EXECUTION OF OCEAN
C
C  NOTE: LQMSUM IS THE SUM OF (IM-1)/2 FOR IM=1,IMTP1
C        LHSUM IS THE SUM OF IM-1 FOR IM=1,IMTP1
C
      IMSAVE=IM
C
C  ASSEMBLE INDEX ARRAY
C
      DO 10 I=1,IMTX8
      IND(I)=I
 10   CONTINUE
C
C  CALCULATE AND SAVE ALL COSINES WHICH WILL BE NEEDED
C
      IBASE=0
      JBASE=0
      DO 50 IM=1,IMTP1
      FIMR=1.0/FLOAT(IM)
C
      IMM1=IM-1
      IF (IMM1.EQ.0) GOTO 25
      DO 20 I=1,IMM1
      DENMSV(IBASE+I)=1.0/(1.0-COS(PI*FLOAT(I)*FIMR))
 20   CONTINUE
 25   IDBASE(IM)=IBASE
      IBASE=IBASE+IMM1
C
      IMQC=(IM-1)/2
      IF (IMQC.EQ.0) GOTO 35
      DO 30 I=1,IMQC
      COSSAV(JBASE+I)=COS(PI*FLOAT(I)*FIMR)
 30   CONTINUE
 35   ICBASE(IM)=JBASE
      JBASE=JBASE+IMQC
C
 50   CONTINUE
C
C  CALCULATE ADJUSTMENTS FOR GENERAL FOURIER CASE IF IM=2*N
C
      DO 60 IM=1,IMT
      COSNPI(IM)=CIRCLE(MOD(IM-1,4)+1)
 60   CONTINUE
      INITDN=.TRUE.
      IM=IMSAVE
C
C
C  CALCULATE SOME USEFUL CONSTANTS
C
 90   IF(MM.EQ.2 .AND. N.EQ.0) THEN
        DO 92 I=1,IM
          S(I)=0.0
 92     CONTINUE
        GO TO 3950
      ENDIF
      IF (MM.EQ.1) THEN
        NMAX=N-1
      ELSE
        NMAX=N
      ENDIF
      NMAXP1=NMAX+1
      C1=0.5*FLOAT(NMAX)+0.25
      C2=FLOAT(NMAX)+0.5
C
      IF (MM.EQ.2) THEN
        LCY=2*(IM+1)
        FNORM=2.0/FLOAT(IM+1)
      ELSE
        LCY=2*IM
        FNORM=2.0/FLOAT(IM)
      ENDIF
      LH=LCY/2
      LHM1=LH-1
      LQM=(LH-1)/2
      LCYP1=LCY+1
C
      IMX4=IM*4
      IMX8=IM*8
C
C  AVERAGE INCOMING ARRAY
C
      SSUM=0.0
      DO 100 I=1,IM
 100  SSUM=SSUM + S(I)
C
C     MM = 1  DERIVATIVE MUST BE ZERO AT BOUNDARIES (COSINE)
C     MM = 2  VALUE MUST BE ZERO AT BOUNDARIES (SINE)
C     MM = 3  CYCLIC BOUNDARY CONDITIONS (GENERAL FOURIER SERIES)
C
      FIM=FLOAT(IM)
      FIMR=1.0/FIM
      STEMP=SSUM*FIMR
      IF (N.GT.1 .OR. MM.NE.1) GO TO 400
      DO 300 I=1,IM
 300  S(I)=STEMP
      GO TO 3950
 400  CONTINUE
      IF(MM.NE.2) THEN
        DO 450 I=1,IM
        S(I)=S(I)-STEMP
 450    CONTINUE
      ENDIF
      IF (ISS.GT.0) GO TO 3000
C
C
C  ASSEMBLE APPROPRIATE 1-CYCLE (2*PI) COSINE ARRAY
C
C  USE STORED 1/4 CYCLE TO CALCULATE FIRST 1/2 CYCLE
      JBASE=ICBASE(LH)
      DO 700 I=1,LQM
      COSINE(I)=COSSAV(JBASE+I)
 700  CONTINUE
      DO 701 I=1,LQM
 701  COSINE(LH-I)=-COSSAV(JBASE+I)
C  FILL IN COS(PI/2) IF LH IS EVEN
      IF (2*(LQM+1).EQ.LH) COSINE(LQM+1)=0.0
C  FILL IN COS(PI) IN ANY CASE
      COSINE(LH)=-1.0
C  FILL IN REST OF CYCLE
      DO 710 I=1,LH
      COSINE(LH+I)=-COSINE(I)
 710  CONTINUE
C
C  ASSEMBLE DENOMINATOR ARRAY
C
      IBASE=IDBASE(LH)
      FXA=0.25
      DO 720 I=1,LHM1
      DENOM(I)=FXA*DENMSV(IBASE+I)
 720  CONTINUE
      DENOM(LH)=0.125
      DO 721 I=1,LHM1
      TEMP(I)=DENOM(LH-I)
 721  CONTINUE
      DO 722 I=1,LHM1
 722  DENOM(LH+I)=TEMP(I)
      DENOM(LCY)=0.0
      DO 730 I=LCYP1,IMX4
      DENOM(I)=DENOM(I-LCY)
 730  CONTINUE
C
C  ASSEMBLE APPROPRIATE SUBSCRIPT ARRAYS
C
C  CALCULATE NEEDED INDICES
C
      IF (MM.EQ.3) THEN
        FACT1=2*NMAX
        FACT2=2*NMAXP1
      ELSE
        FACT1=NMAX
        FACT2=NMAXP1
      ENDIF
      DO 740 I=1,IMX4
      INDX(I)=NINT(IND(I)*FACT1)
 740  CONTINUE
      DO 741 I=1,IMX4
 741  INDX(IMX4+I)=NINT(IND(I)*FACT2)
C  CALCULATE PARAMETERS FOR REDUCING INDICES
      MAXIND=NINT(IMX4*FACT2)
      NCYC=(MAXIND-1)/LCY + 1
      MAXNDX=LCY
      IF (MAXNDX.GE.MAXIND) GOTO 790
      DO 750 NPWR=1,NCYC+2
      MAXNDX=2*MAXNDX
      IF (MAXNDX.GE.MAXIND) GOTO 760
 750  CONTINUE
      STOP 'ERROR 1F'
 760  DO 770 NP=1,NPWR
      MAXNDX=MAXNDX/2
      DO 765 I=1,IMX8
        IF(INDX(I).GT.MAXNDX) INDX(I)=INDX(I)-MAXNDX
 765  CONTINUE
 770  CONTINUE
 790  CONTINUE
C
C  GATHER COEFFICIENTS
C
      DO 810 J=1,IMX8
      COF(J)=COSINE(INDX(J))
 810  CONTINUE
C
C
C  ASSEMBLE TRANSFORMATION ARRAY WHICH WILL FILTER S
C
C
      IF(MM.EQ.1) THEN
C  COSINE TRANSFORM
      IOFF1=LCY
      IOFF2=LCY+IMX4
      FXA=0.5
      DO 1200 J=1,IM
      JOFF=(J-1)*IMT
      DO 1100 I=1,IM
      FTARR(JOFF+I)=(COF(I-J+IOFF1)-COF(I-J+IOFF2))*DENOM(I-J+IOFF1)
     &             +(COF(I+J-1)-COF(IMX4+I+J-1))*DENOM(I+J-1) - FXA
 1100 CONTINUE
 1200 CONTINUE
      DO 1201 J=1,IM
 1201 FTARR(J*IMTP1-IMT)=FTARR(J*IMTP1-IMT)+C1
      ELSE IF(MM.EQ.2) THEN
C
C  SINE TRANSFORM
      IOFF1=LCY
      IOFF2=LCY+IMX4
      DO 1500 J=1,IM
      JOFF=(J-1)*IMT
      DO 1400 I=1,IM
      FTARR(JOFF+I)=(COF(I-J+IOFF1)-COF(I-J+IOFF2))*DENOM(I-J+IOFF1)
     &             -(COF(I+J)-COF(IMX4+I+J))*DENOM(I+J)
 1400 CONTINUE
 1500 CONTINUE
      DO 1501 J=1,IM
 1501 FTARR(J*IMTP1-IMT)=FTARR(J*IMTP1-IMT)+C1
      ELSE IF(MM.EQ.3) THEN
C
C  GENERAL FOURIER TRANSFORM
      IF (2*N.EQ.IM) THEN
        GENADJ=0.5
      ELSE
        GENADJ=0.0
      ENDIF
      IOFF1=LCY
      IOFF2=LCY+IMX4
      FXA=2.0
      FXB=0.5
      DO 1800 J=1,IM
      JOFF=(J-1)*IMT
      DO 1700 I=1,IM
      FTARR(JOFF+I)=
     & (FXA*(COF(I-J+IOFF1)-COF(I-J+IOFF2)))*DENOM(2*I-2*J+IOFF1)
     & -FXB - GENADJ*COSNPI(I)*COSNPI(J)
 1700 CONTINUE
 1800 CONTINUE
      DO 1801 J=1,IM
 1801 FTARR(J*IMTP1-IMT)=FTARR(J*IMTP1-IMT)+C2
      ENDIF
C
C  FILTER S
C
 3000 DO 3010 I=1,IM
      SPRIME(I)=0.
 3010 CONTINUE
      DO 3100 I=1,IM
      IOFF=(I-1)*IMT
      DO 3100 J=1,IM
C  NOTE THAT FTARR(J,I)=FTARR(I,J), SO FOLLOWING IS LEGAL
      SPRIME(J)=SPRIME(J)+S(I)*FTARR(IOFF+J)
 3100 CONTINUE
      DO 3110 I=1,IM
      SPRIME(I)=FNORM*SPRIME(I)
 3110 CONTINUE
      IF(MM.EQ.2) THEN
        DO 3150 I=1,IM
          S(I)=SPRIME(I)
 3150   CONTINUE
        GO TO 3950
      ENDIF
C
      SSM=0.0
      DO 3800 I=1,IM
      SSM=SSM+SPRIME(I)
 3800 CONTINUE
      SSM=(SSUM-SSM)*FIMR
      DO 3900 I=1,IM
      S(I)=SSM+SPRIME(I)
 3900 CONTINUE
 3950 CONTINUE
      RETURN
 4000 PRINT 4001, IM,MM,N,ISS
 4001 FORMAT (' BAD ARGUMENT(S) IN CALL TO FILTER'/' IM,MM,N,ISS = ',
     &  4I10)
      STOP 'ERROR 2F'
C***********************************************************************
      END
