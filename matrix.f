c Major tidy-up of ocean model source code.
c SJP 2007/06/16
c
c Added IMPLICIT NONE statement, plus variable declarations as necessary.
c SJP 2007/05/30
c
c Changed size of ARRAY in order to avoid array bounds violations.
c SJP 2003/05/02
c
c Parameter definitions moved to include file OPARAMS.f.
c SJP 2003/04/29
c
c $Log: matrix.f,v $
c Revision 1.4  1997/12/19 01:25:38  ldr
c Changes from ACH for 21 OGCM levels and optinal GM mixing scheme...
c
c Revision 1.3  1994/03/30  12:34:36  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c
c Revision 1.2  93/12/17  15:33:05  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.1  93/02/05  16:22:20  ldr
c Initial revision
c 
      SUBROUTINE MATRIX(ARRAY,IRDIM,ISTRT,IM,JM,KK,SCALE)
C*********************  RESET VALUES IN SUBROUTINE MATRIX  ************
C
C=======================================================================
C                                                                    ===
C  MATRIX IS A GENERAL TWO-DIMENSIONAL ARRAY PRINTING ROUTINE,       ===
C      WHERE:                                                        ===
C      ARRAY=THE ARRAY TO BE PRINTED                                 ===
C      IRDIM=THE 1ST DIMENSION OF ARRAY                              ===
C      ISTRT=THE 1ST ELEMENT OF THE 1ST DIMENSION TO BE PRINTED      ===
C      IM   =THE LAST ELEMENT OF THE 1ST DIMENSION TO BE PRINTED     ===
C      JM   =THE LAST ELEMENT OF THE 2ND DIMENSION TO BE PRINTED     ===
C            --THE ROWS ARE PRINTED IN DESCENDING ORDER--            ===
C            (IF THIS IS 0, KK IS USED)                              ===
C      KK   =THE LAST ELEMENT OF THE 2ND DIMENSION TO BE PRINTED     ===
C            --THE ROWS ARE PRINTED IN ASCENDING ORDER--             ===
C            (IF THIS IS 0, JM IS USED)                              ===
C      SCALE=A SCALING FACTOR BY WHICH ARRAY IS DIVIDED BEFORE       ===
C            PRINTING.  (IF THIS IS ZERO, NO SCALING IS DONE.)       ===
C            IF SCALE=0, 10 COLUMNS ARE PRINTED ACROSS IN E FORMAT   ===
C            IF SCALE>0, 20 COLUMNS ARE PRINTED ACROSS IN F FORMAT   ===
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
C---------------------------------------------------------------------
C  DEFINE LOCAL DATA
C---------------------------------------------------------------------
C
      integer num, irdim, istrt, im, jm, kk, is, ie, idif, i,
     &        jmorkm, jork, l
      real array, pline, scale, scaler
      DIMENSION ARRAY(IRDIM,JM+KK),NUM(20),PLINE(20)
C
C---------------------------------------------------------------------
C  BEGIN EXECUTABLE CODE
C---------------------------------------------------------------------
C
      IF(SCALE.NE.0) GO TO 500
      DO 251 IS=ISTRT,IM,10
      IE=IS+9
      IF(IE.GT.IM) IE=IM
      IDIF=IE-IS+1
      DO 100 I=1,IDIF
 100  NUM(I)=IS+I-1
      PRINT 9990, (NUM(I),I=1,IDIF)
 9990 FORMAT(10I13,/)
      JMORKM=JM+KK
      DO 252 JORK=1,JMORKM
      IF(JM.NE.0) L=JMORKM+1-JORK
      IF(KK.NE.0) L=JORK
      PRINT 9966, L, (ARRAY(I,L),I=IS,IE)
 252  CONTINUE
      PRINT 9984
 251  CONTINUE
 9966 FORMAT(1X,I2,10(1PE13.5))
 9984 FORMAT(1X)
      RETURN
 500  CONTINUE
      SCALER=1.0/SCALE
      DO 751 IS = ISTRT,IM,12
      IE = IS+11
      IF(IE.GT.IM) IE=IM
      IDIF=IE-IS+1
      DO 600 I=1,IDIF
 600  NUM(I)=IS+I-1
      PRINT 9991, (NUM(I),I=1,IDIF)
 9991 FORMAT(3X,20I6,/)
      JMORKM=JM+KK
      DO 752 JORK=1,JMORKM
      IF(JM.NE.0) L=JMORKM+1-JORK
      IF(KK.NE.0) L=JORK
      DO 753 I=1,IDIF
      PLINE(I)=ARRAY(IS+I-1,L)*SCALER
 753  CONTINUE
      PRINT 9967,L,(PLINE(I),I=1,IDIF)
 752  CONTINUE
      PRINT 9984
 751  CONTINUE
 9967 FORMAT(1X,I3,1X,20F6.2)
      RETURN
      END
