c Major tidy-up of ocean model source code.
c SJP 2007/06/16
c
c Added IMPLICIT NONE statement, plus variable declarations as necessary.
c SJP 2007/05/30
c
c Parameter definitions moved to include file OPARAMS.f.
c SJP 2003/04/29
c
c $Log: findex.f,v $
c Revision 1.4  1997/12/19 01:25:42  ldr
c Changes from ACH for 21 OGCM levels and optinal GM mixing scheme...
c
c Revision 1.3  1994/03/30  12:34:23  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c
c Revision 1.2  93/12/17  15:32:33  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.1  93/02/05  16:22:06  ldr
c Initial revision
c 
      SUBROUTINE FINDEX(FKXX,JJMAX,KMAX,JF1,JF2,IMAX,ISF,IEF)
C
C=======================================================================
C                                                                    ===
C  FINDEX FINDS AND PRINTS STARTING AND ENDING INDICES               ===
C         FOR FILTERING, WHERE:                                      ===
C             FKXX =FIELD OF MAXIMUM LEVELS FOR THE QUANTITY         ===
C                   BEING FILTERED                                   ===
C             JJMAX=NUMBER OF ROWS TO BE FILTERED                    ===
C             KMAX =MAXIMUM NUMBER OF LEVELS TO BE FILTERED          ===
C             JF1  =LAST ROW IN THE SOUTH TO BE FILTERED             ===
C             JF2  =FIRST ROW IN THE NORTH TO BE FILTERED            ===
C             IMAX =MAXIMUM I INDEX TO BE FILTERED                   ===
C             ISF  =RETURNED VALUES OF STARTING INDICES              ===
C             IEF  =RETURNED VALUES OF ENDING INDICES                ===
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
C  DEFINE LOCAL DATA AND DIMENSION ARGUMENT ARRAYS
C---------------------------------------------------------------------
C
      integer isf, ief, iis, iie, jjmax, kmax, jf1, jf2, imax, jj,
     &        j, k, l, i, lm, llast
      real fkxx
      DIMENSION FKXX(IMT,JMT)
     & ,ISF(JJMAX,LSEGF,KMAX),IEF(JJMAX,LSEGF,KMAX)
     & ,IIS(LSEGF+1),IIE(LSEGF+1)
C
C---------------------------------------------------------------------
C  BEGIN EXECUTABLE CODE
C---------------------------------------------------------------------
C
C FIND START AND END INDICES
C
      JJ = 0
      DO 100 J = JFRST,JMTM1
        IF (J.GT.JF1 .AND. J.LT.JF2) GOTO 100
        JJ = JJ+1
        DO 80 K = 1,KMAX
           DO 30 L = 1,LSEGF+1
              IIS(L) = 0
              IIE(L) = 0
 30        CONTINUE
           L = 1
           IF (FKXX(2,J).GE.K) THEN
              IIS(1) = 2
           ENDIF
           DO 50 I = 2,IMAX-1
              IF (FKXX(I-1,J).LT.K .AND. FKXX(I,J).GE.K) THEN
                 IIS(L) = I
              ENDIF
              IF (FKXX(I,J).GE.K .AND. FKXX(I+1,J).LT.K) THEN
                 IF (I.NE.IIS(L) .OR. (I.EQ.2 .AND. FKXX(1,J).GE.K))THEN
                    IIE(L) = I
                    L = L+1
                 ELSE
                    IIS(L) = 0
                 ENDIF
              ENDIF
 50        CONTINUE
           IF (FKXX(IMAX-1,J).GE.K .AND. FKXX(IMAX,J).GE.K) THEN
              IIE(L) = IMAX-1
              L = L+1
           ENDIF
           LM = L-1
           IF (LM.GT.1) THEN
              IF (IIS(1).EQ.2 .AND. IIE(LM).EQ.IMAX-1
     &                        .AND. FKXX(1,J).GE.K) THEN
                 IIS(1) = IIS(LM)
                 IIE(1) = IIE(1) + IMAX-2
                 IIS(LM) = 0
                 IIE(LM) = 0
                 LM = LM-1
              ENDIF
           ENDIF
           IF (LM .GT. LSEGF) THEN
              PRINT 1000, J, K
              STOP 34
           ENDIF
 1000      FORMAT (1X,'LSEGF NOT LARGE ENOUGH.  J=',I4,'    K=',I3)
           DO 70 L = 1,LSEGF
              ISF(JJ,L,K) = IIS(L)
              IEF(JJ,L,K) = IIE(L)
 70        CONTINUE
 80     CONTINUE
 100  CONTINUE
C
C PRINT THEM
C
      LLAST=LSEGF
      IF (LLAST .GT. 11) LLAST=11
      JJ=JJ+1
      DO 200 J=JMTM1,JFRST,-1
         IF (J.GT.JF1 .AND. J.LT.JF2) GO TO 200
         JJ = JJ-1
         IF (KMAX .GT. 1) THEN
            PRINT 1010,J
            DO 150 K=1,KMAX
               PRINT 1011,K,(ISF(JJ,L,K),IEF(JJ,L,K),L=1,LLAST)
 150        CONTINUE
         ELSE
            PRINT 1011,J,(ISF(JJ,L,1),IEF(JJ,L,1),L=1,LLAST)
         ENDIF
 200  CONTINUE
 1010 FORMAT (/1X,'INDICES FOR ROW ',I3,':')
 1011 FORMAT (1X,I9,3X,11(I5,I4))
      RETURN
      END
