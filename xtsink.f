c Commented out carbonaceous tracers as not used in Mk3L for the moment.
c BP  2010/08/06
c
c $Log: xtsink.f,v $
c Revision 1.5  2000/12/08 03:58:54  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.4  2000/11/14 06:47:45  rot032
c Changes from LDR for aerosol model
c
c Revision 1.3  2000/08/16 02:59:25  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.2  2000/06/20 02:27:11  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.1  2000/02/01 04:56:12  rot032
c Initial revision
c
      SUBROUTINE XTSINK (PTMST,  PXTM1,
     &                   PXTE)
C
C   *XTSINK*  CALCULATES THE DECREASE OF TRACER CONCENTRATION
C             FOR  A GIVEN HALF-LIFE-TIME
C
C   JOHANN FEICHTER               UNI-HAMBURG    08-91
C
C   PURPOSE
C  ---------------
C   THE MASS MIXING-RATIO OF TRACERS IS MULTIPLIED WITH
C   EXP(ALOG(0.5)*TIME-STEP/HALF-LIFE-TIME).
C   THIS ROUTINE COULD ALSO BE USED FOR EMISSION OR SINK
C   ABOVE THE SURFACE
C
C   INTERFACE
C  -----------
C   *XTSINK* IS CALLED FROM *PHYSC*
C
C   NO EXTERNALS
C

C Global parameters
      include 'PARAMS.f'
      include 'ECPARM.f'

C Argument list
      REAL PTMST
      REAL PXTM1(KLON,KLEV,NTRAC)
      REAL PXTE(KLON,KLEV,NTRAC)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables

C Local data, functions etc
      real sxtsink(ntrac)
c      data sxtsink / 3*0., 99000.,0.,99000.,0. / !Do this below

C Start code : ----------------------------------------------------------

      pxte(:,:,:)=0. !Very important!
      sxtsink(:)=0.

C BP commented out this portion for debugging Mk3L because ntrac=3 only
C      if(ntrac.gt.3)then
C        sxtsink(ITRACBC)=99000.
C        sxtsink(ITRACOC)=99000.
C      endif

      PQTMST=1./PTMST
      ZFAC=ALOG(0.5)*PTMST
C
      DO 150 JT=1,NTRAC
       IF (SXTSINK(JT).NE.0.) THEN
        ZDECAY=EXP(ZFAC/SXTSINK(JT))
C
        DO 140 JK=1,KLEV
         DO 120 JL=1,KLON
          ZXTP1=PXTM1(JL,JK,JT)+PXTE(JL,JK,JT)*PTMST
          ZXTP1=ZXTP1*ZDECAY
          ZDXTDT=(ZXTP1-PXTM1(JL,JK,JT))*PQTMST-PXTE(JL,JK,JT)
          PXTE(JL,JK,JT)=PXTE(JL,JK,JT)+ZDXTDT
          PXTE(JL,JK,JT+1)=PXTE(JL,JK,JT+1)-ZDXTDT 
  120    CONTINUE
C
  140   CONTINUE
       ENDIF
  150 CONTINUE
C
      RETURN
      END
