c Implementing a Robert time filter.
c SJP 2008/12/17
c
c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/16
c
c Defines the COMMON block /FULLWD/, which was previously defined within the
c ocean model source code.
c SJP 2007/05/29

      integer ndiskb, ndisk, ndiska, mix, mxp, nergy, mscan, kar
      real engint, engext, ttdtot, buoy, plicin, plicex, ektot, dtabs,
     &     tvar, spsin, spcos, pnu2m
      LOGICAL EB
      COMMON /FULLWD/ ENGINT(8),ENGEXT(8),TTDTOT(6,NT),BUOY,PLICIN,
     & PLICEX,EKTOT,DTABS(NT),TVAR(NT),NDISKB,NDISK,NDISKA,
     & MIX,MXP,NERGY,MSCAN,KAR(KM),EB,SPSIN(IMT),SPCOS(IMT),
     & pnu2m
