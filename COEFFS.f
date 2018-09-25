c Enabling user control over the ocean model statistics that are saved to file.
c SJP 2009/04/19
c
c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/16
c
c Defines the COMMON block /COEFFS/, which was previously defined within the
c ocean model source code.
c SJP 2007/05/29

      real rhocw, amfac, ctop, gamma
      COMMON/COEFFS/RHOCW,AMFAC(JMT),CTOP,GAMMA
