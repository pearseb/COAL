c Major tidy-up of ocean model source code.
c SJP 2007/06/16
c
c Defines the COMMON block /FLXBLOK/, which was previously defined within the
c ocean model source code.
c SJP 2007/05/29

      real flux
      COMMON/FLXBLOK/FLUX(IMT,JMT,2)
