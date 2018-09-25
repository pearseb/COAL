c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Defines the COMMON block /SCALAR/, which was previously defined within the
c ocean model source code.
c SJP 2007/05/29

      real c2dtts, c2dtuv, c2dtsf, omega, radius, grav, radian, pi,
     &     swldeg
      COMMON /SCALAR/ C2DTTS,C2DTUV,C2DTSF,OMEGA,RADIUS,GRAV,RADIAN,PI,
     & SWLDEG
