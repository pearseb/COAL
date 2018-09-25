c Purpose
c -------
c Defines the COMMON block /AOGCM/, which contains parameters relevant to the
c coupling of the atmosphere to the ocean.
c
c Contains
c --------
c ISYNC		0 = concurrent coupling (SSTs passed to atmosphere are those
c                                        from the previous ocean timestep)
c               1 = sequential coupling (SSTs passed to atmosphere are
c                                        interpolated in time between ocean
c                                        timesteps)
c OVOLUME	Volume of ocean (m^3)
c
c History
c -------
c 2009 Apr 21	Steven Phipps	Original version

      integer isync
      real ovolume
      common /aogcm/ isync, ovolume
