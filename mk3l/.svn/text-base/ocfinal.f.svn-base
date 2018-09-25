c Modified for the conversion of the ocean model and coupled model restart
c files to netCDF.
c SJP 2008/03/08
c
c Rationalise the NAMELIST input for the ocean model.
c SJP 2007/11/24
c
c Modified to remove the redundant manipulation of the array CDEPTH.
c SJP 2007/06/19
c
c Major tidy-up of ocean model source code.
c SJP 2007/06/18
c
c Added the line "include 'PARAMS.f'", as this line is no longer included via
c the header file OPARAMS.f.
c SJP 2007/05/31
c
c Added IMPLICIT NONE statement, plus variable declarations as necessary.
c SJP 2007/05/30
c
c Transferred COMMON block declarations to separate header files.
c SJP 2007/05/29
c
c Redundant code removed.
c SJP 2004/01/06
c
c Modified to make use of /orestart/.
c SJP 2003/09/03
c
c Call to OCLOSE removed, as they did nothing. Parameter definitions moved
c to include file OPARAMS.f.
c SJP 2003/04/29
c
c Writes to fort.93 commented out, as this data is not required.
c SJP 2002/06/27
c
c $Log: ocfinal.f,v $
c Revision 1.5  2000/06/20 02:08:31  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.4  1997/12/19 01:25:31  ldr
c Changes from ACH for 21 OGCM levels and optinal GM mixing scheme...
c
c Include array for storing depth of convective mixed layer
c Change to 21 levels in ocean model, and insertion of
c eddy-induced transport (major changes delineated)
c
c Revision 1.3  1994/03/30  12:36:05  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c
c Revision 1.2  93/12/17  15:34:35  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.1  93/08/10  15:03:50  ldr
c Initial revision
c 
c Revision 1.1  93/02/05  16:22:22  ldr
c Initial revision
c 
      subroutine ocfinal(lcouple)
C
C=======================================================================
C         Closing routine for ocean model                            ===
C=======================================================================
C
      implicit none

      include 'OPARAMS.f'
      include 'OCEAN_NML.f'

      logical lcouple
C
C=======================================================================
C  END RUN FOR OCEAN MODEL  ======================================
C=======================================================================
C

c...  Write the data to the ocean model restart file, if required
      if (na .eq. 1) then
        write (*, *)
        write (*, *) "Writing data to ocean model restart file..."
        write (*, *)
        call orest_write
      end if

c...  If running in coupled mode, write the surface fluxes to the coupled model
c...  restart file
      if (lcouple) then
        write (*, *)
        write (*, *) "Writing data to coupled model restart file..."
        write (*, *)
        call oflux_write
      end if

      return
      END
