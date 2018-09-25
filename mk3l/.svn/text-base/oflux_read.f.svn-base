c Purpose
c -------
c Reads the surface fluxes from a coupled model restart file.
c
c Inputs
c ------
c None
c
c Outputs
c -------
c itt		Timestep counter
c
c History
c -------
c 2008 Mar 7	Steven Phipps	Original version

      subroutine oflux_read(itt)

      implicit none

      include 'netcdf.inc'

C Global parameters
      include 'OPARAMS.f'

C Argument list
      integer itt

C Global data blocks
      include 'A2O.f'

C Local work arrays and variables
      character*8 file
      integer ncid, status, varid
      parameter (file = "oflux.nc")

C Start code : ------------------------------------------------------------

c...  Open netCDF file
      status = nf_open(trim(file), nf_nowrite, ncid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not open file ", file
        write (*, *)
        stop
      end if

c...  Get the timestep counter [ITT]
      status = nf_inq_varid(ncid, "itt", varid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not get variable ID"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  itt"
        write (*, *)
        stop
      end if
      status = nf_get_var_int(ncid, varid, itt)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not read data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  itt"
        write (*, *)
        stop
      end if

c...  Get the zonal wind stress [OTX]
      status = nf_inq_varid(ncid, "otx", varid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not get variable ID"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  otx"
        write (*, *)
        stop
      end if
      status = nf_get_var_double(ncid, varid, otx)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not read data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  otx"
        write (*, *)
        stop
      end if

c...  Get the meridional wind stress [OTY]
      status = nf_inq_varid(ncid, "oty", varid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not get variable ID"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  oty"
        write (*, *)
        stop
      end if
      status = nf_get_var_double(ncid, varid, oty)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not read data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  oty"
        write (*, *)
        stop
      end if

c...  Get the surface salinity tendency [OSALF]
      status = nf_inq_varid(ncid, "osalf", varid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not get variable ID"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osalf"
        write (*, *)
        stop
      end if
      status = nf_get_var_double(ncid, varid, osalf)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not read data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osalf"
        write (*, *)
        stop
      end if

c...  Get the surface heat flux [OSURH]
      status = nf_inq_varid(ncid, "osurh", varid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not get variable ID"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osurh"
        write (*, *)
        stop
      end if
      status = nf_get_var_double(ncid, varid, osurh)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not read data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osurh"
        write (*, *)
        stop
      end if
c rjm ... get additional variables for obgc
      if (nt.ge.3) then
      status = nf_inq_varid(ncid, "osrad", varid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not get variable ID"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osrad"
        write (*, *)
        stop
      end if
      status = nf_get_var_double(ncid, varid, osrad)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not read data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osrad"
        write (*, *)
        stop
      end if
      status = nf_inq_varid(ncid, "oswnd", varid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not get variable ID"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  oswnd"
        write (*, *)
        stop
      end if
      status = nf_get_var_double(ncid, varid, oswnd)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not read data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  oswnd"
        write (*, *)
        stop
      end if
      status = nf_inq_varid(ncid, "osice", varid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not get variable ID"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osice"
        write (*, *)
        stop
      end if
      status = nf_get_var_double(ncid, varid, osice)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not read data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osice"
        write (*, *)
        stop
      end if
      endif
c rjm 

c...  Close netCDF file
      status = nf_close(ncid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not close file ", file
        write (*, *)
        stop
      end if

      return
      end subroutine
