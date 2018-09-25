c Purpose
c -------
c Writes the surface fluxes to a coupled model restart file.
c
c Inputs
c ------
c None
c
c Outputs
c -------
c None
c
c History
c -------
c 2008 Mar 7	Steven Phipps	Original version
c 2009 Apr 22	Steven Phipps	Modified for five-character experiment names

      subroutine oflux_write()

      implicit none

      include 'netcdf.inc'

C Global parameters
      include 'OPARAMS.f'

C Argument list

C Global data blocks
      include 'A2O.f'
      include 'FILES.f'
      include 'ORESTART.f'
      include 'VERSION.f'

C Local work arrays and variables
      character*8 file
      integer ittid, ncid, osalfid, osurhid, otxid, otyid, status, xid,
     &        yid, osradid, oswndid, osiceid
      parameter (file = "oflux.nc")

C Start code : ------------------------------------------------------------

c...  Create netCDF file
      status = nf_create(trim(file), nf_write, ncid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not create file ", file
        write (*, *)
        stop
      end if

c...  Define dimensions
      status = nf_def_dim(ncid, "x", imt-2, xid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define dimension"
        write (*, *)
        write (*, *) "***  File      =  ", file
        write (*, *) "***  Dimension =  x"
        write (*, *)
        stop
      end if
      status = nf_def_dim(ncid, "y", jmt-2, yid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define dimension"
        write (*, *)
        write (*, *) "***  File      =  ", file
        write (*, *) "***  Dimension =  y"
        write (*, *)
        stop
      end if

c...  Define variables
      status = nf_def_var(ncid, "itt", nf_int, 0, 0, ittid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define variable"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  itt"
        write (*, *)
        stop
      end if
      status = nf_def_var(ncid, "otx", nf_double, 2, (/ xid, yid /),
     &                    otxid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define variable"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  otx"
        write (*, *)
        stop
      end if
      status = nf_def_var(ncid, "oty", nf_double, 2, (/ xid, yid /),
     &                    otyid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define variable"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  oty"
        write (*, *)
        stop
      end if
      status = nf_def_var(ncid, "osalf", nf_double, 2, (/ xid, yid /),
     &                    osalfid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define variable"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osalf"
        write (*, *)
        stop
      end if
      status = nf_def_var(ncid, "osurh", nf_double, 2, (/ xid, yid /),
     &                    osurhid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define variable"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osurh"
        write (*, *)
        stop
      end if
c rjm
      if (nt.ge.3) then
      status = nf_def_var(ncid, "osrad", nf_double, 2, (/ xid, yid /),
     &                    osradid)
      if (status .ne. nf_noerr) then
        write (*, *) "***  netCDF error: Could not define variable"
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osrad"
        stop
      end if
      status = nf_def_var(ncid, "oswnd", nf_double, 2, (/ xid, yid /),
     &                    oswndid)
      if (status .ne. nf_noerr) then
        write (*, *) "***  netCDF error: Could not define variable"
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  oswnd"
        stop
      end if
      status = nf_def_var(ncid, "osice", nf_double, 2, (/ xid, yid /),
     &                    osiceid)
      if (status .ne. nf_noerr) then
        write (*, *) "***  netCDF error: Could not define variable"
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osice"
        stop
      end if

      endif
c rjm

c...  Put global attributes
      status = nf_put_att_text(ncid, nf_global, "version",
     &                         len(trim(version)), trim(version))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not put attribute"
        write (*, *)
        write (*, *) "***  File      =  ", file
        write (*, *) "***  Attribute =  version"
        write (*, *)
        stop
      end if
      status = nf_put_att_text(ncid, nf_global, "runtype", 5, runtype)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not put attribute"
        write (*, *)
        write (*, *) "***  File      =  ", file
        write (*, *) "***  Attribute =  runtype"
        write (*, *)
        stop
      end if

c...  Exit define mode
      status = nf_enddef(ncid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not exit define mode"
        write (*, *)
        write (*, *) "***  File =  ", file
        write (*, *)
        stop
      end if

c...  Write the timestep counter [ITT]
      status = nf_put_var_int(ncid, ittid, itt)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  itt"
        write (*, *)
        stop
      end if

c...  Write the zonal wind stress [OTX]
      status = nf_put_var_double(ncid, otxid, otx)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  otx"
        write (*, *)
        stop
      end if

c...  Write the meridional wind stress [OTY]
      status = nf_put_var_double(ncid, otyid, oty)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  oty"
        write (*, *) 
        stop
      end if

c...  Write the surface salinity tendency [OSALF]
      status = nf_put_var_double(ncid, osalfid, osalf)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osalf"
        write (*, *)
        stop
      end if

c...  Write the surface heat flux [OSURH]
      status = nf_put_var_double(ncid, osurhid, osurh)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osurh"
        write (*, *)
        stop
      end if
c rjm
	if (nt.ge.3) then
c...  Write the surface heat flux [OSRAD]
      status = nf_put_var_double(ncid, osradid, osrad)
      if (status .ne. nf_noerr) then
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osrad"
        stop
      end if
c...  Write the surface heat flux [OSWND]
      status = nf_put_var_double(ncid, oswndid, oswnd)
      if (status .ne. nf_noerr) then
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  oswnd"
        stop
      end if
c...  Write the surface heat flux [OSICE]
      status = nf_put_var_double(ncid, osiceid, osice)
      if (status .ne. nf_noerr) then
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  osice"
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
