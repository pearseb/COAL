c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c Fix the call to NINT, resolving a warning issued by the g95 Fortran compiler.
c SJP 2009/04/15
c
c Specify the precision of all integer arrays and real variables that are
c used as arguments to netCDF routines. This ensures portability of the code
c across 32- and 64-bit platforms.
c SJP 2001/12/13
c
c Updated from v2.4 to v3 of netCDF.
c SJP 2001/11/21
c
c Removed declarations of integers ncopn, ncsfil and ncvid, which conflicted
c declarations in 'netcdf.inc'.
c SJP 2001/11/15
c
c $Log: ncput.f,v $
c Revision 1.5  2000/06/20 02:08:32  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.4  1995/05/02 06:10:36  ldr
c Remove spurious header message.
c
c Revision 1.3  1995/03/14  04:58:52  ldr
c Changes for new version of ncput which handles reduced T63 arrays.
c
c Revision 1.2  1993/07/20  10:35:40  ldr
c Replace exclamation marks in full line comments with C's so as not
c to confuse fix_do on the VP.
c
c Revision 1.1  93/02/05  16:46:25  ldr
c Initial revision
c 

      subroutine nc_put(nlon, nlat, path, varname, year, month, matin)

c Write data to netCDF file

c Author: Harvey Davies

      implicit none

c arguments (all input)
      integer       nlon ! number of longitudes
      integer       nlat ! number of latitudes
      character*(*) path ! netCDF filename
      character*(*) varname ! netCDF variable name
      integer       year ! subscript of netCDF variable = 1, 2, 3, ...
      integer       month ! subscript of netCDF variable = 1, 2, ..., 12
      real          matin(nlon,nlat) ! input data matrix

c parameters
      real misval ! missing value
      parameter (misval = -7777000.0)

      include 'netcdf.inc'

c local variables
      integer     cdfid ! netCDF file id
      integer     count(4) ! hyperslab edge vector
      integer     i ! subscript
      integer     j ! subscript
      integer     rcode ! return code
      integer     start(4) ! hyperslab corner vector
      integer     varid ! netCDF variable id of matrix
      integer     varidy ! netCDF variable id of year
      integer     year1 ! 1st year in file
      integer     old_mode
      real*4      matout(nlon,nlat) ! output data matrix

c open file with read-write access
      rcode = nf_open(path, nf_write, cdfid)
      if (rcode .ne. nf_noerr) stop "***  netCDF error in nc_put"

c set the fill mode to NF_NOFILL
      rcode = nf_set_fill(cdfid, nf_nofill, old_mode)
      if (rcode .ne. nf_noerr) stop "***  netCDF error in nc_put"

c get variable ID for year
      rcode = nf_inq_varid(cdfid, 'year', varidy)
      if (rcode .ne. nf_noerr) stop "***  netCDF error in nc_put"

c get variable ID for data
      rcode = nf_inq_varid(cdfid, varname, varid)
      if (rcode .ne. nf_noerr) stop "***  netCDF error in nc_put"

c get 1st year in file
      rcode = nf_get_var1_int(cdfid, varidy, 1, year1)
      if (rcode .ne. nf_noerr) stop "***  netCDF error in nc_put"

c write data to file
      start(1) = 1
      start(2) = 1
      start(3) = month
      start(4) = 1 + year - year1
      count(1) = nlon
      count(2) = nlat
      count(3) = 1
      count(4) = 1
      do j = 1, nlat
	do i = 1, nlon
	  if ( matin(i,j) .gt. misval ) then
	    matout(i,j) = real(matin(i, j), 4)
	  else
	    matout(i,j) = real(nf_fill_float, 4)
	  end if
	end do
      end do
      rcode = nf_put_vara_real(cdfid, varid, start, count, matout)
      if (rcode .ne. nf_noerr) stop "***  netCDF error in nc_put"

c write year to file
      rcode = nf_put_var1_int(cdfid, varidy, start(4), year)
      if (rcode .ne. nf_noerr) stop "***  netCDF error in nc_put"

c close file
      rcode = nf_close(cdfid)
      if (rcode .ne. nf_noerr) stop "***  netCDF error in nc_put"

      end
