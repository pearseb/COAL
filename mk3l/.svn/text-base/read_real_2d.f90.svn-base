subroutine read_real_2d(file, variable, nx, ny, data, title)

! Purpose
! -------
! Generic routine to read a two-dimensional real array from a netCDF file.
!
! Arguments
! --------
! In:		file		name of netCDF file
!		variable	variable name
!		nx		size of x-dimension
!		ny		size of y-dimension
!
! Out:		data		data
!		title		title of dataset
!
! History
! -------
! 2008 Feb 5	Steven Phipps	Original version

  implicit none

  include 'netcdf.inc'

  ! Declare arguments
  character(len=*), intent(in) :: file, variable
  character(len=60), intent(out) :: title
  integer, intent(in) :: nx, ny
  real, dimension(nx, ny), intent(out) :: data

  ! Declare local variables
  integer :: ncid, status, tlen, varid
  logical :: is_title

  ! ---------------------
  ! BEGIN EXECUTABLE CODE
  ! ---------------------

  ! Open netCDF file
  status = nf_open(file, nf_nowrite, ncid)
  if (status /= nf_noerr) then
    write (*, *)
    write (*, *) "***  netCDF error: Could not open file ", file
    write (*, *)
    stop
  end if

  ! Get variable ID
  status = nf_inq_varid(ncid, variable, varid)
  if (status /= nf_noerr) then
    write (*, *)
    write (*, *) "***  netCDF error: Could not get variable ID"
    write (*, *)
    write (*, *) "***  File     =  ", file
    write (*, *) "***  Variable =  ", variable
    write (*, *)
    stop
  end if

  ! Get data
  status = nf_get_var_double(ncid, varid, data)
  if (status /= nf_noerr) then
    write (*, *)
    write (*, *) "***  netCDF error: Could not read data"
    write (*, *)
    write (*, *) "***  File     =  ", file
    write (*, *) "***  Variable =  ", variable
    write (*, *)
    stop
  end if

  ! Get length of title - if there is an error, assume that the dataset has no
  !                       title
  status = nf_inq_attlen(ncid, nf_global, "title", tlen)
  if (status /= nf_noerr) then
    is_title = .false.
  else
    is_title = .true.
  end if

  ! Set title to a blank character string
  title = "                                        " // &
          "                                        "

  ! Get title of dataset, if there is one
  if (is_title) then

    ! Get title of dataset
    status = nf_get_att_text(ncid, nf_global, "title", title)
    if (status /= nf_noerr) then
      write (*, *)
      write (*, *) "***  netCDF error: Could not read title from file ", file
      write (*, *)
      stop
    end if

  end if

  ! Close netCDF file
  status = nf_close(ncid)
  if (status /= nf_noerr) then
    write (*, *)
    write (*, *) "***  netCDF error: Could not close file ", file
    write (*, *)
    stop
  end if

  return

end subroutine read_real_2d
