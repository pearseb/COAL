subroutine hstring(hist)

! Purpose
! -------
! Returns a string containing the user name, host name, date and time.
!
! Arguments
! ---------
! In:	 	None
!
! Out:		hist		output string
!
! History
! -------
! 2008 Nov 20	Steven Phipps	Original version, replacing hstring.c
! 2009 Apr 14   Steven Phipps   Fixing syntax error identified by the g95
!                               Fortran compiler

  implicit none

  ! Declare arguments
  character(len=82), intent(out) :: hist

  ! Declare local variables
  character(len=5) :: zone
  character(len=8) :: date
  character(len=10) :: logname, time
  character(len=30) :: hostname
  integer, dimension(8) :: date_time

  ! ---------------------
  ! BEGIN EXECUTABLE CODE
  ! ---------------------

  ! Get user name
  call get_environment_variable("LOGNAME", logname)

  ! Get host name
  call get_environment_variable("HOSTNAME", hostname)

  ! Get current date and time
  call date_and_time(date, time, zone, date_time)

  ! Initialise output string
  hist = " "

  ! Construct output string
  write (hist, "('Created by ', a, ' on ', a, ', ', i4.4, '/', i2.2, '/', &
               & i2.2, ' ', i2.2, ':', i2.2, ':', i2.2, ' ', a)") &
    trim(logname), trim(hostname), date_time(1), date_time(2), date_time(3), &
    date_time(5), date_time(6), date_time(7), zone

end subroutine hstring
