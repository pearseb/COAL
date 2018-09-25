c $Log: errcheck.f,v $
c Revision 1.1  1998/05/27 02:03:44  ldr
c Initial revision
c
	subroutine errcheck(ierr, filename, progname)
c
c	Simple subprogram to check error status of data files and print a
c	useful message if the data file is missing.
c

	character*17 filename
	character*10 progname

	if(ierr.ne.0) write(0,9)
     &	'**************************************************',
     &  'Cannot open file ',filename ,' in ', progname,
     &  '**************************************************'

 9	format(/,1x,a50,/,1x,a17,a17,a4,a10,/,1x,a50)

	return
	end
