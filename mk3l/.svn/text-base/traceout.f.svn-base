c $Log: traceout.f,v $
c Revision 1.4  2001/02/22 05:34:42  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.3  1993/12/17 15:34:14  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.2  91/03/01  17:45:33  ldr
c Changes unit numbers of output files. (To free range 21-50 for impacts data)
c 
c Revision 1.1  91/02/22  16:38:19  ldr
c Initial release V3-0
c 
      subroutine traceout(iflag,nsteps)

c Little routine to output positions of particles in tracer experiment
c IFLAG = 0 is used to write diagnostic every 8 hours
c IFLAG = 1 is used to write restart file as end of run
      
      implicit none
C Global parameters
      integer nps
      parameter (nps=63*48)

C Argument list
      integer iflag
      integer nsteps

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      real lonpos,latpos,sigpos
      integer npts
      common/positn/lonpos(nps),latpos(nps),sigpos(nps),npts

C Local work arrays and variables
      character*50 sfn

      integer kk

C Local data, functions etc

C Start code : ----------------------------------------------------------

      if(iflag.eq.0)then

        if(mod(nsteps-1,8).ne.0) return
        write(sfn,556)nsteps-1
 556    format('lonlatsig',i4.4)
        open(unit=51,file=sfn,form='formatted',status='unknown')
        write(51,4444)(lonpos(kk),kk=1,npts)
        write(51,4444)(latpos(kk),kk=1,npts)
        write(51,4444)(sigpos(kk),kk=1,npts)
 4444   format(12f6.3)
        close(unit=51)

      else

        open(unit=51,file='partposout',form='unformatted',
     &  status='unknown')
        write(51)npts
        write(51)lonpos
        write(51)latpos
        write(51)sigpos
        close(unit=51)

      endif

      return
      end
