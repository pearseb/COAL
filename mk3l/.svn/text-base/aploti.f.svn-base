c $Log: aploti.f,v $
c Revision 1.3  2001/10/12 02:06:57  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.2  1992/11/20 12:11:44  ldr
c Changed for extra line plots. (hbg)
c
c Revision 1.1  91/02/22  16:36:45  ldr
c Initial release V3-0
c 
      subroutine aploti(xdata,chi,aplot,lg)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      real xdata
      character*1 chi
      character*2 aplot(lat2,41)
      integer lg

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      integer ix

      character*2 ch2
      character*1 chx(2)
      equivalence (chx(1),ch2)

C Local data, functions etc

C Start code : ----------------------------------------------------------

C**** ROUTINE TO CHECK PAGE PLOTTING LIMITS -
C****  FLAG WITH 'X' IF OUT OF BOUNDS.

      chx(1)=chi
      chx(2)=' '
      ix=nint(0.1*xdata+21.0)
      if(ix.le.41)go to 10
        ix=41
        chx(2)='X'
        aplot(lg,ix)=ch2
        return
   10 if(ix.ge.1)go to 12
        ix=1
        chx(2)='X'
        aplot(lg,ix)=ch2
        return
   12 ch2=aplot(lg,ix)
      if((chx(1).ne.' ').and.(chx(1).ne.'-'))then
c.... if character exists, then swap
        chx(2)=chx(1)
        chx(1)=chi
      else
        chx(1)=chi
        chx(2)=' '
      end if
      aplot(lg,ix)=ch2

      return
      end
