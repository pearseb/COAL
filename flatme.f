c $Log: flatme.f,v $
c Revision 1.14  2001/02/22 05:34:39  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.13  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.12  1998/12/10  00:55:37  ldr
c HBG changes to V5-1-21
c
c Revision 1.11  1996/06/13  02:06:34  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.10  1996/04/17  01:32:52  ldr
c Changes from TIE to allow initialization of monthly mean netcdf files
c from within model.
c
c Revision 1.9  1996/03/21  03:18:42  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.8  1995/05/02  07:25:58  ldr
c Merge HBG V4-5-43h2 stuff into latest version.
c
c Revision 1.7  95/05/02  07:14:37  ldr
c Other changes from HBG that went into coupled run but not into V4-5-43h.
c 
c Revision 1.6  94/03/30  12:34:24  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.5  93/11/29  14:49:23  ldr
c Use icu_sflg and icv_sflg for treatment of ice U,V diagnostics, so they
c are consistent with other "stats".
c 
c Revision 1.4  93/11/29  12:33:12  ldr
c SPO's changes to rationalize printing of icuv statistics.
c 
c Revision 1.3  93/08/10  15:27:12  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c 
c Revision 1.2  93/02/03  12:44:43  ldr
c SPO's minor changes to coupled model and sea-ice scheme of V4-1.
c 
c Revision 1.1  92/06/16  11:54:47  ldr
c Initial revision
c 
      subroutine flatme(ijdyn)

c dumps ice model velocities averaged over number of steps carried out
c for the ice model (the ice model step has a factor of ijdyn per atmos step)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer ijdyn

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'COMDICE.f'
      include 'FILES.f'
      include 'STFLAGS.f'
      include 'TIMEX.f'

      real redice,groice,hav
      common/red/redice(ln2,lat),groice(ln2,lat),hav(ln2,lat)

C Local work arrays and variables
      integer i
      integer j

      real avg

C Local data, functions etc

C Start code : ----------------------------------------------------------

      avg=float(ijdyn)/nsteps
      do j=1,plat
         do i=1,plon
            fav(i,j)=fav(i,j)*avg
            uav(i,j)=uav(i,j)*avg
            vav(i,j)=vav(i,j)*avg
            wl(i,j)=wl(i,j)*avg*1.e8
            wd(i,j)=wd(i,j)*avg*1.e8
         end do
      end do

      if (icu_sflg) call colldyu('icu',1.0,str,uav)
      if (icv_sflg) call colldyu('icv',1.0,str,vav)
      if (ich_sflg) call collph ('ich',1.0,str,hav)
      if (ire_sflg) call collph ('ire',1.0,str,redice)
      if (gro_sflg) call collph ('gro',1.0,str,groice)

      if(div_sflg)then
       call colldyu('wls',1.0,str,wl(1,1))
       call colldyu('wdf',1.0,str,wd(1,1))
      endif

      return
      end
