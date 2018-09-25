c $Log: radh_l.f,v $
c Revision 1.2  2001/10/12 02:06:57  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.1  2001/02/22 05:36:47  rot032
c Initial revision
c
      subroutine radh_copy(lg,ipass,sg,sga,tg,rg,htk,
     &  cl,cll,clm,clh,rgdn,blg,sgdn,als)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'

C Argument list
      integer lg
      integer ipass
      real sg(ln2)
      real sga(ln2)
      real tg(ln2)
      real rg(ln2)
      real htk(ln2,nl)
      real cl(ln2)
      real cll(ln2)
      real clm(ln2)
      real clh(ln2)
      real rgdn(ln2)
      real blg(ln2)
      real sgdn(ln2)
      real als(ln2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'MASIV3.f'
      include 'PMASIV3.f'

C Local work arrays and variables
      integer k
      integer mg

      real xxx

C Start code : ----------------------------------------------------------

c Copy radiation output fields into global arrays 
c  to save for next NRAD timesteps

      if(ipass.eq.1)then ! non-leads data

        do mg=1,ln2
          sgsav(mg,lg)=sg(mg)
          sgamp(mg,lg)=sga(mg)
c     Save the value excluding Ts^4 part.  This is allowed to change.
          xxx=stefbo*tg(mg)**4
          rgsav(mg,lg)=rg(mg)-xxx
          htk(mg,1)=htk(mg,1)-xxx
        enddo

        do k=1,nl
        do mg=1,ln2
          htksav(mg,k,lg)=htk(mg,k)
        enddo
        enddo

      elseif(ipass.eq.2)then ! leads data

        do mg=1,ln2
          psgsav(mg,lg)=sg(mg)
          psgamp(mg,lg)=sga(mg)
c     Save the value excluding Ts^4 part.  This is allowed to change.
          xxx=stefbo*tg(mg)**4
          prgsav(mg,lg)=rg(mg)-xxx
          htk(mg,1)=htk(mg,1)-xxx
        enddo

        do k=1,nl
        do mg=1,ln2
          phtksav(mg,k,lg)=htk(mg,k)
        enddo
        enddo

      endif

      do mg=1,ln2
c Total cloud 
        cl(mg)=(cll(mg)+clm(mg)-cll(mg)*clm(mg))
     &      *(1.0-clh(mg))+clh(mg)
        rgdn(mg)=blg(mg)-rg(mg)
        sgdn(mg)=sg(mg)/(1-als(mg))
        htk(mg,1)=htk(mg,1)+stefbo*tg(mg)**4
      enddo

      return
      end
C---------------------------------------------------------------------
      subroutine radh_old(lg,ipass,htk,rg,ifroz)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer lg
      integer ipass
      real htk(ln2,nl)
      real rg(ln2)
      integer ifroz(2)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'MASIV3.f'
      include 'PMASIV3.f'

C Local work arrays and variables
      integer k
      integer ma
      integer mg
      integer ns

C Start code : ----------------------------------------------------------

c Copy radiation fields from global to working arrays.
c No need to do this on radiation time steps

      if(ipass.eq.1)then ! non-leads data

        do k=1,nl
        do mg=1,ln2
           htk(mg,k)=htksav(mg,k,lg)
        enddo
        enddo

c... (rg(mg)-rgsav(mg,lg)) now represents the updated stefbo*tg(mg)**4
        do mg=1,ln2
          htk(mg,1)=htk(mg,1)+(rg(mg)-rgsav(mg,lg))
        enddo

      elseif(ipass.eq.2)then ! leads data

c.....  phtksav may not be set when a new latitude
c.....  row has just frozen a point
           do ns=1,2
            if(ifroz(ns).eq.1)then
             do  k=1,nl
             do  mg=1,lon
                ma=mg+(ns-1)*lon
                phtksav(ma,k,lg)=htksav(ma,k,lg)
             enddo
             enddo
            endif
           enddo

        do k=1,nl
        do mg=1,ln2
           htk(mg,k)=phtksav(mg,k,lg)
        enddo
        enddo

        do mg=1,ln2
          htk(mg,1)=htk(mg,1)+(rg(mg)-prgsav(mg,lg))
        enddo

      endif

      return
      end
