c $Log: initice.f,v $
c Revision 1.4  2001/02/22 05:34:39  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.3  1999/05/20 06:23:50  rot032
c HBG changes to V5-2
c
c Revision 1.2  1996/10/24  01:03:32  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.1  1996/02/19  03:40:08  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:   from common/freeze in FREEZE.f
c                  tfi - temperature at bottom of ice
c
c     Output:  from common/masiv5 in this subroutine
c                  statsice - global statistics relating to sea-ice scheme
c
c     In/Out:  from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c              from common/masiv4 in MASIV4.f
c                  savegrid - global surface variables
c

      subroutine initice

c Initialize sea ice scheme if starting from restart file that 
c contains no sea ice data. Assumes January conditions.

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FREEZE.f'
      include 'LSMI.f'
      include 'MASIV4.f'

      real statsice
      common/masiv5/statsice(ln2,14,lat)

C Local work arrays and variables
      real hs(lon,lat2),hi(lon,lat2),arlon(lon),arlat(lat2)
      real hs1(lon,lat2),hi1(lon,lat2)
      character*80 header

      integer ierr
      integer ix
      integer iy
      integer j
      integer lg
      integer lgns
      integer ma
      integer maxlg
      integer minlg
      integer mg
      integer ns

      real tf
      real ts

C Local data, functions etc

C Start code : ----------------------------------------------------------
      
      tf=tfi   
      do lg=1,lat
          do mg=1,ln2
            statsice(mg,1,lg)=savegrid(mg,12,lg)*0.01
            statsice(mg,2,lg)=savegrid(mg,13,lg)
            ts=savegrid(mg,3,lg)
            statsice(mg,3,lg)=ts+(tf-ts)/6
            statsice(mg,4,lg)=ts+(tf-ts)/2
            statsice(mg,5,lg)=ts+5*(tf-ts)/6
            statsice(mg,6,lg)=ts+(tf-ts)/3
            statsice(mg,7,lg)=0.0
            statsice(mg,8,lg)=0.0
            statsice(mg,9,lg)=tf
          enddo
      enddo

      do lg=1,lat ! Fix SH leads
        do mg=1,lon
          ma=mg+lon
          if(imsl(ma,lg).eq.1)statsice(ma,8,lg)=0.300
        enddo
      enddo
      do lg=1,lat ! Fix NH leads
        do mg=1,lon
          if(imsl(mg,lg).eq.1)then
            if( lg.le.lat/4)then
              statsice(mg,8,lg)=0.020
            else
              statsice(mg,8,lg)=0.2
            endif
          endif
        enddo
      enddo

      if(lw.eq.22)then          !R21 initialization of ice and snow depths

        open(7,file='icesn.dat',form='formatted',iostat=ierr)
         call errcheck(ierr,'icesn.dat        ','initice   ')
         read(7,697)header
         read(7,699)(arlon(j),j=1,lon)
         read(7,693)
         do 53 iy=1,lat2
           read(7,698)arlat(iy),(hs(ix,iy),ix=1,lon)
 53      continue
         read(7,697)header
         read(7,699)(arlon(j),j=1,lon)
         read(7,693)
         do 54 iy=1,lat2
           read(7,698)arlat(iy),(hi(ix,iy),ix=1,lon)
 54      continue
        close(7)
 693    format(1x)
 697    format(80a,/)
 698    format(f6.1,8f8.2,/,7(6x,8f8.2,/))
 699    format(6x,8f8.2)

        do lg=1,12
          do mg=1,lon
            statsice(mg+lon,1,lg)=hs(mg,lg)*0.01       
            if(imsl(mg+lon,lg).ne.4)then
              savegrid(mg+lon,12,lg)=hs(mg,lg)
              savegrid(mg+lon,13,lg)=hi(mg,lg)
            end if
            statsice(mg+lon,2,lg)=hi(mg,lg)
          enddo
        enddo

        open(7,file='arctic.snice',form='formatted',iostat=ierr)
         call errcheck(ierr,'arctic.snice     ','initice   ')
         do 63 iy=1,14
           read(7,678)(hs1(ix,iy),ix=1,lon)
 63      continue
         do 64 iy=1,14
           read(7,678)(hi1(ix,iy),ix=1,lon)
 64      continue
        close(7)
 678    format(8(6x,8f8.2,/))

        do lg=1,14
          do mg=1,lon
            statsice(mg,1,lg)=hs1(mg,lg)*0.01       
            if(imsl(mg,lg).ne.4)then
              savegrid(mg,12,lg)=hs1(mg,lg)
              savegrid(mg,13,lg)=hi1(mg,lg)
            end if
            statsice(mg,2,lg)=hi1(mg,lg)
          enddo
        enddo
        write(6,*)'read in new Antarctic data'
c     correct imsl values for southern ocean ice initialized from  ty2 run
        do lg=1,12
          do mg=1,lon
            if(hi(mg,lg).gt.0.0)then
              if(imsl(mg+lon,lg).ne.1)then
                imsl(mg+lon,lg)=1
                if(imsl(mg+lon,lg+1).eq.3)imsl(mg+lon,lg+1)=-2
              endif
            endif
          enddo
        enddo
        do lg=1,14
          do mg=1,lon
            if(hi1(mg,lg).gt.0.0)then
              if(imsl(mg,lg).ne.1)then
                imsl(mg,lg)=1
                if(imsl(mg,lg+1).eq.3)imsl(mg,lg+1)=-2
              endif
            endif
          enddo
        enddo
        
      elseif(lw.eq.43)then      !R42
        
c R42 initialization of ice and snow depth (in fact this code is general).
c In each hemisphere do linear interpolation for hi b/w lg=minlg and lg=maxlg 
c In SH ice depth is hi = 0 at lg >= maxlg increasing to hi = 1m at lg <= minlg
c In NH ice depth is hi = 0 at lg >= maxlg increasing to hi = 3m at lg <= minlg
c Just initialize snow depth over ice (hs) to 20cm everywhere
        
        minlg=nint(0.167*lat)   !For SH
        maxlg=nint(0.32*lat)    !For SH
        do lg=1,lat
          do mg=1,lon
            hi(mg,lg)=max(min(float(maxlg-lg)/(maxlg-minlg),1.), 0.)
            hs(mg,lg)=20.
          enddo
        enddo
        
        minlg=nint(0.11*lat)    !For NH
        maxlg=nint(0.39*lat)    !For NH
        do lg=1,lat
          do mg=1,lon
            hi(mg,lat2p-lg)=
     &           3.*max(min(float(maxlg-lg)/(maxlg-minlg),1.), 0.)
            hs(mg,lat2p-lg)=20.
          enddo
        enddo
        
c     Now overwrite values in savegrid and statsice arrays at ice points only
        
        do lgns=1,lat2
          ns=2-(lgns-1)/lat
          lg=(ns-1)*lgns+(lat2+1-lgns)*(2-ns)
          do mg=1,lon
            ma=mg+(ns-1)*lon
            if(imsl(ma,lg).eq.1)then
              statsice(ma,1,lg)=hs(mg,lgns)*0.01
              statsice(ma,2,lg)=hi(mg,lgns)
              savegrid(ma,12,lg)=hs(mg,lgns)
              savegrid(ma,13,lg)=hi(mg,lgns)
            endif
          enddo
        enddo


      endif                     ! lw

      return
      end
