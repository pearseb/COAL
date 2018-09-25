c $Log: ecread.f,v $
c Revision 1.6  2001/10/12 02:13:45  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.5  2001/06/04 02:26:58  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.4  2000/12/08 03:58:53  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.3  2000/11/16 01:19:01  rot032
c Code to read corrected volcano data.
c
c Revision 1.2  2000/11/14 06:47:45  rot032
c Changes from LDR for aerosol model
c
c Revision 1.1  2000/06/20 03:38:02  rot032
c Initial revision
c
      subroutine ecread

c This routine reads in fields needed by the interactive aerosol treatment
c Note that vertical levels are counted from the TOA as in the ECHAM model.

C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'ECPARM.f'

C Argument list

C Global data blocks

      include 'TIMEX.f' !Input month
      include 'ECFIELD.f' !Output source and oxidant fields
      include 'ECAMEAN.f' !Output annual-mean source fields
      include 'FEWFLAGS.f' !Input PI_emissions

C Local work arrays and variables

      character cres*3, cmon*3, header*80, fname*30
      real dum2d(lon,lat2)

C Local data, functions etc


C Start code : ----------------------------------------------------------

      write(cres,'(a1,i2)')trunc,mw-1
      print*,'Reading ECHAM data files for resolution ',cres

c Read annual-mean fields and reorder into arrays with double-length latitudes

      open(3,file='echam-anndat2.'//cres,status='old',iostat=ierr)
      call errcheck(ierr, 'echam-anndat2.'//cres, 'ecread    ')

      read(3,'(a)')header
      write(6,*)'Header is ',header
      read(3,*)dum2d

      do j=1,lat
        do i=1,lon
          vso2(i,j)=dum2d(i,lat2p-j) !NH
          vso2(i+lon,j)=dum2d(i,j) !SH
        enddo
      enddo

      read(3,'(a)')header
      write(6,*)'Header is ',header
      read(3,*)dum2d

      do j=1,lat
        do i=1,lon
          hvolc(i,j)=dum2d(i,lat2p-j) !NH
          hvolc(i+lon,j)=dum2d(i,j) !SH
        enddo
      enddo

      read(3,'(a)')header
      write(6,*)'Header is ',header
      read(3,*)dum2d

      do j=1,lat
        do i=1,lon
          bcsrc(i,j)=dum2d(i,lat2p-j) !NH
          bcsrc(i+lon,j)=dum2d(i,j) !SH
        enddo
      enddo

      read(3,'(a)')header
      write(6,*)'Header is ',header
      read(3,*)dum2d

      do j=1,lat
        do i=1,lon
          ocsrc(i,j)=dum2d(i,lat2p-j) !NH
          ocsrc(i+lon,j)=dum2d(i,j) !SH
        enddo
      enddo

      close (3)



c Read monthly mean fields
c mondat2 denotes reduced data files (no dust, seasalt), so natorg becomes 80.
c mondat3 denotes as mondat2, but with corrected oxidant fields in top layer.
c mondat4 denotes as mondat2, but with updated ocean DMS data set (raw Kettle et al '99 data).
c Then split mondat4 into multi-level (ml) and single-level (sl) fields to save disk in
c transient runs.

      write(cmon,'(a1,i2.2)')'m',month
      fname='echam-mondat4_ml.'//cres//'.'//cmon
      print*,'fname = ',fname
      open(3,file=fname,status='old',iostat=ierr)
      call errcheck(ierr, fname, 'ecread    ')

c First four are 3d oxidant fields (oh, h2o2, o3, no2) 1-72 for 18-level model.
c Next few are 2d source fields (dmsterr, so2biom, so2anth1,
c so2anth2, dmssea, bcbiom, ocbiom) 73-79 for 18-level model.
c Then 3d dust (80-97) and seasalt (98-115), followed by 2d natorg (116).

      do indx=1,4*nl
        read(3,'(a)')header
        write(6,*)'Header is ',header
        read(3,*)dum2d

        do j=1,lat
          do i=1,lon
            field(i,indx,j)=dum2d(i,lat2p-j) !NH
            field(i+lon,indx,j)=dum2d(i,j) !SH
          enddo
        enddo
        
      enddo

c Now do the single-level fields

      fname='echam-mondat4_sl.'//cres//'.'//cmon
      print*,'fname = ',fname
      open(3,file=fname,status='old',iostat=ierr)
      call errcheck(ierr, fname, 'ecread    ')

      do indx=4*nl+1,numfl2
        read(3,'(a)')header
        write(6,*)'Header is ',header
        read(3,*)dum2d

        do j=1,lat
          do i=1,lon
            field(i,indx,j)=dum2d(i,lat2p-j) !NH
            field(i+lon,indx,j)=dum2d(i,j) !SH
          enddo
        enddo
        
      enddo

      if ( PI_emissions ) then
        print*,'Warning: Setting aerosol emissions to pre-industrial...'
        field(:,74,:)=0.1*field(:,74,:)
        field(:,75,:)=0.
        field(:,76,:)=0.
        field(:,78,:)=0.1*field(:,78,:)
        field(:,79,:)=0.1*field(:,79,:)
      endif

      return
      end
