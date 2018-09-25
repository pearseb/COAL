c Purpose
c -------
c Writes the OBGC data to a file.
c
c Inputs
c ------
c None
c
c Outputs
c -------
c None
c
c History
c -------

      subroutine obgc_write()

      implicit none
      include 'netcdf.inc'
C Global parameters
      include 'OPARAMS.f'

C Argument list

C Global data blocks
      include 'FILES.f'
      include 'ORESTART.f'
      include 'VERSION.f'

      include 'ONEDIM.f'
      include 'bio.h'
      include 'extra.h'

C Local work arrays and variables
      character*4 file
      character*6 qnumber
      character*13 qfile
      integer ncid, status, imtid, jmtid, kmid, lsegid, lsegfid,
     &        nisleid, njtbftid, njtbfuid, nstep2id, ittid, ttsecid,
     &        areaid, volumeid, pbid, pid, hrid, ptd2id, iszid, iezid,
     &        isisid, ieisid, jsisid, jeisid, istfid, ietfid, isufid,
     &        ieufid, iszfid, iezfid, tid, sid, uid, vid, kmtid, kmuid,
     &        lontsid, lattsid, ztsid, lmtid, timid
      integer i, j, k, l

      real*4 xt(imt), lonts(imt), latts(jmt), zts(km), tim0(1)
      real*8 missing_value, pi, r2d, radius
      parameter (file = "obgc")
      parameter (missing_value=-9999.9)

* rjm
      character*6 obgc_tr
      integer obgc_id(nt),ofbgc_id(nbio2d),denID,sedDID,sedSID,fixID,
     +        cpID,npID,o2remID,no3remID,MCbID,opsedID,icsedID,opsedfID,
     +        o2negID,popID,optotID,opremID,opdenID,expID,
     +        udifID,vdifID,wdifID,gar84ID,d15ID
      integer m,inum
      save inum 
      data inum /100000/
* rjm

C Start code : ------------------------------------------------------------

c Set up values for spatial axes.  
c  mac, feb09

      pi = 3.14159265
      r2d = 180.0/pi
      radius = 6.37122e8

c Code for determining longitude taken from code in ocend.f used to
c  write the ocean output file (fort.40).  
      xt(1) = 0.0 - 1.5 * 360.0/real(imt-2)
      do i = 2, imt
        xt(i) = xt(i-1) + (dxt(i) * r2d / radius)
      enddo
      do i=1,imt
c Code copied from convert_averages/construct_axes.F90 to calc longitude
c  Gives the same answer for full resolution.  mac, feb09
c       lonts(i)=(real(i-2, 4)-0.5) * 360.0 / real(imt-2, 4)
       lonts(i)=xt(i)
      enddo

      do j=1,jmt
       latts(j)=phit(j)*r2d
      enddo

      do k=1,km
c convert cm to m for depth coordinate
       zts(k)=zdzz(k)*0.01
      enddo

      tim0(1)=ttsec/86400

c...  Create netCDF file
	inum=inum+1
	write(qfile,'(a4,i6,a3)') file,inum,".nc"
        write (*, *) "***  netCDF output ", qfile
      status = nf_create(trim(qfile), nf_write, ncid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not create file ", qfile
        write (*, *)
        stop
      end if
      	
c...  Define dimensions
      status = nf_def_dim(ncid, "lonts", imt-2, imtid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define dimension"
        write (*, *)
        write (*, *) "***  File      =  ", file
        write (*, *) "***  Dimension =  lonts"
        write (*, *)
        stop
      end if
      status = nf_def_dim(ncid, "latts", jmt-2, jmtid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define dimension"
        write (*, *)
        write (*, *) "***  File      =  ", file
        write (*, *) "***  Dimension =  jmt"
        write (*, *)
        stop
      end if
      status = nf_def_dim(ncid, "zts", km, kmid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define dimension"
        write (*, *)
        write (*, *) "***  File      =  ", file
        write (*, *) "***  Dimension =  km"
        write (*, *)
        stop
      end if
c      status = nf_def_dim(ncid, "lseg", 12, lsegid)
c      if (status .ne. nf_noerr) then
c        write (*, *)
c        write (*, *) "***  netCDF error: Could not define dimension"
c        write (*, *)
c        write (*, *) "***  File      =  ", file
c        write (*, *) "***  Dimension =  lseg"
c        write (*, *)
c        stop
c      end if
c      status = nf_def_dim(ncid, "time", NF_UNLIMITED, lmtid)
      status = nf_def_dim(ncid, "time", nf_unlimited, lmtid)
      if (status /= nf_noerr) stop "netCDF error is defining lmtid"

c...  Define variables
      status = nf_def_var(ncid, "lonts", nf_float, 1, imtid, lontsid)
      if (status /= nf_noerr) stop "netCDF error in write_averages"
      status = nf_put_att_text(ncid, lontsid, "units", 12, 
     & "degrees_east")
      if (status /= nf_noerr) stop "netCDF error in write_averages"
      status = nf_put_att_text(ncid, lontsid, "long_name", 27, 
     &                      "Longitude of TS grid points")
      if (status /= nf_noerr) stop "netCDF error in write_averages"
      status = nf_put_att_text(ncid, lontsid, "point_spacing", 4, 
     & "even")
      if (status /= nf_noerr) stop "netCDF error in write_averages"
      status = nf_put_att_text(ncid, lontsid, "modulo", 1, " ")
      if (status /= nf_noerr) stop "netCDF error in write_averages"

      status = nf_def_var(ncid, "latts", nf_float, 1, jmtid, lattsid)
      if (status /= nf_noerr) stop "netCDF error in write_averages"
      status = nf_put_att_text(ncid, lattsid, "units", 13, 
     & "degrees_north")
      if (status /= nf_noerr) stop "netCDF error in write_averages"
      status = nf_put_att_text(ncid, lattsid, "long_name", 26, 
     &                      "Latitude of TS grid points")
      if (status /= nf_noerr) stop "netCDF error in write_averages"
      status = nf_put_att_text(ncid, lattsid, "point_spacing", 4, 
     & "even")
      if (status /= nf_noerr) stop "netCDF error in write_averages"

      status = nf_def_var(ncid, "zts", nf_float, 1, kmid, ztsid)
      if (status /= nf_noerr) stop "netCDF error in write_averages1"
      status = nf_put_att_text(ncid, ztsid, "units", 1, 
     & "m")
      if (status /= nf_noerr) stop "netCDF error in write_averages2"
      status = nf_put_att_text(ncid, ztsid, "long_name", 23, 
     &                      "Depth of TS grid points")
      if (status /= nf_noerr) stop "netCDF error in write_averages3"
      status = nf_put_att_text(ncid, ztsid, "point_spacing", 6, 
     & "uneven")
      if (status /= nf_noerr) stop "netCDF error in write_averages4"
      status = nf_put_att_text(ncid, ztsid, "positive", 4, 
     & "down")
      if (status /= nf_noerr) stop "netCDF error in write_averages5"
      
      status = nf_def_var(ncid, "time", nf_float, 1, lmtid, timid)
      if (status /= nf_noerr) stop "netCDF error defining time"
      status = nf_put_att_text(ncid, timid, "long_name", 25, 
     &                      "Time since start of model")
      if (status /= nf_noerr) stop "netCDF error defining time"
      status = nf_put_att_text(ncid, timid, "units", 5, 
     &                      "days")
      if (status /= nf_noerr) stop "netCDF error defining time"

      status = nf_def_var(ncid, "itt", nf_int, 0, 0, ittid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define variable"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  itt"
        write (*, *)
        stop
      end if
c      status = nf_def_var(ncid, "EP", nf_double, 2, (/ imtid, jmtid /),
c     &                    pbid)
c      if (status .ne. nf_noerr) then
c        write (*, *)
c        write (*, *) "***  netCDF error: Could not define variable"
c        write (*, *)
c        write (*, *) "***  File     =  ", file
c        write (*, *) "***  Variable =  EP"
c        write (*, *)
c        stop
c      end if
c      status = nf_def_var(ncid,"pco2", nf_double, 2, (/ imtid, jmtid /),
c     &                    pid)
c      if (status .ne. nf_noerr) then
c        write (*, *)
c        write (*, *) "***  netCDF error: Could not define variable"
c        write (*, *)
c        write (*, *) "***  File     =  ", file
c        write (*, *) "***  Variable =  pco2"
c        write (*, *)
c        stop
c      end if

* rjm obgc
* Loop through the biogeochem. tracers.  
      do m=1,nt
       write(obgc_tr,'(a5,i1)') "obgc0",m
       if (m.gt. 9) write(obgc_tr,'(a4,i2)') "obgc",m
       print*,obgc_tr

       status = nf_def_var(ncid, obgc_tr, nf_double, 4,
     &         (/ imtid, jmtid, kmid, lmtid/), obgc_id(m) )
       if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
       end if
       write(obgc_tr,'(a1,i1)') "0",m
       if (m.gt. 9) write(obgc_tr,'(i2)') m
       status = nf_put_att_text(ncid, obgc_id(m), "long_name", 15, 
     &     trname(m) // obgc_tr )
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, obgc_id(m), "missing_value", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, obgc_id(m), "_FillValue", 
     &     nf_double, 1, missing_value )
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
      enddo

* Loop through the biogeochem. fluxes.  
      do m=1,nbio2d
       write(obgc_tr,'(a5,i1)') "bflx0",m
       if (m.gt. 9) write(obgc_tr,'(a4,i2)') "bflx",m
       print*,obgc_tr

       status = nf_def_var(ncid, obgc_tr, nf_double, 3,
     &         (/ imtid, jmtid, lmtid /), ofbgc_id(m) )
       if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
       end if
       write(obgc_tr,'(a1,i1)') "0",m
       if (m.gt. 9) write(obgc_tr,'(i2)') m
       status = nf_put_att_text(ncid, ofbgc_id(m), "long_name", 20, 
     &     trname(m) //"flux " // obgc_tr )
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, ofbgc_id(m), "missing_value", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, ofbgc_id(m), "_FillValue", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
      enddo
      
      if (sedfluxes) then
c... define POCsed variable  
       write(obgc_tr,'(a6)') "OPsedG"
       print*,obgc_tr
       status = nf_def_var(ncid, obgc_tr, nf_double, 4,
     &         (/ imtid, jmtid, kmid, lmtid /), opsedID)
       if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
       end if
       status = nf_put_att_text(ncid, opsedID, "long_name", 20, 
     &     "Sed flux of gen P" )
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, opsedID, "missing_value", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, opsedID, "_FillValue", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"

c... define PONsed variable  
       if (fix) then
       write(obgc_tr,'(a6)') "OPsedF"
       print*,obgc_tr
       status = nf_def_var(ncid, obgc_tr, nf_double, 4,
     &         (/ imtid, jmtid, kmid, lmtid /), opsedfID)
       if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
       end if
       status = nf_put_att_text(ncid, opsedfID, "long_name", 20, 
     &     "Sed flux of fix P" )
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, opsedfID, "missing_value", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, opsedfID, "_FillValue", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       endif

c... define PICsed variable  
       write(obgc_tr,'(a6)') "PICsed"
       print*,obgc_tr
       status = nf_def_var(ncid, obgc_tr, nf_double, 4,
     &         (/ imtid, jmtid, kmid, lmtid /), icsedID)
       if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
       end if
       status = nf_put_att_text(ncid, icsedID, "long_name", 20, 
     &     "Sed flux of C-inorg" )
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, icsedID, "missing_value", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, icsedID, "_FillValue", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"

c... define sedimentary denitrification variable
      if (den) then
      write(obgc_tr,'(a6)') "SedDen"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), sedDID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, sedDID, "long_name", 20, 
     &            "Sediment Den" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, sedDID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, sedDID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
      endif

c... define sedimentary sulfate reduction variable      
      write(obgc_tr,'(a6)') "SedSUL"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), sedSID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, sedSID, "long_name", 20, 
     &            "Sediment Sul" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, sedSID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, sedSID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"

      endif  ! --> sedfluxes 

      if (n_n15.gt.0) then
c... define d15Norg variable  
       write(obgc_tr,'(a6)') "d15org"
       print*,obgc_tr
       status = nf_def_var(ncid, obgc_tr, nf_double, 4,
     &         (/ imtid, jmtid, kmid, lmtid /), d15ID)
       if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
       end if
       status = nf_put_att_text(ncid, d15ID, "long_name", 20, 
     &     "Organic d15N" )
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, d15ID, "missing_value", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, d15ID, "_FillValue", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
      endif  ! --> n_n15.gt.0 

       
      if (vary_stoich) then
c... define carb2P variable  
       write(obgc_tr,'(a4)') "CtoP"
       print*,obgc_tr
       status = nf_def_var(ncid, obgc_tr, nf_double, 3,
     &         (/ imtid, jmtid, lmtid /), cpID)
       if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
       end if
       status = nf_put_att_text(ncid, cpID, "long_name", 20, 
     &     "C:P Ratio of OM " )
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, cpID, "missing_value", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, cpID, "_FillValue", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"

c... define NtoP variable
       if (n_no3.gt.0) then
       write(obgc_tr,'(a4)') "NtoP"
       print*,obgc_tr
       status = nf_def_var(ncid, obgc_tr, nf_double, 3,
     &         (/ imtid, jmtid, lmtid /), npID)
       if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
       end if
       status = nf_put_att_text(ncid, npID, "long_name", 20, 
     &     "N:P Ratio of OM " )
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, npID, "missing_value", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, npID, "_FillValue", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"

c... define -NO3toP variable
       if (den) then
       write(obgc_tr,'(a6)') "NO3toP"
       print*,obgc_tr
       status = nf_def_var(ncid, obgc_tr, nf_double, 3,
     &         (/ imtid, jmtid, lmtid /), no3remID)
       if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
       end if
       status = nf_put_att_text(ncid, no3remID, "long_name", 20, 
     &     "-NO3:P Ratio of OM " )
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, no3remID, "missing_value", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, no3remID, "_FillValue", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       endif !--> if den is on
       endif !--> if n_no3 tracer is on

c... define -O2toP variable  
       write(obgc_tr,'(a5)') "O2toP"
       print*,obgc_tr
       status = nf_def_var(ncid, obgc_tr, nf_double, 3,
     &         (/ imtid, jmtid, lmtid /), o2remID)
       if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
       end if
       status = nf_put_att_text(ncid, o2remID, "long_name", 20, 
     &     "-O2:P Ratio of OM " )
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, o2remID, "missing_value", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
       status = nf_put_att_double(ncid, o2remID, "_FillValue", 
     &     nf_double, 1, missing_value)
       if (status /= nf_noerr) stop "netCDF error in obgc_write"
      
      endif !--> vary_stoich

      if (den) then
c... define denitrification variable      
      write(obgc_tr,'(a5)') "WCDen"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), denID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, denID, "long_name", 20, 
     &            "Water column Den" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, denID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, denID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
      endif  
         
         
      if (fix) then
c... define N fixation variable      
      write(obgc_tr,'(a4)') "NFIX"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), fixID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, fixID, "long_name", 20, 
     &            "N Fixation" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, fixID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, fixID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
      endif

      if (diff_out) then
c... define N exp variable      
      write(obgc_tr,'(a4)') "Nexp"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), expID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, expID, "long_name", 20, 
     &            "N export" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, expID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, expID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
      endif


c... define negative oxygen variable      
      write(obgc_tr,'(a6)') "OXYNEG"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), o2negID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, o2negID, "long_name", 20, 
     &            "O2 Negatives" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, o2negID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, o2negID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"

      if (ReminDiagnostics) then
c... define Particulate Organic Phosphate variable      
      write(obgc_tr,'(a3)') "POP"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), popID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, popID, "long_name", 20, 
     &            "Initial Org Phos" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, popID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, popID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"

c... define Organic Phosphate after conservation variable      
      write(obgc_tr,'(a5)') "OPTOT"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), optotID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, optotID, "long_name", 20, 
     &            "Conserved Org Phos" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, optotID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, optotID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"

c... define Organic Phosphate remineralised by oxygen variable      
      write(obgc_tr,'(a5)') "OPREM"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), opremID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, opremID, "long_name", 20, 
     &            "Org P remin O2" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, opremID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, opremID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"

c... define Organic Phosphate remineralised by nitrate variable      
      write(obgc_tr,'(a5)') "OPDEN"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), opdenID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, opdenID, "long_name", 20, 
     &            "Org P remin NO3" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, opdenID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, opdenID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
      endif

      if (ReminTemp_Marsay .or. ReminTemp_Q10 .or. ReminPico) then
c... define Martin curve exponent      
      write(obgc_tr,'(a5)') "MCexp"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 3, 
     &         (/ imtid, jmtid, lmtid /), MCbID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, MCbID, "long_name", 20, 
     &            "Martin exponent" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, MCbID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, MCbID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
      endif
      
      if (diff_out) then
c... define zonal diffusion variable      
      write(obgc_tr,'(a4)') "Udif"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), udifID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, udifID, "long_name", 20, 
     &            "Zonal diffuision" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, udifID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, udifID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"

c... define meridional diffusion variable      
      write(obgc_tr,'(a4)') "Vdif"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), vdifID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, vdifID, "long_name", 20, 
     &            "Meridional diffusion" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, vdifID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, vdifID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"

c... define vertical diffusion variable      
      write(obgc_tr,'(a4)') "Wdif"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), wdifID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, wdifID, "long_name", 20, 
     &            "Vertical diffusion" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, wdifID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, wdifID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"

c... define vertical diffusion according to Gargett      
      write(obgc_tr,'(a5)') "Gar84"
      print*, obgc_tr
      status = nf_def_var(ncid, obgc_tr, nf_double, 4, 
     &         (/ imtid, jmtid, kmid, lmtid /), gar84ID)
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not define variable"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable =  ", obgc_tr
         write (*, *)
         stop
      end if
         status = nf_put_att_text(ncid, gar84ID, "long_name", 20, 
     &            "Gargett vert diffusion" )
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, gar84ID, "missing_value", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
         status = nf_put_att_double(ncid, gar84ID, "_FillValue", 
     &            nf_double, 1, missing_value)
         if (status /= nf_noerr) stop "netCDF error in obgc_write"
      endif


c... define basin mask (depth) variable         
      status = nf_def_var(ncid, "kmt", nf_int, 2, (/ imtid, jmtid /),
     &                    kmtid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not define variable"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  kmt"
        write (*, *)
        stop
      end if


c...  Put global attributes
      status = nf_put_att_text(ncid, nf_global, "version", 10, version)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not put attribute"
        write (*, *)
        write (*, *) "***  File      =  ", file
        write (*, *) "***  Attribute =  version"
        write (*, *)
        stop
      end if
      status = nf_put_att_text(ncid, nf_global, "runtype", 3, runtype)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not put attribute"
        write (*, *)
        write (*, *) "***  File      =  ", file
        write (*, *) "***  Attribute =  runtype"
        write (*, *)
        stop
      end if



c...  Exit define mode
	print*,itt,missing_value
      status = nf_enddef(ncid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not exit define mode"
        write (*, *)
        write (*, *) "***  File =  ", file
        write (*, *)
*        stop
      end if

c...  Write the timestep counter [ITT]
c      status = nf_put_var_int(ncid, ittid, itt)
      status = nf_put_var_int(ncid, ittid, int(ttsec/86400))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  itt"
        write (*, *)
        stop
      end if

c      status = nf_put_var_real(ncid, timid, tim0(1))
       status = nf_put_vara_real(ncid, timid, 
     &  (/1/), (/1/), tim0(1))
       if (status /= nf_noerr) stop "netCDF error in writing time"

      status = nf_put_var_real(ncid, lontsid, lonts(2:imt-1))
      if (status /= nf_noerr) stop "netCDF error in write_averages"

      status = nf_put_var_real(ncid, lattsid, latts(2:jmt-1))
      if (status /= nf_noerr) stop "netCDF error in write_averages"

      status = nf_put_var_real(ncid, ztsid, zts)
      if (status /= nf_noerr) stop "netCDF error in write_averages"

* rjm obgc
c...  Write the obgc tracers [obgc01]
      do m=1,nt
	write(obgc_tr,'(a5,i1)') "obgc0",m
	if (m.gt. 9) write(obgc_tr,'(a4,i2)') "obgc",m
	
      status =nf_put_vara_double(ncid, obgc_id(m),
     & (/1,1,1,1/), (/imt-2, jmt-2, km, 1/), 
     & ave_tr(2:imt-1, 2:jmt-1, :, m))
*      status =nf_put_var_double(ncid, 
*	1  obgc_id(m),ave_tr(2:imtm1,2:jmtm1,:,m))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable = ",obgc_tr
        write (*, *)
        stop
      end if

      enddo
      
c...  Write the obgc fluxes [bflx01]
	do m=1,nbio2d
	write(obgc_tr,'(a5,i1)') "bflx0",m
	if (m.gt. 9) write(obgc_tr,'(a4,i2)') "bflx",m
      status =nf_put_vara_double(ncid, ofbgc_id(m), 
     & (/1,1,1/), (/imt-2, jmt-2, 1/), 
     & ave_flux(2:imt-1, 2:jmt-1, m))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable = ",obgc_tr
        write (*, *)
        stop
      end if

      enddo
* rjm
*pjb

      if (sedfluxes) then
c...  Write the opsed variable
      write(obgc_tr,'(a6)') "OPsedG"
      status =nf_put_vara_double(ncid, opsedID, 
     & (/1,1,1,1/), (/imt-2, jmt-2, km, 1/), 
     & ave_opsed(2:imt-1, 2:jmt-1, :, 1))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable = ",obgc_tr
        write (*, *)
        stop
      end if
      if (fix) then
c...  Write the opsed_f variable
      write(obgc_tr,'(a6)') "OPsedF"
      status =nf_put_vara_double(ncid, opsedfID, 
     & (/1,1,1,1/), (/imt-2, jmt-2, km, 1/), 
     & ave_opsed_f(2:imt-1, 2:jmt-1, :, 1))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable = ",obgc_tr
        write (*, *)
        stop
      end if
      endif ! --> fix
c...  Write the icsed variable
      write(obgc_tr,'(a6)') "PICsed"
      status =nf_put_vara_double(ncid, icsedID, 
     & (/1,1,1,1/), (/imt-2, jmt-2, km, 1/), 
     & ave_icsed(2:imt-1, 2:jmt-1, :, 1))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable = ",obgc_tr
        write (*, *)
        stop
      end if
c... write sedimentary denitrification variable to file
      if (den) then
      write(obgc_tr,'(a6)') "SedDen"
      status = nf_put_vara_double( ncid, sedDID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_Sden(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if
      endif
c... write sedimentary sulfate reduction variable to file
      write(obgc_tr,'(a6)') "SedSul"
      status = nf_put_vara_double( ncid, sedSID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_Ssul(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if
      endif ! --> sedfluxes

      if (n_n15.gt.0) then
c...  Write the d15Norg variable
      write(obgc_tr,'(a6)') "d15org"
      status =nf_put_vara_double(ncid, d15ID, 
     & (/1,1,1,1/), (/imt-2, jmt-2, km, 1/), 
     & ave_d15N(2:imt-1, 2:jmt-1, :, 1))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable = ",obgc_tr
        write (*, *)
        stop
      end if
      endif ! --> n_n15.gt.0
      
      if (vary_stoich) then
c...  Write the CtoP variable
      write(obgc_tr,'(a4)') "CtoP"
      status =nf_put_vara_double(ncid, cpID, 
     & (/1,1,1/), (/imt-2, jmt-2, 1/), 
     & ave_carb2P(2:imt-1, 2:jmt-1, 1))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable = ",obgc_tr
        write (*, *)
        stop
      end if

      if (n_no3.gt.0) then
c...  Write the NtoP variable
      write(obgc_tr,'(a4)') "NtoP"
      status =nf_put_vara_double(ncid, npID, 
     & (/1,1,1/), (/imt-2, jmt-2, 1/), 
     & ave_NtoP(2:imt-1, 2:jmt-1, 1))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable = ",obgc_tr
        write (*, *)
        stop
      end if
c...  Write the -NO3toP variable
      if (den) then
      write(obgc_tr,'(a6)') "NO3toP"
      status =nf_put_vara_double(ncid, no3remID, 
     & (/1,1,1/), (/imt-2, jmt-2, 1/), 
     & ave_no3rem(2:imt-1, 2:jmt-1, 1))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable = ",obgc_tr
        write (*, *)
        stop
      endif
      endif !--> den
      endif !--> n_no3

c...  Write the -O2toP variable
      write(obgc_tr,'(a5)') "O2toP"
      status =nf_put_vara_double(ncid, o2remID, 
     & (/1,1,1/), (/imt-2, jmt-2, 1/), 
     & ave_o2rem(2:imt-1, 2:jmt-1, 1))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable = ",obgc_tr
        write (*, *)
        stop
      end if
      endif !--> vary_stoich

      if (den) then
c... write the denitrification variable to file
      write(obgc_tr,'(a5)') "WCDen"
      status = nf_put_vara_double( ncid, denID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_Pden(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if
      endif

      if (diff_out) then
c... write the N export variable to file
      write(obgc_tr,'(a4)') "Nexp"
      status = nf_put_vara_double( ncid, expID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_exp(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if
      endif

      if (fix) then
c... write the N fixation variable to file
      write(obgc_tr,'(a4)') "NFIX"
      status = nf_put_vara_double( ncid, fixID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_fix(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if
      endif
      
c... write the oxygen negatives variable to file
      write(obgc_tr,'(a6)') "OXYNEG"
      status = nf_put_vara_double( ncid, o2negID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_oxyneg(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if

      if (ReminDiagnostics) then
c... write the POP variable to file
      write(obgc_tr,'(a3)') "POP"
      status = nf_put_vara_double( ncid, popID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_pop(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if

c... write the OPTOT variable to file
      write(obgc_tr,'(a5)') "OPTOT"
      status = nf_put_vara_double( ncid, optotID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_optot(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if

c... write the oprem variable to file
      write(obgc_tr,'(a5)') "OPREM"
      status = nf_put_vara_double( ncid, opremID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_oprem(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if

c... write the OPTOT variable to file
      write(obgc_tr,'(a5)') "OPDEN"
      status = nf_put_vara_double( ncid, opdenID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_opden(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if
      endif

      if (ReminTemp_Marsay .or. ReminTemp_Q10 .or. ReminPico) then
c... write the N fixation variable to file
      write(obgc_tr,'(a5)') "MCexp"
      status = nf_put_vara_double( ncid, MCbID,
     &         (/1,1,1/), (/imt-2, jmt-2, 1/),
     &         ave_MCb(2:imt-1, 2:jmt-1, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if
      endif
      
      if (diff_out) then
c... write the zonal diffusion variable to file
      write(obgc_tr,'(a4)') "Udif"
      status = nf_put_vara_double( ncid, udifID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_udif(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if

c... write the meridional diffusion variable to file
      write(obgc_tr,'(a4)') "Vdif"
      status = nf_put_vara_double( ncid, vdifID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_vdif(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if

c... write the vertical diffusion variable to file
      write(obgc_tr,'(a4)') "Wdif"
      status = nf_put_vara_double( ncid, wdifID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_wdif(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if

c... write the OPTOT variable to file
      write(obgc_tr,'(a5)') "Gar84"
      status = nf_put_vara_double( ncid, gar84ID,
     &         (/1,1,1,1/), (/imt-2, jmt-2, km, 1/),
     &         ave_Gar84(2:imt-1, 2:jmt-1, :, 1) )
      if (status .ne. nf_noerr) then
         write (*, *)
         write (*, *) "***  netCDF error: Could not write data"
         write (*, *)
         write (*, *) "***  File     =  ", file
         write (*, *) "***  Variable = ",obgc_tr
         write (*, *)
         stop
      end if
      endif
*pjb      

c...  Write the number of vertical levels on the tracer grid [KMT]
      status = nf_put_var_int(ncid, kmtid, kmt(2:imt-1,2:jmt-1))
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not write data"
        write (*, *)
        write (*, *) "***  File     =  ", file
        write (*, *) "***  Variable =  kmt"
        write (*, *)
        stop
      end if

c...  Close netCDF file
      status = nf_close(ncid)
      if (status .ne. nf_noerr) then
        write (*, *)
        write (*, *) "***  netCDF error: Could not close file ", file
        write (*, *)
        stop
      end if

      return
      end subroutine
