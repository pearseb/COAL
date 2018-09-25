c nsib comparison is modified to lsm_type = "nsib "
c AJA 2009/01/22
c
c Modified to enable 100,000-year runs.
c SJP 2008/02/28
c
c Removed the line "include 'OCEANPARS.f'", as none of the parameters defined
c in this header file are referenced.
c SJP 2007/05/31
c
c Modified for v6.2-sjp1.4.
c SJP 2004/08/09
c
c $Log: finterp.f,v $
c Revision 1.12  2001/02/22 05:34:45  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.11  1999/05/20 06:23:57  rot032
c HBG changes to V5-2
c
c Revision 1.10  1997/12/23  00:23:39  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.9  1997/12/17  23:23:06  ldr
c Changes from MRD for parallel execution on NEC.
c
c Revision 1.8.1.1  1997/12/19  02:03:19  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.8  1997/01/10  06:28:08  ldr
c Add new NSIB fields.
c
c Revision 1.7  1997/01/10  06:06:36  ldr
c Replace pr and pi by prr and pri respectively, so that these routines work
c again.
c
c Revision 1.6  1996/10/24  01:44:21  ldr
c Don't set ndays to 0 when writing interpolated RS file.
c
c Revision 1.5  1996/10/24  01:03:33  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.4  1996/06/13  02:08:59  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.3  1996/03/21  03:41:39  ldr
c Merge of TIE and LDR changes.
c
c Revision 1.2  1996/03/21  03:19:17  ldr
c Changes from TIE to remove nsemilag.eq.0 and tidy physical constants
c
c Revision 1.1.1.1  1996/03/21  03:29:07  ldr
c A couple of housekeeping tidy-ups from LDR.
c
c Revision 1.1  1996/02/19  04:06:11  ldr
c Initial revision
c

c Interpolate restart file data to NLW vertical levels and write in
c new unified Rhomboidal/Triangular format.
c
c     INPUT/OUTPUT
c     Input:   from common/fewflags in FEWFLAGS.f
c                  qcloud - T if using prognostic cloud scheme
c                  qflux  - set to T for qflux run, default F
c                  semice - if T, use the Semtner sea-ice model
c
c              from common/files in FILES.f
c                  orfilename - name of output restart file
c
c              from common/lsmi in LSMI.f
c                  imsl - land surface type indicator
c
c              from common/masiv1 in this subroutine
c                  ronmx - global array of measure of east-west surface stress
c                  sonmx - global array of measure of north-sth surface stress
c                  iphys ) - current phys time step and previous phys time
c                  iphysm)   step, these indices swap around
c
c              from common/masiv3 in MASIV3.f
c                  rgsav - incoming long-wave radiation
c                  sgamp - amplitude of net solar radiation at the ground
c                  sgsav - net solar radiation at ground
c
c              from common/masiv4 in MASIV4.f
c                  savegrid - global surface variables
c
c              from common/masiv4 in this subroutine
c                  statsice - global statistics relating to sea-ice scheme
c
c              from common/pmasiv3 in PMASIV3.f
c                  prgsav - rgsav over leads(see above)
c                  psgamp - sgamp over leads  "    "
c                  psgsav - sgsav over leads  "    "
c
c              from common/qflxdat in QFLXDAT.f
c                  ochfa - ocean heat flux average
c
c              from common/surf in this subroutine 
c                  aftfh  - reciprocal of air resistance for ground 
c                           covered with foliage
c                  aftfhg - reciprocal of air resistance for bare ground
c                  ocondx - precipitation at previous timestep
c                  opreci - frozen precipitation at previous timestep
c
c              from common/surf1 in SURF1.f
c                  mc   - moisture depth on canopy
c                  tgf  - canopy temperature
c                  tgg  - bare ground temperature
c
c              from common/timex in TIMEX.f
c                  mstep - timestep in minutes
c
c              from argument
c                  exptyp - restart file header
c
c     In/Out:  from common/cnsta in CNSTA.f
c                  sig - sigma values
c
c              from common/fldmri in FLDMRI.f
c                  pmi  - spectral surface pressure, imaginary part
c                  pmr  - spectral surface pressure, real part
c                  psimi - stream function, imaginary part
c                  psimr - stream function, real part 
c                  temi - spectral temp, imaginary part at t-1
c                  temr - spectral temp, real part at t-1
c                  xhimi - velocity potential at t-1, imaginary part
c                  xhimr - velocity potential at t-1, real part
c
c              from common/fldri in FLDRI.f
c                  pi - spectral pressure, imaginary part
c                  pr - spectral pressure, real part
c                  psii - spectral stream function, imaginary part
c                  psir - spectral stream function, real part
c                  tei - spectral temp, imaginary part
c                  ter - spectral temp, real part
c                  xhii - velocity potental, imaginary part
c                  xhir - velocity potental, real part
c
c              from common/mountn in this subroutine
c                  phisi- spectral surface geopotential height,imaginary part
c                  phisr- spectral surface geopotential height, real part
c
c              from common/physical in PHYSICAL.f
c                  ron - measure of east-west stress at surface
c                  son - measure of north-south stress at surface
c
c              from common/qcloud1 in QCLOUD1.f
c                  qfb  - cloud ice at t     qlb  - cloud liquid water at t
c                  qfbm - cloud ice at t-1   qlbm - cloud liquid water at t-1
c                  cfb  - cloud fraction
c
c              from common/rmgrid in RMGRID.f
c                  rgt  - moisture tendency
c                  rmg  - pressure weighted moistures at current timestep
c                  rmmg - pressure weighted mositures at previous timestep
c
c              from common/rshead in this subroutine
c                  header2 - second header on restart file
c
c              from common/timex in TIMEX.f
c                  ndays - day counter for model run
c
c              from common/uvpgd in UVPGD.f
c                  sdot - d(sigma)/d(t)
c

      subroutine finterp(exptyp)

      implicit none
C Global parameters
      include 'PARAMS.f'
      include 'PHYSPARAMS.f'
      include 'VERSION.f'
      integer nln
      parameter (nln=24) !New no. of levels

C Argument list
      character*50 exptyp

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CNSTA.f'
      include 'DIFM.f'
      include 'FEWFLAGS.f'
      include 'FILES.f'
      include 'FLDRI.f'
      include 'FLDMRI.f'
      include 'LSMI.f'
      include 'MASIV3.f'
      include 'MASIV4.f'
      include 'PMASIV3.f'
      include 'QCLOUD1.f'
      include 'QFLXDAT.f'
      include 'RMGRID.f'
      include 'SURF1.f'
      include 'TIMEX.f'
      include 'UVPGD.f'

      real ronmx,sonmx
      integer iphysm,iphys
      common/masiv1/ronmx(ln2,nl,lat,2),sonmx(ln2,nl,lat,2)
     &             ,iphysm,iphys

      real statsice
      common/masiv5/statsice(ln2,14,lat)

      real phisr,phisi
      common/mountn/phisr(lw,mw),phisi(lw,mw)

      character*50 header2
      common/rshead/header2

      real aftfh,aftfhg,ocondx,opreci,pmc
      integer isoilm
      common/surf/aftfh(ln2,lat),aftfhg(ln2,lat),ocondx(ln2,lat)
     &          ,isoilm(ln2,lat),opreci(ln2,lat)
     &          ,pmc(ln2,lat)

C Local work arrays and variables
      real qout(ln2,lat)
      real stice(ln2,9)
      real sign(nln),sighn(0:nln)
      real fieldn1(lw,mw,nln),fieldn2(lw,mw,nln)
      real fieldng(ln2,nln,lat)
      real sdotng(lon,nln,lat2)
      real fieldnl(ln2,nln)
      real griddata(ln2,ngrid),ron(ln2,nl),son(ln2,nl)
      complex wfld(nw)

      integer irestf
      integer k
      integer lg
      integer mg

      real alf
      real den
      real rk
      real sigt

C Local data, functions etc

C Start code : ----------------------------------------------------------

c---- To vertically interpolate model data to new levels
c---- and write a restart file for the Mk3 model.

c Compute new sigma levels (HBG formula)

      sigt=0.004458
      den=.75*( 1./nln-.5 - .5*(2./nln-1.)**3 )
      alf=(sigt-.25 -.25*(2./nln-1.)**3 )/den
c set up general cubic half-sigma levels (normal formula bottom up)
      do k=0,nln
        rk=(.5*(nln+1)-k-.5)/nln
        sighn(k)=.5 + 1.5*alf*rk - 2.*(3.*alf-2.)*rk**3
      enddo
      do k=1,nln
        sign(k)=0.5*(sighn(k-1)+sighn(k))
      enddo
      print *,'New sigh values: ',(sighn(k),k=0,nln)
      print *,'New sig values: ',(sign(k),k=1,nln)


      irestf=18

c****  OPEN THE RESTART FILES

      open(unit=irestf,file=orfilename,form='unformatted')

C**** READ/WRITE HEADER

      write(6,93)exptyp
  93  format(2x,'Header on RS file is : ',a50)
      write(irestf)exptyp

c Second header

#ifdef ALPH
      header2(1:4)="ALPH"
#else
      header2(1:4)="LINU"
#endif

      header2(5:5)=' '
      header2(6:6)=trunc        !'R' or 'T' truncation type
      write(header2(7:8),'(i2)')lw-1
      write(header2(9:11),'(i3)')nln
      header2(12:14)='L  '
      header2(15:24)=version

C****
      write(6,*)'WARNING: NOT setting ndays=0 in finterp.'
      write(0,*)'WARNING: NOT setting ndays=0 in finterp.'
c      ndays=0
      write (irestf) ndays,mstep
      write(irestf)header2

c**** output spectral geopotential
      call setw (phisr,phisi,wfld)
      write (irestf) wfld

C**** READ/WRITE PSI(flux)
      call fint(psir,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn1)                     !Output
      call fint(psii,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn2)                     !Output
      do 10 k=1,nln
         call setw (fieldn1(1,1,k),fieldn2(1,1,k),wfld)
 10      write (irestf) wfld
      call fint(psimr,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn1)                     !Output
      call fint(psimi,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn2)                     !Output
      do 12 k=1,nln
         call setw (fieldn1(1,1,k),fieldn2(1,1,k),wfld)
 12      write (irestf) wfld

C**** READ/WRITE XHI(flux)
      call fint(xhir,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn1)                     !Output
      call fint(xhii,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn2)                     !Output
      do 14 k=1,nln
         call setw (fieldn1(1,1,k),fieldn2(1,1,k),wfld)
 14      write (irestf) wfld
      call fint(xhimr,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn1)                     !Output
      call fint(xhimi,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn2)                     !Output
      do 16 k=1,nln
         call setw (fieldn1(1,1,k),fieldn2(1,1,k),wfld)
 16      write (irestf) wfld

C**** NOW SET UP TO WRITE (P*T') NOT (P*T)
      call fint(ter,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn1)                     !Output
      call fint(tei,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn2)                     !Output
      do 18 k=1,nln
         call setw (fieldn1(1,1,k),fieldn2(1,1,k),wfld)
 18      write (irestf) wfld
      call fint(temr,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn1)                     !Output
      call fint(temi,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn2)                     !Output
      do 20 k=1,nln
         call setw (fieldn1(1,1,k),fieldn2(1,1,k),wfld)
 20      write (irestf) wfld

C**** READ/WRITE P*
      call setw (prr,pri,wfld)
      write (irestf) wfld
      call setw (pmr,pmi,wfld)
      write (irestf) wfld

C**** OUTPUT GRID MOISTURE + LAST TENDENCY
      call fintg(rmg,ln2,lat,nl,nln,sig,sign, !Inputs
     &          fieldng)                     !Output
      do 30 k=1,nln
      do 31 lg=1,lat
      do 31 mg=1,ln2
   31 qout(mg,lg)=fieldng(mg,k,lg)
   30 write (irestf)qout

      call fintg(rmmg,ln2,lat,nl,nln,sig,sign, !Inputs
     &          fieldng)                     !Output
      do 32 k=1,nln
      do 33 lg=1,lat
      do 33 mg=1,ln2
   33 qout(mg,lg)=fieldng(mg,k,lg)
   32 write (irestf)qout

      call fintg(rgt,ln2,lat,nl,nln,sig,sign, !Inputs
     &          fieldng)                     !Output
      do 34 k=1,nln
      do 35 lg=1,lat
      do 35 mg=1,ln2
   35 qout(mg,lg)=fieldng(mg,k,lg)
   34 write (irestf)qout
C**** OUTPUT SURFACE TYPE INDICATOR ARRAY
      write (irestf) imsl
C****
C****     OUTPUT PHYSICAL STATISTICS
C****
      do 80 lg=1,lat
      do 855 k=1,ngrid
      do 855 mg=1,ln2
 855  griddata(mg,k)=savegrid(mg,k,lg)
C**** OUTPUT PHYSICAL RESTART DATA TO RESTART FILE.
      write (irestf) griddata
      if(nln.eq.nl)then
        do k=1,nl
          do mg=1,ln2
            ron(mg,k)=ronmx(mg,k,lg,iphysm)
            son(mg,k)=sonmx(mg,k,lg,iphysm)
          enddo
        enddo
        write(irestf)ron,son
      else
        do k=1,nln
          do mg=1,ln2
            fieldnl(mg,k)=0. !Setting ron,son to 0
          enddo
        enddo
        write(irestf)fieldnl,fieldnl
      endif
 80   continue

c Extra stuff for SPO's sea-ice model

      if(semice)then
        do lg=1,lat
          do k=1,9
              do mg=1,ln2
                stice(mg,k)=statsice(mg,k,lg)
              enddo
          enddo
          write(irestf)stice
        enddo
      endif

c**** output dpsi,dxhi (spectral diffusion energy
c****    loss conversion to temp tendancy terms)
c**** Make split vertical levels standard for all resolutions

      call fint(dpsir,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn1)                     !Output
      call fint(dpsii,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn2)                     !Output
      do k=1,nln
        call setw (fieldn1(1,1,k),fieldn2(1,1,k),wfld)
        write (irestf) wfld
      enddo

      call fint(dxhir,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn1)                     !Output
      call fint(dxhii,lw,mw,nl,nln,sig,sign, !Inputs
     &          fieldn2)                     !Output
      do k=1,nln
        call setw (fieldn1(1,1,k),fieldn2(1,1,k),wfld)
        write (irestf) wfld
      enddo

c  save vertical velocity (sdot)

      call fsdot(sdot,lon,lat2,nl,nln,sig,sign, !Inputs
     &           sdotng)                        !Output
      write(irestf)sdotng


      if ( qflux) write(irestf) ochfa

c Write nsib fields
      if(lsm_type .eq. "nsib ")then
c.... Mk3 restart data for sib
        write (irestf) wb
        write (irestf) wbice
        write (irestf) tggsl
        write (irestf) tggsn
        write (irestf) ssdnn
        write (irestf) ssdn3
        write (irestf) smass
        write (irestf) gflux
        write (irestf) sgflux
        write (irestf) isflag
        write (irestf) snage
        write (irestf) osnowd
        write (irestf) tgf
        write (irestf) mc
        write (irestf) aftfh
        write (irestf) aftfhg
        write (irestf) ocondx
        write (irestf) pmc
      endif

c Add radiation arrays to restart file.

      write(irestf)sgsav
      write(irestf)rgsav
      write(irestf)sgamp
      write(irestf)psgsav
      write(irestf)prgsav
      write(irestf)psgamp

c Write cloud water and cloud ice

      if(qcloud)then
        call fintg(qlb,ln2,lat,nl,nln,sig,sign, !Inputs
     &             fieldng)                     !Output
        write(irestf)fieldng

        call fintg(qlbm,ln2,lat,nl,nln,sig,sign, !Inputs
     &             fieldng)                      !Output
        write(irestf)fieldng

        call fintg(qfb,ln2,lat,nl,nln,sig,sign, !Inputs
     &             fieldng)                     !Output
        write(irestf)fieldng

        call fintg(qfbm,ln2,lat,nl,nln,sig,sign, !Inputs
     &             fieldng)                      !Output
        write(irestf)fieldng

        call fintg(cfb,ln2,lat,nl,nln,sig,sign, !Inputs
     &             fieldng)                     !Output
        write(irestf)fieldng

        write(irestf)opreci

      endif
C****
C**** --- END OF REWRITE OF RESTART FILE ---
C****
      endfile irestf
      close(unit=irestf)
      write(6,140)irestf,ndays
  140 format(t10,'restart file(',i2,') written at day',i9)

      return
      end
c******************************************************************************
      subroutine fint(field,nx,ny,nl,nln,sig,sign, !Inputs
     &                 fieldn)                     !Output

      implicit none
C Global parameters

C Argument list
      integer nl    !No. of vertical levels (old,new)
      integer nx,ny !Horizontal array size (x,y)
      real field(nx,ny,nl)
      real sig(nl)
      integer nln   !No. of vertical levels (new)
      real sign(nln)
      real fieldn(nx,ny,nln)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      integer i
      integer j
      integer k
      integer kk
      integer k1
      integer k2

      real fsig

C Local data, functions etc

C Start code : ----------------------------------------------------------

      if(nln.eq.nl)then
        do k=1,nl
          do j=1,ny
            do i=1,nx
              fieldn(i,j,k)=field(i,j,k)
            enddo
          enddo
        enddo
      else
        do k=1,nln
          if(sign(k).ge.sig(1))then !just use old k=1 field
            do j=1,ny
              do i=1,nx
                fieldn(i,j,k)=field(i,j,1)
              enddo
            enddo
          elseif(sign(k).le.sig(nl))then !just use old k=nl field
            do j=1,ny
              do i=1,nx
                fieldn(i,j,k)=field(i,j,nl)
              enddo
            enddo
          else                  !do linear interpolation
            kk=1
            do while (sig(kk).gt.sign(k))
              kk=kk+1
            enddo
            k1=kk-1
            k2=kk               !New sign(k) lies b/w old levels k1 and k2
            fsig=(sig(k1)-sign(k))/(sig(k1)-sig(k2))
            do j=1,ny
              do i=1,nx
                fieldn(i,j,k)=fsig*field(i,j,k2)+(1-fsig)*field(i,j,k1)
              enddo
            enddo
          endif
        enddo
      endif

      return
      end
c******************************************************************************
      subroutine fintg(field,nx,ny,nl,nln,sig,sign, !Inputs
     &                 fieldn)                !Output

      implicit none
C Global parameters

C Argument list
      integer nl    !No. of vertical levels (old)
      integer nx,ny !Horizontal array size (x,y)
      real field(nx,nl,ny)
      real sig(nl)
      integer nln   !No. of vertical levels (new)
      real sign(nln)
      real fieldn(nx,nln,ny)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables

C Local data, functions etc
      integer i
      integer j
      integer k
      integer kk
      integer k1
      integer k2

      real fsig

C Start code : ----------------------------------------------------------

      if(nln.eq.nl)then
        do k=1,nl
          do j=1,ny
              do i=1,nx
                fieldn(i,k,j)=field(i,k,j)
              enddo
          enddo
        enddo
      else
        do k=1,nln
          if(sign(k).ge.sig(1))then !just use old k=1 field
            do j=1,ny
                do i=1,nx
                  fieldn(i,k,j)=field(i,1,j)
                enddo
            enddo
          elseif(sign(k).le.sig(nl))then !just use old k=nl field
            do j=1,ny
                do i=1,nx
                  fieldn(i,k,j)=field(i,nl,j)
                enddo
            enddo
          else                  !do linear interpolation
            kk=1
            do while (sig(kk).gt.sign(k))
              kk=kk+1
            enddo
            k1=kk-1
            k2=kk               !New sign(k) lies b/w old levels k1 and k2
            fsig=(sig(k1)-sign(k))/(sig(k1)-sig(k2))
            do j=1,ny
                do i=1,nx
                  fieldn(i,k,j)=fsig*field(i,k2,j)
     &                 +(1-fsig)*field(i,k1,j)
                enddo
            enddo
          endif
        enddo
      endif      

      return
      end
c******************************************************************************
      subroutine fsdot(field,nx,ny,nl,nln,sig,sign, !Inputs
     &                 fieldn)                !Output

      implicit none
C Global parameters

C Argument list
      integer nl    !No. of vertical levels (old)
      integer nx,ny !Horizontal array size (x,y)
      real field(nx,nl,ny)
      real sig(nl)
      integer nln   !No. of vertical levels (new)
      real sign(nln)
      real fieldn(nx,nln,ny)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables
      integer i
      integer j
      integer k
      integer kk
      integer k1
      integer k2

      real fsig

C Local data, functions etc

C Start code : ----------------------------------------------------------

      if(nln.eq.nl)then
        do k=1,nl
            do j=1,ny
              do i=1,nx
                fieldn(i,k,j)=field(i,k,j)
              enddo
            enddo
        enddo
      else
        do k=1,nln
          if(sign(k).ge.sig(1))then !just use old k=1 field
              do j=1,ny
                do i=1,nx
                  fieldn(i,k,j)=field(i,1,j)
                enddo
              enddo
          elseif(sign(k).le.sig(nl))then !just use old k=nl field
              do j=1,ny
                do i=1,nx
                  fieldn(i,k,j)=field(i,nl,j)
                enddo
            enddo
          else                  !do linear interpolation
            kk=1
            do while (sig(kk).gt.sign(k))
              kk=kk+1
            enddo
            k1=kk-1
            k2=kk               !New sign(k) lies b/w old levels k1 and k2
            fsig=(sig(k1)-sign(k))/(sig(k1)-sig(k2))
              do j=1,ny
                do i=1,nx
                  fieldn(i,k,j)=fsig*field(i,k2,j)
     &                 +(1-fsig)*field(i,k1,j)
                enddo
              enddo
          endif
        enddo
      endif      

      return
      end
