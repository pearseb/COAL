c $Log: landrun.f,v $
c Revision 1.13  2001/10/12 02:06:56  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.12  2001/02/22 05:34:36  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.11  2001/02/12 05:39:40  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.10  1998/05/27 02:07:34  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.9  1997/12/19  02:03:09  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.8  1997/03/06  23:33:59  ldr
c Corrections from HBG to the recent tidy-ups to ocean routines.
c
c Revision 1.7  1996/10/24  01:03:29  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.6  1996/06/13  02:08:52  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.5  1994/08/08  16:22:46  ldr
c Fix up RCS header at top of file.
c
c     INPUT/OUTPUT
c
c     Output:  from common/rungrid in this subroutine
c                  landlg, landmg, minlgx  ) runoff relocation
c                  minmgx, spread, topland ) data
c
      subroutine landrun

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'FEWFLAGS.f'
      include 'TIMEX.f'

      real topland,spread
      integer landlg,landmg,minlgx,minmgx
      common/rungrid/topland(nlanp),spread(nlanp,8)
     & ,landlg(nlanp),landmg(nlanp),minlgx(ndips),minmgx(ndips)

C Global data for T63 river routing scheme with delay
      real gradwe(lon,lat2),gradsn(lon,0:lat2)
      real topminT63(lon,lat2)
C at T63, 309 points are troughs
      integer ntrpnt(0:309)
      integer ntrmg(0:309),ntrlg(0:309)
      integer ntre(0:309),ntrn(0:309),ntrw(0:309),ntrs(0:309)
      logical basin(0:309)
      real distf(0:309)
      common/rivflow/gradwe,gradsn
     & ,ntrpnt,ntrmg,ntrlg,ntre,ntrn,ntrw,ntrs,basin,distf

C Local work arrays and variables
      integer ierr
      integer kk
      integer lg
      integer lgn
      integer lgs
      integer mg
      integer mge
      integer mgw

      real tdl10
      real topdiff
      real topij
      real topx
      real ugradsn
      real ugradwe

      logical teste
      logical testn
      logical testw
      logical tests

      character*20 lfn

C Local data, functions etc

C Start code : ----------------------------------------------------------

C Read data from files for routing runoff to the oceans.

      IF(newriver)THEN

c For T63 can use new river routing scheme
c Get files to drive this system. Topography gradients first.

      write(lfn,110)mw-1
  110 format('landgrads1T',i2)
      open(unit=7,file=lfn,form='formatted',status='old',iostat=ierr)
        call errcheck(ierr,'landgrads???     ','landrun   ')
        read(7,*)
        read(7,*)gradwe ! (lon,lat2)
        read(7,*)
        read(7,*)gradsn ! (lon,0:lat2)
C Get the (mg,lg) index for each T63 trough
        read(7,*)
        do kk=0,309
          read(7,115)ntrmg(kk),ntrlg(kk)
        enddo
  115 format(6x,i3,2x,i3)
C Get the T63 trough flow directional data
        read(7,*)
        do kk=0,309
          read(7,120)ntre(kk),ntrn(kk),ntrw(kk),ntrs(kk)
     &       ,basin(kk)
  120 format(3x,4i3,4x,l1)
        enddo
C Get the topography minimum per T63 box data set
        read(7,*)
        read(7,*)topminT63
      close(unit=7)

      do lg=1,96
      do mg=1,192
c The input data gradwe and gradsn is the topography gradient (m/m)
c Convert to flow velocity according to Miller et al(1994), J Climate, 914-928
        ugradwe=0.35*sqrt(abs(gradwe(mg,lg))/0.00005)
        ugradsn=0.35*sqrt(abs(gradsn(mg,lg))/0.00005)
c Limit these velocities to between 0.15m/sec and 5m/sec
        if(ugradwe.ne.0.0) ugradwe=max(0.15,min(5.0,ugradwe))
        if(ugradsn.ne.0.0) ugradsn=max(0.15,min(5.0,ugradsn))
c Scale these velocities to suit T63 model : now
c  fraction of river/reservoir removed per model timestep
        gradwe(mg,lg)=(mstep/15.0)*0.0256*sign(ugradwe,gradwe(mg,lg))
        gradsn(mg,lg)=(mstep/15.0)*0.0256*sign(ugradsn,gradsn(mg,lg))
      enddo
      enddo

C Set up a pointer to show direction of flow from trough point X
c
c     *************
c     * 4 * 3 * 2 *
c     *************
c     * 5 * X * 1 *  1=Easterly, 2=North-Easterly, 3=Northerly etc
c     *************
c     * 6 * 7 * 8 *
c     *************

      do kk=0,309
        ntrpnt(kk)=0 ! if "Basin"
        teste=ntre(kk).gt.0
        testn=ntrn(kk).gt.0
        testw=ntrw(kk).gt.0
        tests=ntrs(kk).gt.0
        if(teste.and.
     &  (.not.testn).and.(.not.testw).and.(.not.tests))
     &  ntrpnt(kk)=1
        if(teste.and.testn)ntrpnt(kk)=2
        if(testn.and.
     &  (.not.teste).and.(.not.testw).and.(.not.tests))
     &  ntrpnt(kk)=3
        if(testn.and.testw)ntrpnt(kk)=4
        if(testw.and.
     &  (.not.testn).and.(.not.teste).and.(.not.tests))
     &  ntrpnt(kk)=5
        if(testw.and.tests)ntrpnt(kk)=6
        if(tests.and.
     &  (.not.testn).and.(.not.testw).and.(.not.teste))
     &  ntrpnt(kk)=7
        if(tests.and.teste)ntrpnt(kk)=8
      enddo

C Set up the distance factor for troughs
      do kk=0,309
        mg=ntrmg(kk)
        lg=ntrlg(kk)
        lgn=lg+ntrn(kk)
        lgs=lg-ntrs(kk)
        mge=mg+ntre(kk)
        mgw=mg-ntrw(kk)
c       distf(kk)=min(2.0,sqrt(float(mge-mgw)**2+float(lgn-lgs)**2))
        distf(kk)=sqrt(float(mge-mgw)**2+float(lgn-lgs)**2)
c Make Niagara Falls and downstream flow slow
        if((kk.eq.166).or.(kk.eq.176))distf(kk)=10.0

        if(ntrpnt(kk).eq.1)topx=topminT63(mge,lg)  ! E flow
        if(ntrpnt(kk).eq.2)topx=topminT63(mge,lgn) ! NE flow
        if(ntrpnt(kk).eq.3)topx=topminT63(mg,lgn)  ! N flow
        if(ntrpnt(kk).eq.4)topx=topminT63(mgw,lgn) ! NW flow
        if(ntrpnt(kk).eq.5)topx=topminT63(mgw,lg)  ! W flow
        if(ntrpnt(kk).eq.6)topx=topminT63(mgw,lgs) ! SW flow
        if(ntrpnt(kk).eq.7)topx=topminT63(mg,lgs)  ! S flow
        if(ntrpnt(kk).eq.8)topx=topminT63(mge,lgs) ! SE flow

        topij=topminT63(mg,lg)
c topdiff is difference between lowest pixel height (topij) of
c  the trough (mg,lg) and the downstream grid point (topx) to which
c  the river is directed to flow. A few downstream points have
c  a greater height (caused by Basins etc, or special requirements).
c  Note : set the topdiff value to have a minimum of 10 m
c The flowrate downstream will be made (weakly) dependent upon the
c  topdiff value.
        topdiff=max(10.0,topij-topx)
        tdl10=log10(topdiff) ! Values from 1 to just over 3

c Reduce the distance factor for steep topography
c (i.e.increase flowrate)
        distf(kk)=distf(kk)/tdl10

      enddo

      ELSE

c---- get the runoff relocation data from file
      write(lfn,100)mw-1
  100 format('landrun',i2)
      write(6,1616)lfn
 1616 format('Reading land to ocean runoff file=',a20)
      open(unit=7,file=lfn,form='formatted',status='old',iostat=ierr)
       call errcheck(ierr,'landrun??        ','landrun   ')
       read(7,*)topland
       read(7,*)spread
       read(7,*)landlg
       read(7,*)landmg
       read(7,*)minmgx,minlgx
      close(unit=7)

      ENDIF

c----
c---- see ocforce.f for use of this data
c----

      return
      end
