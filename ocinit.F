c Enhancing user control over the coupling between the atmosphere and ocean.
c SJP 2009/04/21
c
c Minor fixes to resolve warnings issued by the g95 Fortran compiler.
c SJP 2009/04/14
c
c Modified for coupling of the new, high-resolution version of the ocean model
c to the atmosphere.
c SJP 2007/12/20
c
c Major tidy-up of the ocean model source code.
c SJP 2007/06/16
c
c Replaced the line "include 'OCEANPARS.f'" with "include 'OPARAMS.f'",
c enabling the header file OCEANPARS.f to be removed from the model source
c code.
c SJP 2007/05/31
c
c Added IMPLICIT NONE statement, plus variable declarations as necessary.
c SJP 2007/05/30
c
c Transferred COMMON block declarations to separate header files.
c SJP 2007/05/29
c
c Modified to apply Delta-S corrections to the OGCM SSSs, as well as applying
c Delta-T corrections to the OGCM SSTs.
c SJP 2004/09/10
c
c$Log: ocinit.f,v $
cRevision 1.20  2001/10/12 02:13:45  rot032
cLDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
cRevision 1.19  2000/06/21 05:38:04  rot032
cFixips from HBG for NEC.
c
cRevision 1.17  1999/05/20 06:23:54  rot032
cHBG changes to V5-2
c
c Revision 1.16  1998/12/10  00:55:45  ldr
c HBG changes to V5-1-21
c
c Revision 1.15  1994/08/08  17:21:53  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.14  93/12/23  15:31:26  ldr
c Changes to V4-4-54l from HBG for coupled model of V4-5
c 
c Revision 1.13  93/12/20  16:22:06  ldr
c Minor changes to V4-4-53l from HBG for coupled model
c 
c Revision 1.12  93/12/17  11:51:40  ldr
c Changes to V4-4-45l from HBG for coupled model
c 
c Revision 1.11  93/11/29  11:38:34  ldr
c Changes to V4-4-32l from HBG for coupled model
c 
c Revision 1.10  93/11/03  11:44:23  ldr
c Changes to V4-4-24l from HBG for coupled model
c 
c Revision 1.9  93/08/10  15:27:41  ldr
c Changes made by HBG to V4-3-28l to rationalize control of coupled model.
c 
c Revision 1.8  93/07/12  14:14:24  ldr
c Minor changes from SPO for coupled model V4-4.
c 
c Revision 1.7  93/02/05  16:56:05  ldr
c Put /bcogcm in include file (and correct its subscripts!)
c 
c Revision 1.6  93/02/03  12:44:50  ldr
c SPO's minor changes to coupled model and sea-ice scheme of V4-1.
c 
      subroutine ocinit(nato,icstp)

c Set up SST's and zero arrays holding fluxes for coupled runs

      implicit none

C Global parameters
      include 'PARAMS.f'
      include 'OPARAMS.f'

C Argument list
      integer nato
      integer icstp

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'BCOGCM.f'
      include 'A2O.f'
      include 'O2A.f'
      include 'FEWFLAGS.f'
      include 'MASIV4.f'
      include 'TIMEX.f'
      include 'AMM.f'
      include 'TM.f'
      include 'ONEDIM.f'
      include 'SFC1.f'
      include 'AOGCM.f'

C Local work arrays and variables
      integer i, io, j, jo, jon, jos, nv
      real tadjun, tadjus

#ifdef OCEAN_LOW

      real osst2(ln2, lat)
      integer is
      real grad, sum

#endif

C Local data, functions etc

C Start code : ----------------------------------------------------------

      if(isync.eq.0)then
        ammratlm=0.0
      else
        ammratlm=float(icstp+1)/float(nato)
      endif

c Set up SST's : Note that the Atmos model requires that they are 
c NEGATIVE. Now 'osst' is already negative (see step.f).
c The Delta-T temp correction is derived from (Levitus - Ocean model) sst
c (i.e still postive). Thus subtract the Delta-T correction.

#ifdef OCEAN_LOW

C**** INTERPOLATE SST TO LAND POINTS (L) WHICH THE AGCM THINKS ARE SEA
C**** POINTS.
C
C    2345678901234567890123456789012345678901234567890123456789012345
C57           .         .         .         .         .         .
C56  ................................................................
C55           .         .         .         .         .   ********
C54    S      .         .         .         .         . **********
C53           .      *****        .         .         . L  *******
C52           .  ***********  **  .         .  * ********   *****
C51  ....****.************************..***************.**..****.....
C50     *******************************L**************. **  *** . S
C49    ** ***************************   ************* . LL   *  .
C48    ** **********************  *      *  .******** . **      .
C47      *********************   *.         . ********. ***     .   *
C46  **************************................************.........*
C45  *************************    .         .  *********** *    .
C44  ****** **.***************    .         .  ***********      .
C43  *   *   *.**************     .         .  **********       .   *
C42      *****.*************  *   .         .  *********        .   *
C41  .*....****************..*..................********............L
C40  ***    ***************L*     .         .   *******         .   *
C39  **********************       .         .    ******         .   *
C38  *********. ***********       .         .     **  .         .  **
C37  ***********  ********        .         .     **  .         .  **
C36  ***********..***.***..........................*...............**
C35  **********   **  ***         .         .       **.         .  **
C34  *********.    *   **         .         .        **         .  **
C33  **********    *   L*         .         .         * **      .  **
C32  *********.        L.  S      .         .         .****     .  **
C31  ..*******.........L................................****.........
C30    ****** .        *LL*       .         .         .******   .
C29    ****** .        **LSS  SS  .         .         .*******  .
C28    L***** .         S     *S  .         .         .******** .
C27     ***** .         .     **  .         .         .******** .
C26  ...*****..................L........................*******......
C25     ****  .         .   ** *  .         .         . ******* .
C24     ****L*.         .  *****  .         .         . ******  .
C23     ****L*.         . ******  .         .         . ******  .
C22     ****L*.         . ******* .         .         . *****   .
C21  ...***...............*******.......................*****........
C20     ***   .         . ******* .         .         . ****    .
C19      *    .         . *  ***  .         .         . ****    .
C18           .         .     **  .         .         . ***     .
C17           .         .      L  . *       .         . **      .
C16  ..........................*...*....................**...........
C15           .         .         .         .         . **      .
C14           .         .         .         .         . *       .
C13           .         .         .         .         . S       .
C12           .         .         .         .         .         .
C11  ................................................................
C10           .         .         .         .         .         .
C 9           .         .         .         .         .   S     .
C 8          *****  ************* .         .         .  *      .
C 7  *******************************        .         .  *      .   *
C 6  *****************************.............**********************
C 5  ****************************************************************
C 4  ****************************************************************
C 3  ****************************************************************
C 2  ****************************************************************
C    2345678901234567890123456789012345678901234567890123456789012345
C
C**** There are 18 points that the atmospheric model requires SST data
C**** which the ocean model has as land points. Generate pseudo SSTs.
C****  TWO TIME LEVELS IS=1,2.
C
      DO 9150 IS=1,2

      OSST(28,17,IS)=0.5*(OSST(27,17,IS)+OSST(29,17,IS))

      GRAD=(OSST(9,25,IS)-OSST(9,21,IS))/(DYU(21)+DYU(22)+DYU(23)+
     &      DYU(24))
      SUM=0.0
      DO 9002 J=22,24
      SUM=SUM+DYU(J-1)
 9002 OSST(9,J,IS)=OSST(9,21,IS)+SUM*GRAD

      OSST(28,26,IS)=0.5*(OSST(27,26,IS)+OSST(29,26,IS))

      OSST( 4,28,IS)=OSST( 3,28,IS)

      GRAD=(OSST(22,28,IS)-OSST(22,31,IS))/(DYU(28)+DYU(29)+DYU(30))
      SUM=0.0
      DO 9004 J=29,30
      SUM=SUM+DYU(J-1)
 9004 OSST(22,J,IS)=OSST(22,28,IS)+SUM*GRAD

      OSST(21,30,IS)=OSST(21,31,IS)

      DO 9006 J=31,32
 9006 OSST(20,J,IS)=0.5*(OSST(19,J,IS)+OSST(21,J,IS))

      OSST(20,33,IS)=OSST(19,33,IS)

      OSST(24,40,IS)=OSST(24,39,IS)+DYU(39)*
     &       (OSST(24,41,IS)-OSST(24,39,IS))/
     &    (DYU(39)+DYU(40))

      OSST(65,41,IS)=0.5*(OSST(64,41,IS)+OSST(2,41,IS))

      OSST(53,49,IS)=(2.0*OSST(52,49,IS)+OSST(55,49,IS))/3.0
      OSST(54,49,IS)=(OSST(52,49,IS)+2.0*OSST(55,49,IS))/3.0

      OSST(36,50,IS)=OSST(36,49,IS)+DYU(49)*
     &          (OSST(36,51,IS)-OSST(36,49,IS))/
     &    (DYU(49)+DYU(50))

      OSST(53,53,IS)=0.5*(OSST(52,53,IS)+OSST(54,53,IS))

 9150     CONTINUE

c...  Loop over the two time levels
      do is = 1, 2

c... Copy SSTs into temporary array OSST2
        do j = 2, jmt/2
          do i = 2, imt-1
            osst2(i-1, j-1) = osst(i, jmt+1-j, is)
            osst2(i-1+lon, j-1) = osst(i, j, is)
          end do
        end do

c...  Ensure conservation of global-mean SST
        call conserve(osst2, osst2, 3)

c...  Copy corrected SSTs back into OSST
        do j = 2, jmt/2
          do i = 2, imt-1
            osst(i, jmt+1-j, is) = osst2(i-1, j-1)
            osst(i, j, is) = osst2(i-1+lon, j-1)
          end do
        end do

      end do

      do nv = 1, 2
        do j = 1, lat
          do i = 1, lon
            savegrid(i, nv+4, j) = osst(i+1, jmt-j, nv)
            savegrid(i+lon, nv+4, j) = osst(i+1, j+1, nv)
          end do
        end do
      end do

#else

c...  When using the high-resolution version of the ocean model, the SST
c...  for each AGCM gridbox is obtained by calculating the area-weighted
c...  average of the SST for each of the four corresponding OGCM gridboxes.

      do nv = 1, 2
        do j = 1, lat
          jon = jmt - 2 * j
          jos = 2 * j
          do i = 1, lon
            io = 2 * i
            savegrid(i, nv+4, j) =
     &        (cst(jon) * osst(io, jon, nv) +
     &         cst(jon) * osst(io+1, jon, nv) +
     &         cst(jon+1) * osst(io, jon+1, nv) +
     &         cst(jon+1) * osst(io+1, jon+1, nv)) /
     &        (2.0 * (cst(jon) + cst(jon+1)))
            savegrid(i+lon, nv+4, j) =
     &        (cst(jos) * osst(io, jos, nv) +
     &         cst(jos) * osst(io+1, jos, nv) +
     &         cst(jos+1) * osst(io, jos+1, nv) +
     &         cst(jos+1) * osst(io+1, jos+1, nv)) /
     &        (2.0 * (cst(jos) + cst(jos+1)))
          end do
        end do
      end do

#endif

      if (fluxadj) then
        do nv = 1, 2
          do j = 1, lat
            do i = 1, lon
              tadjun = tmo(i, j, 1, 1) * (1.0 - ratlm) +
     &                 tmo(i, j, 1, 2) * ratlm
              tadjus = tmo(i, j, 2, 1) * (1.0 - ratlm) +
     &                 tmo(i, j, 2, 2) * ratlm
              savegrid(i, nv+4, j) = savegrid(i, nv+4, j) - tadjun
              savegrid(i+lon, nv+4, j) = savegrid(i+lon, nv+4, j) -
     &                                   tadjus
            end do
          end do
        end do
      end if

c...  Set up the SSSs for the atmosphere model. Note that, as the SSSs are
c...  POSITIVE, the Delta-S corrections are ADDED.

#ifdef OCEAN_LOW

      do nv = 1, 2
        do j = 1, lat2
          do i = 1, lon
            asal(i, j, nv) = osal(i+1, j+1, nv)
          end do
        end do
      end do

#else

c...  When using the high-resolution version of the ocean model, the SSS
c...  for each AGCM gridbox is obtained by calculating the area-weighted
c...  average of the SSS for each of the four corresponding OGCM gridboxes.

      do nv = 1, 2
        do j = 1, lat2
          jo = 2 * j
          do i = 1, lon
            io = 2 * i
            asal(i, j, nv) = (cst(jo) * osal(io, jo, nv) +
     &                        cst(jo) * osal(io+1, jo, nv) +
     &                        cst(jo+1) * osal(io, jo+1, nv) +
     &                        cst(jo+1) * osal(io+1, jo+1, nv)) /
     &                       (2.0 * (cst(jo) + cst(jo+1)))
          end do
        end do
      end do

#endif

      if (fluxadj) then
        do nv = 1, 2
          do j = 1, lat
            do i = 1, lon
              tadjun = smo(i, j, 1, 1) * (1.0 - ratlm) +
     &                 smo(i, j, 1, 2) * ratlm
              tadjus = smo(i, j, 2, 1) * (1.0 - ratlm) +
     &                 smo(i, j, 2, 2) * ratlm
              asal(i, lat2+1-j, nv) = asal(i, lat2+1-j, nv) + tadjun
              asal(i, j, nv) = asal(i, j, nv) + tadjus
            end do
          end do
        end do
      end if

c Zero arrays holding fluxes to pass back to OGCM

      do j = 1, jmt-2
        do i = 1, imt-2
          otx(i, j) = 0.0
          oty(i, j) = 0.0
          osalf(i, j) = 0.0
          osurh(i, j) = 0.0
c rjm
	   oswnd(i,j) = 0.0
	   osrad(i,j) = 0.0
	   osice(i,j) = 0.0
c rjm
        end do
      end do

      return
      end
