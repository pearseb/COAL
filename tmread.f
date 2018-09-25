c Modified for the conversion of the coupled model auxiliary files to netCDF.
c In the process, these files have been renamed as follows:
c   dtm1av      -> dtm.nc
c   hfcor.dat12 -> hfcor.nc
c   s1coravth   -> ssscor.nc
c   sfcor.dat12 -> sfcor.nc
c   t1coravth   -> sstcor.nc
c   txcor.dat12 -> txcor.nc
c   tycor.dat12 -> tycor.nc
c SJP 2008/02/05
c
c Transferred COMMON blocks to separate header files, as follows:
c /TM/  ->  TM.f
c SJP 2007/05/29
c
c Re-inserted read from s1coravth, as we are now using this data after all!
c SJP 2004/09/10
c
c Commented out read from s1coravth, as this data is not used.
c SJP 2002/02/15
c
c $Log: tmread.f,v $
c Revision 1.11  2001/02/22 05:34:37  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.10  1999/05/20 06:23:49  rot032
c HBG changes to V5-2
c
c Revision 1.9  1998/05/27  02:31:37  ldr
c Tidy up last change.
c
c Revision 1.8  1998/05/27  02:10:17  ldr
c Merge TIE and ACH changes.
c
c Revision 1.7  1998/05/26  04:48:55  ldr
c Final changes (inc. ACH salinity stuff) for 10k run.
c
c Revision 1.6.1.1  1998/05/27  02:07:36  ldr
c Tidy up "file not found" error messages (from TIE).
c
c Revision 1.6  1995/05/02  07:49:15  ldr
c Merge of HBG stuff in V4-5-43h3 into latest version.
c
c Revision 1.5  1995/05/02  07:31:32  ldr
c Extra changes from HBG to files which already went in as V4-5-43h.
c
c Revision 1.4  1994/08/04  16:56:45  ldr
c Changes to V4-5-30l from HBG (clouds,coupled), IGW (tracers) and SPO (coupled)
c
c Revision 1.3  93/12/17  15:34:13  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.2  93/11/29  11:38:43  ldr
c Changes to V4-4-32l from HBG for coupled model
c 
c Revision 1.1  93/07/12  14:15:29  ldr
c Initial revision
c  
      subroutine tmread(month) 

c subroutine to read in temperature/salinity correction terms for
c coupled model.
c The correction terms are at the 1st of each month.

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer month

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'TM.f'
      real tsst1,tsst2,tdtm1,tdtm2
      common/compsst/tsst1(ln2,lat),tsst2(ln2,lat)
     &              ,tdtm1(ln2,lat),tdtm2(ln2,lat)

C Local work arrays and variables
      character title*60
      integer i, j, lg, mg, monthp1, ns

C Local data, functions etc
      logical got_data
      real dtm_save(lon, lat2, 12), ssscor_save(lon, lat2, 12),
     &     sstcor_save(lon, lat2, 12)

      data got_data /.false./
      save got_data, dtm_save, ssscor_save, sstcor_save

C Start code : ----------------------------------------------------------

c...  ***********************************
c...  ***  READ FLUX ADJUSTMENT DATA  ***
c...  ***********************************

c...  Read the flux adjustment data, if this has not already been done
      if (.not. got_data) then

c...  Sea surface temperature adjustments
        write (*, *)
        write (*, *) "Reading sea surface temperature adjustment data :"
        write (*, *)
        write (*, *) "File     =  sstcor.nc"
        write (*, *) "Variable =  sstcor"

        call read_real_3d("sstcor.nc", "sstcor", lon, lat2, 12,
     &                    sstcor_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

c...  Sea surface salinity adjustments
        write (*, *)
        write (*, *) "Reading sea surface salinity adjustment data :"
        write (*, *)
        write (*, *) "File     =  ssscor.nc"
        write (*, *) "Variable =  ssscor"

        call read_real_3d("ssscor.nc", "ssscor", lon, lat2, 12,
     &                    ssscor_save, title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

c...  Mixed-layer ocean temperature anomalies
        write (*, *)
        write (*, *)
     &    "Reading mixed-layer ocean temperature anomaly data :"
        write (*, *)
        write (*, *) "File     =  dtm.nc"
        write (*, *) "Variable =  dtm"

        call read_real_3d("dtm.nc", "dtm", lon, lat2, 12, dtm_save,
     &                    title)

        write (*, *) "Title    =  ", trim(title)
        write (*, *)

        got_data = .true.

      end if

c...  Transfer the data to the arrays TMO, SMO, TDTM1 and TDTM2
      monthp1 = month + 1
      if (monthp1 .eq. 13) monthp1 = 1
      do ns = 1, 2
        do lg = 1, lat
          if (ns .eq. 1) then
            j = lat2 + 1 - lg
          else
            j = lg
          end if
          do i = 1, lon
            if (ns .eq. 1) then
              mg = i
            else
              mg = i + lon
            end if
            tmo(i, lg, ns, 1) = sstcor_save(i, j, month)
            tmo(i, lg, ns, 2) = sstcor_save(i, j, monthp1)
            smo(i, lg, ns, 1) = ssscor_save(i, j, month)
            smo(i, lg, ns, 2) = ssscor_save(i, j, monthp1)
            tdtm1(mg, lg) = dtm_save(i, j, month)
            tdtm2(mg, lg) = dtm_save(i, j, monthp1)
          end do
        end do
      end do

      return
      end
