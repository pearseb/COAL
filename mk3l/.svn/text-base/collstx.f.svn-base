c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c Modified for five-character experiment names.
c SJP 2009/04/22
c
c Modified for five-digit year numbers.
c SJP 2004/09/29
c
c $Log: collstx.f,v $
c Revision 1.15  2001/02/28 04:36:38  rot032
c Further tidy ups from HBG
c
c Revision 1.14  2001/02/22 05:34:41  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.13  1999/06/22 04:11:12  rot032
c Add geopotential height diagnostic and tidy datard.f.
c
c Revision 1.12  1999/05/20 06:23:52  rot032
c HBG changes to V5-2
c
c Revision 1.11  1996/10/24  01:02:34  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.10  1995/05/02  07:25:58  ldr
c Merge HBG V4-5-43h2 stuff into latest version.
c
c Revision 1.9  95/03/14  05:02:17  ldr
c Merge of INS changes for new ncput into LDR's latest version.
c 
c Revision 1.8  1995/03/14  04:58:52  ldr
c Changes for new version of ncput which handles reduced T63 arrays.
c
c Revision 1.7.1.1  95/05/02  07:14:37  ldr
c Other changes from HBG that went into coupled run but not into V4-5-43h.
c 
c Revision 1.7  93/12/17  15:31:56  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.6  93/10/05  13:05:34  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c 
c Revision 1.5  93/06/18  14:43:15  ldr
c Add u*3 diagnostic.
c 
c Revision 1.4  93/03/15  17:08:09  ldr
c Corrected NetCDF filenames.
c 
c Revision 1.3  93/02/05  16:45:28  ldr
c Introduce NetCDF option for stats files and rationalize filenames.
c 
c Revision 1.2  92/12/09  12:20:12  ldr
c Some changes to get coupled common blocks into include files.
c 
c Revision 1.1  92/07/08  12:34:02  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:   from common/masiv5 in this subroutine
c                  statsice - global statistics relating to sea-ice scheme
c
c              from arguments
c                  con - factoring constant
c                  str - character string YYYYYMM.XXXXX where XXXXX is runtype
c                  kvar- variable number
c
c     In/Out:  from common/timex in TIMEX.f
c                  iyear - year counter for model run
c                  month - month counter
c
c
      subroutine collstx(kvar,con,str)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer kvar
      real con
      character str*13

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'TIMEX.f'    !Input iyear, month

      real statsice
      common/masiv5/statsice(ln2,14,lat)

C Local work arrays and variables
      real tx(lon,lat2)
      character*13 filename

      integer lg
      integer mg

C Local data, functions etc
      character*3 varname(14)
      data varname /'xxx','xxx','xxx','xxx','xxx','xxx','xxx','xxx',
     & 'xxx','ico','icd','itf','isf','u@3' /

C Start code : ----------------------------------------------------------

c**** ROUTINE TO COLLECT ICE MODEL STATS FROM THE N1,S1,N2,S2...
c**** METHOD OF STORAGE IN /MASIV5/STATSICE INTO ARRAYS SUITABLE FOR
c**** PLOTTING,ETC : EG. (LON,LAT*2) WHERE 1,1 IS (SOUTH POLE,GM).
c----
c---- Note : Data is collected from arrays arranged (ln2,kvar,lat)
c----

c Work out filename:
c Character string STR*13 contains YYYYYMM.XXXXX where XXXXX is runtype
      filename='s'//varname(kvar)//'_'//str(9:13)//'.nc' !Harvey's suggestion

c Transfer data to plotting array : (1,1) is SP,GM
c and Write data
c---- write full field
         do 20 lg=1,lat
            do 20 mg=1,lon
            tx(mg,lg)=statsice(mg+lon,kvar,lg)*con
 20         tx(mg,lat2p-lg)=statsice(mg,kvar,lg)*con

      call nc_put(lon, lat2, filename, varname(kvar), iyear, month, tx)

      return
      end
