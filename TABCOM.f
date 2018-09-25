c $Log: TABCOM.f,v $
c Revision 1.5  1994/12/07 23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.4  1993/12/17  15:31:33  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.3  93/08/19  15:07:37  ldr
c Minor cosmetic changes.
c 
c Revision 1.2  93/08/03  10:12:41  ldr
c Rearrange common  blocks to avoid misaligned data on SGI.
c 
c Revision 1.1  92/04/15  11:13:39  mrd
c Initial revision
c 
c     common block tabcom contains quantities precomputed in subroutine 
c     table for use in the longwave radiation program:  

      real em1(28,180)    ! e1 function, evaluated over the 0-560 and 
                          ! 1200-2200 cm-1 intervals
      real em1wde(28,180) ! e1 function, evaluated over the 160-560 cm-1
                          ! interval
      real table1(28,180) ! e2 function, evaluated over the 0-560 and 
                          !1200-2200 cm-1 intervals
      real table2 (28,180)! Temperature derivative of table1
      real table3(28,180) ! Mass derivative of table1 
      real em3(28,180)    ! e3 function, evaluated over the 0-560 and 
                          ! 1200-2200 cm-1 intervals
      real source(28,nbly)! Planck function, evaluated at specified temps. for
                          ! bands used in CTS calculations
      real dsrce(28,nbly) ! Temperature derivative of source
      integer ind(imax)   ! Index, with value ind(i)=i. Used in fst88 
      integer indx2(lp1v) ! Index values used in obtaining "lower triangle" 
                          ! elements of avephi,etc.,in fst88
      integer kmaxv(lp1)  ! Index values used in obtaining "upper triangle" 
                          ! elements of avephi,etc.,in fst88
      integer  kmaxvm     ! kmaxv(l), used for do loop indices 
 
      common /tabcom/ em1, em1wde, table1, table2, table3, em3,
     &                source, dsrce, ind, indx2, kmaxv, kmaxvm
