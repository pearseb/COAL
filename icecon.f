c Added the flag SUBICE, which controls whether or not the sub-ice heat input
c is used.
c SJP 2004/01/05
c
c HFACA removed from /FICECON/, as this array is never used.
c SJP 2003/06/20
c
c $Log: icecon.f,v $
c Revision 1.14  2001/10/12 02:13:44  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.13  2000/11/21 01:18:39  rot032
c Minor re-tuning from SPO.
c
c Revision 1.12  1999/06/16 06:21:53  rot032
c HBG changes to V5-3
c
c Revision 1.11  1999/05/20 06:23:47  rot032
c HBG changes to V5-2
c
c Revision 1.10  1998/12/10  00:55:27  ldr
c HBG changes to V5-1-21
c
c Revision 1.9  1997/12/23  00:23:33  ldr
c Merge of the HBG/SPO/EAK changes with other changes since V5-1.
c
c Revision 1.8  1997/12/03  05:35:42  ldr
c Put IGW's qflux code with variable MLO depth into main version of model.
c
c Revision 1.7.1.1  1997/12/19  02:03:08  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.7  1996/10/24  01:02:53  ldr
c Comments, Legendre transforms and tidy-ups from TIE
c
c Revision 1.6  1996/06/13  02:06:48  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.5  1994/04/29  15:05:59  ldr
c Correct flxia to make consistent with SPO's latest version.
c
c Revision 1.4  94/03/30  12:34:31  ldr
c Changes to V4-5 from HBG and SPO for coupled and qflux runs
c 
c Revision 1.3  93/12/17  15:32:48  ldr
c Hack V4-4-52l to change all continuation chars to &
c 
c Revision 1.2  93/06/16  14:17:28  ldr
c Reduced flxia from 4 to 2.
c 
c Revision 1.1  93/03/22  11:19:34  ldr
c Initial revision
c
c     INPUT/OUTPUT
c     Input:   from common/fewflags in FEWFLAGS.f
c                  semice - if T, use the Semtner sea-ice model
c
c     Output:  from common/ficecon in this subroutine
c                  flxia - subice heat input
c  
      subroutine icecon
      include 'PARAMS.f'
      common /ficecon/ flxia(lat,2)
      include 'FEWFLAGS.f'

c**** Set flxi values (Subice heat input)

      if (subice) then

        if (semice) then

          if ( qflux ) then
            do lg=1,lat
              flxia(lg,1)=0.0
              flxia(lg,2)=0.0
            enddo
          elseif (.not.qflux) then
c--    flxi is now set for NH and (larger) for SH
           if(lw.eq.22)then
             do lg=1,lat
               flxia(lg,1)=2.0
               flxia(lg,2)=15.0
             enddo
           else !T63
             do lg=1,lat
               flxia(lg,1)=4.0
               flxia(lg,2)=15.0
             enddo
           endif
          endif

        else

c---- If not SEMICE, then sub iceheat input is 2W/M**2 everywhere
           do 60 lg=1,lat
              flxia(lg,1)=2.0
 60           flxia(lg,2)=2.0

        end if

      else

        do lg = 1, lat
          flxia(lg, 1) = 0.0
          flxia(lg, 2) = 0.0 
        end do

      end if

      return
      end
