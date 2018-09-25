c $Log: aplota.f,v $
c Revision 1.14  2001/10/12 02:06:58  rot032
c HBG changes, chiefly for new river-routing scheme and conservation
c
c Revision 1.13  1996/06/13 02:05:33  ldr
c Tidy-ups from TIE using cflint to remove redundant code
c
c Revision 1.12  1995/08/14  05:29:54  ldr
c HBG's new improved printing routines.
c
c Revision 1.11  1993/12/17  15:31:40  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.10  93/10/19  11:22:24  ldr
c Add extra diagnostics for qcloud scheme.
c 
c Revision 1.9  93/08/19  15:07:42  ldr
c Minor cosmetic changes.
c 
c Revision 1.8  93/03/25  16:24:31  ldr
c Correction for R42 (HBG).
c 
c Revision 1.7  93/03/12  10:00:40  ldr
c HBG change to indicate which part of plot is NP.
c 
c Revision 1.6  93/01/08  17:16:08  ldr
c Tidy up Hal's character arrays so that it works on SGI.
c 
c Revision 1.5  92/12/09  14:42:46  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.4  92/12/07  15:24:59  ldr
c Changes by HBG to add new diagnostics and make general for R42.
c 
c Revision 1.3  92/11/20  12:11:19  ldr
c Changed for extra line plots. (hbg)
c 
c Revision 1.2  91/03/13  12:54:32  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:36:43  ldr
c Initial release V3-0
c 
      subroutine aplota(ind)

      implicit none
C Global parameters
      include 'PARAMS.f'

C Argument list
      integer ind

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'RADAV.f'

      real gev,grn,ghf,grh,ght,gcl,gcll,gclm,gclh,grg
     & ,grt,gsg,gsin,gsou,galb
     & ,gta,gts,gtsl,gtss,gsbl,ggw,grns,grnc,grtclr,gsouclr
     & ,gclc,gasbal,gatbal,gsfbal
     & ,gochfx,gflxic,gflxml,gflxse,gflxmj,grebal
     & ,gsbl_l,gsbl_i,gsbl_s,gcolwv
     & ,guz,gvz,gtz,gww
      common/maty/gev,grn
     & ,ghf,grh(nl),ght(nl),gcl,gcll,gclm,gclh,grg
     & ,grt,gsg,gsin,gsou,galb
     & ,gta,gts,gtsl,gtss,gsbl,ggw,grns,grnc,grtclr,gsouclr
     & ,gclc,gasbal,gatbal,gsfbal
     & ,gochfx,gflxic,gflxml,gflxse,gflxmj,grebal
     & ,gsbl_l,gsbl_i,gsbl_s,gcolwv
     & ,guz(nl),gvz(nl),gtz(nl),gww(nl)

C Local work arrays and variables
      integer latpr(lat2)

      character*130 linech
      character*4 ch4

      character*2 aplot(lat2,41)
      character*1 aplot1(2*lat2,41)
      equivalence (aplot1(1,1),aplot(1,1))

      integer i
      integer latx
      integer lg
      integer lga
      integer lgb
      integer ll
      integer loops
      integer lx
      integer ly
      integer mga
      integer nchrs

C Local data, functions etc
      character*4 cher(21)
      data cher/'10.0',' 9.5',' 9.0',' 8.5',' 8.0',' 7.5',' 7.0',
     &   ' 6.5',' 6.0',' 5.5',' 5.0',' 4.5',' 4.0',' 3.5',' 3.0',
     &   ' 2.5',' 2.0',' 1.5',' 1.0',' 0.5',' 0.0'/

C Start code : ----------------------------------------------------------

C**** ROUTINE TO PLOT ZONAL ATMOSPHERIC STATISTICS
C**** AS A PAGE PLOT (LINE PRINTER OUTPUT).

c****
c**** IF IND.EQ.1
c**** TO PAGE PLOT ASBAL,ATBAL,SFBAL (R,A,O) (+/- 200 W/M**2)
c**** IF IND.EQ.2
c**** TO PAGE PLOT LOW,MIDDLE,HIGH CL (L,M,H) (0 - 100 %)
c**** IF IND.EQ.3
c**** TO PAGE PLOT NETT TOP OF ATMOS RADIATION BALANCE
c****   NETT SOLAR=S, NETT OLR=R (0 - 400 W/M**2)
c**** IF IND.EQ.4
c**** TO PAGE PLOT EVAP/RAIN (E,R) (0 - 10 mms/day)
c****
      do 10 ly=1,41
         do 10 lx=1,lat2
 10      aplot(lx,ly)='  '
      if (ind.ne.3) then
         do 20 lx=1,lat2
 20         aplot(lx,21)='--'
      end if
      do 30 lx=1,lat
         latpr(lx)=lx-(lx/10)*10
 30      latpr(lat2+1-lx)=latpr(lx)
c---- decide on the number of parts to a plot - the latitudes are
c---- being shown from North Pole at left running across the page
c---- The maximum number of lats across the page with the current
c---- format is 58.
      loops=lat2/58+1
      latx=lat2/loops
      nchrs=4+3+latx*2+3+4
c---- put ' * ' in the character array for printing
      linech(5:5)=' '
      linech(6:6)='*'
      linech(7:7)=' '
      lga=8+latx*2
      lgb=lga+3
      linech(lga:lga+2)=' * '
      if (ind.eq.1) then
         do 40 lx=1,lat
            call aploti (asbalz(lx,1),'R',aplot,lx)
            call aploti (atbalz(lx,1),'A',aplot,lx)
 40         call aploti (sfbalz(lx,1),'O',aplot,lx)
         do 50 lx=lat,1,-1
            lg=lat2+1-lx
            call aploti (asbalz(lx,2),'R',aplot,lg)
            call aploti (atbalz(lx,2),'A',aplot,lg)
 50         call aploti (sfbalz(lx,2),'O',aplot,lg)
         do 100 ll=1,loops
            write (6,60) ll,loops
 60      format (//,1x,'plot of radiation balance(R), atmos heating(A),'
     &,' ocean balance(O) : part ',i1,' of ',i1)
            mga=(ll-1)*latx
            if(ll.eq.1)write (6,70) (latpr(lx),lx=1,latx)
            if(ll.gt.1)write (6,71) (latpr(lx),lx=mga+1,mga+latx)
 70         format (6x,'NP',58i2)
 71         format (8x,58i2)
            do 90 ly=41,1,-1
               i=10*ly-210
               write (ch4,'(i4)') i
               linech(1:4)=ch4(1:4)
               linech(lgb:lgb+3)=ch4(1:4)
               do 80 lx=1,latx*2
 80               linech(lx+7:lx+7)=aplot1(2*mga+lx,ly)
 90            write (6,110) (linech(lx:lx),lx=1,nchrs)
            if(ll.eq.1)write (6,70) (latpr(lx),lx=1,latx)
            if(ll.gt.1)write (6,71) (latpr(lx),lx=mga+1,mga+latx)
 100     continue
 110     format (1x,130a1)
         write (6,120) gasbal,gatbal,gsfbal
 120     format (1h ,'global mean heating rates(R,A,O):',3(1x,f6.1))
         write (6,130)
 130     format (1h ,'R=(solin-solout-rt) ; A=(solin-solout-sg)+(rg-rt)'
     &,'+hf+rn ; O=sg-rg-hf-ev-icf')
      end if
      if (ind.eq.2) then
c**** scale the clouds to be from 0-200 to match atmos/ocean balance
c**** plot above
         do 140 lx=1,lat
            call aploti (200.0*cllz(lx,1),'L',aplot,lx)
            call aploti (200.0*clmz(lx,1),'M',aplot,lx)
            call aploti (200.0*clhz(lx,1),'H',aplot,lx)
 140        call aploti (200.0*clz(lx,1)-200.0,'T',aplot,lx)
         do 150 lx=lat,1,-1
            lg=lat2+1-lx
            call aploti (200.0*cllz(lx,2),'L',aplot,lg)
            call aploti (200.0*clmz(lx,2),'M',aplot,lg)
            call aploti (200.0*clhz(lx,2),'H',aplot,lg)
 150        call aploti (200.0*clz(lx,2)-200.0,'T',aplot,lg)
         do 210 ll=1,loops
            write (6,160) ll,loops
 160     format (//,1x,'plot of cloud : low (L), middle (M), high (H) ',
     &'and total (T) below : part ',i1,' of ',i1)
            mga=(ll-1)*latx
            if(ll.eq.1)write (6,70) (latpr(lx),lx=1,latx)
            if(ll.gt.1)write (6,71) (latpr(lx),lx=mga+1,mga+latx)
            do 180 ly=41,21,-1
               i=5*ly-105
               write (ch4,'(i4)') i
               linech(1:4)=ch4(1:4)
               linech(lgb:lgb+3)=ch4(1:4)
               do 170 lx=1,latx*2
 170              linech(lx+7:lx+7)=aplot1(2*mga+lx,ly)
 180           write (6,110) (linech(lx:lx),lx=1,nchrs)
            do 200 ly=20,1,-1
               i=5*ly-5
               write (ch4,'(i4)') i
               linech(1:4)=ch4(1:4)
               linech(lgb:lgb+3)=ch4(1:4)
               do 190 lx=1,latx*2
 190              linech(lx+7:lx+7)=aplot1(2*mga+lx,ly)
 200           write (6,110) (linech(lx:lx),lx=1,nchrs)
            if(ll.eq.1)write (6,70) (latpr(lx),lx=1,latx)
            if(ll.gt.1)write (6,71) (latpr(lx),lx=mga+1,mga+latx)
 210     continue
      end if
      if (ind.eq.3) then
c**** Plot the top of atmosphere radiation balance
c**** Scale with -200
         do 220 lx=1,lat
            call aploti (sinz(lx,1)-souz(lx,1)-200.0,'S',aplot,lx)
 220        call aploti (rtz(lx,1)-200.0,'R',aplot,lx)
         do 230 lx=lat,1,-1
            lg=lat2+1-lx
            call aploti (sinz(lx,2)-souz(lx,2)-200.0,'S',aplot,lg)
 230        call aploti (rtz(lx,2)-200.0,'R',aplot,lg)
         do 270 ll=1,loops
            write (6,240) ll,loops
 240        format (//,1x,'plot of Top of Atmos Solar(in - out) (S), O',
     &'LR (R) : part ',i1,' of ',i1)
            mga=(ll-1)*latx
            if(ll.eq.1)write (6,70) (latpr(lx),lx=1,latx)
            if(ll.gt.1)write (6,71) (latpr(lx),lx=mga+1,mga+latx)
            do 260 ly=41,1,-1
               i=10*ly-10
               write (ch4,'(i4)') i
               linech(1:4)=ch4(1:4)
               linech(lgb:lgb+3)=ch4(1:4)
               do 250 lx=1,latx*2
 250              linech(lx+7:lx+7)=aplot1(2*mga+lx,ly)
 260           write (6,110) (linech(lx:lx),lx=1,nchrs)
            if(ll.eq.1)write (6,70) (latpr(lx),lx=1,latx)
            if(ll.gt.1)write (6,71) (latpr(lx),lx=mga+1,mga+latx)
 270     continue
      end if
      if (ind.eq.4) then
c**** Plot the evaporation(E) and rainfall(R) together
c**** Scale to be 0 to 200 (for 0-10 mms/day)
c**** Plot the supersaturation rainfall(S) and convection rainfall
c**** (C) below. Scale to be -200 to 0 (for 0-10mms/day)
         do 280 lx=1,lat
            call aploti (evz(lx,1)*20.0,'E',aplot,lx)
            call aploti (rnz(lx,1)*20.0,'R',aplot,lx)
            call aploti (rnsz(lx,1)*20.0-200.0,'S',aplot,lx)
 280        call aploti (rncz(lx,1)*20.0-200.0,'C',aplot,lx)
         do 290 lx=lat,1,-1
            lg=lat2+1-lx
            call aploti (evz(lx,2)*20.0,'E',aplot,lg)
            call aploti (rnz(lx,2)*20.0,'R',aplot,lg)
            call aploti (rnsz(lx,2)*20.0-200.0,'S',aplot,lg)
 290        call aploti (rncz(lx,2)*20.0-200.0,'C',aplot,lg)
         do 350 ll=1,loops
            write (6,300) ll,loops
 300        format (//,1x,'plot of Evap(E) and Rain(R), with rain split'
     &,' into saturation(S) and convective (C) below : part ',i1,' of ',
     &i1)
            mga=(ll-1)*latx
            if(ll.eq.1)write (6,70) (latpr(lx),lx=1,latx)
            if(ll.gt.1)write (6,71) (latpr(lx),lx=mga+1,mga+latx)
            do 320 ly=41,21,-1
               ch4=cher(42-ly)
               linech(1:4)=ch4(1:4)
               linech(lgb:lgb+3)=ch4(1:4)
               do 310 lx=1,latx*2
 310              linech(lx+7:lx+7)=aplot1(2*mga+lx,ly)
 320           write (6,110) (linech(lx:lx),lx=1,nchrs)
            do 340 ly=20,1,-1
               ch4=cher(22-ly)
               linech(1:4)=ch4(1:4)
               linech(lgb:lgb+3)=ch4(1:4)
               do 330 lx=1,latx*2
 330              linech(lx+7:lx+7)=aplot1(2*mga+lx,ly)
 340           write (6,110) (linech(lx:lx),lx=1,nchrs)
            if(ll.eq.1)write (6,70) (latpr(lx),lx=1,latx)
            if(ll.gt.1)write (6,71) (latpr(lx),lx=mga+1,mga+latx)
 350     continue
      end if
      return
      end
