c $Log: prtcl.f,v $
c Revision 1.15  2001/02/28 04:36:39  rot032
c Further tidy ups from HBG
c
c Revision 1.14  2001/02/22 05:34:43  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.13  1999/06/16 06:21:54  rot032
c HBG changes to V5-3
c
c Revision 1.12  1999/05/20 06:23:55  rot032
c HBG changes to V5-2
c
c Revision 1.11  1998/12/10  00:55:51  ldr
c HBG changes to V5-1-21
c
c Revision 1.10  1997/12/19  02:03:17  ldr
c Changes from HBG, SPO and EAK for T63 CGCM, ice fixes and soil.
c
c Revision 1.9  1995/08/14  05:29:54  ldr
c HBG's new improved printing routines.
c
c Revision 1.8  1993/12/17  15:33:27  ldr
c Hack V4-4-52l to change all continuation chars to &
c
c Revision 1.7  93/07/06  16:08:54  ldr
c      New map for ice depth averaged over leads area. This requires extra array
c      in common block cloudm1. The mapping in surfset.f reqires pl(mg).
c 
c Revision 1.6  93/03/12  10:04:17  ldr
c Minor changes (HBG).
c 
c Revision 1.5  93/02/03  12:44:53  ldr
c SPO's minor changes to coupled model and sea-ice scheme of V4-1.
c 
c Revision 1.4  92/12/09  14:44:07  ldr
c Replaced all n's with nl's for compatibility with ogcm.
c 
c Revision 1.3  92/04/22  11:52:19  ldr
c Generalized for R21/R42.
c 
c Revision 1.2  91/03/13  12:59:36  ldr
c Common blocks put into include file form.
c Any changes made since V3-0 and this V3-1 merged.
c 
c Revision 1.1  91/02/22  16:37:50  ldr
c Initial release V3-0
c 
      subroutine prtcl

C Global parameters
      include 'PARAMS.f'

C Argument list

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks
      include 'CHMAP.f'
c     include 'FEWFLAGS.f' ! leads
      include 'PRINTT.f'
      character*1 cnmlp(ln2,lat),rainp(ln2,lat)
     & ,snmap(ln2,lat),snwice(ln2,lat),snicel(ln2,lat)
      common/cloud1/cnmlp,rainp,snmap,snwice,snicel
      common /negqb/ qbneg(lat,2,nl),qmean(lat,2,nl)

C Local work arrays and variables
      integer i1num(lon)
      character*5 clch(nl)

C Local data, functions etc

C Start code : ----------------------------------------------------------

c---- map print of istantaneous model indicators

      if (.not.gmap1) return
      do 10 mg=1,lon
 10      i1num(mg)=mg-(mg/10)*10
c---- maximum number of integers/characters per line is 128
      maxnum=128
      loops=(lon-1)/maxnum+1
      maxnum=lon/loops

c---- optional print of rain and convection points
      if (.not.cvrnm) go to 60
      do 30 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         if(lx.eq.loops)mgb=lon
         write (6,20) mga,mgb
 20      format (1x,' mid,low,pen conv points: mid=m,low=l,pen=p : mg='
     &    ,i3,' to ',i3)
 30      call prtclz (cnmlp,mga,mgb,i1num)
      do 50 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         if(lx.eq.loops)mgb=lon
         write (6,40) mga,mgb
 40      format (1x,' rain 1-2=l, 3-5=m, 6+=h  and combs :l+m(4),l+h(5)'
     &,',m+h(6),l+m+h(7) : mg=',i3,' to ',i3)
 50      call prtclz (rainp,mga,mgb,i1num)

c---- optional print of snow and ice points, and 2 ground wetness
 60   if (.not.gwicm) go to 150
      do 80 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         if(lx.eq.loops)mgb=lon
         write (6,70) mga,mgb
 70      format (10x,'  snowd (in 2**n cms) : mg=',i3,' to ',i3)
 80      call prtclz (snmap,mga,mgb,i1num)
      do 100 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         if(lx.eq.loops)mgb=lon
         write (6,90) mga,mgb
 90      format (10x,'  snow(*),ice depth (feet) : mg=',i3,' to ',i3)
 100     call prtclz (snwice,mga,mgb,i1num)

c     if(leads)then
c     do 102 lx=1,loops
c        mgb=lx*maxnum
c        mga=mgb-maxnum+1
c        if(lx.eq.loops)mgb=lon
c        write (6,92) mga,mgb
c92      format (10x,'  snow(*),ice depth * leads area (feet) : mg='
c    & ,i3,' to ',i3)
c102     call prtclz (snicel,mga,mgb,i1num)
c     end if

c     do 120 lx=1,loops
c        mgb=lx*maxnum
c        mga=mgb-maxnum+1
c        if(lx.eq.loops)mgb=lon
c        write (6,110) mga,mgb
c110     format (10x,'wb ground wetness(0-8) available for evapotr.:
c    &    mg=',i3,' to ',i3)
c120     call prtclz (chwb,mga,mgb,i1num)
      do 140 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         if(lx.eq.loops)mgb=lon
         write (6,130) mga,mgb
 130     format (10x,'wg ground wetness(0-9,t) available for evapotr.:
     &    mg=',i3,' to ',i3)
 140     call prtclz (chwg,mga,mgb,i1num)

c---- optional print of relative humidity stats
 150  if (.not.rhnmm) go to 270
c---- max number of model levels is 18 in the next print outs
c---- skip levels if nl>18
      ksk=1+(nl-1)/18
      do 152 k=1,nl
 152     write (clch(k),"(' Rm',i2.2)") k

c     write (6,190)
c190  format (10x,'instantaneous mean negative P*.q before Physics')
c     write (6,235) (clch(k),k=ksk,nl,ksk)
c     do 200 lg=1,lat
c200     write (6,210) lg,(qbneg(lg,1,k),k=1,nl,ksk)
c210  format (1x,i2,2x,18f7.4)
c     do 220 lg=lat,1,-1
c220     write (6,210) lg,(qbneg(lg,2,k),k=1,nl,ksk)
c     write (6,235) (clch(k),k=ksk,nl,ksk)

      write (6,230)
 230  format (10x,'instantaneous mean  P*.q after Physics')
      write (6,235) (clch(k),k=ksk,nl,ksk)
 235  format(5x,18(2x,a5))
      do 240 lg=1,lat
 240     write (6,250) lg,(qmean(lg,1,k),k=1,nl,ksk)
 250  format (1x,i2,2x,18f7.3)
      do 260 lg=lat,1,-1
 260     write (6,250) lg,(qmean(lg,2,k),k=1,nl,ksk)
      write (6,235) (clch(k),k=ksk,nl,ksk)

c---- optional print of cloud maps
 270  if (.not.cldm) return
      do 290 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         if(lx.eq.loops)mgb=lon
         write (6,280) mga,mgb
 280     format (1x,'     low level cloud : mg=',i3,' to ',i3)
 290     call prtclz (chlowc,mga,mgb,i1num)
      do 310 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         if(lx.eq.loops)mgb=lon
         write (6,300) mga,mgb
 300     format (1x,'     mid level cloud : mg=',i3,' to ',i3)
 310     call prtclz (chmidc,mga,mgb,i1num)
      do 330 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         if(lx.eq.loops)mgb=lon
         write (6,320) mga,mgb
 320     format (1x,'      high level cloud : mg=',i3,' to ',i3)
 330     call prtclz (chhic,mga,mgb,i1num)
      do 350 lx=1,loops
         mgb=lx*maxnum
         mga=mgb-maxnum+1
         if(lx.eq.loops)mgb=lon
         write (6,340) mga,mgb
 340     format (1x,'      total cloud : mg=',i3,' to ',i3)
 350     call prtclz (chtotc,mga,mgb,i1num)
      return
      end
      subroutine prtclz (charr,mga,mgb,i1num)

C Global parameters
      include 'PARAMS.f'

C Argument list
      character*1 charr(ln2,lat)
      integer mga,mgb
      integer i1num(lon)

C Local shared common blocks (see TASK COMMON/TASKLOCAL above)

C Global data blocks

C Local work arrays and variables

C Local data, functions etc

C Start code : ----------------------------------------------------------

c---- maximum number of integers/characters per line is 128
c     (maxnum=128)
      write (6,10) (i1num(mg),mg=mga,mgb)
 10   format (3x,128i1)
      do 20 lg=1,lat
         l=lg-(lg/10)*10
 20      write (6,40) l,(charr(mg,lg),mg=mga,mgb)
      do 30 lg=lat,1,-1
         l=lg-(lg/10)*10
 30      write (6,40) l,(charr(mg+lon,lg),mg=mga,mgb)
 40   format (1x,i1,1x,128a1)

      return
      end
