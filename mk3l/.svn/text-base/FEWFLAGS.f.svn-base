c Converting output of atmosphere model from 16-bit INTEGER to 32-bit REAL.
c SJP 2009/04/24
c
c variable  nsib is taken out and replaced with a character(5)
c variable lsm_type
c AJA 2009/01/22
c
c Added the flag SUBICE, which controls whether or not the sub-ice heat input
c is used.
c SJP 2004/01/05
c
c $Log: FEWFLAGS.f,v $
c Revision 1.46  2001/10/12 02:13:46  rot032
c LDR changes, to bring sulfur cycle to run P76, and Mk2 OGCM fixes
c
c Revision 1.45  2001/02/22 05:34:45  rot032
c Changes from HBG to complete concatenation of NH/SH latitudes
c
c Revision 1.44  2001/02/12 05:39:57  rot032
c New river routing scheme and tidy-ups from HBG
c
c Revision 1.43  2000/11/14 06:55:45  rot032
c New PI_emissions flag.
c
c Revision 1.42  2000/11/14 03:11:39  rot032
c Energy/water balance for land-surface and sea-ice, plus tidy ups from HBG
c
c Revision 1.41  2000/08/16 02:59:26  rot032
c NCEP initialization option and CGCM changes from HBG
c
c Revision 1.40  2000/06/20 02:54:07  rot032
c Merge of HBG and LDR changes
c
c Revision 1.39  2000/06/20 02:08:35  rot032
c HBG changes to V5-3-3, mainly for coupled model
c
c Revision 1.38.1.1  2000/06/20 02:27:12  rot032
c Initial coupling of ECHAM chemical transport model (LDR)
c
c Revision 1.38  1999/05/20 06:23:58  rot032
c HBG changes to V5-2
c
c Revision 1.37  1998/12/10  00:56:01  ldr
c HBG changes to V5-1-21
c
c Revision 1.36  1997/06/13  06:03:53  ldr
c Changes from HBG to include UKMO convection (if ukconv=T).
c
c Revision 1.35  1997/06/11  02:21:33  ldr
c Include indirect sulfate aerosol effects.
c
c Revision 1.34  1996/12/23  03:58:02  mrd
c Add new gravity wave drag scheme as an option.
c
c Revision 1.33  1995/08/18  05:52:51  ldr
c HBG's changes to V4-7-12l for hybrid coordinate.
c
c Revision 1.32  1994/12/22  23:51:06  ldr
c Mods to cloud scheme of V4-6-13l to produce faster version for HBG and
c to allow vertically subgrid cloud below inversions in BL. Set slowice=T
c for accurate (but slow) ice fall calculation. Also, Ri now includes hlfusionn
c
c Revision 1.31  1994/12/07  23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.30  1994/08/08  17:20:44  ldr
c Strip off excessive RCS comments at top of file.
c
c Revision 1.29  94/03/30  10:17:12  ldr
c Add new debug diagnostic.
c 
c Revision 1.28  93/11/03  13:08:12  ldr
c Make locean a namelist flag and Replace silly flag zfl_sflg with savezflux.
c 
c Revision 1.27  93/09/30  12:40:26  ldr
c Changes to add option for automatic naming of output restart files.
c 
c Revision 1.26  93/09/15  17:09:04  ldr
c Changes to get inital version of cloud water scheme going.
c 
c Revision 1.25  93/08/19  15:07:19  ldr
c Minor cosmetic changes.
c 
c Revision 1.24  93/08/19  11:50:06  ldr
c Replace INCD with NDSTOP in namelist and add required housekeeping code.
c 
c Revision 1.23  93/08/17  11:57:50  ldr
c Add extra flag saveicuv to control dumping of SPO's U,V diagnostics.
c 
c Revision 1.22  93/08/06  11:42:32  ldr
c Introduce new flag savefcor to control saving of flux correction data.
c 
c Revision 1.21  93/08/03  10:12:43  ldr
c Rearrange common  blocks to avoid misaligned data on SGI.
c 
c Revision 1.20  93/07/14  14:51:02  ldr
c  ECMWF implicit treatment of vorticity eqn is now an option (HBG).
c 
      logical ifwds,iener,ispec,cdmap,mlomap,filewrflag
     &    ,ltrace,statsflag,rainflag,tempflag,uvtflag,qflux
     &    ,lcouple,locean,sdiagflag,saveqflux,savegbrf,saveglmean
     &    ,semice,leads,idyn,plotheat,plotclds,plotnetr,plotevrn
     &    ,clforflag,sltrace,chenflag,laust,impvor,savefcor
     &    ,saveicuv,qcloud,autoname,savezflux,debug,hybrid
     &    ,ukconv,kuoflag,ncarpbl,jsfdiff,amipo3,coupled_aero,fluxadj
     &    ,ncepagcm,PI_emissions,newriver,subice

      integer incd,nsstop,months,nsemilag,lastmonth,nthreads,idayp
     &    ,ksc,glmean_interval,ndstop,insdebug,lgdebug,mgdebug,ngwd
     &    ,naerosol_d,naerosol_i(2)

      character*5 lsm_type

      common/fewflags/ifwds,iener,ispec,cdmap,mlomap,filewrflag
     &    ,ltrace,statsflag,rainflag,tempflag,uvtflag,qflux
     &    ,lcouple,locean,sdiagflag,saveqflux,savegbrf,saveglmean
     &    ,semice,leads,idyn,plotheat,plotclds,plotnetr,plotevrn
     &    ,clforflag,sltrace,chenflag,laust,impvor,savefcor
     &    ,saveicuv,qcloud,autoname,savezflux,debug,hybrid
     &    ,ukconv,kuoflag,ncarpbl,jsfdiff,amipo3,coupled_aero,fluxadj
     &    ,ncepagcm,PI_emissions,newriver
     &    ,incd,nsstop,months,nsemilag,lastmonth,nthreads,idayp
     &    ,ksc,glmean_interval,ndstop,insdebug,lgdebug,mgdebug,ngwd
     &    ,naerosol_d,naerosol_i,subice,lsm_type
