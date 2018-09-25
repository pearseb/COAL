c $Log: ECPARM.f,v $
c Revision 1.7  2001/06/04 02:26:52  rot032
c Changes from LDR to bring sulfur cycle to run P41
c
c Revision 1.6  2001/03/07 04:28:56  rot032
c Changes from LDR to bring sulfur cycle to run N63
c
c Revision 1.5  2000/12/11 04:12:06  rot032
c Set ne back to 0 for default version.
c
c Revision 1.4  2000/12/08 03:58:51  rot032
c Changes to aerosol scheme (and small newrain fix) from LDR
c
c Revision 1.3  2000/06/20 07:50:13  rot032
c Change order of type and parameter declarations.
c
c Revision 1.2  2000/06/20 04:01:49  rot032
c Make NE=0 the default.
c
c Revision 1.1  2000/06/20 02:22:22  rot032
c Initial revision
c
 
      integer ne,late,lone,lat2e,ln2e,klon,klev,nlev,klp2,klevp1,
     &     numfl2,itracso2,itracbc,itracoc,ntrac

      parameter (ne=1)    !Memory saving trick: NE=1 to use; NE=0 to save memory.
      parameter (late=ne*(lat-1)+1) ! = lat if ne = 1
      parameter (lone=ne*(lon-1)+1) ! = lon if ne = 1
      parameter (lat2e=ne*(lat2-1)+1) ! = lat2 if ne = 1
      parameter (ln2e=ne*(ln2-1)+1) ! = ln2 if ne = 1

      parameter (KLON=LN2E, KLEV=NL, NLEV=NL)
      parameter (KLP2=KLON+2, KLEVP1=KLEV+1)
c      parameter (NUMFL2=6*nl+8)
      parameter (NUMFL2=4*nl+8) !Reduced data file (no dust, seasalt)
      parameter (ITRACSO2=2) !Index for SO2 tracer
      parameter (ITRACBC=4)  !Index for BC    "
      parameter (ITRACOC=6)  !Index for OC    "
c Following must be changed several times in ukall.f as to be consistent with here.
      parameter (NTRAC=3) !No. of tracers: DMS, SO2, SO4, BCO, BCI, OCO, OCI

c Following array determines which tracers experience wet deposition

      logical lwetdep(7) !Hard-coded for maximum of 7 tracers 
      data lwetdep /.false., 2*.true., .false.,.true.,.false.,.true./
      save lwetdep
      
