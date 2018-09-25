c $Log: SI.f,v $
c Revision 1.2  1994/12/07 23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.1  1993/10/06  11:25:27  ldr
c Initial revision
c
      real um(nl,nl)
      real am(nl,nl)
      real amh(nl,nl)
      common /si/ um, am, amh
