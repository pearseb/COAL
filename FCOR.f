c Wind stress adjustment arrays re-inserted.
c SJP 2004/11/18
c
c Wind stress adjustment arrays commented out.
c SJP 2003/09/04
c
c $Log: FCOR.f,v $
c Revision 1.3  1994/12/07 23:53:36  mrd
c Explicitly declare all variables in include files.
c
c Revision 1.2  1993/10/05  13:07:53  ldr
c Changes to V4-4-15l from HBG for T63 and coupled model
c
c Revision 1.1  93/08/10  15:03:25  ldr
c Initial revision
c 
c Revision 1.6  93/03/23  12:44:53  ldr
c 
      real hfcor(imt,jmt,2),
     &     sfcor(imt,jmt,2),
     &     txcor(imt,jmt,2),
     &     tycor(imt,jmt,2)
      common /fcor/ hfcor, sfcor, txcor, tycor
