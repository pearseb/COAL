c Defines the COMMON block /ETRANS3/, which was previously defined within the
c ocean model source code.
c SJP 2007/05/29

      real arxrz, aryrz, fvnedd, fuwedd, fvstedd
      common/etrans3/arxrz(imt,kmp1),aryrz(imt,kmp1),
     &     fvnedd(imt,km),fuwedd(imt,km),fvstedd(imt,km)
