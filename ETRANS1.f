c Defines the COMMON block /ETRANS1/, which was previously defined within the
c ocean model source code.
c SJP 2007/05/29

      real uedd, vedd, wedd
      common/etrans1/uedd(imt,km),vedd(imt,km),wedd(imt,kmp1)
