c Defines the COMMON block /TM/, which was previously defined within the
c ocean model source code.
c SJP 2007/05/29

      real tmo, smo
      common/tm/tmo(lon,lat,2,2),smo(lon,lat,2,2)
