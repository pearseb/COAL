c Defines common block ORBRAD, which contains the epoch, the solar
c constant and the Earth's three orbital parameters.
c
c SJP 2001/12/23

c Common block ORBRAD contains:
c
c BPYEAR   The epoch, in years before 1950
c
c CSOLAR   The solar constant, in W/m^2
c
c EC       The eccentricity of the Earth's orbit
c
c PERIL    The longitude of the Earth's perihelion, in degrees
c
c OBLQTY   The obliquity of the Earth's axis, in degrees

       integer bpyear
       real csolar, ec, peril, oblqty
       common /orbrad/ bpyear, csolar, ec, peril, oblqty
