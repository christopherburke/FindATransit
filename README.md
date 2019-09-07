# FindATransit
Identify transiting exoplanets well placed for an observatory

Builds a database of confirmed transiting exoplanet ephemerides from the NASA exoplanet archive, and identifies transits that are well placed and timed in the night sky for an observatory to observe.
'python findatransit.py -h' - lists help
Definitely recommend 
'python findatransit.py --full_plus' - to identify transits with half a duration of observability before ingress and after egress.

First time running will build the ephemeris database and takes a while.  It saves database to a local file and if that file exists next time it is run, it will skip building and run much faster.
