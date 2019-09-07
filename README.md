# FindATransit
Identify transiting exoplanets well placed for an observatory

Builds a database of confirmed transiting exoplanet ephemerides from the NASA exoplanet archive, and identifies transits that are well placed and timed in the night sky for an observatory to observe.
'python findatransit.py -h' - lists help
Definitely recommend 
'python findatransit.py --full_plus' - to identify transits with half a duration of observability before ingress and after egress.  Command line options exist to filter the transits to min depth, durations, mag of target, etc.. to find the most promising for an observatories setup.

First time running will build the ephemeris database and takes a while.  It saves database to a local file and if that file exists next time it is run, it will skip building and run much faster.

This will definitely only run on python 3 and recent numpy, astropy, and astroplan.  You've been warned.  It may run into the IERS file download issue with astropy that can be fixed by 'clearing the astropy cache' . Google how to do that if you run into an error along the lines of ValueError: Column year failed to convert: invalid literal for int() with base 10: '<!'
