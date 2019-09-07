#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 19:31:53 2019

@author: cjburke
"""

import numpy as np
import matplotlib.pyplot as plt
import argparse
import sys
import datetime
import os
from urllib.parse import urlencode as urlencode

import urllib
#import urllib2

from astropy.time import Time
import astropy.units as u
from astroplan import EclipsingSystem
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord
from astropy.coordinates import (get_sun, AltAz, get_moon)

def idx_filter(idx, *array_list):
    new_array_list = []
    for array in array_list:
        new_array_list.append(array[idx])
    return new_array_list


def read_ephem_db(fileName, minDep=None, maxDur=None, obsLat=None):
    dtypeseq = ['U40']
    dtypeseq.extend(['f8'] * 7)
    dtypeseq.extend(['i4'] * 2)
    dataBlock = np.genfromtxt(fileName, delimiter='|', dtype=dtypeseq)
    gtName = dataBlock['f0']
    gtPer = dataBlock['f1']
    gtEpc = dataBlock['f2']
    gtDur = dataBlock['f3']
    gtDep = dataBlock['f4']
    gtRa = dataBlock['f5']
    gtDec = dataBlock['f6']
    gtMag = dataBlock['f7']
    gtDepGd = dataBlock['f8']
    gtDurGd = dataBlock['f9']
    print('Read in {0:d} targets'.format(len(gtPer)))
    # filter out targets too shallow
    if minDep is not None:
        idx = np.where((gtDep > minDep))[0]
        gtName = gtName[idx]
        gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, gtDepGd, gtDurGd = idx_filter(idx, \
                    gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, gtDepGd, gtDurGd)        
    # filter out targets with duration too long
    if maxDur is not None:
        idx = np.where((gtDur < maxDur))[0]
        gtName = gtName[idx]
        gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, gtDepGd, gtDurGd = idx_filter(idx, \
                    gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, gtDepGd, gtDurGd)        

    # filter out targets too far latitude away from observatory
    if obsLat is not None:
        maxDec = np.min([90.0, obsLat+90.0])
        minDec = np.max([-90.0, obsLat-90.0])
        idx = np.where((gtDec <= maxDec) & (gtDec >= minDec))[0]
        gtName = gtName[idx]
        gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, gtDepGd, gtDurGd = idx_filter(idx, \
                    gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, gtDepGd, gtDurGd)        
        
    print('After Observatory Latitutde, depth, and duration filter {0:d} targets left'.format(len(gtPer)))   

    # Make list of astroplan ephemeris objects from remaining targets
    ephemList = []
    auxList = []
    for i, curName in enumerate(gtName):
        curEpc = Time(gtEpc[i], format='jd')
        curPer = gtPer[i] * u.day
        curDur = gtDur[i] * u.hr
        curEphemObj = EclipsingSystem(primary_eclipse_time=curEpc, \
                                      orbital_period = curPer, \
                                      duration = curDur, \
                                      name = curName)
        curAuxTup = (gtRa[i], gtDec[i], gtMag[i], gtDep[i], gtDepGd[i], gtDurGd[i])
        ephemList.append(curEphemObj)
        auxList.append(curAuxTup)
    return ephemList, auxList
        

def make_ephem_db(fileName, minDep=None, maxDur=None, maxMag=None):
    # Load known transiting planet table from NEXSCI
    # build query string
    whereString = 'pl_tranflag = 1'
    if maxMag is not None:
        whereString = whereString + ' and st_optmag < {0:6.3f}'.format(maxMag)
    url = 'https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?'
    data = {'table':'planets', \
            'select':'pl_name,pl_orbper,pl_tranmid,pl_trandur,ra,dec,st_optmag,pl_trandep,pl_ratror,pl_rads,st_rad,pl_ratdor', \
            'format':'csv', \
            'where':whereString}
    url_values = urlencode(data)
    #print url_values
    queryData = urllib.request.urlopen(url + url_values)


    returnPage = queryData.read()
    dtypeseq = ['U40']
    dtypeseq.extend(['f8'] * 11)
    dataBlock = np.genfromtxt(returnPage.splitlines(), delimiter=',', skip_header=1, \
                        dtype=dtypeseq)
    gtName = dataBlock['f0']
    gtPer = dataBlock['f1']
    gtEpc = dataBlock['f2']
    gtDur = dataBlock['f3']*24.0 # Is in days convert to hours
    gtRa = dataBlock['f4']
    gtDec = dataBlock['f5']
    gtMag = dataBlock['f6']
    gtDep = dataBlock['f7']
    gtRor = dataBlock['f8']
    gtRpSol = dataBlock['f9']
    gtRs = dataBlock['f10']
    gtAoRs = dataBlock['f11']
    nT = len(gtPer)
    gtEphemGd = np.ones((nT,), dtype=np.int)
    gtDepGd = np.ones((nT,), dtype=np.int)
    gtDurGd = np.ones((nT,), dtype=np.int)
    # Check for missing ephemeris values
    idxBd = np.where((np.logical_not(np.isfinite(gtPer))) | \
                     (np.logical_not(np.isfinite(gtEpc))) | \
                     (np.logical_not(np.isfinite(gtDur))) | \
                     (np.logical_not(np.isfinite(gtDep))))[0]

    if len(idxBd) > 0:
        for curIdxBd in idxBd:
            gtEphemGd[curIdxBd] = 0
            gtDepGd[curIdxBd] = 0
            gtDurGd[curIdxBd] = 0
#            print('Bad Ephemeris for Known Planet')
            print(gtName[curIdxBd])
            # Load known multi data planet table from NEXSCI   
            whereString = 'mpl_name like \'{0}\''.format(gtName[curIdxBd])
            url = 'https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?'
            data = {'table':'exomultpars', \
                    'select':'mpl_name,mpl_orbper,mpl_tranmid,mpl_trandur,mpl_trandep,mpl_ratror,mst_rad,mpl_rads,mpl_ratdor', \
                    'format':'csv', \
                    'where':whereString}
            url_values = urlencode(data)
            #print url_values
            queryData = urllib.request.urlopen(url + url_values)
            returnPage = queryData.read()
            dtypeseq = ['U40']
            dtypeseq.extend(['f8'] * 8)
            dataBlock = np.genfromtxt(returnPage.splitlines(), delimiter=',', skip_header=1, \
                        dtype=dtypeseq)
            curName = dataBlock['f0']
            curPer = np.atleast_1d(dataBlock['f1'])
            curEpc = np.atleast_1d(dataBlock['f2'])
            curDur = np.atleast_1d(dataBlock['f3'])*24.0 # Is in days convert to hours
            curDep = np.atleast_1d(dataBlock['f4'])
            curRor = np.atleast_1d(dataBlock['f5'])
            curRs = np.atleast_1d(dataBlock['f6'])
            curRpSol = np.atleast_1d(dataBlock['f7'])
            curAoRs = np.atleast_1d(dataBlock['f8'])
                
            if curPer.size > 0:
                gotEphemValues = False
                gotDurValues = False
                gotDepValues = False
                gotRorValues = False
                gotRpSolValues = False
                gotAoRsValues = False
                for ip in range(len(curPer)):
                    if np.isfinite(curPer[ip]) and np.isfinite(curEpc[ip]) and not gotEphemValues:
                        gotEphemValues = True
                        gtPer[curIdxBd] = curPer[ip]
                        gtEpc[curIdxBd] = curEpc[ip]
                        gtEphemGd[curIdxBd] = 1
                        #print('Found Ephem {0} {1}'.format(gtPer[curIdxBd], gtEpc[curIdxBd]))
                    if np.isfinite(curDur[ip]) and not gotDurValues:
                        gotDurValues = True
                        gtDur[curIdxBd] = curDur[ip]
                        gtDurGd[curIdxBd] = 1
                        #print('Found Dur {0}'.format(gtDur[curIdxBd]))
                    if np.isfinite(curDep[ip]) and not gotDepValues:
                        gotDepValues = True
                        gtDep[curIdxBd] = curDep[ip]
                        gtDepGd[curIdxBd] = 1
                    if np.isfinite(curRor[ip]) and not gotRorValues:
                        gotRorValues = True
                        gtRor[curIdxBd] = curRor[ip]
                        #print('Found Dep {0}'.format(gtDep[curIdxBd]))
                    if np.isfinite(curRs[ip]) and np.isfinite(curRpSol[ip]) and not gotRpSolValues:
                        gotRpSolValues = True
                        gtRs[curIdxBd] = curRs[ip]
                        gtRpSol[curIdxBd] = curRpSol[ip]
                    if np.isfinite(curAoRs[ip]) and not gotAoRsValues:
                        gotAoRsValues = True
                        gtAoRs[curIdxBd] = curAoRs[ip]

                # Do a backup depth estimate from ror if dep not availble
                if not gotDepValues and gotRorValues:
                    gtDep[curIdxBd] = gtRor[curIdxBd] * gtRor[curIdxBd] * 100.0
                    gotDepValues = True
                    gtDepGd[curIdxBd] = 2
                if not gotDepValues and gotRpSolValues:
                    gtRor[curIdxBd] = gtRs[curIdxBd] / gtRs[curIdxBd]
                    gtDep[curIdxBd] = gtRor[curIdxBd] * gtRor[curIdxBd] * 100.0
                    gotDepValues = True
                    gtDepGd[curIdxBd] = 2
                # Do a backup duration estimate from a/Rs
                if not gotDurValues and gotAoRsValues and gotEphemValues:
                    gtDur[curIdxBd] = (gtPer[curIdxBd]*24.0) / 4.0 / gtAoRs[curIdxBd]
                    gotDurValues = True
                    gtDurGd[curIdxBd] = 2
                if not (gotEphemValues and gotDurValues and gotDepValues):
                    print('Not Found Part to {0} {1} {2} {3}'.format(gtName[curIdxBd], gotEphemValues, gotDurValues, gotDepValues))
            else:
                print('Not Matching to {0}'.format(gtName[curIdxBd]))
                print(curPer, curEpc)
           
    idx = np.where((np.isfinite(gtPer)) & (gtEphemGd>0) & (gtDurGd>0) & (gtDepGd>0))[0]
    gtName = gtName[idx]
    gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, gtEphemGd, gtDepGd, gtDurGd = idx_filter(idx, \
        gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, \
        gtEphemGd, gtDepGd, gtDurGd)

#    if minDep is not None:
#        idx = np.where((gtDep>minDep))[0]
#        gtName = gtName[idx]
#        gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, gtEphemGd, gtDepGd, gtDurGd = idx_filter(idx, \
#            gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, \
#            gtEphemGd, gtDepGd, gtDurGd)
#        
#    if maxDur is not None:
#        idx = np.where((gtDur<maxDur))[0]
#        gtName = gtName[idx]
#        gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, gtEphemGd, gtDepGd, gtDurGd = idx_filter(idx, \
#            gtPer, gtEpc, gtDur, gtRa, gtDec, gtMag, gtDep, \
#            gtEphemGd, gtDepGd, gtDurGd)

    print("Tot # PCs: {0:d}".format(len(gtName)))
    fp = open(fileName, 'w')
    for i, curName in enumerate(gtName):
        strout = '{0} | {1:12.7f} | {2:14.6f} | {3:7.4f} | {4:10.6f} | {5:10.5f} | {6:10.5f} | {7:6.3f} | {8:1d} | {9:1d}'.format(\
            curName, gtPer[i], gtEpc[i], gtDur[i], gtDep[i], gtRa[i], gtDec[i], \
                           gtMag[i], gtDepGd[i], gtDurGd[i])
        fp.write('{0}\n'.format(strout))
    fp.close()


if __name__ == '__main__':
    
    # Parse the command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mag", type=float, default=12.0,\
                        help="V Mag limit for target")
    parser.add_argument("-n", "--nday", type=int, default=7, \
                        help="Number of days results to show")
    parser.add_argument("-d", "--maxdur", type=float, default=8.0,\
                        help="Maximum duration [hr] to consider")
    parser.add_argument("-f", "--full", action='store_true', \
                        help="Full Transit visible above airmass constraint")
    parser.add_argument("-a", "--maxairmass", type=float, default=2.0,\
                        help="Maximum airmass to consider.  If -f/-full set then both ingress and egress must be above the constraint  If not set then only mid transit be above airmass limit")
    parser.add_argument("-p", "--depth", type=float, default=0.5, \
                        help="Minimum depth (%) for consideration")
    parser.add_argument("-g", "--longitude", type=float, default=-71.43784, \
                        help="Longitude in Degrees of observatory")
    parser.add_argument("-t", "--latitude", type=float, default=42.5792583, \
                        help="Latitude in Degrees of observatory")  
    parser.add_argument("-r", "--reload", action='store_true', \
                        help="Reload Planet Ephemerides from NASA Exoplanet Archive")  
    parser.add_argument("-fp", "--full_plus", action='store_true', \
                        help="Full Transit PLUS half duration before and after ingress/egress visible above airmass constraint")    
    args = parser.parse_args()

    args.full_plus = True
    # First check to see if the file storing Planet Ephemerides exists
    #  Retrieve if it does not or if reload argument is given.
    #  This takes a long time to complete for reasons of the odd way exoplanet
    #  archive is storing planets
    defaultFilename = 'ExoArchiveData.txt'
    if  (not os.path.isfile(defaultFilename)) or args.reload:
        # Read table of data 
        print('Start making database of known planet ephemerides')
        make_ephem_db(defaultFilename, maxMag=args.mag)
    ephemList, auxList = read_ephem_db(defaultFilename, minDep=args.depth, maxDur=args.maxdur, obsLat=args.latitude)

    # Setup the observatory parameters
    location = EarthLocation.from_geodetic(lon=args.longitude*u.deg, lat=args.latitude*u.deg, height=0.0*u.m)

    # full_plus implies full
    if args.full_plus:
        args.full = True
    # Now lets iterate through the list of valid targets and see which are observable
    #  in the next N days
    nTarg = len(ephemList)
    nDay = args.nday
    curTime = Time.now()
    print('Name | RA [deg] | Dec [deg] | RA [hms] | Dec [hms] | Vmag | Depth [%] | Duration [hr] | Period [day] | Ingress[UTC] | MidTransit[UTC] | Egress[UTC] | Ingress[SecZ] | MidTransit[SecZ] | Egress[SecZ] | TransitMidpoint [JD] | MoonIllum | MoonSep[deg]')
    for i, curEphem in enumerate(ephemList):
        curAux = auxList[i]
        # If requesting the full plus more transit duration visible then double the durations now
        if args.full_plus:
            curEphem.duration = curEphem.duration * 2.0
        # Calculate the maximum number of events that can occur in the requested window nDay
        curPer = curEphem.period.value
        nEvents = int(np.ceil(nDay/curPer))
        # Get the next midtransit times
        mid_times = curEphem.next_primary_eclipse_time(curTime, n_eclipses=nEvents)
        ineg_times = curEphem.next_primary_ingress_egress_time(curTime, n_eclipses=nEvents)
        # Setup coordinate object
        coords = SkyCoord(ra=curAux[0]*u.deg, dec=curAux[1]*u.deg)
        # Calculate sun position at these times
        sunPos = get_sun(mid_times)
        sunPosIn = get_sun(ineg_times[:,0])
        sunPosEg = get_sun(ineg_times[:,1])
        # Calculate moon position
        moonPos = get_moon(mid_times)
        # From astroplan source code calc moon illumination
        elon = sunPos.separation(moonPos)
        moonAngles = np.arctan2(sunPos.distance*np.sin(elon), \
                                moonPos.distance - sunPos.distance*np.cos(elon))
        moonIllum = (1.0 + np.cos(moonAngles))/2.0

        # Define the observers altaz system at these times
        obsAltAz = AltAz(location=location, obstime=mid_times)
        obsAltAzIn = AltAz(location=location, obstime=ineg_times[:,0])
        obsAltAzEg = AltAz(location=location, obstime=ineg_times[:,1])

        moonSep = []
        for j in range(nEvents):
            moonSep.append(coords.separation(moonPos[j].transform_to('icrs')))

        # Transform suns position to observers system
        for j in range(len(mid_times)):
            sunObsPos = sunPos[j].transform_to(obsAltAz[j])
            # See if night
            if sunObsPos.alt.degree < -6.0:
                # see if target meets airmass constraint
                targObsPos = coords.transform_to(obsAltAz[j])
                if (targObsPos.secz.value < args.maxairmass) and (targObsPos.secz.value > 0.99):
# Mid transit time meets dark sky and airmass constraints if wanting full transit check ingress
# and egress now
                    if args.full:
                        
                        sunObsPosIn = sunPosIn[j].transform_to(obsAltAzIn[j])
                        sunObsPosEg = sunPosEg[j].transform_to(obsAltAzEg[j])
                        if sunObsPosIn.alt.degree < -6.0 and sunObsPosEg.alt.degree < -6.0: # dark at both ingress and egress
                            targObsPosIn = coords.transform_to(obsAltAzIn[j])
                            targObsPosEg = coords.transform_to(obsAltAzEg[j])
                            inOK = (targObsPosIn.secz.value < args.maxairmass) and (targObsPosIn.secz.value > 0.99)
                            egOK = (targObsPosEg.secz.value < args.maxairmass) and (targObsPosEg.secz.value > 0.99)
                            if inOK and egOK:
                                durat = curEphem.duration.value
                                if args.full_plus:
                                    durat = durat / 2.0
                                strout = '{0} | {1:9.5f} | {2:9.5f} | {3} | {4} | {5:5.2f} | {6:7.4f} | {7:6.3f} | {8:9.5f} | {9} | {10} | {11} | {12:6.3f} | {13:6.3f} | {14:6.3f} | {15:16.8f} | {16:4.2f} | {17:5.1f}'.format(\
                                          curEphem.name, coords.ra.degree, coords.dec.degree,\
                                          coords.ra.to_string(u.hour), coords.dec.to_string(), \
                                          curAux[2], curAux[3], durat, curEphem.period.value, \
                                          ineg_times[j,0].isot, mid_times[j].isot, ineg_times[j,1].isot, \
                                          targObsPosIn.secz.value, targObsPos.secz.value, targObsPosEg.secz.value, \
                                          curEphem.epoch.value, moonIllum[j], moonSep[j].value)
                                print(strout)
                    else:
                        print(curEphem.name, sunObsPos.alt.degree, targObsPos.secz.value, mid_times[j])
   
        #print('help')
                
        
        

    
    print('done')
    
