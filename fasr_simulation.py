"""
Module for doing simulations for the Frequency Agile Solar Radiotelescope
"""
# History
#   2020-Dec-16 BC
#       Initialized the module with routines for reading KSC site file and making simple simulation of a
#       type III radio burst observed by VLA at ~1.2 GHz

import numpy as np
from simutil import simutil as su
import matplotlib.pyplot as plt
from astropy import units as u
import os, re
from taskinit import casalog, vptool, tb
from simobserve_cli import simobserve_cli as simobserve

def kml2coords(kmlfile='sites/KSC_Antenna_Sites.kml'):
    """
    Read Google Earth kml file and return the site information
    :param sitefile: Input kml file from Google Earth
    :return:
        antnames: names of the antennas, string array
        lons: longitudes of the antennnas in degrees, float array
        lats: latitudes of the antennnas in degrees, float array
    """
    from pykml import parser
    with open(kmlfile) as f:
        doc = parser.parse(f)
    antnames = []
    lons = []
    lats = []
    for e in doc.findall('.//{http://www.opengis.net/kml/2.2}Placemark'):
        try:
            name = e.name.name.text
        except:
            name = e.name.text

        try:
            coord = e.Point.coordinates.text.split(',')
            antnames.append(name)
            lons.append(float(coord[0]))
            lats.append(float(coord[1]))
        except:
            pass

    return lons, lats, antnames


def coords2kml(lons, lats, antnames=None, kmlfile='sites/KSC_Antenna_Sites_7m.kml'):
    """
    Use a list of longitudes and latitudes to create a kml file for Google Earth to visualize
    :param lons: longitudes of the antennnas in degrees, float array
    :param lats: latitudes of the antennnas in degrees, float array
    :param antnames: (optional) names of the antennas, string array
    :return sitefile: output Google Earth kml file
    """
    from pykml.factory import KML_ElementMaker as KML
    from lxml import etree
    if not antnames:
        antnames = [str(i+1) for i in range(len(lons))]
    fld = KML.Folder()
    for (antname, lon, lat) in zip(antnames, lons, lats):
        pm = KML.Placemark(
                KML.name(antname),
                KML.Point(KML.coordinates('{0:.13f}, {1:.13f}'.format(lon, lat)))
        )
        fld.append(pm)

    with open(kmlfile, 'w') as f:
        f.write('<kml xmlns="http://www.opengis.net/kml/2.2"> \n')
        f.write('<Document> \n')
        f.write(etree.tostring(fld, pretty_print=True))
        f.write('</Document> \n')
        f.write('</kml> \n')
        #print(etree.tostring(pm, pretty_print=True))



def site2coords(sitefile='sites/KSC_Antenna_Sites.txt', ants=[0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 12],
                 cfgfile='ksc-all.cfg', dishdiam=7, write2casa=True, verbose=True):
    """
    %%%%% This functionality has been depreciated. Preferred to use the Google Earth kml method %%%%%
    Read txt file of KSC antenna sites and write antenna configuration file in CASA
    :param sitefile: Input site file that contains longitudes and latitudes of the antenna sites (from KSC team)
    :param ants: Antenna indices to select from the input files
    :param dishdiam: diameter of each dish, in meters
    :param overwrite: Overwrite the cfg file in CASA? Default is True
    :return:
        antnames: names of the antennas, string array
        lons: longitudes of the antennnas in degrees, float array
        lats: latitudes of the antennnas in degrees, float array
    """
    # read longitude & latitude from KSC array
    rad2arcsec = (180. * 3600.) / np.pi
    f = open(sitefile, 'r')
    lines = f.readlines()
    # 3 lines for each antenna site
    antnames = []
    antidxs = []
    lons = []
    lats = []
    xs = []
    ys = []
    zs = []
    idx0 = 3
    casalog.post('Generated antenna configuration file from {0:s} as {1:s}'.format(sitefile, cfgfile))
    if verbose:
        print('Generated antenna configuration file from {0:s} as {1:s}'.format(sitefile, cfgfile))
    for i in range(0, len(lines), 3):
        line = lines[i]
        if line == '\r\n':
            break
        else:
            (antname, antidx) = line.rstrip().split(' (')
            antnames.append(antname)
            antidxs.append(antidx[:-1])
            lon = float(re.findall(r'[-+]?\d*\.\d+|\d+', lines[i + 1].rstrip().strip())[0])
            lons.append(lon)
            lat = float(re.findall(r'[-+]?\d*\.\d+|\d+', lines[i + 2].rstrip().strip())[0])
            lats.append(lat)
    return lons, lats, antnames


def coords2cfg(lons, lats, antnames=None, antidxs=None, cfgfile='ksc-7m.cfg', arrayname='KSC',
               dishdiam=7, write2casa=True, verbose=True):
    util = su()
    xs = []
    ys = []
    zs = []
    for (lon, lat) in zip(lons, lats):
        x, y, z = util.long2xyz(np.radians(lon), np.radians(lat), 3., datum='WGS84')
        xs.append(x)
        ys.append(y)
        zs.append(z)

    if antidxs:
        xs = [xs[i] for i in antidxs]
        ys = [ys[i] for i in antidxs]
        zs = [zs[i] for i in antidxs]

    xs = np.array(xs)
    ys = np.array(ys)
    zs = np.array(zs)

    f = open(cfgfile, 'w')
    f.write('#observatory=KSC\n')
    # f.write('#COFA={0:.5f}, {1:.5f}\n'.format(lats[idx0], lons[idx0]))
    f.write('#coordsys=XYZ\n')
    f.write('# x y z diam pad\n')
    for i in range(len(xs)):
        f.write('{0} {1} {2} {3} K{4:02d}\n'.format(xs[i], ys[i], zs[i], dishdiam, i))
    f.close()

    if write2casa:
        repodir = os.getenv("CASAPATH").split(' ')[0] + "/data/alma/simmos/"
        os.system('cp ' + cfgfile + ' ' + repodir)
        casalog.post('Antenna configuration file {0:s} written to CASA'.format(cfgfile))
        if verbose:
            print('Antenna configuration file {0:s} written to CASA'.format(cfgfile))

    # set voltage patterns and primary beams for KSC 7 m. This is a placeholder for now (but required for PB correction)
    vp = vptool()
    if len(vp.getvp(telescope='KSC').keys()) == 0:
        vprec = vp.setpbairy(telescope='KSC', dishdiam='{0:.1f}m'.format(dishdiam),
                             blockagediam='0.75m', maxrad='1.784deg',
                             reffreq='1.0GHz', dopb=True)

    obstable = os.getenv("CASAPATH").split(' ')[0] + "/data/geodetic/Observatories"
    tb.open(obstable, nomodify=True)
    if arrayname not in tb.getcol('Name'):
        print('{0:s} not in Observatories. Adding a record to {1:s}'.format(arrayname, obstable))
        obsdict_dsc = {'MJD': 59147.0, 'Name': 'KSC', 'Type': 'WGS84', 'Long': -80.693, 'Lat': 28.508, 'Height': 3.00,
                       'X': 0.0, 'Y': 0.0, 'Z': 0.0, 'Source': 'EOVSA Team'}
        tb.open(obstable, nomodify=False)
        nrows = tb.nrows()
        tb.addrows(1)
        for i in obsdict_dsc.keys():
            tb.putcell(i, nrows, obsdict_dsc[i])
        tb.close()

def kml2cfg(kmlfile='sites/KSC_Antenna_Sites_7m.kml', dishdiam=7, cfgfile='ksc-7m.cfg', arrayname='KSC',
                write2casa=True, verbose=True):
    """
    Read txt file of KSC antenna sites and write antenna configuration file in CASA
    :param kmlfile: Input site file that contains longitudes and latitudes of the antenna sites (from KSC team)
    :param dishdiam: diameter of each dish, in meters
    :param overwrite: Overwrite the cfg file in CASA? Default is True
    :return:
    """
    (antnames, lons, lats) = kml2coords(kmlfile)
    coords2cfg(lons, lats, antnames=antnames, antidxs=None, cfgfile=cfgfile, arrayname=arrayname,
               dishdiam=dishdiam, write2casa=write2casa, verbose=verbose)


def ksc_sim_type3(projname='sim_type3', skymodel='skymodels/sun_t191028.75_0.05s_spw0-3_mfs_I.image',
                  incell='1arcsec', incenter='1.2GHz', inwidth='1MHz', hourangle='transit',
                  refdate='2014/11/01', totaltime='120s', antennalist='ksc-all.cfg'):
    simobserve(project=projname, skymodel=skymodel, indirection='', incell=incell, incenter=incenter,
               inwidth=inwidth, obsmode='int', hourangle=hourangle, refdate=refdate, totaltime=totaltime,
               antennalist=antennalist)
