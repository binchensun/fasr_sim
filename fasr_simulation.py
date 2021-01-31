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
from taskinit import casalog, vptool, tb, cltool, iatool, qa
from simobserve_cli import simobserve_cli as simobserve
from tclean_cli import tclean_cli as tclean
from viewer_cli import viewer_cli as viewer

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
    (lons, lats, antnames) = kml2coords(kmlfile)
    coords2cfg(lons, lats, antnames=antnames, antidxs=None, cfgfile=cfgfile, arrayname=arrayname,
               dishdiam=dishdiam, write2casa=write2casa, verbose=verbose)


def ksc_sim_gauss(projname='sim_7m_array_1GHz_gaus', antennalist='ksc-7m.cfg', dishdiam=7, imagename=None,
                  indirection='J2000 14h26m46.0s -14d31m22.0s', incell='1arcsec', frequency='1.0GHz', inwidth='1MHz',
                  radec_offset=[100., 100.], flux=50., majoraxis='40arcsec', minoraxis='27arcsec',
                  positionangle='45.0deg', imsize=[512, 512]):
    """
    Simulate observation of an input Gaussian model
    :param projname: project name for simobserve
    :param antennalist: cfg file of the array configuration
    :param dishdiam: diameter of each dish, in meters
    :param imagename: name (and path) for output clean image, psf, etc.
    :param indirection: phase center of the observation
    :param incell: pixel scale of the model/simulated image
    :param frequency: central frequency
    :param inwidth: frequency bandwidth
    :param radec_offset: offset of the Gaussian source from the phasecenter, in arcsec
    :param flux: total flux of the Gaussian source, in solar flux unit (sfu)
    :param majoraxis: FWHM size of the Gaussian source along the major axis
    :param minoraxis: FWHM size of the Gaussian source along the minor axis
    :param positionangle: position angle of the Gaussian source
    :param imsize: size of the model/simulated image, in pixels (x and y)
    :return:
    """

    # set voltage patterns and primary beams for KSC 7 m. This is a placeholder for now (but required for PB correction)
    vp = vptool()
    if len(vp.getvp(telescope='KSC').keys()) == 0:
        vprec = vp.setpbairy(telescope='KSC', dishdiam='{0:.1f}m'.format(dishdiam),
                             blockagediam='0.75m', maxrad='1.784deg',
                             reffreq='1.0GHz', dopb=True)

    # make a Gaussian source
    cl = cltool()
    ia = iatool()
    cl.addcomponent(dir=indirection, flux=flux*1e4, fluxunit='Jy', freq=frequency, shape="Gaussian",
                    majoraxis=majoraxis, minoraxis=minoraxis, positionangle=positionangle)
    ia.fromshape("Gaussian.im", imsize + [1, 1], overwrite=True)
    cs = ia.coordsys()
    cs.setunits(['rad', 'rad', '', 'Hz'])
    cell_rad = qa.convert(qa.quantity(incell), "rad")['value']
    cs.setincrement([-cell_rad, cell_rad], 'direction')
    ra_ref = qa.toangle(indirection.split(' ')[1])
    dec_ref = qa.toangle(indirection.split(' ')[2])
    ra = ra_ref['value'] - radec_offset[0] / 3600.
    dec = dec_ref['value'] - radec_offset[1] / 3600.
    cs.setreferencevalue([qa.convert('{0:.4f}deg'.format(ra), 'rad')['value'],
                          qa.convert('{0:.4f}deg'.format(dec), 'rad')['value']], type="direction")
    cs.setreferencevalue("1.0GHz", 'spectral')
    cs.setincrement('10MHz', 'spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(), subtract=False)
    ia.close()

    simobserve(project=projname, skymodel='Gaussian.im', indirection=indirection,
               incell=incell, incenter=frequency, inwidth=inwidth, hourangle='transit', refdate='2014/11/01',
               totaltime='120s', antennalist=antennalist, obsmode='int', overwrite=True)

    if not imagename:
        imagename = projname + '/tst'
    tclean(vis=projname+'/'+projname + '.' + antennalist.split('.')[0] + '.ms', imagename=imagename,
           imsize=imsize, cell=incell, phasecenter=indirection, niter=200, interactive=False)

    viewer(projname+'/'+projname + '.' + antennalist.split('.')[0] + '.skymodel')
    viewer(imagename+'.psf')
    viewer(imagename+'.image')


