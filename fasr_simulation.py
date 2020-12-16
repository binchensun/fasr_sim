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
from taskinit import casalog
from taskinit import vptool
from simobserve_cli import simobserve_cli as simobserve


def ksc_site2cfg(sitefile='sites/KSC_Antenna_Sites.txt', ants=[0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 12],
                 cfgfile='ksc-all.cfg', dishdiam=7, write2casa=True, verbose=True):
    """
    Read txt file of KSC antenna sites and write antenna configuration file in CASA
    :param sitefile: Input site file that contains longitudes and latitudes of the antenna sites (from KSC team)
    :param ants: Antenna indices to select from the input files
    :param dishdiam: diameter of each dish, in meters
    :param overwrite: Overwrite the cfg file in CASA? Default is True
    :return:
    """
    rad2arcsec = (180. * 3600.) / np.pi
    util = su()
    # read longitude & latitude from KSC array
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
            x, y, z = util.long2xyz(np.radians(lon), np.radians(lat), 3., datum='WGS84')
            xs.append(x)
            ys.append(y)
            zs.append(z)
            # use antenna site B as the reference
            if i == idx0:
                (x0, y0, z0) = (x, y, z)

    xs_all = np.array(xs)[ants]
    ys_all = np.array(ys)[ants]
    zs_all = np.array(zs)[ants]
    # antnames = np.array(antnames)[ants]
    # antidxs = np.array(antidxs)[ants]

    f = open(cfgfile, 'w')
    f.write('#observatory=KSC\n')
    # f.write('#COFA={0:.5f}, {1:.5f}\n'.format(lats[idx0], lons[idx0]))
    f.write('#coordsys=XYZ\n')
    f.write('# x y z diam pad\n')
    for i in range(len(xs_all)):
        f.write('{0} {1} {2} {3} K{4:02d}\n'.format(xs_all[i], ys_all[i], zs_all[i], dishdiam, i))
    f.close()
    casalog.post('Generated antenna configuration file from {0:s} as {1:s}'.format(sitefile, cfgfile))
    if verbose:
        print('Generated antenna configuration file from {0:s} as {1:s}'.format(sitefile, cfgfile))

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


def ksc_sim_type3(projname='sim_type3', skymodel='skymodels/sun_t191028.75_0.05s_spw0-3_mfs_I.image',
                  incell='1arcsec', incenter='1.2GHz', inwidth='1MHz', hourangle='transit',
                  refdate='2014/11/01', totaltime='120s', antennalist='ksc-all.cfg'):
    simobserve(project=projname, skymodel=skymodel, indirection='', incell=incell, incenter=incenter,
               inwidth=inwidth, obsmode='int', hourangle=hourangle, refdate=refdate, totaltime=totaltime,
               antennalist=antennalist)
