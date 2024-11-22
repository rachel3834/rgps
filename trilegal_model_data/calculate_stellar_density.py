from os import path
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import healpy as hp
import argparse
import glob
import json
import skyproj
import matplotlib.pyplot as plt

NSIDE = 64
NPIX = hp.nside2npix(NSIDE)
HPRADIUS = hp.max_pixrad(NSIDE, degrees=True)
HPAREA = hp.nside2pixarea(NSIDE, degrees = True)
def build_stellar_density_map(args):
    """
    Function to build a HEALpixel map of stellar density based on input from the Trilegal Galactic model v1.6.
    This function reads in a set of Trilegal output files generated for a grid of points across the region
    of interest.  It applies magnitude cuts to simulate those stars that would fall within Roman's
    magnitude range in the filterset offered by its Wide Field Imager.

    The code expects model data output from the following simulation:
    Galactic Model data generated from Trilegal v1.6
    Using the webform at http://stev.oapd.inaf.it/cgi-bin/trilegal

    Default settings used throughout except:
    Limiting magnitude in 4th filter =30mag
    No dust extinction
    Field area = 0.001sq.deg to keep star counts to a level where the output files are manageable.

    The star counts for each pointing are therefore scaled up by x1000.

    :param args: Commandline arguments
    :return: output HEALpixel maps as FITS binary tables
    """

    # Load the grid of Trilegal model output files.
    # Since the header information provided in these files does not include
    # the details of the simulation used to create them, the files are expected to
    # follow a systematic naming convention:
    # l<gal_longitude_degrees>_b<gal_latitude_degrees>.dat with n used to prefix negative values
    model_file_list = glob.glob(path.join(args.model_data_dir, 'l*_b*.dat'))
    model_data = [parse_trilegal_output(f) for f in model_file_list]

    # Calculate the density of stars for each point in the grid, and store the results in
    # a HEALpixel map.
    # Stars from the model output are selected based on Roman WFI's saturation and limiting magnitudes
    # These were estimated from https://roman.gsfc.nasa.gov/science/anticipated_performance_tables.html
    # The number of stars is factored because the models could only be generated for a very
    # small stamp, just 0.001sq.deg due to the excessive number of stars.  This factor scales the
    # star count to estimate that for a 1sq.deg. stamp.
    density_map = np.zeros(NPIX)
    factor = HPAREA / 0.001
    x = []
    y = []
    z = []
    for s, data in model_data:
        pixels, s = calc_healpixels_skycoord(s)
        result = np.logical_and(data['F213'] >= 16.0, data['F213'] <= 25.0)
        sidx = (np.where(result)[0]).astype(int)
        nstars = len(sidx) * 1000
        density_map[pixels] = nstars
        x.append(s.ra.deg)
        y.append(s.dec.deg)
        z.append(np.log10(nstars))
    z = np.array(z)

    # For plotting purposes, skyproj takes RA, Dec pairs of points
    fig, ax = plt.subplots(figsize=(8, 5))
    sp = skyproj.HammerSkyproj(ax=ax, galactic=True, longitude_ticks='symmetric', celestial=True)
    cmap = plt.get_cmap('viridis')
    z = z / z.max()
    for i in range(0,len(x),1):
        sp.plot(x[i], y[i], c=cmap(float(z[i])), marker='s')
    plt.show()

def calc_healpixels_skycoord(s):
    """Function converts a SkyCoord in the Galactic frame to the corrsponding
    set of HEALpixel indices.
    If the radius of the region is smaller than half that of the HEALpixel map
    resolution, then a minimum radius of 1 HEALpixel is imposed

    Input:
    :s:  SkyCoord in Galactic frame of the center of a pointing
    Output:
    :pixels:  list of HEALpixel indices covered by a circle of HPRADIUS
    """

    skycoord = s.transform_to('icrs')
    phi = np.deg2rad(skycoord.ra.deg)
    theta = (np.pi / 2.0) - np.deg2rad(skycoord.dec.deg)
    xyz = hp.ang2vec(theta, phi)
    pixels = hp.query_disc(NSIDE, xyz, HPRADIUS)

    return pixels, skycoord

def parse_trilegal_output(file_path):
    """
    Function to load the output ASCII file of the Trilegal model for a single pointing.

    Since the header information provided in these files does not include
    the details of the simulation used to create them, the files are expected to
    follow a systematic naming convention:
    l<gal_longitude_degrees>_b<gal_latitude_degrees>.dat with n used to prefix negative values

    :param file_path:
    :return: Table of galactic model output, SkyCoord of center of region used to generate the model data
    """

    # Expected columns in trilegal output for Roman filters:
    columns = [
        'Gc',
        'logAge',
        '[M/H]',
        'm_ini',
        'logL',
        'logTe',
        'logg',
        'm-M0',
        'Av',
        'm2/m1',
        'mbol',
        'F062',
        'F087',
        'F106',
        'F129',
        'F158',
        'F184',
        'F146',
        'F213',
        'SNprism',
        'Grism_1stOrder',
        'Grism_0thOrder'
        'Mact'
    ]

    # Parse the filename to get the galactic coordinates of the center of the region
    # simulated to generate the model data
    entries = str(path.basename(file_path.replace('.dat',''))).split('_')
    coords = [
        int(x.replace('l','').replace('b','').replace('n','-')) for x in entries
    ]
    field_center = SkyCoord(coords[0], coords[1], frame='galactic', unit=(u.deg, u.deg))

    if not path.isfile(file_path):
        raise IOError('Cannot find Trilegal model file ' + file_path)

    with open(file_path, 'r') as f:
        data = np.loadtxt(f)

    column_list = [Column(name=columns[i], data=data[:,i]) for i in range(0,len(columns),1)]

    return field_center, Table(column_list)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('model_data_dir', help='Path to input Trilegal model data directory')
    args = parser.parse_args()

    build_stellar_density_map(args)
