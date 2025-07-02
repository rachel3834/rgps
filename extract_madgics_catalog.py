# Code to extract the required positional information from the full MADGICS catalog

# Catalog of targets from MADGICS kindly provided by Andrew Saydjari on the condition
# that only positional information is made public since the full results are not
# yet published.

from os import path, getcwd
from astropy.io import fits
from astropy_healpix import HEALPix
from astropy.coordinates import SkyCoord, ICRS
from astropy import units as u
import numpy as np
import argparse
import healpy as hp
import json
import config_utils

def extract_madgics(args):

    # Load the broader configuration parameters of the whole simulation
    sim_config = config_utils.read_config(path.join(getcwd(), 'config', 'sim_config.json'))

    # Load full catalog from its FITS binary table
    with fits.open(args.input_file) as hdul:
        header = hdul[1].header
        data = hdul[1].data

    # Extract the RA, DEC, GLAT, GLONG columns with target galactic coordinates
    col_ra = -1
    col_dec = -1
    col_lat = -1
    col_long = -1
    for i,col in enumerate(hdul[1].columns):
        if col.name == 'RA': col_ra = i
        if col.name == 'DEC': col_dec = i
        if col.name == 'GLAT': col_lat = i
        if col.name == 'GLON': col_long = i

    dataset = []
    for row in data:
        dataset.append([row[col_ra][0], row[col_dec][0], row[col_lat][0], row[col_long][0]])
    dataset = np.array(dataset)

    # The full MADGICS catalog has >2.6 million objects, so it's unwealdy for our purposes.
    # Here we convert it to count the number of objects per healpixel.
    NSIDE = sim_config['NSIDE']
    NPIX = hp.nside2npix(NSIDE)
    proj = HEALPix(nside=NSIDE, order='ring', frame='icrs')
    coords = SkyCoord(dataset[:,0], dataset[:,1], frame='icrs', unit=(u.deg, u.deg))
    pixels = proj.skycoord_to_healpix(coords)

    skymap = np.zeros(NPIX)
    for i in range(0,NPIX,1):
        idx = np.where(pixels == i)[0]
        skymap[i] = len(idx)

    # Output the skymap of object counts
    jdata = {
        "nside": NSIDE,
        "healpix_resolution_deg": 0.8392936452111668,
        "n_healpix": NPIX,
        "healpix_map": skymap.tolist()
    }

    output_file = '/Users/rstreet/software/rgps/config/MADGICS_object_map.json'
    with open(output_file, 'w') as f:
        json.dump(jdata, f, indent=4)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='Path to input catalog file')
    parser.add_argument('output_file', help='Path to output catalog')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    extract_madgics(args)