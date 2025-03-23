import json
from os import path, getcwd
import numpy as np
import healpy as hp
import csv
from astropy.io import fits
import config_utils
import regions

NSIDE = 64

def load_survey_footprints(sim_config, root_dir):
    """
    Function to load a set of pre-defined survey footprints in the form of HEALpixel maps

    :param root_dir:
    :return: survey_footprints dict
    """

    survey_footprints = {}

    # Load Rubin Galactic Plane survey footprint
    survey_footprints['rubin_galactic_plane'] = load_rubin_galplane_footprint(root_dir)

    # Load the DECaPS2 survey footprint
    survey_footprints['DECaPS2'] = load_DECaPS2_footprint(sim_config)

    return survey_footprints

def load_catalog(root_dir, catalog_name):
    """
    Function to load a catalog of regions, defined as centroid locations plus a radial extent from that
    centroid, which is assumed to be circular.  Centroids can be defined in RA, Dec or (l,b), but
    all quantities should be in decimal degrees.  Since the catalogs have been shared in CSV format
    from different authors, this function handles the formatting specific to each catalog.

    :param root_dir: Path to the config directory where the catalog file can be found
    :param catalog_name: Name of the catalog file
    :return:
    pointing_set: List of dictionaries in the form {"pointing": [l_center, b_center, radius]} in decimal degrees
    """

    catalog_file = path.join(root_dir, catalog_name)
    if not path.isfile(catalog_file):
        raise IOError('Cannot find requested catalog ' + catalog_file)

    pointing_set = []
    pixscale = hp.max_pixrad(NSIDE,degrees=True)

    if catalog_name == 'SFRs_for_Roman.csv':

        with open(catalog_file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for i,row in enumerate(csv_reader):
                if i > 0:
                    pointing_set.append({"pointing": [float(row[4]), float(row[5]), float(row[7])]})

    elif catalog_name == 'table-4LAC-DR3-l.fits':

        with fits.open(catalog_file) as hdul:
            data = hdul[1].data
            for row in data:
                pointing_set.append({"pointing": [float(row[4]), float(row[5]), 0.3]})

    elif catalog_name == 'baumgardt_harris_GCs.csv':

        with open(catalog_file, newline='') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ', quotechar='|')
            for i,row in enumerate(csv_reader):
                if i >= 1:
                    entries = row[0].split(',')

                    # The cluster radius is set to 1.5*half-light radius
                    # The half-light radius is r_hl(pc) / R_sun(kpc) converted to deg
                    # However, the radius needs to have a minimum of at least one HEALpixel to register
                    # on the map
                    try:
                        radius = max((1.5 * (0.001 * float(entries[16]) / float(entries[5])) * 2.0 * (180.0 / np.pi)),
                                 pixscale)
                        pointing_set.append({"pointing": [float(entries[3]), float(entries[4]), radius]})

                    # Skip malformed catalog entries
                    except ValueError:
                        pass

    elif catalog_name == 'hunt_openclusters.fits':

        with fits.open(catalog_file) as hdul:
            data = hdul[1].data
            for row in data:
                pointing_set.append({"pointing": [float(row[8]), float(row[9]), row[11]]})

    elif catalog_name == 'villasenor_combined_pointings_table.csv':

        with open(catalog_file, newline='') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ', quotechar='|')
            for i, row in enumerate(csv_reader):
                if i >= 1:
                    entries = ''.join(row).split(',')
                    radius = np.arctan((float(entries[3]) / 2.0) / (float(entries[2]) * 1000.0)) * 180.0 / np.pi
                    pointing_set.append({"pointing": [float(entries[4]), float(entries[5]), radius],
                                         "priority": float(entries[6])})

    elif catalog_name == 'blackcat_catalog.csv':

        with open(catalog_file, newline='') as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=' ', quotechar='|')
            for i, row in enumerate(csv_reader):
                if i >= 1:
                    entries = ''.join(row).split(',')
                    pointing_set.append({"pointing": [float(entries[3]), float(entries[4]), 0.1]})

    return pointing_set

def load_rubin_galplane_footprint(root_dir):
    """
    Function to load the survey footprint map for Rubin in the Galactic Plane
    :return: HEALpixel array
    """

    file_path = path.join(root_dir, 'config', 'rubin_galplane_survey_footprint.json')
    with open(file_path, 'r') as f:
        spec = json.load(f)

    return np.array(spec['healpix_map'])

def load_DECaPS2_footprint(sim_config):
    """
    Function to load the survey footprint for the DECaPS2 survey of the Galactic Plane by DECam.
    The footprint data was taken from Saydjari et al.(2023), ApJSS, 264, 28, figure 1.

    Returns HEALpixel map for consistency with other survey regions
    """

    # Since no exact table of pointings was included, the survey boundaries have been estimated
    # based on Figure 1.
    survey_config = {
        "DECaPS2": {
            "F213":[{   # Nominal filter used to fit required structure
                "l": [
                    -120.0,
                    5.0
                ],
                "b": [
                    -10.0,
                    10.0
                ],
                "nvisits": 1,
                "duration": 730.0,
                "visit_interval": [None],
                "name": "DECaPS2"
                }],
        "comment": "None",
        "ready_for_use": "True"
        }
    }

    survey_regions = regions.build_region_maps(sim_config, survey_config)

    return survey_regions['DECaPS2']['F213'][0].region_map

if __name__ == '__main__':
    root_dir = '/'
    #survey_footprints = load_survey_footprints(root_dir)
    #print(survey_footprints)

    sim_config = config_utils.read_config(path.join(getcwd(), 'config', 'sim_config.json'))

    survey_regions = load_DECaPS2_footprint(sim_config)

    print(survey_regions)