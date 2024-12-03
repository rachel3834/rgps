import json
from os import path
import numpy as np
import csv
from astropy.io import fits

def load_survey_footprints(root_dir):
    """
    Function to load a set of pre-defined survey footprints in the form of HEALpixel maps

    :param root_dir:
    :return: survey_footprints dict
    """

    survey_footprints = {}

    # Load Rubin Galactic Plane survey footprint
    survey_footprints['rubin_galactic_plane'] = load_rubin_galplane_footprint(root_dir)

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


if __name__ == '__main__':
    root_dir = '/Users/rstreet/software/rgps'
    survey_footprints = load_survey_footprints(root_dir)
    print(survey_footprints)