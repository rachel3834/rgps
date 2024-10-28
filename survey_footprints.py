import json
from os import path
import numpy as np

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