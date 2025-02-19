import pytest
from os import path, getcwd
from sys import path as pythonpath
pythonpath.append(path.join(getcwd(), '..'))
import config_utils
import healpy as hp
import numpy as np
import regions

SIM_CONFIG = config_utils.read_config(path.join(getcwd(), '..', 'config', 'sim_config.json'))
NPIX = hp.nside2npix(SIM_CONFIG['NSIDE'])

@pytest.mark.parametrize(
    "test_input, expected",
    [
        (
            {
                "label": "Option 1a (Wide field=65%)",
                "optic": "F106",
                "l_center": 18.65,
                "b_center": 0.0,
                "l_width": 17.3,
                "b_height": 4.0,
                "radius": None,
                "predefined_pixels": True,
                "NSIDE": 64,
                "NPIX": 49152,
                "pixels": [33346, 33602, 33603, 33858, 33090, 33345, 33346, 33602],
                "pixel_priority": [0.0]*NPIX
            },
            regions.CelestialRegion()
        )
    ])
def test_create_region_from_json(test_input, expected):
    """Test the function that creates a CelestialObject from a dictionary"""

    from regions import create_region_from_json

    test_result = create_region_from_json(test_input)

    # Check this returns a CelestialObject
    assert(type(test_result) == type(expected))

    # Check the parameters are as expected
    keys = ['optic', 'l_center', 'b_center', 'l_width', 'b_height']
    for key in keys:
        assert(getattr(test_result, key) == test_input[key])
@pytest.mark.parametrize(
    "test_input, expected",
    [
        (
            path.join(getcwd(), 'data', 'test_survey_definition.json'),
            {
                'test_survey_1': {
                    'F213': [
                        regions.CelestialRegion(
                            {
                                'label': "test_survey1",
                                'optic': 'F213',
                                'l_center': 18.65,
                                'b_center': 0.0,
                                'l_width': 17.3,
                                'b_height': 4.0,
                                'radius': None,
                                'predefined_pixels': True,
                                'NSIDE': 64,
                                'NPIX': 49152,
                                'pixels': np.array([]),
                                'pixel_priority': np.zeros(NPIX)
                            }
                        )
                    ]
                }
            }
        )
    ])
def test_load_regions_from_file(test_input, expected):
    """Test the loading of region definitions from file"""

    from regions import load_regions_from_file

    survey_regions = load_regions_from_file(SIM_CONFIG, test_input)

    # This should return a dictionary of CelestialRegions for each survey description given
    assert(type(survey_regions) == type({}))

    # Check the parameters match for each region in the set
    keys = ['optic', 'l_center', 'b_center', 'l_width', 'b_height']
    array_keys = ['pixels', 'pixel_priority']
    for name, survey_def in survey_regions.items():
        for optic in expected[name].keys():
            for i,r in enumerate(survey_def[optic]):
                for key in keys:
                    assert (getattr(r, key) == getattr(expected[name][optic][i], key))
                for akey in array_keys:
                    assert (getattr(r, akey).all() == getattr(expected[name][optic][i], akey).all())
