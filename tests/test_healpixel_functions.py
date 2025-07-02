import pytest
from os import path, getcwd
from sys import path as pythonpath
pythonpath.append(path.join(getcwd(), '..'))
import config_utils
import healpy as hp
import numpy as np
import regions
import copy
from astropy_healpix import HEALPix

SIM_CONFIG = config_utils.read_config(path.join(getcwd(), '..', 'config', 'sim_config.json'))
NPIX = hp.nside2npix(SIM_CONFIG['NSIDE'])

@pytest.mark.parametrize(
    "nside, test_region",
    [
        (
                SIM_CONFIG['NSIDE'],
                regions.CelestialRegion(
                        {
                            'label': "test_survey1",
                            'optic': 'F213',
                            'l_center': 18.65,
                            'b_center': 0.0,
                            'radius': 0.3,
                            'predefined_pixels': False,
                            'NSIDE': SIM_CONFIG['NSIDE'],
                            'NPIX': NPIX,
                            'pixels': np.arange(0, NPIX, 1, dtype='int'),
                            'pixel_priority': np.ones(NPIX)
                        }
                    )
        )
    ])
def test_region_pixels(nside, test_region):
    """
    Test the conversion of l,b coordinates to HEALpixels
    """

    # Use Astropy functions to calculate the HEALpixels that lie within the circular region
    # RESULTS IN GALACTIC COORDINATES!
    ahp = HEALPix(nside=nside, order='ring')
    test_region.calc_ap_healpixels_for_circular_region(ahp)
    ap_pixels = list(copy.deepcopy(test_region.pixels))
    ap_pixels.sort()

    # Use HEALpy functions to calculate the HEALpixels within the circular region
    # RESULTS IN EQUATORIAL COORDINATES!
    test_region.calc_hp_healpixels_for_circular_region()
    test_region.rot_pixels()    # Convert to Galactic coordinates
    hp_pixels = list(copy.deepcopy(test_region.pixels))
    hp_pixels.sort()

    print('Astropy: ', ap_pixels, len(ap_pixels))
    print('HEALpy: ', hp_pixels, len(hp_pixels))

    assert(ap_pixels == hp_pixels)