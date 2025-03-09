from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np
def get_region_tiles(region_list):
    """
    Function to return a set of scatter points representing a region.
    Regions can be defined as a list of boxes in the form: {'l': [l_min, l_max], 'b': [b_min, b_max]}
    all units in degrees
    Function returns a corresponding list of the scatter points for each region
    """
    points = []

    for i, region in enumerate(region_list):
        # Convert the region into a meshgrid of points.  These need to be in RA,Dec
        nl = int(region['l'][1] - region['l'][0])
        nb = int(region['b'][1] - region['b'][0])
        l = np.linspace(region['l'][0], region['l'][1], nl)
        b = np.linspace(region['b'][0], region['b'][1], nb)
        x_1, y_1 = np.meshgrid(l, b)

        # Convert to RA, Dec coordinates
        s = SkyCoord(l=x_1, b=y_1, frame='galactic', unit=(u.deg, u.deg))
        s = s.transform_to('icrs')
        points.append(s)

    return points