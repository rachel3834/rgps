import healpy as hp
from os import getcwd
import survey_footprints
from mw_plot import MWSkyMap, MWSkyMapBokeh
from astropy_healpix import HEALPix
from astropy import units as u
from astropy.coordinates import Galactic, TETE, SkyCoord
import numpy as np
import matplotlib.pyplot as plt

class CelestialRegion:
    """Class to describe a region on sky, including its position and
    extend in an on-sky visualization

    Note: NSIDE = 64
    """

    def __init__(self, params={}):
        self.label = None
        self.optic = None
        self.l_center = None
        self.b_center = None
        self.l_width = None
        self.b_height = None
        self.predefined_pixels = False
        self.pixel_priority = None
        self.NSIDE = 64
        self.NPIX = hp.nside2npix(self.NSIDE)
        self.pixels = None

        for key, value in params.items():
            if key in dir(self):
                setattr(self, key, value)

        if self.l_width:
            self.halfwidth = self.l_width * u.deg / 2.0
        if self.b_height:
            self.halfheight = self.b_height * u.deg / 2.0

    def calc_ap_healpixels_for_region(self, ahp):

        if not self.predefined_pixels:
            self.skycoord = SkyCoord(self.l_center * u.deg,
                                     self.b_center * u.deg,
                                     frame=Galactic())
            self.pixels = ahp.cone_search_skycoord(self.skycoord,
                                                   self.halfwidth)

    def calc_hp_healpixels_for_circular_region(self):
        """Method calculates the HEALpixels included within the region.
        If the radius of the region is smaller than half that of the HEALpixel map
        resolution, then a minimum radius of 1 HEALpixel is imposed"""

        if not self.predefined_pixels:
            self.skycoord = SkyCoord(self.l_center * u.deg,
                                     self.b_center * u.deg,
                                     frame=Galactic())
            self.skycoord = self.skycoord.transform_to('icrs')
            phi = np.deg2rad(self.skycoord.ra.deg)
            theta = (np.pi / 2.0) - np.deg2rad(self.skycoord.dec.deg)
            radius = max(np.deg2rad(self.halfwidth.data),
                         np.deg2rad(hp.max_pixrad(self.NSIDE, degrees=True) / 2.0))
            xyz = hp.ang2vec(theta, phi)
            self.pixels = hp.query_disc(self.NSIDE, xyz, radius)

    def calc_hp_healpixels_for_box_region(self):
        """Method calculates the HEALpixels included within a box-shaped region.
        If the radius of the region is smaller than half that of the HEALpixel map
        resolution, then a minimum radius of 1 HEALpixel is imposed"""

        if not self.predefined_pixels:
            pixels = []

            l_min = self.l_center - self.halfwidth
            l_max = self.l_center + self.halfwidth
            b_min = self.b_center - self.halfheight
            b_max = self.b_center + self.halfheight

            for l in np.arange(l_min.value, l_max.value, 1.0):
                for b in np.arange(b_min.value, b_max.value, 1.0):
                    self.skycoord = SkyCoord(l, b, frame=Galactic(), unit=(u.deg, u.deg))

                    self.skycoord = self.skycoord.transform_to('icrs')
                    phi = np.deg2rad(self.skycoord.ra.deg)
                    theta = (np.pi / 2.0) - np.deg2rad(self.skycoord.dec.deg)
                    radius = np.deg2rad(1.0)
                    xyz = hp.ang2vec(theta, phi)
                    new_pixels = hp.query_disc(self.NSIDE, xyz, radius)
                    pixels += new_pixels.tolist()

            self.pixels = pixels

    def pixels_to_skycoords(self):
        """
        Method to convert pixels in a CelestialRegion back to a set of SkyCoords for the center of
        each HEALpixel

        Returns:
            SkyCoord object with multiple pointings in ICRS coordinates
        """

        hp = HEALPix(nside=self.NSIDE, order='ring', frame='icrs')
        pixels = self.pixels
        if type(pixels) == type(np.array([])):
            pixels = self.pixels.tolist()
        s = hp.healpix_to_skycoord(pixels)

        return s

    def make_map(self):
        """
        Method to create a HEALpixel array with the included pixels given a value of one and
        zero elsewhere
        """

        self.region_map = np.zeros(self.NPIX)
        self.region_map += self.pixel_priority

    def summary(self):
        return self.label + ': l_center=' + str(self.l_center) + ', b_center=' \
            + str(self.b_center) + ' n_pixels=' + str(len(self.pixels)) + ' ' + str(self.optic)

    def sky_plot(self):
        mw1 = MWSkyMap(projection='aitoff', grayscale=False, grid='galactic', background='infrared', figsize=(16, 10))
        mw1.title = self.label + ' ' + self.optic
        s = self.pixels_to_skycoords()
        mw1.scatter(s.ra.deg * u.deg, s.dec.deg * u.deg, c="r", s=5, alpha=0.4)
        plt.rcParams.update({'font.size': 22})

def create_region(params):
    """
    Function to generate a CelestialRegion object from a dictionary describing the boundaries of the region.

    Regions can be specified in the following ways:
    - Dictionary of galactic longitude, latitude [min,max] ranges: {"l": [10.0, 60.0], "b": [-2.5, 2.5]}
    - Predefined survey footprint, loaded from a HEALpixel map: {"survey_footprint": "rubin_galactic_plane"}
        The name of the map given must correspond to a map recognized by the survey_footprints module
    - Pointing given as galactic longitude, latitude and radius (deg): {"pointing": [53.37088107, -35.76976315, 0.43]}

    All parameter dictionaries should also contain keywords:
        label: String giving the name used to refer to the region
        optic: String indicating the optical element that this region should be observed with
                (These two are normally concatenated to create a distinct label)

    :input:
        params: dict  Single dictionary in one of the formats given above
    :return:
        r: CelestialRegion object
    """

    # Box regions defined as min, max ranges in (l,b):
    if 'l' in params.keys() and 'b' in params.keys():

        lspan = params['l'][1] - params['l'][0]
        bspan = params['b'][1] - params['b'][0]
        params = {
            'l_center': (params['l'][0] + lspan / 2.0) * u.deg,
            'b_center': (params['b'][0] + bspan / 2.0) * u.deg,
            'l_width': lspan,
            'b_height': bspan,
            'label': params['label'],
            'optic': params['optic']
        }
        r = CelestialRegion(params)

        # Calculate which HEALpixels belong to this region.
        # This method creates the pixels list attribute
        r.calc_hp_healpixels_for_box_region()

    # Irregular regions defined as arrays of HEALpixels which are loaded from
    # a pre-existing configuration file.
    elif 'survey_footprint' in params.keys():
        survey_regions = survey_footprints.load_survey_footprints(getcwd())
        r = CelestialRegion()
        r.label = params['label']
        r.optic = params['optic']
        r.region_map = survey_regions[params['survey_footprint']]
        r.pixels = (np.where(r.region_map > 0.0)[0]).tolist()

    # Circular regions defined by a central galactic longitude, latitude and radial extent,
    # all in units of degrees.
    elif 'pointing' in params.keys():
        params = {
            'l_center': (params['pointing'][0]),
            'b_center': (params['pointing'][1]),
            'l_width': params['pointing'][2],
            'b_height': params['pointing'][2],
            'label': params['label'],
            'optic': params['optic']
        }
        r = CelestialRegion(params)

        # Calculate which HEALpixels belong to this region.
        r.calc_hp_healpixels_for_circular_region()

    # Generate pixel map array for the region
    # If the region is valid, the list of included pixels will be non-zero.
    # Each pixel within a region is given a value of 1 - essentially being a 'vote' for that pixel,
    # for each science case.
    if len(r.pixels) > 0:
        r.pixel_priority = np.zeros(r.NPIX)
        r.pixel_priority[r.pixels] = 1.0
        r.predefined_pixels = True
        r.make_map()

    return r

def combine_regions(region_list):
    """
    Function to combine a list of CelestialRegions into a single object.
    This is designed to merge multiple regions that are defined separately for whatever reason.

    :param region_list: list of CelestialRegions with valid pixel maps

    :return: r_merge: CelestialRegion of the combined footprint
    """

    def merge_string_param(par, r_merge, r):
        val = getattr(r_merge, par)
        new_val = getattr(r, par)

        if len(val) == 0:
            val = getattr(r, par)
        elif new_val not in val:
            val = val + '_' + new_val

        setattr(r_merge, par, val)

        return r_merge

    r_merge = CelestialRegion()
    r_merge.label = ''
    r_merge.optic = ''

    map = np.zeros(r_merge.NPIX)
    pixels = np.array([], dtype='int')

    for r in region_list:
        r_merge = merge_string_param('label', r_merge, r)
        r_merge = merge_string_param('optic', r_merge, r)
        map += r.region_map
        pixels = np.concatenate((pixels, r.pixels))

    r_merge.map = map
    uniq = set()
    uniq_pixels = [int(x) for x in pixels if x not in uniq and (uniq.add(x) or True)]
    r_merge.pixels = uniq_pixels

    return r_merge
