from os import getcwd, path
import survey_footprints
import regions
from mw_plot import MWSkyMap
from astropy_healpix import HEALPix
from astropy import units as u
from astropy.coordinates import Galactic, SkyCoord
import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import json

OPTICAL_COMPONENTS = ['F087', 'F106', 'F129', 'F158', 'F184', 'F213', 'F146', 'G150', 'P127']

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
        self.radius = None
        self.predefined_pixels = False
        self.NSIDE = 64
        self.NPIX = hp.nside2npix(self.NSIDE)
        self.pixels = np.array([], dtype='int')
        self.pixel_priority = np.zeros(self.NPIX)
        self.array_keys = ['pixels', 'pixel_priority']

        for key, value in params.items():
            if key not in self.array_keys and key in dir(self):
                setattr(self, key, value)

        if 'pixels' in params.keys():
            self.pixels = np.array(params['pixels'], dtype='int')
        if 'pixel_priority' in params.keys():
            self.pixel_priority = np.array(params['pixel_priority'], dtype='float')

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
            radius = max(np.deg2rad(self.radius),
                         np.deg2rad(hp.max_pixrad(self.NSIDE, degrees=True)))
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
        if len(self.pixel_priority) == self.NPIX:
            self.region_map += self.pixel_priority
        else:
            raise Warning('make_map: pixel_priority array has inconsistent number of pixels')

    def summary(self):
        return self.label + ': l_center=' + str(self.l_center) + ', b_center=' \
            + str(self.b_center) + ' n_pixels=' + str(len(self.pixels)) + ' ' + str(self.optic)

    def sky_plot(self):
        mw1 = MWSkyMap(projection='aitoff', grayscale=False, grid='galactic', background='infrared', figsize=(16, 10))
        mw1.title = self.label + ' ' + self.optic
        s = self.pixels_to_skycoords()
        mw1.scatter(s.ra.deg * u.deg, s.dec.deg * u.deg, c="r", s=5, alpha=0.4)
        plt.rcParams.update({'font.size': 22})

    def to_json(self):
        region = {
            "label": self.label,
            "optic": self.optic,
            "l_width": self.l_width,
            "b_height": self.b_height,
            "radius": self.radius,
            "predefined_pixels": self.predefined_pixels,
            "NSIDE": self.NSIDE,
            "NPIX": self.NPIX
        }
        for key in ['l_center', 'b_center']:
            datum = getattr(self, key)
            if datum is not None:
                try:
                    region[key] = datum.value
                except AttributeError:
                    region[key] = datum
            else:
                region[key] = None

        try:
            region['pixels'] = self.pixels.tolist()
        except AttributeError:
            region['pixels'] = self.pixels
        try:
            region['pixel_priority'] = self.pixel_priority.tolist()
        except:
            region['pixel_priority'] = self.pixel_priority

        return region

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
        rparams = {
            'l_center': (params['l'][0] + lspan / 2.0) * u.deg,
            'b_center': (params['b'][0] + bspan / 2.0) * u.deg,
            'l_width': lspan,
            'b_height': bspan,
            'label': params['label'],
            'optic': params['optic']
        }
        r = CelestialRegion(rparams)

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
        rparams = {
            'l_center': (params['pointing'][0]),
            'b_center': (params['pointing'][1]),
            'radius': params['pointing'][2],
            'label': params['label'],
            'optic': params['optic']
        }
        r = CelestialRegion(rparams)

        # Calculate which HEALpixels belong to this region.
        r.calc_hp_healpixels_for_circular_region()

    # Generate pixel map array for the region
    # If the region is valid, the list of included pixels will be non-zero.
    # Each pixel within a region is given a value of 1 - essentially being a 'vote' for that pixel,
    # for each science case.
    r.pixel_priority = np.zeros(r.NPIX)
    if len(r.pixels) > 0:
        r.pixel_priority[r.pixels] = 1.0
        r.predefined_pixels = True
    else:
        raise IOError('Warning: region ' + r.label + ', ' + r.optic
              + ' has zero valid pixels, parameters: ' + repr(params))
    r.make_map()

    return r

def create_region_from_json(params):

    r = CelestialRegion(params)
    r.pixels = np.array(r.pixels, dtype='int')
    r.pixel_priority = np.array(r.pixel_priority, dtype='float')
    r.make_map()
    r.predefined_pixels = np.array(r.predefined_pixels)

    return r

def create_region_set(params):
    """
    Function to create a set of regions from a list of region dictionaries.

    :param params: list of region specification dictionaries
    :return: list of regions
    """

    region_list = []

    cat_dir = path.join(getcwd(), 'config')
    pointing_set = survey_footprints.load_catalog(cat_dir, params['catalog'])
    for pointing in pointing_set:
        pointing['label'] = params['label']
        pointing['optic'] = params['optic']
        r = create_region(pointing)
        region_list.append(r)

    return region_list

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
    pix_priority = np.zeros(r_merge.NPIX)
    pixels = np.array([], dtype='int')

    for r in region_list:
        r_merge = merge_string_param('label', r_merge, r)
        r_merge = merge_string_param('optic', r_merge, r)
        map += r.region_map
        pixels = np.concatenate((pixels, r.pixels))
        pix_priority += r.pixel_priority

    r_merge.region_map = map
    r_merge.pixel_priority = pix_priority
    uniq = set()
    uniq_pixels = [int(x) for x in pixels if x not in uniq and (uniq.add(x) or True)]
    r_merge.pixels = uniq_pixels

    return r_merge

def load_regions_from_file(sim_config, file_path):
    """
    Function to load a set of Celestial Regions from file, where the region maps have been
    pre-computed for efficient handling.
    """

    region_set = {}

    with open(file_path, 'r') as f:
        content = json.load(f)

    survey_regions = {name: {} for name in content.keys()}

    for name, survey_params in content.items():

        for optic in sim_config['OPTICAL_COMPONENTS']:
            if optic in survey_params.keys():
                survey_regions[name][optic] = []

                for params in survey_params[optic]:
                    r = create_region_from_json(params)
                    survey_regions[name][optic].append(r)

    return survey_regions

def extract_requested_regions(science_cases):
    """
    Function to extract the on-sky regions requested for a selected set of science cases,
    organized according to filter.

    Parameters:
        science_cases: dict     A subset of the science cases from the rgps_science_cases.json file

    Returns:
        requested_regions dict
    """

    # Extracting the set of regions from the set of science cases
    optical_components = ['F087', 'F106', 'F129', 'F158', 'F184', 'F213', 'F146', 'G150', 'P127']
    requested_regions = {optic: [] for optic in optical_components}

    for author, info in science_cases.items():
        if info['ready_for_use']:
            for optic in optical_components:
                if optic in info.keys():
                    for region in info[optic]:
                        region['label'] = author
                        region['optic'] = optic
                        requested_regions[optic].append(region)

    return requested_regions

def calc_healpixel_regions(requested_regions):
    """
    Function to generate the CelestialRegion objects for a set of science cases, and
    calculate the HEALpixel maps.

    Parameters:
        requested_regions   dict    A subset of science cases, organized by filter

    Returns:
        desired_regions     dict of CelestialRegions, organized by filter
    """

    # Repackaging the full set of regions according to the regions requested for each filter,
    # and computing the HEALpixels included in each requested region
    desired_regions = {}

    for optic, region_list in requested_regions.items():
        regions_for_optic = []
        for box in region_list:
            if 'catalog' in box.keys():
                region_set = regions.create_region_set(box)
            else:
                region_set = [regions.create_region(box)]

            for r in region_set:
                # If the region is valid, the list of included pixels will be non-zero.
                # Each pixel within a region is given a value of 1 - essentially being a 'vote' for that pixel,
                # for each science case.
                if len(r.pixels) > 0:
                    r.pixel_priority = np.zeros(r.NPIX)
                    r.pixel_priority[r.pixels] = 1.0
                    r.predefined_pixels = True
                    r.make_map()

                    regions_for_optic.append(r)

        desired_regions[optic] = regions_for_optic

    return desired_regions

def combine_regions_per_filter(desired_regions):
    """
    Convenience function to combine the HEALpixel regions for a set of science cases organized
    per filter.

    Parameters:
        science_cases: dict     A subset of the science cases from the rgps_science_cases.json file

    Returns:
        combined_regions dict
    """

    # Building combined maps of all HEALpixels requested by all science cases in each filter
    combined_regions = {}

    # In order to use the plotting method of the CelestialRegion object, we can create separate regions for the combined maps
    for optic, region_list in desired_regions.items():
        if len(region_list) > 0:
            r_merge = regions.combine_regions(region_list)
            r_merge.optic = optic
            r_merge.label = 'Combined survey footprint'

            combined_regions[optic] = r_merge

    return combined_regions

def build_region_maps(sim_config, survey_definitions):
    """
    Function to calculate the region maps for a dictionary of CelestialRegions
    index by author or name.
    """

    requested_regions = {}

    for name, info in survey_definitions.items():
        if info['ready_for_use']:
            requested_regions[name] = {f: [] for f in sim_config['OPTICAL_COMPONENTS']}

            for optic in sim_config['OPTICAL_COMPONENTS']:
                if optic in info.keys():
                    for region in info[optic]:
                        region['label'] = name
                        region['optic'] = optic

                        if 'catalog' in region.keys():
                            region_set = create_region_set(region)

                        else:
                            region_set = [create_region(region)]

                        for r in region_set:
                            # If the region is valid, the list of included pixels will be non-zero.
                            # Each pixel within a region is given a value of 1 - essentially being a 'vote' for that pixel,
                            # for each science case.
                            if len(r.pixels) > 0:
                                r.pixel_priority = np.zeros(r.NPIX)
                                r.pixel_priority[r.pixels] = 1.0
                                r.predefined_pixels = True
                                r.make_map()

                                requested_regions[name][optic].append(r)

    return requested_regions

