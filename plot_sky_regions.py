from os import path
import config_utils
import regions
from regions import CelestialRegion
from astropy import units as u
from mw_plot import MWSkyMap
import matplotlib.pyplot as plt
import numpy as np
from astropy_healpix import HEALPix
from astropy.coordinates import SkyCoord

def plot_all_regions():
    """
    Function to generate all-sky plots of all of the sky regions requested in community
    White Papers and Science Pitches
    """

    # Load the config file containing the definitions of the requested survey regions
    config = config_utils.read_config(path.join('./config', 'rgps_survey_regions.json'))

    # Full list of optical elements (filters, prism and grism) available for Roman WFI
    optical_components = ['F087', 'F106', 'F129', 'F158', 'F184', 'F213', 'F146', 'G150', 'P127']

    # The configurations are grouped on a per-author basis, and may include a list of
    # regions for each filter.
    # Working through each configuration, produce a sky plot for each combination of
    # author and optical element, since the regions may be different
    # but if an author requested multiple pointings/regions for a given optical element, these
    # should be combined on the same map.
    for author, info in config.items():
        if info['ready_for_use']:
            for optic in optical_components:
                if optic in info.keys():

                    # Combine all regions this author has requested for the current optic
                    # into a single Celestial region
                    region_list = []
                    for params in info[optic]:
                        params['label'] = author
                        params['optic'] = optic
                        region_list.append(regions.create_region(params))
                    r_merge = regions.combine_regions(region_list)

                    # Use the CelestialRegion's built-in methods to plot a sky map of the region
                    # and save it to file
                    r_merge.sky_plot()
                    plt.tight_layout()
                    plt.savefig(path.join('./survey_maps', r_merge.label + '_' + r_merge.optic + '.png'))

                    print('Plotted region map for ' + r_merge.summary())
                    exit()

if __name__ == '__main__':
    plot_all_regions()