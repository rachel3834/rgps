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
import argparse

def plot_all_regions(args):
    """
    Function to generate all-sky plots of all of the sky regions requested in community
    White Papers and Science Pitches
    """

    # Load the config file containing the definitions of the requested survey regions
    config = config_utils.read_config(path.join('./config', 'rgps_science_cases.json'))

    # Full list of optical elements (filters, prism and grism) available for Roman WFI
    optical_components = ['F087', 'F106', 'F129', 'F158', 'F184', 'F213', 'F146', 'G150', 'P127']

    # Apply user selection of regions to plot, based on the author's name
    if str(args.region).lower() == 'all':
        author_list = config.keys()
    else:
        if args.region in config.keys():
            author_list = [args.region]
        else:
            raise IOError('No configuration available for region ' + args.region)

    # The configurations are grouped on a per-author basis, and may include a list of
    # regions for each filter.
    # Working through each configuration, produce a sky plot for each combination of
    # author and optical element, since the regions may be different
    # but if an author requested multiple pointings/regions for a given optical element, these
    # should be combined on the same map.
    for author in author_list:
        info = config[author]
        if info['ready_for_use']:
            print('Plotting regions for ' + author)
            for optic in optical_components:
                if optic in info.keys():
                    print(' -> optic ' + optic)

                    # Combine all regions this author has requested for the current optic
                    # into a single Celestial region
                    region_list = []
                    for params in info[optic]:
                        params['label'] = author
                        params['optic'] = optic
                        if 'catalog' in params.keys():
                            region_set = regions.create_region_set(params)
                        else:
                            region_set = [regions.create_region(params)]

                        region_list += region_set

                    r_merge = regions.combine_regions(region_list)

                    # Use the CelestialRegion's built-in methods to plot a sky map of the region
                    # and save it to file
                    r_merge.sky_plot()
                    plt.tight_layout()
                    plt.savefig(path.join('./survey_maps', r_merge.label + '_' + r_merge.optic + '.png'))
                    plt.close()

                    print('--> Plotted region map for ' + r_merge.summary())

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('region', help='Name of region to plot or all')
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = get_args()
    plot_all_regions(args)