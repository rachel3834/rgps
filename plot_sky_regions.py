from os import path, getcwd
import config_utils
import regions
import matplotlib.pyplot as plt
import argparse
import numpy as np
from astropy.coordinates import Galactic, TETE, SkyCoord, ICRS
from astropy import units as u

def plot_all_regions(args):
    """
    Function to generate all-sky plots of all of the sky regions requested in community
    White Papers and Science Pitches
    """

    # Load the config file containing the definitions of the requested survey regions
    sim_config = config_utils.read_config(path.join(getcwd(), 'config', 'sim_config.json'))
    science_cases = config_utils.read_config(path.join('config', 'rgps_science_cases.json'))

    # Apply user selection of regions to plot, based on the author's name
    if str(args.region).lower() == 'all':
        author_list = science_cases.keys()
    else:
        if args.region in science_cases.keys():
            author_list = [args.region]
        else:
            raise IOError('No configuration available for region ' + args.region)

    # The configurations are grouped on a science topic basis.
    # Working through each configuration, produce a sky plot for each combination of
    # author and optical element, since the regions may be different
    # but if an author requested multiple pointings/regions for a given optical element, these
    # should be combined on the same map.
    if str(args.region).lower() != 'all':
        info = science_cases[author]

        if info['ready_for_use']:
            # Load the regions for this author from the corresponding regions JSON file
            science_regions = regions.load_regions_from_file(sim_config,
                                                             path.join(getcwd(), 'region_data',
                                                                       'rgps_science_regions_' + info['category'] + '.json'))

            print('Plotting regions for ' + author)
            for optic in sim_config['OPTICAL_COMPONENTS']:
                if optic in info.keys():
                    print(' -> optic ' + optic)

                    # Combine all regions this author has requested for the current optic
                    # into a single Celestial region
                    region_list = science_regions[optic]

                    if len(region_list) > 1:
                        r_merge = regions.combine_regions(region_list)
                    elif len(region_list) == 1:
                        r_merge = region_list[0]
                    else:
                        r_merge = None

                    # Use the CelestialRegion's built-in methods to plot a sky map of the region
                    # and save it to file
                    if r_merge:
                        r_merge.sky_plot()
                        plt.tight_layout()
                        plt.savefig(path.join('survey_maps', r_merge.label + '_' + r_merge.optic + '.png'))
                        plt.close()

                        print('--> Plotted region map for ' + r_merge.summary())


def plot_outline(skymapplot, survey_region, ssmall=5.0, outline_color='red'):
    """
    Function to plot the outline of a survey footprint, given the survey region boundaries in the form of
    survey_region = { 'l': [lmin, lmax], 'b': [bmin, bmax] }
    """

    l0 = survey_region['l'][0]
    l1 = survey_region['l'][1]
    if l0 > 0.0 and l0 > 180.0 and l1 < 180.0:
        lrangeset = [np.arange(l0, 359.9, 0.1), np.arange(0.0, l1, 0.1)]
        brangeset = [
            np.arange(survey_region['b'][0], survey_region['b'][1], 0.1),
            np.arange(survey_region['b'][0], survey_region['b'][1], 0.1),
        ]
    else:
        lrangeset = [np.arange(l0, l1, 0.1)]
        brangeset = [np.arange(survey_region['b'][0], survey_region['b'][1], 0.1)]

    # Plot a boundary region if lrange or brange has non-zero length:
    for xx in range(0, len(lrangeset), 1):
        lrange = lrangeset[xx]
        brange = brangeset[xx]

        if len(lrange) > 0 and len(brange) > 0:
            # Since the plotting package supports only scatter plots, calculate a range of points
            # to represent the outer boundaries.
            llist = []
            blist = []

            # Right-hand edge of box
            llist += [survey_region['l'][0]] * len(brange)
            blist += brange.tolist()

            # Left-hand edge of box
            llist += [survey_region['l'][1]] * len(brange)
            blist += brange.tolist()

            # Lower edge of box
            llist += lrange.tolist()
            blist += [survey_region['b'][0]] * len(lrange)

            # Upper edge of box
            llist += lrange.tolist()
            blist += [survey_region['b'][1]] * len(lrange)

            alpha = 0.4

        # Otherwise, this is a single-point region, so plot it as such
        else:
            llist = [survey_region['l'][0]]
            blist = [survey_region['b'][0]]
            alpha = 1.0

        # Add this outline to the plot; this has to be done in ICRS coordinates
        s = SkyCoord(llist, blist, frame='galactic', unit=(u.deg, u.deg))
        s = s.transform_to(ICRS)
        skymapplot.scatter(s.ra.deg * u.deg, s.dec.deg * u.deg, c=outline_color, s=ssmall, alpha=alpha)

    return skymapplot

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('region', help='Name of region to plot or all')
    args = parser.parse_args()

    return args

if __name__ == '__main__':
    args = get_args()
    plot_all_regions(args)