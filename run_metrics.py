from os import path, getcwd
import argparse
from astropy.table import Table, Column
import numpy as np
import metrics
import config_utils
import regions

def calculate_metrics(args):
    """
    Function to run some or all of the metric functions for the science cases and survey designs

    Input:
        args    arguments object    User-specified options
    """

    # Load simulation parameters
    sim_config = config_utils.read_config(path.join(getcwd(), 'config', 'sim_config.json'))
    print('Loaded simulation configuration')

    # Load the defined survey strategy options from file
    all_survey_regions = regions.load_regions_from_file(sim_config,
                                                        path.join(getcwd(), 'region_data', 'rgps_survey_regions.json'))
    print('Loaded survey design information')

    # Load the science case definitions; this is used as a reference for their parameters
    science_cases = config_utils.read_config(path.join(getcwd(), 'config', 'rgps_science_cases.json'))

    # Build a list of the topical categories
    science_categories = []
    for author, info in science_cases.items():
        if info['ready_for_use'] and info['category'] not in science_categories:
            science_categories.append(info['category'])

    # Parse the category option
    if str(args.category).lower() == 'all':
        category_list = science_categories
    else:
        category_list = [args.category]

    # Analysing each requested category in turn,
    # load the science regions corresponding to each category from file
    for category in category_list:
        region_file = path.join(getcwd(),'region_data','rgps_science_regions_' + category + '.json')
        if not path.isfile(region_file):
            raise IOError('Missing science region file: ' + region_file)

        science_regions = regions.load_regions_from_file(sim_config, region_file)
        print('Loaded science use cases information for ' + category)

        if 'all' in str(args.survey).lower():
            survey_regions = all_survey_regions
        else:
            if args.survey in all_survey_regions.keys():
                survey_regions = {}
                survey_regions[args.survey] = all_survey_regions[args.survey]
            else:
                raise IOError('Requested survey design (' + args.survey + ') not recognized')
        print('Selected ' + str(len(science_regions)) + ' science cases and '
                 + str(len(survey_regions)) + ' survey designs to analyse')

        # Identify which metrics the user requested to run:
        all_metrics = {
            'M1_survey_footprint': metrics.M1_survey_footprint,
            'M2_star_counts': metrics.M2_star_counts,
            'M3_extended_region_count': metrics.M3_extended_region_count,
            'M4_proper_motion_precision': metrics.M4_proper_motion_precision,
            'M5_sky_area_optical_elements': metrics.M5_sky_area_optical_elements,
            'M6_sky_area_nvisits': metrics.M6_sky_area_nvisits,
        }
        if 'all' in str(args.metric).lower():
            metric_set = all_metrics
        else:
            metric_set = {}
            if args.metric in all_metrics.keys():
                metric_set[args.metric] = all_metrics[args.metric]
            else:
                raise IOError('Metric ' + args.metric + ' not recognised.  Available metrics are '
                              + ', '.join(all_metrics.keys()))
        print('Evaluating metrics: ' + ', '.join(list(metric_set)))

        # Calculate metrics
        for metric_name, metric_func in metric_set.items():
            print('Calculating metric ' + metric_name + ' for all selected science cases...')

            if metric_name in [
                'M1_survey_footprint',
                'M3_extended_region_count',
                'M5_sky_area_optical_elements',
                'M6_sky_area_nvisits',
            ]:
                results = metric_func(sim_config, science_regions, survey_regions)

            if metric_name == 'M2_star_counts':
                # Load the galatic model stellar density data
                galactic_model_file = path.join(getcwd(), 'trilegal_model_data', 'trilegal_nir_stellar_density_extinction.json')
                galactic_model_data = config_utils.read_config(galactic_model_file)
                stellar_density_data = {optic: np.array(galactic_model_data['healpix_map_' + optic])
                                        for optic in sim_config['OPTICAL_COMPONENTS']}

                results = metrics.M2_star_counts(sim_config, survey_regions, stellar_density_data)

            if metric_name == 'M4_proper_motion_precision':
                results = metrics.M4_proper_motion_precision(sim_config, survey_regions)

            # Store results
            if metric_name == 'M2_star_counts' \
                    or metric_name == 'M3_extended_region_count' \
                    or metric_name == 'M4_proper_motion_precision':
                output_file = path.join(args.data_dir, metric_name + '_results.txt')
            else:
                output_file = path.join(args.data_dir, metric_name + '_' + category + '_results.txt')
            results.write(output_file, format='ascii', delimiter=' ', overwrite=True)

    print('Completed metric analysis')

def get_args():
    """Function to gather commandline arguments"""

    parser = argparse.ArgumentParser()
    parser.add_argument('category', help='Category of science case to evaluate for')
    parser.add_argument('survey', help='Name of survey design to evaluate or ALL')
    parser.add_argument('metric', help='Name of metric to evaluate or ALL')
    parser.add_argument('data_dir', help='Path to output directory')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    calculate_metrics(args)