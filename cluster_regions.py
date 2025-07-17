# Program to build region maps for catalogs of star clusters

# For most science cases, the general-purpose rgps_survey_regions.py code can build
# a set of CelestialRegions for all regions specified for a given science case or survey region.
# However, this code can exceed the memory limits of a machine in cases where a large input
# catalog of distinct regions is used, because it generates a full HEALpixel map array for each
# target region.  While the arrays themselves can be generated, outputing the results to a single
# file tends to result in the OS' OOM killer intervening.

# This code replicates the essential functions of rgps_survey_regions.py but does not calculate
# a full HEALpixel map for each CelestialRegion, just the HEALpixel pixels overlapping the target,
# since this is sufficient for some metrics.
# While this technically fails DRY principles, it was done to avoid overloading the general-purpose
# code any further.

from os import path, getcwd
import argparse
import survey_footprints
import config_utils
import regions
import json

def build_regions(args):

    # Load the RGPS survey configurations.  These are listed according to
    science_cases = config_utils.read_config(args.spec_file)

    # Load the broader configuration parameters of the whole simulation
    sim_config = config_utils.read_config(path.join(getcwd(), 'config', 'sim_config.json'))

    # Identify the specific science case configuration requested
    if args.author not in science_cases.keys():
        raise IOError('Cannot find science case configuration for ' + args.author)

    science_strategy = science_cases[args.author]

    # Check ready for use:
    if not science_strategy['ready_for_use']:
        raise IOError('Science strategy for ' + args.author + ' marked not ready for use')

    # Build the regions for all targets without HEALpixel map arrays
    requested_regions = build_region_pixel_sets(sim_config, args, science_strategy)

    # Output survey regions in JSON format
    output_regions(args, sim_config, requested_regions)

def build_region_pixel_sets(sim_config, args, science_strategy):

    requested_regions = {}
    requested_regions[args.author] = {f: [] for f in sim_config['OPTICAL_COMPONENTS']}

    for optic in sim_config['OPTICAL_COMPONENTS']:
        if optic in science_strategy.keys():
            for region in science_strategy[optic]:
                # try:
                region['label'] = args.author + '_' + region['name']
                region['optic'] = optic
                region['extended_object_catalog'] = science_strategy['extended_object_catalog']
                region['time_domain'] = science_strategy['time_domain']
                region['topics'] = science_strategy['topics']
                try:
                    region['code'] = science_strategy['code']
                except:
                    print(args.author, science_strategy)
                    raise IOError()

                if 'category' in science_strategy.keys():
                    region['category'] = science_strategy['category']

                region_set = create_region_set_pixels_only(sim_config, region)

                for r in region_set:
                    if len(r.pixels) > 0:
                        requested_regions[args.author][optic].append(r)

    return requested_regions

def output_regions(args, sim_config, survey_conf):
    """
    Function to output a set of Celestial Regions index by survey name and optical elements
    to a JSON file.
    """

    regions = {}
    namelist = list(survey_conf.keys())
    for name in namelist:
        optic_regions = survey_conf[name]
        regions[name] = {f: [] for f in sim_config['OPTICAL_COMPONENTS']}
        for optic, region_set, in optic_regions.items():
            for r in region_set:
                regions[name][optic].append(r.to_json())

    jstr = json.dumps(regions, indent=4)

    output_path = path.join(args.output_dir, 'rgps_science_regions_' + args.author + '.json')

    with open(output_path, 'w') as f:
        f.write(jstr)

def create_region_set_pixels_only(sim_config, params):
    """
    Function to create a set of regions from a list of region dictionaries WITHOUT HEALPIXEL ARRAYS

    :param params: list of region specification dictionaries
    :return: list of regions
    """

    region_list = []

    cat_dir = path.join(sim_config['root_dir'], 'config')
    pointing_set = survey_footprints.load_catalog(sim_config, cat_dir, params['catalog'])
    ntargets = str(len(pointing_set))

    for i,pointing in enumerate(pointing_set):
        if i%10 == 0:
            print('Computing region for target ' + str(i) + ' out of ' + ntargets)
        rparams = {
            'l_center': (pointing['pointing'][0]),
            'b_center': (pointing['pointing'][1]),
            'radius': pointing['pointing'][2],
            'name': params['name'] + '_' + str(i),
            'label': params['label'],
            'optic': params['optic'],
            'nvisits': params['nvisits'],
            'duration': params['duration'],
            'visit_interval': params['visit_interval'],
            'extended_object_catalog': params['extended_object_catalog'],
            'topics': params['topics'],
            'code': params['code']
        }
        if 'category' in params.keys():
            rparams['category'] = params['category']
        if 'time_domain' in params.keys():
            rparams['time_domain'] = params['time_domain']

        r = regions.CelestialRegion(sim_config['NSIDE'], rparams)

        # Switch off the regions' HEALpixel arrays
        r.pixel_priority = None
        r.region_map = None

        # Calculate which HEALpixels belong to this region.
        r.calc_hp_healpixels_for_circular_region()

        # Convert map HEALpixel coordinates to Galactic
        r.rot_pixels()

        region_list.append(r)

    return region_list


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('spec_file', help='Path to the JSON file containing the region definitions')
    parser.add_argument('author', help='Select specific author')
    parser.add_argument('output_dir', help='Path to the output directory')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    build_regions(args)