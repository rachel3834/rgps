from os import path
import json
import config_utils
import regions
import argparse

# Configure path to local repository
root_dir = '/'

def build_regions(args):
    """
    Function to load the Roman Galactic Plane survey configuration
    from the rgps_survey_definition.json file, convert the
    specifications into CelestialRegion objects and calculate
    the region maps for each configuration, and output these objects
    to a JSON file which can be reloaded easily without
    recalculation.
    """

    # Load the RGPS survey configurations.  These are listed according to
    survey_config = config_utils.read_config(path.join(root_dir, 'config', args.spec_file))

    # Extract the set of regions from the set of science cases
    survey_regions = regions.build_region_maps(survey_config)

    # Output survey regions in JSON format
    output_regions(args, survey_regions)

def output_regions(args, survey_regions):
    """
    Function to output a set of Celestial Regions index by survey name and optical elements
    to a JSON file.
    """

    regions = {}
    for name, region_set in survey_regions.items():
        regions[name] = {f: [] for f in SIM_CONFIG['OPTICAL_COMPONENTS']}
        for optic, r, in region_set.items():
            r.to_json()
            regions[name][optic].append(r)

    json.dumps(args.output_file, regions)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('spec_file', help='Path to JSON file containing the region definitions')
    parser.add_argument('output_file', help='Path to JSON output file')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    build_regions(args)