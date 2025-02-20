from os import path, getcwd
import json
import config_utils
import regions
import argparse

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
    all_survey_configs = config_utils.read_config(args.spec_file)

    # Load the broader configuration parameters of the whole simulation
    sim_config = config_utils.read_config(path.join(getcwd(), 'config', 'sim_config.json'))

    # Select regions for a specific survey strategy or all depending on user input
    if 'all' in str(args.use_case).lower():
        survey_config = all_survey_configs
    else:
        survey_config = {args.use_case: all_survey_configs[args.use_case]}

    # Extract the set of regions from the set of science cases
    survey_regions = regions.build_region_maps(sim_config, survey_config)

    # Output survey regions in JSON format
    output_regions(args, sim_config, survey_regions)

def output_regions(args, sim_config, survey_regions):
    """
    Function to output a set of Celestial Regions index by survey name and optical elements
    to a JSON file.
    """

    regions = {}
    namelist = list(survey_regions.keys())
    for name in namelist:
        optic_regions = survey_regions[name]
        regions[name] = {f: [] for f in sim_config['OPTICAL_COMPONENTS']}
        for optic, region_set, in optic_regions.items():
            for r in region_set:
                regions[name][optic].append(r.to_json())

    jstr = json.dumps(regions, indent=4)

    with open(args.output_file, 'w') as f:
        f.write(jstr)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('spec_file', help='Path to the JSON file containing the region definitions')
    parser.add_argument('use_case', help='Select specific strategy or science case or ALL')
    parser.add_argument('output_file', help='Path to the JSON output file')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    build_regions(args)