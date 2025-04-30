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
    science_cases = config_utils.read_config(args.spec_file)

    # Load the broader configuration parameters of the whole simulation
    sim_config = config_utils.read_config(path.join(getcwd(), 'config', 'sim_config.json'))

    # Make a list of all configured science categories for future reference
    science_categories = []
    for author, info in science_cases.items():
        if info['ready_for_use'] and info['category'] not in science_categories:
            science_categories.append(str(info['category']).lower())

    # Select regions for a specific survey strategy or all depending on user input
    survey_configs = {}
    if 'all' in str(args.use_case).lower():
        survey_configs = {x: {} for x in science_categories}

        for category in science_categories:
            for author, info in science_cases.items():
                if str(info['category']).lower() == category and info['ready_for_use']:
                    survey_configs[category][author] = info

    elif 'time_domain' in str(args.use_case).lower():
        survey_configs['time_domain'] = {name: par for name,par in science_cases.items() if par['time_domain'] and par['ready_for_use']}

    elif 'spectroscopy' in str(args.use_case).lower():
        survey_configs['spectroscopy'] = {}
        for name, par in science_cases.items():
            if 'P127' in par.keys() or 'G150' in par.keys() and par['ready_for_use']:
                survey_configs['spectroscopy'][name] = par

    elif 'extended_object_catalog' in str(args.use_case).lower():
        survey_configs['extended_object_catalog'] = {name: par for name,par in science_cases.items() if par['extended_object_catalog'] and par['ready_for_use']}

    elif str(args.use_case).lower() in science_categories:
        survey_configs[str(args.use_case)] = {name: par for name,par in science_cases.items() if args.use_case in par['category'] and par['ready_for_use']}

    else:
        if science_cases[args.use_case]['ready_for_use']:
            survey_configs['science_case'] = {args.use_case: science_cases[args.use_case]}
        else:
            raise IOError('Science case ' + args.use_case + ' flagged as not ready for use')

    print('Selected survey requirements for the following categories:')
    for category, survey_conf in survey_configs.items():
        print(category + ': includes ' + str(len(survey_conf)) + ' science cases')

    # Extract the set of regions from the set of science cases
    survey_regions = {}
    for category, survey_conf in survey_configs.items():
        print('Extracting regions for ' + category)
        survey_regions[category] = regions.build_region_maps(sim_config, survey_conf)

    # Output survey regions in JSON format
    output_regions(args, sim_config, survey_regions)

def output_regions(args, sim_config, survey_regions):
    """
    Function to output a set of Celestial Regions index by survey name and optical elements
    to a JSON file.
    """

    for category, survey_conf in survey_regions.items():
        regions = {}
        namelist = list(survey_conf.keys())
        for name in namelist:
            optic_regions = survey_conf[name]
            regions[name] = {f: [] for f in sim_config['OPTICAL_COMPONENTS']}
            for optic, region_set, in optic_regions.items():
                for r in region_set:
                    regions[name][optic].append(r.to_json())

        jstr = json.dumps(regions, indent=4)

        output_path = path.join(args.output_dir, 'rgps_science_regions_' + category + '.json')
        if category == 'science_case':
            output_path = path.join(args.output_dir, 'rgps_science_regions_' + namelist[0] + '.json')

        with open(output_path, 'w') as f:
            f.write(jstr)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('spec_file', help='Path to the JSON file containing the region definitions')
    parser.add_argument('use_case', help='Select specific author, category or ALL')
    parser.add_argument('output_dir', help='Path to the output directory')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    build_regions(args)