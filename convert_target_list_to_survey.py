from os import path
import argparse
from astropy.table import Table, Column
import json
import config_utils

def parse_csv_targetlist(args):
    """
    Function to parse a target list in CSV format to JSON
    """

    # Load the metrics config
    sim_config = config_utils.read_config(path.join(getcwd(), 'config', 'sim_config.json'))

    # Load the target list as a Table
    target_list = Table.read(args.target_file, format='ascii.csv')

    # Definition of survey elements
    # The RGPS consists of a number of elements with distinct survey cadences.
    # These are grouped separately in the survey definitions in order to
    # evaluate the metrics more appropriately.
    # This dictionary defines the prefixes used to indicate each element
    # All regions with other prefixes are assumed to be components of the
    # wide-area survey
    survey_elements = {
        'Deep': 'panchromatic',
        'TDS': 'time_domain'
    }

    # Convert all regions defined in the table into JSON-format regions,
    # and list the filters with which they are to be observed
    # Note: skip last line of time which gives the totals
    regions = {'wide_area': {}, 'panchromatic': {}, 'time_domain': {}}
    for entry in target_list:
        if 'TOTAL' not in entry['Coverage_Type']:
            # Identify which survey element this field belongs to
            survey = 'wide_area'
            for prefix, survey_name in survey_elements.items():
                if prefix in str(entry['Fields']):
                    survey = survey_name

            r = {
                'name': str(entry['Fields']),
                'l': [float(entry['Lmin']), float(entry['Lmax'])],
                'b': [float(entry['Bmin']), float(entry['Bmax'])],
                'Filters_first_epoch': str(entry['Filters_first_epoch']),
                'Filters_second_epoch': str(entry['Filters_second_epoch']),
                'Filters_any_epoch': str(entry['Filters_any_epoch']),
                'Filters_TDS': str(entry['Filters_TDS']),
            }

            # Check for duplicate region entries
            if r['name'] in regions[survey].keys():
                raise IOError('Duplicate entry for region ' + r['name'])
            else:
                regions[survey][r['name']] = r

    # The survey regions are grouped according to the filters required for observation
    # for most metric analysis.  So the next stage is to reorganize the dictionary on a
    # per-optic basis for each survey element


    # Output new regions
    json_string = json.dumps(regions, indent=4)
    with open(args.out_path, 'w') as f:
        f.write(json_string)

def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('target_file', help='Path to target list file in CSV format')
    parser.add_argument('out_path', help='Path to output file')
    args = parser.parse_args()

    return args


if __name__ == '__main__':
    args = get_args()
    parse_csv_targetlist(args)