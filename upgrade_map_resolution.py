import copy
from os import path, getcwd
import config_utils
import json
import healpy as hp
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('in_file_path', help='Path to input file')
parser.add_argument('out_file_path', help='Path to output file')
args = parser.parse_args()

sim_config = config_utils.read_config(path.join(getcwd(), 'config', 'sim_config.json'))
NEW_NSIDE = 256
NEW_NPIX = hp.nside2npix(NEW_NSIDE)

gal_model = False
if gal_model:
    galactic_model_file = path.join(getcwd(), 'trilegal_model_data', 'trilegal_nir_stellar_density_extinction.json')
    galactic_model_data = config_utils.read_config(galactic_model_file)
    stellar_density_data = {optic: np.array(galactic_model_data['healpix_map_' + optic])
                                    for optic in sim_config['OPTICAL_COMPONENTS']}

    output = {
            "label": "Trilegal_v1.6_log10_density",
            "nside": NEW_NSIDE,
            "healpix_resolution_deg": hp.nside2pixarea(NEW_NSIDE),
            "n_healpix": hp.nside2npix(NEW_NSIDE),
        }
    for optic, map in stellar_density_data.items():
        new_map = hp.ud_grade(map, NEW_NSIDE)
        output['healpix_map_'+optic] = new_map.tolist()

    # Serializing json
    json_object = json.dumps(output, indent=4)

    # Writing to sample.json
    with open(path.join('trilegal_model_data', 'trilegal_nir_stellar_density_extinction_'+str(NEW_NSIDE)+'.json'), "w") as f:
        f.write(json_object)

else:

    region_config_data = config_utils.read_config(args.in_file_path)
    new_region_config = copy.deepcopy((region_config_data))

    for case_key, case_config in region_config_data.items():
        new_case_config = copy.deepcopy((case_config))

        for f in sim_config['OPTICAL_COMPONENTS']:
            new_case_config[f] = []

            for r_config in case_config[f]:
                new_r_config = copy.deepcopy(r_config)

                old_map = np.zeros(r_config['NPIX'])
                old_map[r_config['pixels']] = 1.0

                old_priority_map = np.array(r_config['pixel_priority'])

                new_map = hp.ud_grade(old_map, NEW_NSIDE)
                new_priority_map = hp.ud_grade(old_priority_map, NEW_NSIDE)

                new_r_config['NSIDE'] = NEW_NSIDE
                new_r_config['NPIX'] = NEW_NPIX
                new_r_config['pixel_priority'] = new_priority_map.tolist()
                pixels = np.where(new_map > 0.0)[0]
                new_r_config['pixels'] = pixels.tolist()

                new_case_config[f].append(new_r_config)

        new_region_config[case_key] = new_case_config

    json_object = json.dumps(new_region_config, indent=4)

    with open(args.out_file_path, "w") as f:
        f.write(json_object)
