from os import path, getcwd
import config_utils
import json
import healpy as hp
import numpy as np

sim_config = config_utils.read_config(path.join(getcwd(), 'config', 'sim_config.json'))
galactic_model_file = path.join(getcwd(), 'trilegal_model_data', 'trilegal_nir_stellar_density_extinction.json')
galactic_model_data = config_utils.read_config(galactic_model_file)
stellar_density_data = {optic: np.array(galactic_model_data['healpix_map_' + optic])
                                for optic in sim_config['OPTICAL_COMPONENTS']}

NEW_NSIDE = 256

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