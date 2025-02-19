import sys
from os import getcwd, path
pypath = sys.path.append(path.join(getcwd(),'..'))
import config_utils
import regions

# Load simulation-wide parameters
root_dir = getcwd()
SIM_CONFIG = config_utils.read_config(path.join(root_dir, 'config', 'sim_config.json'))
# Establish 
desired_regions = {
  "survey1": {
    "F213": [
      {"l": [10.0, 60.0], "b": [-2.5, 2.5]},
      {"l": [-60.0, -10.0], "b": [-2.5, 2.5]},
      {"l": [-10.0, 30.0], "b": [-8.0,8.0]}
    ],
    "comment": "Galactic Plane Wide Area Survey",
    "proper_motions": "True",
    "time_domain": "False",
    "cadence": "multi-epoch",
    "n_visits_per_field": 8,
    "category": "wide-area",
    "ready_for_use": "True"
  }
}

requested_regions = regions.extract_requested_regions(desired_regions)
test_cases = regions.calc_healpixel_regions(requested_regions)


#@pytest.mark.parametrize(
#    "test_cases, test_config, expected_results",
#    [
#        ()
#    ])
#def test_M1_survey_footprint(test_cases, test_config, expected_results):
#    M1_survey_footprint(science_cases, survey_config)