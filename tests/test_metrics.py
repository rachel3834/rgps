import pytest
from os import path, getcwd
from sys import path as pythonpath
pythonpath.append(path.join(getcwd(), '..'))
import config_utils
import regions
from astropy.table import Table

@pytest.mark.parametrize(
    "test_survey_regions, test_cases, expected_results",
    [
        (
                path.join(getcwd(), 'data', 'test_survey_definition_regions.json'),
                path.join(getcwd(), 'data', 'test_science_regions.json'),
                None
        )
    ])
def test_M1_survey_footprint(test_survey_regions, test_cases, expected_results):

    from metrics import M1_survey_footprint

    # Load simulation parameters
    sim_config = config_utils.read_config(path.join(getcwd(), '..', 'config', 'sim_config.json'))

    # Load the defined survey strategy options from file
    survey_regions = regions.load_regions_from_file(sim_config, test_survey_regions)

    # Load the science cases from file
    science_regions = regions.load_regions_from_file(sim_config, test_cases)

    results = M1_survey_footprint(sim_config, science_regions, survey_regions)
    print(results)

    # Test that the metric returns a table of five columns and non-zero rows
    assert (type(results) == type(Table([])))
    assert (len(results) > 0)
    assert (len(results.colnames) == 5)

    # Test metric values returned are valid percentages
    assert (results['M1_%pix'].data >= 0.0).all() & (results['M1_%pix'].data <= 100.0).all()
    assert (results['M1_%priority'].data >= 0.0).all() & (results['M1_%priority'].data <= 100.0).all()
