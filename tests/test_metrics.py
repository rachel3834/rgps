import pytest
from os import path, getcwd
from sys import path as pythonpath
pythonpath.append(path.join(getcwd(), '..'))
import config_utils
import regions
from astropy.table import Table
import numpy as np

@pytest.mark.parametrize(
    "test_survey_regions_file, test_science_regions_file",
    [
        (
                path.join(getcwd(), 'data', 'test_m1_survey_regions1.json'),
                path.join(getcwd(), 'data', 'test_m1_science_regions1.json')
        ),
        (
                path.join(getcwd(), 'data', 'test_m1_survey_regions2.json'),
                path.join(getcwd(), 'data', 'test_m1_science_regions2.json')
        )
    ])
def test_M1_survey_footprint(test_survey_regions_file, test_science_regions_file):
    """
    Unittest calculates the percentage of HEALpixels included in the survey
    footprints and their priorities for each science case
    """

    from metrics import M1_survey_footprint

    # Load simulation parameters
    sim_config = config_utils.read_config(path.join(getcwd(), '..', 'config', 'sim_config.json'))

    # Load the defined survey strategy options from file
    survey_regions = regions.load_regions_from_file(sim_config, test_survey_regions_file)

    # Load the science cases from file
    science_regions = regions.load_regions_from_file(sim_config, test_science_regions_file)

    results = M1_survey_footprint(sim_config, science_regions, survey_regions)

    # Test that the metric returns a table of five columns and non-zero rows
    assert (type(results) == type(Table([])))
    assert (len(results) > 0)
    assert (len(results.colnames) == 5)

    # Test metric values returned are valid percentages
    assert (results['M1_%pix'].data >= 0.0).all() & (results['M1_%pix'].data <= 100.0).all()
    assert (results['M1_%priority'].data >= 0.0).all() & (results['M1_%priority'].data <= 100.0).all()

    # Second test case is designed to verify output in the case of 100% overlap
    if 'test_m1_survey_regions2.json' in path.basename(test_survey_regions_file):
        assert (results['M1_%pix'][0] == 100.0)
        assert (results['M1_%priority'][0] == 100.0)

@pytest.mark.parametrize(
    "test_survey_regions, galactic_model_file",
    [
        (
                path.join(getcwd(), 'data', 'test_survey_definition_regions.json'),
                path.join(getcwd(), '..', 'trilegal_model_data', 'trilegal_nir_stellar_density_extinction.json'),
        )
    ])
def test_M2_star_counts(test_survey_regions, galactic_model_file):
    """
    Unittest calculates the metric to evaluate the total star counts within each survey footprint
    """

    from metrics import M2_star_counts

    # Load simulation parameters
    sim_config = config_utils.read_config(path.join(getcwd(), '..', 'config', 'sim_config.json'))

    # Load the defined survey strategy options from file
    survey_regions = regions.load_regions_from_file(sim_config, test_survey_regions)

    # Load the galatic model stellar density data
    galactic_model_data = config_utils.read_config(galactic_model_file)
    stellar_density_data = {optic: np.array(galactic_model_data['healpix_map_'+optic])
                                    for optic in sim_config['OPTICAL_COMPONENTS']}
    # Calculate metrics
    results = M2_star_counts(sim_config, survey_regions, stellar_density_data)

    # Test that the metric returns a table of five columns and non-zero rows
    assert (type(results) == type(Table([])))
    assert (len(results) > 0)
    assert (len(results.colnames) == 3)

    # Test metric values returned are credible star counts
    assert ((results['M2_nstars'].data >= 0.0).all())


@pytest.mark.parametrize(
    "test_survey_regions, test_cases",
    [
        (
                path.join(getcwd(), 'data', 'test_m3_survey_regions1.json'),
                path.join(getcwd(), 'data', 'test_m3_science_regions1.json')
        )
    ])
def test_M3_extended_regions(test_survey_regions, test_cases):
    """
    Unittest calculates the metric to evaluate how many entries from a catalog of extended regions
    are included within the survey footprints
    """

    from metrics import M3_extended_region_count

    # Load simulation parameters
    sim_config = config_utils.read_config(path.join(getcwd(), '..', 'config', 'sim_config.json'))

    # Load the defined survey strategy options from file
    survey_regions = regions.load_regions_from_file(sim_config, test_survey_regions)

    # Load the science cases from file
    science_regions = regions.load_regions_from_file(sim_config, test_cases)

    # Calculate metric
    results = M3_extended_region_count(sim_config, science_regions, survey_regions)
    print(results)

    # Test that the metric returns a table of five columns and non-zero rows
    assert (type(results) == type(Table([])))
    assert (len(results) > 0)
    assert (len(results.colnames) == 5)

    # Test metric values returned are valid percentages and the known result that the pixels for
    # one of the requested regions was included in the F129 survey definition.
    assert (results['M3_%regions'][0] > 0.0)
    assert ((results['M3_%regions'].data >= 0.0).all())

@pytest.mark.parametrize(
    "test_survey_regions, expected",
    [
        (
                path.join(getcwd(), 'data', 'test_survey_regions2.json'),
                [100.0, 0.0, 0.0]
        )
    ])
def test_M4_proper_motion_precision(test_survey_regions, expected):
    """
    Unittest for metric to evaluate the area of sky to receive observations to measure proper motions
    to the required precision.
    """

    from metrics import M4_proper_motion_precision

    # Load simulation parameters
    sim_config = config_utils.read_config(path.join(getcwd(), '..', 'config', 'sim_config.json'))

    # Load the defined survey strategy options from file
    survey_regions = regions.load_regions_from_file(sim_config, test_survey_regions)

    # Compute metric
    results = M4_proper_motion_precision(sim_config, survey_regions)

    # Test that the metric returns a table of five columns and non-zero rows
    assert (type(results) == type(Table([])))
    assert (len(results) > 0)
    assert (len(results.colnames) == 4)

    # Test for expected metric results
    for i, expected_value in enumerate(expected):
        assert(expected_value == results['M4_proper_motion_precision'][i])

@pytest.mark.parametrize(
    "test_survey_regions, test_science_cases",
    [
        (
                path.join(getcwd(), 'data', 'test_m5_survey_regions1.json'),
                path.join(getcwd(), 'data', 'test_m5_science_regions1.json')
        )
    ])
def test_M5_sky_area_optical_elements(test_survey_regions, test_science_cases):
    """
    Unittest for metric to evaluate the area of sky to receive observations in one or more filters,
    with the filtersets parameter providing the tuples of filters combinations to check for.
    """

    from metrics import M5_sky_area_optical_elements

    # Load simulation parameters
    sim_config = config_utils.read_config(path.join(getcwd(), '..', 'config', 'sim_config.json'))

    # Load the science cases from file
    science_regions = regions.load_regions_from_file(sim_config, test_science_cases)

    # Load the defined survey strategy options from file
    survey_regions = regions.load_regions_from_file(sim_config, test_survey_regions)

    # Calculate metrics
    results = M5_sky_area_optical_elements(sim_config, science_regions, survey_regions)

    # Test that the metric returns a table of five columns and non-zero rows
    assert (type(results) == type(Table([])))
    assert (len(results) > 0)
    assert (len(results.colnames) == 4)

    # Test metric values returned are valid percentages and the known result that the pixels for
    # one of the requested regions was included in the F129 survey definition.
    assert (len(np.where(results['M5_sky_area_single_filter'].data >= 0.0)[0]) > 0)
    assert (len(np.where(results['M5_sky_area_filter_combo'].data >= 0.0)[0]) > 0)

@pytest.mark.parametrize(
    "test_survey_regions, test_cases",
    [
        (
                path.join(getcwd(), 'data', 'test_m6_survey_regions1.json'),
                path.join(getcwd(), 'data', 'test_m6_science_regions1.json')
        )
    ])
def test_M6_sky_area_nvisits(test_survey_regions, test_cases):
    """
    Unittest for metric to evaluate the percentage of the desired survey region for each science
    case to receive the desired number of visits in each filter.
    """

    from metrics import M6_sky_area_nvisits

    # Load simulation parameters
    sim_config = config_utils.read_config(path.join(getcwd(), '..', 'config', 'sim_config.json'))

    # Load the defined survey strategy options from file
    survey_regions = regions.load_regions_from_file(sim_config, test_survey_regions)

    # Load the science cases from file
    science_regions = regions.load_regions_from_file(sim_config, test_cases)

    # Compute metric
    results = M6_sky_area_nvisits(sim_config, science_regions, survey_regions)
    print(results)

    # Test that the metric returns a table of five columns and non-zero rows
    assert (type(results) == type(Table([])))
    assert (len(results) > 0)
    assert (len(results.colnames) == 7)

    # Test metric value equals 100% for one of the tests
    assert(len(np.where(results['M6_%sky_area_nvisits'].data == 100.0)[0] > 0))

@pytest.mark.parametrize(
    "test_cases",
    [
        (
                path.join(getcwd(), 'data', 'test_science_regions_defurio.json')
        )
    ])
def test_extract_multiband_science(test_cases):
    """
    Test for the function to identify multi-filter testcases
    """

    from metrics import extract_multiband_science

    # Load simulation parameters
    sim_config = config_utils.read_config(path.join(getcwd(), '..', 'config', 'sim_config.json'))

    # Load the science cases from file
    science_regions = regions.load_regions_from_file(sim_config, test_cases)

    # Extract the set of multi-band science cases
    multiband_cases = extract_multiband_science(sim_config, science_regions)

    # Check the results are consistent with the input test data
    input_fset = []
    for f in sim_config['OPTICAL_COMPONENTS']:
        if len(science_regions['De_Furio'][f]) > 0:
            input_fset.append(f)
    input_fset.sort()

    assert(multiband_cases['De_Furio']['filterset'] == input_fset)
