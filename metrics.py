

def M1_survey_footprint(science_cases, survey_config):
    """
    Metric to calculate how well the survey fields defined cover the survey footprint recommended in White
    Papers/Science Pitches.

    This metric makes use of the CelestialRegion.region_map pre-computed for all science cases.
    This map of HEALpixel priority for the science in question is compared with the region_map's of the
    survey definition for each optical element requested by the science case.

    The metric value represents the % of region included in survey design’s field per filter.

    Parameters:
        science_cases dict   Descriptions of the community-proposed science cases
        survey_config dict   Description of the proposed survey configuration

    Returns:
        metric        dict   Metric value calculated for all science cases

        Output format:
            results = {
                'survey_concept': {
                    'optical_element': metric_value
                }
            }
    """

    #