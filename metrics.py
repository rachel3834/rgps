
NSIDE = 64
OPTICAL_COMPONENTS = ['F087', 'F106', 'F129', 'F158', 'F184', 'F213', 'F146', 'G150', 'P127']

def M1_survey_footprint(science_cases, survey_config):
    """
    Metric to calculate how well the survey fields defined cover the survey footprint recommended in White
    Papers/Science Pitches.

    This metric makes use of the CelestialRegion.region_map pre-computed for all science cases.
    This map of HEALpixel priority for the science in question is compared with the region_maps of the
    survey definition for each optical element requested by the science case.

    The metric values returned are:
        percent_map: % of pixels in the desired region included in survey design’s footprint per filter
        priority_map: % of the summed pixel priority values in the desired region included
                        in survey design’s footprint per filter

    Parameters:
        science_cases dict   CelestialRegions representing the community-proposed science cases
        survey_config dict   Description of the proposed survey configuration

    Returns:
        metric        dict   Metric value calculated for all science cases

        Output format:
            results = {
                'survey_concept': {
                    'optical_element': {
                        'percent_map': array of metric_values for each science case,
                        'percent_priority': array of metric_values for each science case,
                        'science_case': list of names of the science cases
                    }
                }
            }
    """

    results = {}

    # For each survey description, loop over all optical components and compare the HEALpixel
    # maps of the survey footprint with the region map for each science case
    for survey_name, region_set in survey_config.items():

        # Loop over all optical components since the requested footprints can be different
        for f in OPTICAL_COMPONENTS:
            rsurvey = region_set[f]

            # Compare the survey definition with each science case
            m1 = []
            m2 = []
            cases = []

            # This assumes one region per science case per filter
            for name, info in science_cases.items():
                cases.append(name)

                # Calculate the number of overlapping pixels, and metric values
                common_pixels = list(set(info[f].pixels).intersection(set(rsurvey.pixels)))

                m1.append((len(common_pixels) / len(info[f].pixels)) * 100.0)
                m2.append((info[f][common_pixels].sum() / info[f].sum()) * 100.0)

            results[survey_name][f] = {
                'percent_map': np.array(m1),
                'percent_priority': np.array(m2),
                'science_cases': cases
            }

    return results

def M2_star_counts(survey_config, galactic_model):
    """
    Metric calculates the number of stars are included in the survey.

    This metric uses star count data derived from the Trilegal model to estimate the
    total number of stars included within each survey footprint for each optical element.

    The metric value return is the total number of stars included.

    Parameters:
        survey_config   dict   Description of the proposed survey configuration
        galactic_model  dict   HEALpixel maps of the Trilegal stellar density data per filter

    Returns:
        metric        dict   Metric value calculated for all science cases

        Output format:
            results = {
                'survey_concept': {
                    'star_count': array of metric_values per optical element
                    }
                }
            }
    """

    results = {}

    PIXAREA = hp.nside2pixarea(NSIDE, degrees=True)

    for survey_name, region_set in survey_config.items():
        metric = np.array(len(OPTICAL_COMPONENTS))

        # Loop over all optical components since the requested footprints can be different
        for k,f in enumerate(OPTICAL_COMPONENTS):
            rsurvey = region_set[f]

            metric[k] = galactic_model[rsurvey.pixels] * PIXAREA

        results[survey_name] = {'star_counts': metric}

    return results

