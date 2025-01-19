
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
        results        dict   Metric value calculated for all science cases

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
        results        dict   Metric value calculated for all science cases

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

def M3_extended_region_count(survey_config, science_cases):
    """
    Metric to evaluate the number from a set of extended-region targets that lie fully
    within the survey footprint.

    This metric is intended for targets that span an extended area on sky, such as star
    clusters and Star Forming Regions, as opposed to point sources.  The desired set of
    targets are parsed into a set of CelestialRegions with HEALpixel maps
    giving the HEALpixels included in the target region.

    The metric values calculated are:
        ntarget: number of targets fully included within the survey footprint

    Parameters:
        survey_config   dict   Description of the proposed survey configuration
        science_cases  dict   Catalog of target regions represented as science cases indexed
                                by catalog name, including (l,b) coordinates and
                                radial extent

    Returns:
        results        dict   Metric value calculated for all science cases

        Output format:
            results = {
                'survey_concept': {
                    'extended_region_count': array of metric_values per optical element
                    }
                }
            }
    """

    # Build region dictionaries for all targets from a catalog for the science cases given
    requested_regions = regions.extract_requested_regions(science_cases)

    # Generate CelestialRegions for each target and calculate the HEALpixel regions within them,
    # now indexed by optical element
    desired_regions = regions.calc_healpixel_regions(requested_regions)

    results = {}

    for survey_name, region_set in survey_config.items():
        metric = np.array(len(OPTICAL_COMPONENTS))

        # Loop over all optical components since the requested footprints can be different
        for k, f in enumerate(OPTICAL_COMPONENTS):
            rsurvey = region_set[f]

            # Calculate the number of HEALpixels that are both within the target region and
            # the survey footprint for each target region
            in_pixels = np.array(
                [ len(list(set(r.pixels).intersection(set(rsurvey.pixels)))) for r in desired_regions[f] ]
            )

            # Calculate the expected number of pixels for each target region
            r_pixels = np.array(
                [len(list(r.pixels)) for r in desired_regions[f]]
            )

            # Metric value is the number of regions where in_pixel >= r_pixels
            # (due to the HEALpixels providing irregular coverage of the regions)
            metric[k] = len(np.where(in_pixels >= r_pixels)[0])

        results[survey_name] = {'extended_region_count': metric}

    return results

def M5_proper_motions(survey_config, science_cases, req_interval=730.0):
    """
    Metric to evaluate the sky area that receives at least 2 observations
    separated by AT LEAST the required interval, as a percentage of the desired sky region
    from different science cases in different filters.  This metric is designed to evaluate
    how well the survey performs for the measurement of proper motions.

    Parameters:
        survey_config   dict   Description of the proposed survey configuration
        science_cases  dict   Catalog of science cases that require multiple epochs
                                of visits per field
        req_interval    float   Desired number of days between sequential observations
                                of the same field
    Returns:
        results        dict   Metric value calculated for all science cases

        Output format:
            results = {
                'survey_concept': {
                    'optical_element': {
                        'percent_area_at_interval': array of metric_values for each science case,
                        'science_case': list of names of the science cases
                    }
                }
            }
    """

    results = {}

    for survey_name, region_set in survey_config.items():

        # Loop over all optical components since the requested footprints can be different
        for k, f in enumerate(OPTICAL_COMPONENTS):
            rsurvey = region_set[f]

            # Create a pixel map of the survey footprint, filling the pixel values with
            # the interval between sequential visits.
            # XXX FOR THIS TO WORK CELESTIALREGIONS need the cadence info
            if rsurvey.n_visits_per_field >= 2:
                idx = np.where(rsurvey.region_map > 0.0)[0]
                cadence_footprint = np.zeros(rsurvey.NPIX)
                cadence_footprint[idx].fill(rsurvey.visit_interval)
            else:
                cadence_footprint = np.zeros(rsurvey.NPIX)

            cases = []
            metric = []

            # This assumes one region per science case per filter
            for name, info in science_cases.items():
                if 'multi-epoch' in info['cadence']:
                    cases.append(name)

                    # Similarly for the science cases, create a desired cadence_footprint
                    if rsurvey.n_visits_per_field >= 2:
                        idx = np.where(info[f].region_map > 0.0)[0]
                        science_cadence_footprint = np.zeros(info[f].NPIX)
                        science_cadence_footprint[idx].fill(info[f].visit_interval)
                    else:
                        science_cadence_footprint = np.zeros(info[f].NPIX)

                    # Calculate the percentage of pixels that receive observations
                    # at at least the required interval
                    jdx1 = np.where(science_cadence_footprint >= req_interval)[0]
                    jdx2 = np.where(cadence_footprint >= req_interval)[0]
                    common_pixels = list(set(jdx1).intersection(set(jdx2)))

                    metric.append((len(common_pixels)/len(jdx1))*100.0)

            results[survey_name][f] = {
                'percent_area_at_interval': np.array(metric),
                'science_cases': cases
            }

    return results

def M7_sky_area_nvisits(survey_config, science_cases):
    """
    Metric to calculate the percentage of the desired survey region to receive the desired number
    of visits in each filter.

    This metric is similar to the M5_proper_motion metric, but it examines the number of visits
    to ensure at least the minimum number of visits are accomplished rather than at minimum intervals.
    This is a proxy for the cadence.

    POSSIBLE EXTENSION: add parameter to describe over what period the observations are achieved to
                give a more direct measure of cadence

    Parameters:
        survey_config   dict   Description of the proposed survey configuration
        science_cases  dict   Catalog of science cases that require multiple epochs
                                of visits per field

    Returns:
        results        dict   Metric value calculated for all science cases

        Output format:
            results = {
                'survey_concept': {
                    'optical_element': {
                        'percent_area_at_cadence': array of metric_values for each science case,
                        'science_case': list of names of the science cases
                    }
                }
            }
    """

    results = {}

    for survey_name, region_set in survey_config.items():

        # Loop over all optical components since the requested footprints can be different
        for k, f in enumerate(OPTICAL_COMPONENTS):
            rsurvey = region_set[f]

            # Create a pixel map of the survey footprint, filling the pixel values with
            # the number of visits per field
            # XXX FOR THIS TO WORK CELESTIALREGIONS need the cadence info; wide-area surveys
            # should have n_visits_per_field = 1
            nvisits_footprint = np.zeros(rsurvey.NPIX)
            idx = np.where(rsurvey.region_map > 0.0)[0]
            if rsurvey.n_visits_per_field >= 2:
                nvisits_footprint[idx].fill(rsurvey.n_visits_per_field)
            else:
                nvisits_footprint[idx].fill(1.0)

            cases = []
            metric = []

            # This assumes one region per science case per filter
            for name, info in science_cases.items():
                if 'multi-epoch' in info['cadence']:
                    cases.append(name)

                    # Similarly for the science cases, create a desired cadence_footprint
                    idx = np.where(info[f].region_map > 0.0)[0]
                    science_nvisits_footprint = np.zeros(info[f].NPIX)
                    if rsurvey.n_visits_per_field >= 2:
                        science_nvisits_footprint[idx].fill(info[f].n_visits_per_field)
                    else:
                        science_nvisits_footprint[idx].fill(1.0)

                    # Calculate the percentage of pixels that receive observations
                    # at at least the required interval.  Here zero or negative values
                    # of the difference map within the desired survey region indicate
                    # that the required number of visits have been achieved
                    diff_map = science_nvisits_footprint - nvisits
                    jdx1 = np.where(diff_map <= 0.0)[0]
                    jdx2 = np.where(science_nvisits_footprint > 0.0)[0]
                    common_pixels = list(set(jdx1).intersection(set(jdx2)))

                    metric.append((len(common_pixels)*PIXAREA/len(jdx2)*PIXAREA)*100.0)

            results[survey_name][f] = {
                'percent_area_at_interval': np.array(metric),
                'science_cases': cases
            }

    return results
