from astropy.table import Table, Column
import numpy as np
import healpy as hp

def M1_survey_footprint(sim_config, science_cases, survey_config):
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
        sim_config    dict   General configuration parameters common to the whole simulation
        science_cases dict   CelestialRegions representing the community-proposed science cases
        survey_config dict   Description of the proposed survey configuration

    Returns:
        results        astropy.table   Metric value calculated for all science cases
    """

    data = []

    # For each science case, loop over all optical components choices requested
    # and compare the HEALpixel maps of the survey footprint with the region map
    # for each science case
    for author, science_strategy in science_cases.items():

        # Loop over all optical components since the requested footprints can be different
        for optic in sim_config['OPTICAL_COMPONENTS']:
            science_regions = science_strategy[optic]

            if len(science_regions) > 0:
                # Loop over all survey strategies calculating the overlap if any in this optic
                for survey_name, survey_definition in survey_config.items():

                    if optic in survey_definition.keys():

                        # Surveys have multiple regions, so we have to calculate the area
                        # summed over all of them
                        survey_regions = survey_definition[optic]

                        # Calculate the number of overlapping pixels between all science regions
                        # and all survey regions for this strategy and optic
                        common_pixels = []
                        total_science_pixels = []
                        pixel_priorities = 0.0
                        total_pixel_priorities = 0.0
                        for rsurvey in survey_regions:
                            for rscience in science_regions:

                                # Calculate the number of overlapping pixels
                                pix_list = list(set(rscience.pixels).intersection(set(rsurvey.pixels)))
                                common_pixels += pix_list
                                total_science_pixels += list(rscience.pixels)

                                # Sum the priority of the overlapping pixels
                                pixel_priorities += rscience.pixel_priority[pix_list].sum()
                                total_pixel_priorities += rscience.pixel_priority.sum()

                        # Calculate metric values
                        if len(common_pixels) > 0.0:
                            m1 = (len(common_pixels) / len(total_science_pixels) * 100.0)
                            m2 = (pixel_priorities / total_pixel_priorities) * 100.0
                        else:
                            m1 = 0.0
                            m2 = 0.0
                        data.append([survey_name, optic, author, m1, m2])

                    # Otherwise this survey strategy doesn't provide any coverage in this
                    # optic, so return a metric value of zero
                    else:
                        data.append([survey_name, optic, author, 0.0, 0.0])

    data = np.array(data)

    # Return a table of the metric results
    results = Table([
        Column(name='Survey_strategy', data=data[:,0], dtype='S30'),
        Column(name='Optic', data=data[:,1], dtype='S5'),
        Column(name='Science_case', data=data[:,2], dtype='S30'),
        Column(name='M1_%pix', data=data[:,3], dtype='f8'),
        Column(name='M1_%priority', data=data[:,4], dtype='f8')
    ])

    return results

def M2_star_counts(sim_config, survey_config, stellar_density_data):
    """
    Metric calculates the number of stars are included in the survey.

    This metric uses star count data derived from the Trilegal model to estimate the
    total number of stars included within each survey footprint for each optical element.

    The metric value return is the total number of stars included.

    Parameters:
        sim_config    dict   General configuration parameters common to the whole simulation
        survey_config   dict   Description of the proposed survey configuration
        stellar_density_data  dict   HEALpixel arrays of the Trilegal stellar density data per filter

    Returns:
        results        astropy.table   Metric value calculated for all science cases
    """

    data = []

    PIXAREA = hp.nside2pixarea(sim_config['NSIDE'], degrees=True)

    for survey_name, survey_definition in survey_config.items():

        for optic in sim_config['OPTICAL_COMPONENTS']:

            # Surveys have multiple regions, so we have to calculate the area
            # summed over all of them
            if optic in survey_definition.keys():
                survey_regions = survey_definition[optic]

                # Sum star counts over all survey regions in the strategy for this optic
                metric = 0.0
                for r in survey_regions:
                    metric += (10**(stellar_density_data[optic][r.pixels]) * PIXAREA).sum()

                data.append([survey_name, optic, metric])
    data = np.array(data)

    # Return a table of the metric results
    results = Table([
        Column(name='Survey_strategy', data=data[:, 0], dtype='S30'),
        Column(name='Optic', data=data[:, 1], dtype='S5'),
        Column(name='M2_nstars', data=data[:, 2], dtype='f8'),
    ])

    return results

def M3_extended_region_count(sim_config, science_cases, survey_config):
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
        sim_config    dict   General configuration parameters common to the whole simulation
        survey_config   dict   Description of the proposed survey configuration
        science_cases  dict   Catalog of target regions represented as science cases indexed
                                by catalog name, including (l,b) coordinates and
                                radial extent

    Returns:
        results        astropy.table   Metric value calculated for all science cases
    """

    data = []

    # Both the survey definition and the science case can include multiple regions for each
    # optical component.  So we need to check for intersections of the HEALpixels for all cases
    for author, science_strategy in science_cases.items():

        for optic in sim_config['OPTICAL_COMPONENTS']:

            if optic in science_strategy.keys() and len(science_strategy[optic]) > 0:

                for survey_name, survey_definition in survey_config.items():

                    if optic in survey_definition.keys():

                        # Calculate the number of HEALpixels that are both within the target region and
                        # the survey footprint for each target region
                        nregions = 0.0
                        for rscience in science_strategy[optic]:
                            in_pixels = []
                            for rsurvey in survey_definition[optic]:
                                in_pixels += list(set(rscience.pixels).intersection(set(rsurvey.pixels)))
                            in_pixels = list(set(in_pixels))

                            if len(in_pixels) >= len(set(rscience.pixels)):
                                nregions += 1.0

                        # Metric value is the percentage of regions where in_pixel >= r_pixels
                        # (due to the HEALpixels providing irregular coverage of the regions)
                        metric = (nregions / float(len(science_strategy[optic])))*100.0
                        data.append([survey_name, optic, author, metric])

                    else:
                        data.append([survey_name, optic, author, 0.0])

    data = np.array(data)

    # Return a table of the metric results
    results = Table([
        Column(name='Survey_strategy', data=data[:, 0], dtype='S30'),
        Column(name='Optic', data=data[:, 1], dtype='S5'),
        Column(name='Science_case', data=data[:, 2], dtype='S30'),
        Column(name='M3_%regions', data=data[:, 3], dtype='f8'),
    ])

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
        for k, f in enumerate(SIM_CONFIG['OPTICAL_COMPONENTS']):
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

def list_pixels_all_regions(region_set):
    """
    Function to return a list of the total unique pixels in a list of regions
    """
    in_pixels = []
    for r in region_set:
        in_pixels += list(set(r.pixels))
    in_pixels = list(set(in_pixels))

    return in_pixels

def M6_sky_area_optical_elements(sim_config, survey_config, filtersets):
    """
    Metric to evaluate the total area of sky to receive observations in each optical element,
    and combinations of the filters, as a proxy for color measurements.
    Note that this does NOT required contemporaneous observations in different filters to derive
    color measurements and is therefore unsuitable to evaluate color measurements of variable objects.

    Parameters:
        sim_config    dict   General configuration parameters common to the whole simulation
        survey_config   dict   Description of the proposed survey configuration
        filtersets      list of tuples  Combinations of filters

    Returns:
        results        astropy.table   Metric values calculated for all survey designs
    """

    PIXAREA = hp.nside2pixarea(sim_config['NSIDE'], degrees=True)

    data = []

    for survey_name, survey_definition in survey_config.items():

        # Calculate the sky area covered in each optical element
        for optic in sim_config['OPTICAL_COMPONENTS']:
            if optic in survey_definition.keys():
                in_pixels = list_pixels_all_regions(survey_definition[optic])
                m1 = len(in_pixels) * PIXAREA
            else:
                m1 = 0.0
            data.append([survey_name, optic, m1, None])

        # For each combination of filters, calculate the sky area covered in both filters
        for filter_combo in filtersets:

            # First check whether all filters in the combination are present in the survey
            # design.  If not, then this metric returns zero
            check = np.array([True if f in survey_definition.keys() and len(survey_definition[f]) > 0 else False for f in filter_combo]).all()

            # If all filters in the set are available, calculate the area of
            # HEALpixels where observations in all filters are present
            if check:
                in_pixels = set(survey_definition[filter_combo[0]][0].pixels)
                in_pixels = [in_pixels.intersection(set(list_pixels_all_regions(survey_definition[f])))
                             for f in filter_combo[1:]]
                m2 = len(in_pixels) * PIXAREA
            else:
                m2 = 0.0

            data.append([survey_name, ','.join(filter_combo), None, m2])

    data = np.array(data)

    # Return a table of the metric results
    results = Table([
        Column(name='Survey_strategy', data=data[:, 0], dtype='S30'),
        Column(name='Optics', data=data[:, 1], dtype='S5'),
        Column(name='M6_%sky_area_single_filter', data=data[:, 2], dtype='f8'),
        Column(name='M6_%sky_area_filter_combo', data=data[:, 3], dtype='f8'),
    ])

    return results

def M7_sky_area_nvisits(sim_config, science_cases, survey_config):
    """
    Metric to calculate the percentage of the desired survey region to receive the desired number
    of visits in each filter.

    This metric is similar to the M5_proper_motion metric, but it examines the number of visits
    to ensure at least the minimum number of visits are accomplished rather than at minimum intervals.
    This is a proxy for the cadence.

    POSSIBLE EXTENSION: add parameter to describe over what period the observations are achieved to
                give a more direct measure of cadence

    Parameters:
        sim_config    dict   General configuration parameters common to the whole simulation
        survey_config   dict   Description of the proposed survey configuration
        science_cases  dict   Catalog of science cases that require multiple epochs
                                of visits per field

    Returns:
        results        astropy.table   Metric value calculated for all science cases
    """

    data = []

    PIXAREA = hp.nside2pixarea(sim_config['NSIDE'], degrees=True)

    # Each science case requests a distinct set of regions for each filter
    for author, science_strategy in science_cases.items():

        # Loop over all optical components since the requested footprints can be different
        optics_requested = [
            optic for optic in sim_config['OPTICAL_COMPONENTS'] if optic in science_strategy.keys()
        ]
        for optic in optics_requested:

            # Note that since a different cadence can be requested for each region and filter,
            # these have to be compared on a region-by-region basis.
            for i, rscience in enumerate(science_strategy[optic]):

                # For all regions in each survey design option, calculate the metric for the
                # current optical component
                for survey_name, survey_definition in survey_config.items():

                    if optic in survey_definition.keys():

                        for j,rsurvey in enumerate(survey_definition[optic]):
                            # Create a pixel map of the overlap between each region requested for the
                            # science and those from the survey footprint.
                            # If the list of common HEALpixels is non-zero,
                            # fill the pixel values with the number of visits per field
                            survey_visits = np.zeros(rsurvey.NPIX)
                            common_pixels = list(set(rscience.pixels).intersection(set(rsurvey.pixels)))

                            if len(common_pixels) > 0:
                                if rsurvey.nvisits is not None and rsurvey.nvisits >= 2:
                                    survey_visits[common_pixels].fill(rsurvey.nvisits)
                                else:
                                    survey_visits[common_pixels].fill(1.0)

                                # Create a similar pixels map of the requested number of visits
                                # over the whole pixel area
                                science_visits = np.zeros(rscience.NPIX)
                                science_visits[rscience.pixels].fill(rscience.nvisits)

                                # Calculate the percentage of pixels that receive observations
                                # at at least the required interval.  Here zero or negative values
                                # of the difference map within the desired survey region indicate
                                # that the required number of visits have been achieved
                                diff_map = science_visits - survey_visits
                                obs_pixels = np.where(diff_map[common_pixels] <= 0.0)[0]

                                metric = (len(obs_pixels) * PIXAREA / len(rscience.pixels) * PIXAREA) * 100.0

                                data.append([survey_name, rsurvey.label, author, rscience.label, optic, metric])

                            # Handle case of no overlap between the science and survey region
                            else:
                                data.append([survey_name, rsurvey.label, author, rscience.label, optic, 0.0])
                    else:
                        data.append([survey_name, rsurvey.label, author, rscience.label, optic, 0.0])
    data = np.array(data)

    # Return a table of the metric results
    results = Table([
        Column(name='Survey_strategy', data=data[:, 0], dtype='S30'),
        Column(name='Survey_region', data=data[:, 1], dtype='S30'),
        Column(name='Science_case', data=data[:, 2], dtype='S30'),
        Column(name='Science_region', data=data[:, 3], dtype='S30'),
        Column(name='Optic', data=data[:, 4], dtype='S5'),
        Column(name='M7_%sky_area_nvisits', data=data[:, 5], dtype='f8'),
    ])

    return results
