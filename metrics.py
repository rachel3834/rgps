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
                        survey_pixels = list_pixels_all_regions(survey_regions)
                        science_pixels = list_pixels_all_regions(science_regions)

                        # Calculate the number of overlapping pixels
                        common_pixels = list(set(science_pixels).intersection(set(survey_pixels)))

                        # Sum the priority of the overlapping pixels
                        pixel_priorities = 0.0
                        total_pixel_priorities = 0.0
                        for rscience in science_regions:
                            pixel_priorities += rscience.pixel_priority[common_pixels].sum()
                            total_pixel_priorities += rscience.pixel_priority.sum()

                        # Calculate metric values
                        if len(common_pixels) > 0.0:
                            m1 = (len(common_pixels) / len(science_pixels) * 100.0)
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
                    metric += (stellar_density_data[optic][r.pixels] * PIXAREA).sum()

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

    # Unlike other metrics, this metric only makes sense when calculated for the set catalogs
    # of known objects.  These are identified in the science_cases with a special flag
    case_list = []
    for author, params in science_cases.items():
        for optic in sim_config['OPTICAL_COMPONENTS']:
            if optic in params.keys() and len(params[optic]) > 0:
                for r in params[optic]:
                    if r.extended_object_catalog and author not in case_list:
                        case_list.append(author)

    # Both the survey definition and the science case can include multiple regions for each
    # optical component.  So we need to check for intersections of the HEALpixels for all cases
    for author in case_list:
        science_strategy = science_cases[author]

        for optic in sim_config['OPTICAL_COMPONENTS']:

            if optic in science_strategy.keys() and len(science_strategy[optic]) > 0:

                for survey_name, survey_definition in survey_config.items():

                    if optic in survey_definition.keys():

                        # Calculate the number of HEALpixels that are both within the target region and
                        # the survey footprint for each target region
                        nregions = 0.0
                        category = 'None'
                        for rscience in science_strategy[optic]:
                            category = rscience.category
                            in_pixels = []
                            for rsurvey in survey_definition[optic]:
                                in_pixels += list(set(rscience.pixels).intersection(set(rsurvey.pixels)))
                            in_pixels = list(set(in_pixels))

                            # Condition only requires some overlap to count the region
                            # it doesn't have to be fully within the survey boundaries
                            #if len(in_pixels) >= len(set(rscience.pixels)):
                            if len(in_pixels) > 0:
                                nregions += 1.0

                        # Metric value is the percentage of regions where in_pixel >= r_pixels
                        # (due to the HEALpixels providing irregular coverage of the regions)
                        metric = (nregions / float(len(science_strategy[optic])))*100.0
                        data.append([survey_name, optic, author, category, metric])
                    else:
                        data.append([survey_name, optic, author, category, 0.0])

    data = np.array(data)

    # Return a table of the metric results
    if len(data) > 0:
        results = Table([
            Column(name='Survey_strategy', data=data[:, 0], dtype='S30'),
            Column(name='Optic', data=data[:, 1], dtype='S5'),
            Column(name='Science_case', data=data[:, 2], dtype='S30'),
            Column(name='Category', data=data[:, 3], dtype='S30'),
            Column(name='M3_%regions', data=data[:, 4], dtype='f8'),
        ])
    else:
        results = Table([
            Column(name='Survey_strategy', data=np.array([]), dtype='S30'),
            Column(name='Optic', data=np.array([]), dtype='S5'),
            Column(name='Science_case', data=np.array([]), dtype='S30'),
            Column(name='Category', data=np.array([]), dtype='S30'),
            Column(name='M3_%regions', data=np.array([]), dtype='f8'),
        ])

    return results

def M4_proper_motion_precision(sim_config, survey_config):
    """
    Proper motions working group proposed the following metric based on the
    measurement uncertainty, sigma, calculated from:
    sigma = 1 mas / T [years] / sqrt(N exposures per epoch)

    The recommended threshold for "detection" is  ~1 mas / yr, since this precision
    makes it possible to measure quantities like Galactic rotation and bulge dynamics.

    This is calculated for all HEALpixels within all region for a given survey design,
    summing visits over all filters chosen. The metrics calculated represent
    the percentage of the survey region that meets the precision threshold

    Possible extension:
    This might be combined with the star counts metric for a figure of merit, e.g.
    N_{bright stars} * (1 / sigma^2)

    Parameters:
        survey_config   dict   Description of the proposed survey configuration
        survey_config   dict   Definition of survey designs

    Returns:
        results     astropy.table   Summary of metric results
    """

    data = []

    NPIX = hp.nside2npix(sim_config['NSIDE'])

    for survey_name, survey_definition in survey_config.items():

        # Some survey designs have return visits only in different filters.
        # So in addition to calculating this metric per filter, we also
        # sum the visits across all optics
        all_visits_map = np.zeros(NPIX)
        all_interval_map = np.zeros(NPIX)
        all_interval_map.fill(730.0)
        all_pixels = np.zeros(NPIX, dtype='int')
        filter_list = []

        for optic in sim_config['OPTICAL_COMPONENTS']:

            if optic in survey_definition.keys():
                for rsurvey in survey_definition[optic]:
                    nvisits_map = np.zeros(NPIX)
                    interval_map = np.zeros(NPIX)
                    interval_map.fill(730.0)    # Start with the maximum possible interval

                    nvisits_map[rsurvey.pixels] += np.array([rsurvey.nvisits] * len(rsurvey.pixels))
                    all_visits_map[rsurvey.pixels] += np.array([rsurvey.nvisits] * len(rsurvey.pixels))
                    all_pixels[rsurvey.pixels] = 1
                    filter_list.append(optic)

                    # Option 1: A single numerical interval of days between observations is
                    # given in the survey definition.
                    if not np.isnan(rsurvey.visit_interval[0]) and len(rsurvey.visit_interval) == 1:
                        interval_map[rsurvey.pixels] = np.minimum(
                            interval_map[rsurvey.pixels],
                            np.array([rsurvey.visit_interval[0]]*len(rsurvey.pixels))
                        )
                        all_interval_map[rsurvey.pixels] = np.minimum(
                            all_interval_map[rsurvey.pixels],
                            np.array([rsurvey.visit_interval[0]]*len(rsurvey.pixels))
                        )

                    # Option 2: A series of numerical intervals are given
                    elif len(rsurvey.visit_interval) > 1:
                        interval_map[rsurvey.pixels] = np.minimum(
                            interval_map[rsurvey.pixels],
                            np.array([np.median(rsurvey.visit_interval)] * len(rsurvey.pixels))
                        )
                        all_interval_map[rsurvey.pixels] = np.minimum(
                            all_interval_map[rsurvey.pixels],
                            np.array([np.median(rsurvey.visit_interval)] * len(rsurvey.pixels))
                        )

                    # Option 3: A list containing a single None value is given,
                    # meaning that there is only a single visit to each field.
                    elif np.isnan(rsurvey.visit_interval[0]) and len(rsurvey.visit_interval) == 1:
                        interval_map[rsurvey.pixels] = np.minimum(
                            interval_map[rsurvey.pixels],
                            np.array([730.0] * len(rsurvey.pixels))
                        )
                        all_interval_map[rsurvey.pixels] = np.minimum(
                            all_interval_map[rsurvey.pixels],
                            np.array([730.0] * len(rsurvey.pixels))
                        )

                    # Substitute NaN for pixels which have the maximal (730.0) day interval between observations
                    jdx = np.where(interval_map >= 730.0)[0]
                    interval_map[jdx] = np.nan

                    # Compute metric as a HEALpixel map, then calculate the percentage of pixels
                    # that meet the 1 mas critiera
                    metric_map = 1.0 / (interval_map / 365.24) / np.sqrt(nvisits_map)   # in mas
                    idx = np.where(metric_map[rsurvey.pixels] <= 1.0)[0]
                    m1 = (len(idx)/len(rsurvey.pixels))*100.0

                    data.append([survey_name, rsurvey.name, optic, m1])

        # Calculate the metric values for all visits summed over all optics and regions
        # Normally if all HEALpixels have visit intervals >730 then effectively there are no revisits,
        # but this is calculated per filter.  If there are multiple filters listed and some pixels
        # have multiple visits but no visit intervals, then it should be assumed that the visits
        # are executed as pairs of filters, one executed at the start of the survey and the other at the end
        # after 2 years.  So we assume an average interval of 1 yr.
        pixel_list = np.where(all_pixels == 1)[0]
        kdx = np.where(all_visits_map > 1.0)[0]
        jdx = np.where(all_interval_map >= 730.0)[0]
        if len(filter_list) > 1 and len(kdx) > 0 and len(jdx) == len(all_interval_map):
            all_interval_map[pixel_list] = 365.24
            jdx = np.where(all_interval_map >= 730.0)[0]
            all_interval_map[jdx] = np.nan
        else:
            all_interval_map[jdx] = np.nan
        all_optic_metric_map = 1.0 / (all_interval_map / 365.24) / np.sqrt(all_visits_map)  # in mas
        idx = np.where(all_optic_metric_map[pixel_list] <= 1.0)[0]
        m1 = (len(idx) / len(pixel_list)) * 100.0
        data.append([survey_name, 'combined_regions', 'all_optics', m1])

    data = np.array(data)

    # Return a table of the metric results
    results = Table([
        Column(name='Survey_strategy', data=data[:, 0], dtype='S30'),
        Column(name='Survey_region', data=data[:, 1], dtype='S30'),
        Column(name='Optic', data=data[:, 2], dtype='S20'),
        Column(name='M4_proper_motion_precision', data=data[:, 3], dtype='f8'),
    ])

    return results

def list_pixels_all_regions(region_set):
    """
    Function to return a list of the total unique pixels in a list of regions
    """
    in_pixels = []
    for r in region_set:
        in_pixels += list(r.pixels)
    in_pixels = list(set(in_pixels))

    return in_pixels

def extract_multiband_science(sim_config, science_cases):
    """
    Function to review the science cases given and identify those cases which request
    multi-filter observations for the same region, and which combinations
    of filters they request.
    """

    multiband_cases = {}
    for author, science_strategy in science_cases.items():
        fset = []
        for f in sim_config['OPTICAL_COMPONENTS']:
            if f in science_strategy.keys() and len(science_strategy[f]) > 0:
                fset.append(f)
        fset = list(set(fset))
        fset.sort()

        if len(fset) > 1:
            multiband_cases[author] = science_strategy
            multiband_cases[author]['filterset'] = fset

    return multiband_cases

def M5_sky_area_optical_elements(sim_config, science_cases, survey_config):
    """
    Metric to evaluate the total area of sky to receive observations in each optical element,
    and combinations of the filters, as a proxy for color measurements.
    Note that this does NOT required contemporaneous observations in different filters to derive
    color measurements and is therefore unsuitable to evaluate color measurements of variable objects.

    Parameters:
        sim_config    dict   General configuration parameters common to the whole simulation
        science_cases dict      Description of the requested science regions
        survey_config   dict   Description of the proposed survey configuration

    Returns:
        results        astropy.table   Metric values calculated for all survey designs
    """

    PIXAREA = hp.nside2pixarea(sim_config['NSIDE'], degrees=True)

    # First establish what combinations of filters were requested by both
    # science cases and the survey design, based on the region data given
    filter_sets = []
    for survey_name, survey_regions in survey_config.items():
        optic_list = [optic for optic in sim_config['OPTICAL_COMPONENTS'] if len(survey_regions[optic]) > 0]
        if not optic_list in filter_sets and len(optic_list) > 1:
            filter_sets.append(optic_list)
    for author, params in science_cases.items():
        optic_list = [optic for optic in sim_config['OPTICAL_COMPONENTS'] if len(params[optic]) > 0]
        if not optic_list in filter_sets and len(optic_list) > 1:
            filter_sets.append(optic_list)

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
        for filter_combo in filter_sets:

            # First check whether all filters in the combination are present in the survey
            # design.  If not, then this metric returns zero
            check = np.array([True if f in survey_definition.keys() and len(survey_definition[f]) > 0 else False for f in filter_combo]).all()

            # If all filters in the set are available, calculate the area of
            # HEALpixels where observations in all filters are present
            if check:
                in_pixels = set(list_pixels_all_regions(survey_definition[filter_combo[0]]))
                for f in filter_combo[1:]:
                    in_pixels = in_pixels.intersection(set(list_pixels_all_regions(survey_definition[f])))
                in_pixels = list(in_pixels)

                m2 = len(in_pixels) * PIXAREA
            else:
                m2 = 0.0

            data.append([survey_name, ','.join(filter_combo), None, m2])

    data = np.array(data)

    # Return a table of the metric results
    results = Table([
        Column(name='Survey_strategy', data=data[:, 0], dtype='S30'),
        Column(name='Optic', data=data[:, 1], dtype='S20'),
        Column(name='M5_sky_area_single_filter', data=data[:, 2], dtype='f8'),
        Column(name='M5_sky_area_filter_combo', data=data[:, 3], dtype='f8'),
    ])

    return results

def M6_sky_area_nvisits(sim_config, science_cases, survey_config):
    """
    Metric to calculate the percentage of the desired survey region to receive the desired number
    of visits in each filter.

    This metric is similar to the M4_proper_motion metric, but it examines the number of visits
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

                    if len(survey_definition[optic]) > 0:

                        # Calculate for time domain regions only
                        region_list = [rsurvey in enumerate(survey_definition[optic] if rsurvey.time_domain]

                        for j,rsurvey in enumerate(region_list):
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
                        data.append([survey_name, survey_name, author, rscience.label, optic, 0.0])

    data = np.array(data)

    # Return a table of the metric results
    results = Table([
        Column(name='Survey_strategy', data=data[:, 0], dtype='S30'),
        Column(name='Survey_region', data=data[:, 1], dtype='S30'),
        Column(name='Science_case', data=data[:, 2], dtype='S40'),
        Column(name='Science_region', data=data[:, 3], dtype='S40'),
        Column(name='Optic', data=data[:, 4], dtype='S5'),
        Column(name='M6_%sky_area_nvisits', data=data[:, 5], dtype='f8'),
    ])

    return results

def M7_multiband_sky_area(sim_config, science_cases, survey_config):
    """
    Metric to evaluate the science regions for which multi-filter observations are requested.
    The sky regions requested are compared with those survey regions covered in multiple filters.

    Parameters:
        sim_config    dict   General configuration parameters common to the whole simulation
        survey_config   dict   Description of the proposed survey configuration
        science_cases dict      Description of the requested science regions
        filtersets      list of tuples  Combinations of filters

    Returns:
        results        astropy.table   Metric values calculated for all survey designs
    """

    PIXAREA = hp.nside2pixarea(sim_config['NSIDE'], degrees=True)

    # Extract from the science cases provided the subset which request observations
    # in multiple passbands.
    # Note that not all of these regions need to overlap.
    multiband_cases = extract_multiband_science(sim_config, science_cases)

    # For each science case which requests more than one filter to observe a given region,
    # calculate the overlap between the requested science region and the area covered by
    # the combined filter coverage.
    data = []
    for author, science_strategy in science_cases.items():

        # Loop over all optical components since the requested footprints can be different
        for optic in sim_config['OPTICAL_COMPONENTS']:
            science_regions = science_strategy[optic]

            if len(science_regions) > 0:
                # Loop over all survey strategies calculating the overlap if any in this optic
                for survey_name, survey_definition in survey_config.items():

                    # For the requested combination of filters, calculate the sky area covered in all filters

                    # First check whether all filters in the combination are present in the survey
                    # design.  If not, then this metric returns zero
                    check = np.array(
                        [True if f in survey_definition.keys() and len(survey_definition[f]) > 0 else False for
                         f in science_strategy['filterset']]).all()

                    # If all filters in the set are available, calculate the area of
                    # HEALpixels where observations in all filters are present
                    if check:
                        in_pixels = set(survey_definition[science_strategy['filterset'][0]][0].pixels)
                        in_pixels = [in_pixels.intersection(set(list_pixels_all_regions(survey_definition[f])))
                                     for f in science_strategy['filterset'][1:]]
                        metric = len(in_pixels) * PIXAREA
                    else:
                        metric = 0.0

                    data.append([survey_name, author, ','.join(science_strategy['filterset']), metric])

    data = np.array(data)

    # Return a table of the metric results
    results = Table([
        Column(name='Survey_strategy', data=data[:, 0], dtype='S30'),
        Column(name='Science_case', data=data[:, 1], dtype='S40'),
        Column(name='Filterset', data=data[:, 2], dtype='S15'),
        Column(name='M7_multiband_sky_area', data=data[:, 3], dtype='f8'),
    ])

    return results
