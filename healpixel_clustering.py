from astropy_healpix import HEALPix
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, Column
import numpy as np

def find_priority_pixels(combined_region, threshold=4.0):
    """
    Function to identify the HEALpixels above the given pixel priority threshold.
    """

    pixels = np.where(combined_region.pixel_priority >= threshold)[0]

    candidate_regions = {'pixel_set': pixels}

    return candidate_regions

def get_lrange(llist):

    idx1 = np.where(llist > 180.0)[0]
    idx2 = np.where(llist <= 180.0)[0]

    l0 = llist[idx1].min()
    l1 = llist[idx2].max()

    return l0, l1


def cluster_pixels(params):
    # Convert the HEALpixels to a set of Skycoords in the galactic frame
    proj = HEALPix(nside=64, order='ring', frame='icrs')
    pixels = params['pixel_set']
    coords = proj.healpix_to_skycoord(pixels)
    coords = coords.transform_to('galactic')

    clusters = {}
    nc = -1
    max_sep = 2.5 * u.deg
    while len(coords) > 0:
        s = coords[0]

        sep = s.separation(coords)

        jdx = np.where(sep <= max_sep)[0]
        kdx = np.where(sep > max_sep)[0]

        l = coords[jdx].l.deg
        lmin = l.min()
        lmax = l.max()
        #print('RANGE: ', l, coords[jdx].b.deg)

        # Handle case if the region straddles the l-coordinate rollover
        if lmin > 0.0 and lmin < 180.0 and lmax > 180.0:
            l0, l1 = get_lrange(l)
        else:
            l0 = l.min()
            l1 = l.max()
        b0 = coords[jdx].b.deg.min()
        b1 = coords[jdx].b.deg.max()
        #print('CORNERS: ', lmin, lmax, l0, l1, b0, b1)

        # Calculate the centroid as the mid-point of the diagonal between
        # the opposing corners of a box encompasing the all HEALpixels in the region
        corner1 = SkyCoord(l0, b0, frame='galactic', unit=(u.deg, u.deg))
        corner2 = SkyCoord(l1, b1, frame='galactic', unit=(u.deg, u.deg))
        pa = corner1.position_angle(corner2)
        sep = corner1.separation(corner2)
        center = corner1.directional_offset_by(pa, sep / 2)

        cluster = {
            'l_center': center.l.deg,
            'b_center': center.b.deg,
            'l0': l0,
            'l1': l1,
            'b0': b0,
            'b1': b1,
            'pixels': pixels[jdx]
        }
        nc += 1

        clusters[nc] = cluster

        # Remove identified pixels from the coords list
        coords = coords[kdx]
        pixels = pixels[kdx]

    params['clusters'] = clusters

    return params

def identify_clusters(candidate_regions, max_sep=2.5):
    """
    Function to identify individual clusters of HEALpixels

    Parameters:
        candidate_regions dict Pixel clusters dictionary
        max_sep float Minimum angular separation between clusters in degrees
    """

    # Pixel clusters are separated by assuming they have to be a maximum angular separation apart
    max_sep = max_sep * u.deg
    data = np.array([])
    ids = []
    centroids = np.array([])
    for cid, cluster in candidate_regions['clusters'].items():
        s = SkyCoord(cluster['l_center'], cluster['b_center'], frame='galactic', unit=(u.deg, u.deg))

        # Check whether we already have a region at this location
        if len(centroids) > 0:
            cluster_centroids = SkyCoord(centroids[:, 0], centroids[:, 1], frame='galactic', unit=(u.deg, u.deg))
            sep = s.separation(cluster_centroids)

            # If not, then add a new cluster
            if (sep >= max_sep).all():
                ids.append(len(ids))
                data = np.vstack((data, [
                    cluster['l_center'], cluster['b_center'],
                    cluster['l0'], cluster['l1'], cluster['b0'], cluster['b1']
                ]))

                centroids = np.vstack((centroids, [cluster['l_center'], cluster['b_center']]))

            # If the cluster is already known, compare the boundaries and extend if need be to
            # form the cluster superset of pixels.
            else:
                # print('Existing clusters: ', cluster_centroids)
                # print('Candiddate cluster ',s)

                # print('Separations: ',sep, max_sep)

                cid = np.where(sep <= max_sep)[0][0]
                if len(data.shape) == 1:
                    cmatch = data
                else:
                    cmatch = data[cid, :]
                superset = False
                if superset:
                    cmatch[2] = min(cmatch[2], cluster['l0'])
                    cmatch[3] = max(cmatch[3], cluster['l1'])
                    cmatch[4] = min(cmatch[4], cluster['b0'])
                    cmatch[5] = max(cmatch[5], cluster['b1'])
                    cmatch[0] = np.median([cmatch[2], cmatch[3]])
                    cmatch[1] = np.median([cmatch[4], cmatch[5]])
                    if len(data.shape) == 1:
                        data = cmatch
                    else:
                        data[cid, :] = cmatch

        # If we have no clusters yet, simply add a new one
        else:
            ids.append(len(ids))
            data = np.array([
                cluster['l_center'], cluster['b_center'],
                cluster['l0'], cluster['l1'], cluster['b0'], cluster['b1']
            ])

            centroids = np.array([[cluster['l_center'], cluster['b_center']]])

    data = np.array(data)

    # Allow for testing with a single filter
    if len(data.shape) == 1:
        regions_table = Table([
            Column(name='ID', data=ids),
            Column(name='l_center', data=[data[0]]),
            Column(name='b_center', data=[data[1]]),
            Column(name='l0', data=[data[2]]),
            Column(name='l1', data=[data[3]]),
            Column(name='b0', data=[data[4]]),
            Column(name='b1', data=[data[5]])
        ])

    else:
        regions_table = Table([
            Column(name='ID', data=ids),
            Column(name='l_center', data=data[:, 0]),
            Column(name='b_center', data=data[:, 1]),
            Column(name='l0', data=data[:, 2]),
            Column(name='l1', data=data[:, 3]),
            Column(name='b0', data=data[:, 4]),
            Column(name='b1', data=data[:, 5])
        ])

    regions_table.pprint_all()

    #regions_table.write(path.join(root_dir, 'time_domain_science', 'all_tda_candidate_tda_fields_table.txt'),
    #                 format='ascii', overwrite=True)

    return regions_table