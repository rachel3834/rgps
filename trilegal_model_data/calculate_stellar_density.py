from os import path
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord, Galactic, TETE
from astropy import units as u
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import healpy as hp
from astropy_healpix import HEALPix
import argparse
import glob
import json
import skyproj
import matplotlib.pyplot as plt

NSIDE = 64
NPIX = hp.nside2npix(NSIDE)
HPRADIUS = hp.max_pixrad(NSIDE, degrees=True)
HPAREA = hp.nside2pixarea(NSIDE, degrees = True)

def build_stellar_density_map(args):
    """
    Function to build a HEALpixel map of stellar density based on input from the Trilegal Galactic model v1.6.
    This function reads in a set of Trilegal output files generated for a grid of points across the region
    of interest.  It applies magnitude cuts to simulate those stars that would fall within Roman's
    magnitude range in the filterset offered by its Wide Field Imager.

    The code expects model data output from the following simulation:
    Galactic Model data generated from Trilegal v1.6
    Using the webform at http://stev.oapd.inaf.it/cgi-bin/trilegal

    Default settings used throughout except:
    Limiting magnitude in 4th filter =30mag
    No dust extinction
    Field area = 0.001sq.deg to keep star counts to a level where the output files are manageable.

    The star counts for each pointing are therefore scaled up by x1000.

    :param args: Commandline arguments
    :return: output HEALpixel maps as FITS binary tables
    """

    # Load the grid of Trilegal model output files.
    # Since the header information provided in these files does not include
    # the details of the simulation used to create them, the files are expected to
    # follow a systematic naming convention:
    # l<gal_longitude_degrees>_b<gal_latitude_degrees>.dat with n used to prefix negative values
    model_file_list = glob.glob(path.join(args.model_data_dir, 'l*_b*.dat'))
    if len(model_file_list) == 0:
        raise IOError('No input model data files found at ' + args.model_data_dir)
    model_data = [parse_trilegal_output(f) for f in model_file_list]

    # Mapping of the Roman optical elements to the column names used in Trilegal files
    optical_components = {
        'F087': 'F087',
        'F106': 'F106',
        'F129': 'F129',
        'F158': 'F158',
        'F184': 'F184',
        'F213': 'F213',
        'F146': 'F146',
        'G150': 'Grism_0thOrder',
        'P127': 'SNprism'
    }

    # Calculate the density of stars for each point in the grid, and store the results in
    # a HEALpixel map, for each filter.
    # Stars from the model output are selected based on Roman WFI's saturation and limiting magnitudes
    # These were estimated from https://roman.gsfc.nasa.gov/science/anticipated_performance_tables.html
    density_maps = {}
    for optic, column_name in optical_components.items():
        density_maps[optic] = calc_density_map(args, model_data, optic, column_name)

    # Output stellar density as log10 HEALpix map to JSON file
    output_density_map(args, density_maps)

def calc_density_map(args, model_data, optic, column_name):
    """
    Function calculates the stellar density for a grid of points on the sky

    The number of stars is factored because the models could only be generated for a very
    small stamp, just 0.001sq.deg due to the excessive number of stars.  This factor scales the
    star count to estimate that for a 1sq.deg. stamp.

    :param args: commandline arguments object
    :param model_data: Model data table for the given optic
    :param optic: string indicating the filter/prism/grism currently under consideration
    :return: HEALpix array of the density map for this optic
    """
    ahp = HEALPix(nside=NSIDE, order='ring', frame=TETE())

    # Calculate the stellar density for all regions of the sky included in the grid
    # Find the minimum stellar density of the available grid - this will be
    # used to fill in the unsampled regions of the map
    factor = HPAREA / 0.001
    z = []
    l = []
    b = []
    lplot = []
    min_density = 1e32
    for sgal, data in model_data:
        # Offset to avoid interpolating over wrap over zero deg
        if sgal.l.deg > 180.0:
            l.append(sgal.l.deg-360.0)
        else:
            l.append(sgal.l.deg)
        lplot.append(sgal.l.deg)
        b.append(sgal.b.deg)
        result = np.logical_and(data[column_name] >= 16.0, data[column_name] <= 25.0)
        sidx = (np.where(result)[0]).astype(int)
        nstars = len(sidx) * factor
        z.append(np.log10(nstars))
        min_density = min(min_density, np.log10(nstars))
    z = np.array(z)
    l = np.array(l)
    b = np.array(b)

    # Plot a map of the original discrete set of model stellar densities
    fig, ax = plt.subplots(figsize=(8, 5))
    sp = skyproj.HammerSkyproj(ax=ax, galactic=True, longitude_ticks='symmetric', celestial=True)
    cmap = plt.get_cmap('viridis')
    z = z / z.max()
    for i in range(0, len(l), 1):
        sp.plot(lplot[i], b[i], c=cmap(float(z[i])), marker='o')
    sp.draw_milky_way()
    plt.savefig(path.join(args.output_dir, 'trilegal_' + optic + '_discrete_map.png'))
    plt.close()

    # Interpolate over the discrete data to obtain a smooth density function
    # Note this only generates values within the grid of available datapoints,
    # so it isn't a full sky map
    # X, Y ranges here are in galactic coordinates.
    #idx1 = np.where(l < 180.0)[0]
    #idx2 = np.where(l > 180.0)[0]
    #x1 = np.linspace(0.0, l[idx1].max(), 50)
    #x2 = np.linspace(l[idx2].min(), 360.0, 50)
    #Xrange = np.concatenate((x2, x1))
    Xrange = np.linspace(l.min(), l.max(), 101)
    Yrange = np.linspace(b.min(), b.max(), 51)
    print('Xrange: ', Xrange.min(), Xrange.max(), l.min(), l.max())
    print('Yrange: ', Yrange.min(), Yrange.max(), b.min(), b.max())

    XX, YY = np.meshgrid(Xrange, Yrange)  # 2D grid for interpolation
    interp = LinearNDInterpolator(list(zip(l, b)), z)
    ZZ = interp(XX, YY)

    #fig, ax = plt.subplots(figsize=(8, 5))
    #plt.imshow(ZZ)
    #plt.show()

    # Now we place the interpolated function onto the HEALpixel sky map
    # Note these coordinates are input in galactic coordinates,
    # but s converts to RA, Dec for plotting
    map = np.zeros(NPIX)
    plotx = []
    ploty = []
    plotz = []
    for ix, x in enumerate(Xrange):
        for iy, y in enumerate(Yrange):
            sgal = SkyCoord(l=x, b=y, frame='galactic', unit=(u.deg, u.deg))
            pixels = calc_healpixels_skycoord(sgal, ahp, coord='galactic')

            if not np.isnan(ZZ[iy,ix]):
                map[pixels] = [max(map[p], ZZ[iy,ix]) for p in pixels]
                plotx.append(sgal.l.deg)
                ploty.append(sgal.b.deg)
                plotz.append(ZZ[iy,ix])

    plotx = np.array(plotx)
    ploty = np.array(ploty)
    plotz = np.array(plotz)
    print('Plotx range: l ', plotx.min(), plotx.max())
    print('Ploty range: b', ploty.min(), ploty.max())


    # Plot a map of the interpolated discrete set of model stellar densities
    fig, ax = plt.subplots(figsize=(8, 5))
    sp = skyproj.HammerSkyproj(ax=ax, galactic=True, longitude_ticks='symmetric', celestial=True)
    cmap = plt.get_cmap('viridis')
    plotz = plotz / plotz.max()
    for i in range(0, len(plotx), 1):
        sp.plot(plotx[i], ploty[i], c=cmap(float(plotz[i])), marker='.', markersize=1)
    sp.draw_milky_way()
    plt.savefig(path.join(args.output_dir, 'trilegal_' + optic + '_map.png'))
    plt.close()

    fig2, ax2 = plt.subplots(2, 1, figsize=(8, 5))
    plt.subplots_adjust(top=0.98, bottom=0.4)

    ax2[0].plot(plotx, plotz, 'r.')
    ax2[0].set_xlabel('l [deg]')
    ax2[0].set_ylabel('Log stellar density')


    ax2[1].plot(ploty, plotz, 'r.')
    ax2[1].set_xlabel('b [deg]')
    ax2[1].set_ylabel('Log stellar density')
    plt.savefig(path.join(args.output_dir, 'trilegal_' + optic + '_xyplot.png'))
    plt.close()

    print('Output interpolated density map for ' + optic)
    exit()

    return map

def output_density_map(args, density_maps):
    output = {
        "label": "Trilegal_v1.6_log10_density",
        "nside": NSIDE,
        "healpix_resolution_deg": HPAREA,
        "n_healpix": NPIX,
    }
    for optic, map in density_maps.items():
        output['healpix_map_'+optic] = map.tolist()

    # Serializing json
    json_object = json.dumps(output, indent=4)

    # Writing to sample.json
    with open(path.join(args.output_dir,'trilegal_nir_stellar_density.json'), "w") as f:
        f.write(json_object)

def calc_healpixels_skycoord(s, ahp, coord='galactic'):
    """Function converts a SkyCoord in the Galactic frame to the corrsponding
    set of HEALpixel indices.
    If the radius of the region is smaller than half that of the HEALpixel map
    resolution, then a minimum radius of 1 HEALpixel is imposed

    Input:
    :s:  SkyCoord in Galactic frame of the center of a pointing
    Output:
    :pixels:  list of HEALpixel indices covered by a circle of HPRADIUS
    """

    if coord == 'galactic':
        skycoord = s.transform_to('icrs')
    else:
        skycoord = s
    #phi = np.deg2rad(skycoord.ra.deg)
    #theta = (np.pi / 2.0) - np.deg2rad(skycoord.dec.deg)
    #xyz = hp.ang2vec(theta, phi)
    #pixels = hp.query_disc(NSIDE, xyz, HPRADIUS)

    pixels = ahp.cone_search_skycoord(skycoord, HPRADIUS*u.deg)

    return pixels

def parse_trilegal_output(file_path):
    """
    Function to load the output ASCII file of the Trilegal model for a single pointing.

    Since the header information provided in these files does not include
    the details of the simulation used to create them, the files are expected to
    follow a systematic naming convention:
    l<gal_longitude_degrees>_b<gal_latitude_degrees>.dat with n used to prefix negative values

    :param file_path:
    :return: Table of galactic model output, SkyCoord of center of region used to generate the model data
    """

    # Expected columns in trilegal output for Roman filters:
    columns = [
        'Gc',
        'logAge',
        '[M/H]',
        'm_ini',
        'logL',
        'logTe',
        'logg',
        'm-M0',
        'Av',
        'm2/m1',
        'mbol',
        'F062',
        'F087',
        'F106',
        'F129',
        'F158',
        'F184',
        'F146',
        'F213',
        'SNprism',
        'Grism_1stOrder',
        'Grism_0thOrder',
        'Mact'
    ]

    # Parse the filename to get the galactic coordinates of the center of the region
    # simulated to generate the model data
    entries = str(path.basename(file_path.replace('.dat',''))).split('_')
    coords = [
        int(x.replace('l','').replace('b','').replace('n','-')) for x in entries
    ]
    field_center = SkyCoord(l=coords[0], b=coords[1], frame='galactic', unit=(u.deg, u.deg))

    if not path.isfile(file_path):
        raise IOError('Cannot find Trilegal model file ' + file_path)

    with open(file_path, 'r') as f:
        data = np.loadtxt(f)

    column_list = [Column(name=columns[i], data=data[:,i]) for i in range(0,len(columns),1)]

    return field_center, Table(column_list)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('model_data_dir', help='Path to input Trilegal model data directory')
    parser.add_argument('output_dir', help='Path to output directory')
    args = parser.parse_args()

    build_stellar_density_map(args)
