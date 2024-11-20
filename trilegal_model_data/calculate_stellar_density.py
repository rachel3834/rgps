from os import path
from astropy.table import Table, Column
import numpy as np
import argparse

def parse_trilegal_output(args):
    """
    Function to load the output ASCII file of the Trilegal model for a single pointing.

    :param file_path:
    :return: Table of galactic model output
    """

    if not path.isfile(args.file_path):
        raise IOError('Cannot find Trilegal model file ' + args.file_path)

    with open(args.file_path, 'r') as f:
        data = np.loadtxt(f)

    print(data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_path', help='Path to input Trilegal file')
    args = parser.parse_args()

    parse_trilegal_output(args)