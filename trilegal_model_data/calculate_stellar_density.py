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
        'Grism_0thOrder'
        'Mact'
    ]

    if not path.isfile(args.file_path):
        raise IOError('Cannot find Trilegal model file ' + args.file_path)

    with open(args.file_path, 'r') as f:
        data = np.loadtxt(f)

    column_list = [Column(name=columns[i], data=data[:,i]) for i in range(0,len(columns),1)]

    return Table(column_list)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file_path', help='Path to input Trilegal file')
    args = parser.parse_args()

    data_table = parse_trilegal_output(args)
    print(data_table)