from os import path, getcwd
import glob
import argparse
from astropy.table import Table, vstack

def combine_results_tables(args):

    # Identify all results files for the selected metric
    file_list = glob.glob(path.join(args.data_dir, args.metric+'*_results.txt'))
    if len(file_list) == 0:
        raise IOError('No ' + args.metric + ' results found at ' + args.data_dir)
    print('Found ' + str(len(file_list)) + ' files of results for ' + args.metric)

    # Stack all the metric results vertically into one table
    for i,results_file in enumerate(file_list):
        if '_combined_results.txt' not in results_file:
            if i == 0:
                metric_table = Table.read(results_file, format='ascii')
            else:
                results = Table.read(results_file, format='ascii')
                if len(results) > 0:
                    metric_table = vstack([metric_table, results])
    # Output
    output_file = path.join(args.data_dir, args.metric+'_combined_results.txt')
    metric_table.write(output_file, format='ascii', delimiter=' ', overwrite=True)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('metric', help='Name of metric')
    parser.add_argument('data_dir', help='Path to metric results data directory')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = get_args()
    combine_results_tables(args)
