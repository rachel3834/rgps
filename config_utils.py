import argparse
import json
import sys
from os import path
import argparse

def read_config(path_to_config_file):

    if path.isfile(path_to_config_file) == False:
        raise IOError("No config file found at given location: "+path_to_config_file)

    config_file = open(path_to_config_file,'r')

    config_dict = json.load(config_file)

    config_file.close()

    for key, value in config_dict.items():
        try:
            if '[' in value and ']' in value and ',' in value:
                entries = value.replace('[','').replace(']','').split(',')
                l = []
                for e in entries:
                    try:
                        l.append(float(e))
                    except:
                        l.append(e)
                config_dict[key] = l
            elif isinstance(value, dict):
                for key2,value2 in value.items():
                    if value2 in ['True', 'False']:
                        if value2 == 'True':
                            value[key2] = True
                        else:
                            value[key2] = False
                config_dict[key] = value
        except TypeError:
            pass

    return config_dict


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config_file', help='Path to the configuration file')
    args = parser.parse_args()

    config = read_config(args.config_file)

    print(config)
