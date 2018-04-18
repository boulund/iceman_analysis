#!/usr/bin/env python3
__doc__ = "Create count matrix for MEGARes mappings"
__author__ = "Fredrik Boulund"
__date__ = "2017"

from sys import argv, exit
import gzip
import os
import argparse
from collections import namedtuple
from itertools import chain

import pandas as pd


def parse_args():

    parser = argparse.ArgumentParser(description=". ".join([__doc__, __author__, __date__]))

    parser.add_argument("COVSTATS", metavar="FILE", nargs="+",
            help="Covstats output from BBMap.")
    parser.add_argument("-a", "--annotations", 
            default="/ctmr/db/MEGARes/latest/megares_annotations_v1.01.csv",
            help="Path to MEGARes annotations CSV file [%(default)s].")
    parser.add_argument("-o", "--output", 
            default="megares_counts.csv",
            help="Output filename [%(default)s].")

    if len(argv) < 2:
        parser.print_help()
        exit(1)

    return parser.parse_args()


def open_file(filename):
    if filename.endswith(".gz"):
        return gzip.open(filename, 'rt')
    else:
        return open(filename)


def parse_covstats(file_handle, index_column=False):
    """Parse mapped read counts from BBMap covstat file.
    """
    line = file_handle.readline()
    if not line.startswith("#"):
        raise Exception("PARSE ERROR: %s: %s", file_handle.filename, line)

    for line in file_handle:
        splitline = line.split()
        if index_column:
            yield splitline[0]
        else:
            yield int(splitline[6]) + int(splitline[7])


def read_megares_annotations(annotations_file):
    """Read MEGARes annotations.
    """
    annotations = pd.read_csv(annotations_file, index_col=0)
    return annotations.to_dict()


def parse_sample_name(filename):
    return filename.split(".")[0]


def merge_columns(covstats_files, index_column, sample_columns):
    """Merge columns into row-wise iterator.
    """

    for row in zip(index_column, *sample_columns):
        yield(row)


def main(covstats_files, megares_annotations_file, output_file):
    """Convert BBMap covstats output to Pandas DataFrame.

    Read and parse counts from BBMap covstats files.
    Concatenate columns into Pandas DataFrame.
    Add MEGARes annotations as DataFrame indexes.
    Write DataFrame to file.
    """

    # Read the index column from the first file,
    # then parse counts from all files.
    index_column_generator = parse_covstats(open_file(covstats_files[0]), index_column=True)
    sample_column_generators = [parse_covstats(open_file(covstats_file)) for covstats_file in covstats_files]
    column_headers = list(chain(["Gene"], [parse_sample_name(covstats_file) for covstats_file in covstats_files]))
    df = pd.DataFrame(merge_columns(covstats_files,
                                    index_column_generator,
                                    sample_column_generators),
                      columns=column_headers)

    # Append the class, group, and mechaninism annotations, 
    # and adjust the indexes
    megares = read_megares_annotations(megares_annotations_file)
    df["class"] = [megares["class"][gene] for gene in df["Gene"]]
    df["group"] = [megares["group"][gene] for gene in df["Gene"]]
    df["mechanism"] = [megares["mechanism"][gene] for gene in df["Gene"]]
    df.set_index(["Gene", "class", "group", "mechanism"], inplace=True)
    print("Dataframe size is:", df.shape)

    # Write DataFrame to file
    df.to_csv(output_file)


if __name__ == "__main__":
    options = parse_args()
    main(options.COVSTATS, options.annotations, options.output)
