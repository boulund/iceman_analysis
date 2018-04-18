#!/usr/bin/env python
# Plot BBCountUnique histograms in Iceman resistome analysis workflow
# (c) Fredrik Boulund 2017

from sys import argv, exit
import os
import argparse

import matplotlib as mpl
mpl.use("Agg")

import pandas as pd
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("histogram_data")
    parser.add_argument("output_file")
    if len(argv) < 2:
        parser.print_help()
        exit(1)
    return parser.parse_args()


def plot_histogram(histogram_data_file, output_image):
    """Plot a histogram for a sample.
    """
    sample_name = os.path.basename(histogram_data_file).split(".", maxsplit=1)[0]
    df = pd.read_table(histogram_data_file, index_col=0)

    fig, ax  = plt.subplots(figsize=(10,7))
    df["first"].plot(ax=ax)
    ax.set_title("{sample}\nKmer saturation (BBCountUnique)".format(sample=sample_name))
    ax.set_xlabel("Sampled reads")
    ax.set_ylabel("Proportion novel kmers observed")
    ax.set_ylim([0,100])
    fig.savefig(output_image)


if __name__ == "__main__":
    options = parse_args()
    plot_histogram(options.histogram_data, options.output_file)

