#!/usr/bin/env python
# Plot PMD scores in Iceman resistome analysis workflow
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
    parser.add_argument("pmd_deamination")
    parser.add_argument("output_file")
    if len(argv) < 2:
        parser.print_help()
        exit(1)
    return parser.parse_args()


def plot_pmd(pmd, output_image):
    """Plot PMD deamination estimates.
    """
    sample_name = os.path.basename(pmd).split(".", maxsplit=1)[0]
    df = pd.read_table(pmd)
    df = df.iloc[:, :-1]
    df.columns = "CT       CA       CG       CC       GA       GT       GC  GG".split()
    print(df)

    fig, ax  = plt.subplots(figsize=(10,7))
    df.plot(ax=ax)
    ax.set_title("{sample}\nPMDtools deamination estimates".format(sample=sample_name))
    ax.set_xlabel("Position(?)")
    ax.set_ylabel("Deamination")
    fig.savefig(output_image)


if __name__ == "__main__":
    options = parse_args()
    plot_pmd(options.pmd_deamination, options.output_file)

