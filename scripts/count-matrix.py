#!/usr/bin/env python3
"""Build the gene count matrix from STAR ReadsPerGene.out.tab files.

Port of count-matrix.py from snakemake-workflows/rna-seq-star-deseq2
(v3.1.1). Reads each tab, picks the column for the unit's strandedness
(none -> 1, yes -> 2, reverse -> 3), names columns by sample, collapses
technical replicates (units of the same sample) by summation, and writes
results/counts/all.tsv.

Usage: count-matrix.py <units.tsv> <output.tsv> <tab...> <log>
"""
import sys
import os

log = sys.argv[-1]
os.makedirs(os.path.dirname(log), exist_ok=True)
sys.stderr = open(log, "w")

import pandas as pd

units_file = sys.argv[1]
output = sys.argv[2]
tabs = sys.argv[3:-1]


def get_column(strandedness):
    if pd.isnull(strandedness) or strandedness == "none":
        return 1  # non stranded protocol
    elif strandedness == "yes":
        return 2  # 3rd column
    elif strandedness == "reverse":
        return 3  # 4th column, usually for Illumina truseq
    else:
        raise ValueError(
            (
                "'strandedness' column should be empty or have the "
                "value 'none', 'yes' or 'reverse', instead has the "
                "value {}"
            ).format(repr(strandedness))
        )


units = pd.read_csv(units_file, sep="\t", dtype=str)
# unit key convention of this port: <sample_name>-<unit_name> (matches the
# upstream results/star/{sample}-{unit}/ directory naming)
units["unit_key"] = units["sample_name"] + "-" + units["unit_name"].astype(str)
key_to_sample = dict(zip(units["unit_key"], units["sample_name"]))
key_to_strand = dict(zip(units["unit_key"], units["strandedness"]))

counts = []
for f in tabs:
    unit_key = os.path.basename(os.path.dirname(f))
    sample = key_to_sample[unit_key]
    strandedness = key_to_strand.get(unit_key, None)
    t = pd.read_table(
        f, index_col=0, usecols=[0, get_column(strandedness)], header=None, skiprows=4
    )
    t.columns = [sample]
    counts.append(t)

matrix = pd.concat(counts, axis=1)
matrix.index.name = "gene"
# collapse technical replicates
matrix = matrix.groupby(matrix.columns, axis=1, sort=False).sum()
matrix.to_csv(output, sep="\t")
