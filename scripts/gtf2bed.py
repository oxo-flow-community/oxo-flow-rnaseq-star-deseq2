#!/usr/bin/env python3
"""Convert a reference GTF to a BED12 transcript annotation.

Port of gtf2bed.py from snakemake-workflows/rna-seq-star-deseq2 (v3.1.1).
Creates a gffutils database and writes one BED12 line per transcript with
the transcript id in the name column.

Usage: gtf2bed.py <genome.gtf> <annotation.bed> <annotation.db>
"""
import sys
import gffutils

gtf = sys.argv[1]
bed_path = sys.argv[2]
db_path = sys.argv[3]

db = gffutils.create_db(
    gtf,
    dbfn=db_path,
    force=True,
    keep_order=True,
    merge_strategy="merge",
    sort_attribute_values=True,
    disable_infer_genes=True,
    disable_infer_transcripts=True,
)

with open(bed_path, "w") as outfileobj:
    for tx in db.features_of_type("transcript", order_by="start"):
        bed = [s.strip() for s in db.bed12(tx).split("\t")]
        bed[3] = tx.id
        outfileobj.write("{}\n".format("\t".join(bed)))
