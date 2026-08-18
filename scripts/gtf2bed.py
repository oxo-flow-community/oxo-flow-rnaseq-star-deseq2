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
        # gffutils emits '-' strand blocks in transcript order
        # (descending genome coordinates); RSeQC computes introns as
        # (block_i.end, block_{i+1}.start) and dies on negative spans
        # (live: 'Count (-501) must be non-negative' in
        # read_distribution's BED.unionBed3). Sort the blocks in
        # ascending genomic order.
        sizes = [int(x) for x in bed[10].rstrip(",").split(",") if x]
        starts = [int(x) for x in bed[11].rstrip(",").split(",") if x]
        order = sorted(range(len(starts)), key=lambda i: starts[i])
        bed[10] = ",".join(str(sizes[i]) for i in order) + ","
        bed[11] = ",".join(str(starts[i]) for i in order) + ","
        outfileobj.write("{}\n".format("\t".join(bed)))
