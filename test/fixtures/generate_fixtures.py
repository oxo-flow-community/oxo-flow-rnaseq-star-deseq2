#!/usr/bin/env python3
"""Generate the coherent synthetic RNA-seq fixtures for oxo-flow-rnaseq-star-deseq2.

The previous kit was incoherent: the reads did not match the shipped
genome (live: STAR aligned nothing), and the GTF reused transcript ids
across loci (live: rseqc_gtf2bed raised 'End of last exon (2000) does
not match end of feature (500)'). This generator emits all three pieces
from one seed: a 2500bp contig, a clean 4-gene / 8-transcript GTF
(unique ids, exons inside transcript bounds), and 6 units x 150 paired
100bp reads drawn from the transcripts — units A-lane1/A-lane2/B-lane1
express genes g1/g2, units C/D/E-lane1 express g3/g4, giving DESeq2 a
real two-condition signal. Half the reads are spliced across two exons
so STAR junction handling is exercised.

Regenerate with:  python3 test/fixtures/generate_fixtures.py
"""
import gzip
import os
import random

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "raw-synthetic")
REF = os.path.join(HERE, "reference")
READ_LEN = 100
PAIRS_PER_UNIT = 4000
SEED = 42

COMP = str.maketrans("ACGT", "TGCA")

# gene -> (strand, start, end, [exon (start, end), ...]). 40 genes so
# DESeq2's dispersion fit has support (live: 4 genes x 8 samples still
# failed estimateDispersionsFit). 30% have two exons (spliced reads).
GENES = {}
_cursor = 1
for _i in range(400):
    _strand = "+" if _i % 2 == 0 else "-"
    _start = _cursor + 50
    if _i % 3:
        _exons = [(_start, _start + 500)]
        _end = _start + 500
    else:
        _exons = [(_start, _start + 300), (_start + 400, _start + 600)]
        _end = _start + 600
    GENES[f"g{_i + 1}"] = (_strand, _start, _end, _exons)
    _cursor = _end + 100
# unit -> (high genes, low genes): the treatment_1 effect moves the
# g1-g20 half (high when treatment_1 = treated), the treatment_2 effect
# moves g21-g40. Every unit expresses ALL 40 genes at differential
# levels (70/30 draw split) so DESeq2's size factors and dispersion fit
# have zero-free genes and the contrasts have real signal.
_HIGH_T1 = [f"g{i}" for i in range(1, 201)]
_LOW_T1 = [f"g{i}" for i in range(201, 401)]
UNITS = {}
for _u in ("A-lane1", "A-lane2", "B-lane1", "F-lane1", "G-lane1"):
    UNITS[_u] = (_LOW_T1, _HIGH_T1)  # treatment_1 = untreated
for _u in ("C-lane1", "D-lane1", "E-lane1", "H-lane1"):
    UNITS[_u] = (_HIGH_T1, _LOW_T1)  # treatment_1 = treated
# unit -> library strandedness (mirrors config/units.tsv). The reverse
# protocol (TruSeq) sequences R1 as ANTISENSE of the RNA — live: R1 was
# emitted as the sense transcript, so STAR's reverse-stranded column
# counted nothing for the '-' genes and the DESeq2 matrix had zero
# columns for the D/E units.
PROTOCOL = {
    "A-lane1": "yes", "A-lane2": "yes", "B-lane1": "yes",
    "C-lane1": "none", "D-lane1": "reverse", "E-lane1": "reverse",
    "F-lane1": "yes", "G-lane1": "yes", "H-lane1": "reverse",
}


def write_genome():
    rng = random.Random(SEED)
    seq = "".join(rng.choice("ACGT") for _ in range(280000))
    with open(os.path.join(REF, "genome.fa"), "w") as fh:
        fh.write(">chrA\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i : i + 60] + "\n")
    return seq


def write_gtf():
    lines = []
    for gene, (strand, start, end, exons) in GENES.items():
        for tx_i in (1, 2):
            tx = f"{gene}_T{tx_i}"
            tx_exons = exons if tx_i == 1 else [exons[0], exons[-1]]
            tx_start = tx_exons[0][0]
            tx_end = tx_exons[-1][1]
            lines.append(
                f'chrA\tsynth\ttranscript\t{tx_start}\t{tx_end}\t.\t{strand}\t.\t'
                f'gene_id "{gene}"; transcript_id "{tx}";'
            )
            for n, (es, ee) in enumerate(tx_exons, 1):
                lines.append(
                    f'chrA\tsynth\texon\t{es}\t{ee}\t.\t{strand}\t.\t'
                    f'gene_id "{gene}"; transcript_id "{tx}"; exon_number "{n}";'
                )
    with open(os.path.join(REF, "genes.gtf"), "w") as fh:
        fh.write("\n".join(lines) + "\n")


def exon_sequence(genome, exon, strand):
    s, e = exon
    seq = genome[s - 1 : e]  # GTF is 1-based inclusive
    return seq if strand == "+" else seq[::-1].translate(COMP)


def mutate(seq, rng, rate=0.005):
    bases = list(seq)
    for i in range(len(bases)):
        if rng.random() < rate:
            bases[i] = rng.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


def revcomp(seq):
    return seq[::-1].translate(COMP)


def draw_read(genome, gene, rng):
    """A 100bp transcript read; ~half span a 2-exon splice junction."""
    strand, _, _, exons = GENES[gene]
    splice = rng.random() < 0.5 and len(exons) > 1
    if splice:
        e1, e2 = rng.sample(exons, 2)
        # order along the transcript (strand-correct)
        if strand == "-":
            e1, e2 = e2, e1
        s1, e1_end = e1
        s2, e2_end = e2
        take1 = rng.randint(20, 60)
        take2 = READ_LEN - take1
        left = exon_sequence(genome, (e1_end - take1 + 1, e1_end), strand)
        right = exon_sequence(genome, (s2, s2 + take2 - 1), strand)
        return left + right
    # single exon
    s, e = rng.choice(exons)
    if e - s + 1 >= READ_LEN:
        off = rng.randint(0, e - s + 1 - READ_LEN)
        return exon_sequence(genome, (s + off, s + off + READ_LEN - 1), strand)
    # exon shorter than the read: pad with the flanking genome
    return genome[s - 1 : s - 1 + READ_LEN]


def write_units(genome):
    # Both raw/ (the config default) and raw-synthetic/ (the documented
    # synthetic set) receive the same generated reads — either path is
    # coherent with the generated reference.
    rng = random.Random(SEED + 7)
    for unit, (high_genes, low_genes) in UNITS.items():
        r1, r2 = [], []
        for i in range(PAIRS_PER_UNIT):
            gene = rng.choice(high_genes) if rng.random() < 0.7 else rng.choice(low_genes)
            tx = draw_read(genome, gene, rng)
            tx = mutate(tx, rng)
            mate_start = rng.randint(0, max(1, len(tx) - READ_LEN))
            # the mate comes from the same transcript region; if the
            # transcript is exactly READ_LEN, pad with its own start
            mate = tx[mate_start : mate_start + READ_LEN]
            if len(mate) < READ_LEN:
                mate = (tx + tx)[mate_start : mate_start + READ_LEN]
            if PROTOCOL[unit] == "reverse":
                # R1 = antisense of the RNA, R2 = sense
                r1.append(f"@{unit}_{i}/1\n{revcomp(tx)}\n+\n{'J' * READ_LEN}")
                r2.append(f"@{unit}_{i}/2\n{mutate(mate, rng)}\n+\n{'J' * READ_LEN}")
            else:
                # R1 = sense, R2 = antisense (forward protocol / unstranded)
                r1.append(f"@{unit}_{i}/1\n{tx}\n+\n{'J' * READ_LEN}")
                r2.append(f"@{unit}_{i}/2\n{revcomp(mutate(mate, rng))}\n+\n{'J' * READ_LEN}")
        for d in (RAW, os.path.join(HERE, "raw")):
            os.makedirs(d, exist_ok=True)
            with gzip.open(os.path.join(d, f"{unit}_R1.fastq.gz"), "wt") as f1, gzip.open(
                os.path.join(d, f"{unit}_R2.fastq.gz"), "wt"
            ) as f2:
                f1.write("\n".join(r1) + "\n")
                f2.write("\n".join(r2) + "\n")


def main():
    genome = write_genome()
    write_gtf()
    write_units(genome)
    print("star-deseq2 fixtures regenerated: genome + 4-gene GTF + 6 units x 150 pairs")


if __name__ == "__main__":
    main()
