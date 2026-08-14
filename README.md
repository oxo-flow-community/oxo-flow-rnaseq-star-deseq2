# RNA-seq with STAR and DESeq2 (oxo-flow port)

[![CI](https://github.com/oxo-flow-community/oxo-flow-rnaseq-star-deseq2/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-rnaseq-star-deseq2/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

RNA-seq differential-expression pipeline: downloads the Ensembl reference
genome + annotation, builds a STAR index, trims reads with fastp, aligns with
STAR (sorted BAM + per-gene counts), runs RSeQC QC + MultiQC, builds the gene
count matrix (collapsing technical replicates), annotates gene symbols via
Ensembl biomaRt, and runs DESeq2 (normalized counts, PCA plots, per-contrast
results with ashr shrinkage and MA plots).

## Source

Ported from **[snakemake-workflows/rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2)**,
version `v3.1.1` (MIT, Copyright (c) 2017 Johannes Köster). This port is
maintained independently and **may lag the upstream** — check the `v3.1.1`
tag above (sha `aa6b17edf3396230165c18709d04cd982bdaaa4c`) and the fidelity
table below for the exact ported state. Ported 2026-08-15.

## Fidelity

Scope: the **default-parameters main execution path** (upstream `rule all`).
Rows cover every upstream rule; "not ported" rows carry a reason. Upstream
rules use snakemake wrappers v7.2.0 (`bio/fastp`, `bio/star/*`,
`bio/multiqc`, `bio/reference/ensembl-*`, `bio/samtools/faidx`,
`bio/bwa/index`) whose conda pins were carried over verbatim.

| Upstream rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| get_genome | `get_genome` | curl (system) + Ensembl FTP/HTTPS | ensembl-sequence wrapper: primary_assembly URL with toplevel fallback; probe/fallback restructured into shell, HTTPS branch only (upstream also probes FTP) |
| get_annotation | `get_annotation` | curl (system) + Ensembl FTP/HTTPS | ensembl-annotation wrapper, identical URL + `gzip -d` logic |
| star_index | `star_index` | STAR 2.7.11b | star/index wrapper verbatim; tmpdir moved to `.oxo-flow/tmp/star_index` |
| fastp_pe | `fastp_pe` | fastp 1.0.1 | fastp wrapper verbatim (extra + adapters + reads + trimmed + json + html ordering); upstream per-unit `fastp_adapters`/`fastp_extra` columns → global `[config] fastp_adapters`/`fastp_extra` (defaults equal upstream defaults) |
| star_align | `star_align` | STAR 2.7.11b | star/align wrapper verbatim: `--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile "<gtf>"` in the upstream extra-string order, `--readFilesCommand gunzip -c`, `--outStd BAM_SortedByCoordinate` to the BAM, `cat` of ReadsPerGene/SJ/Logs out of the tmp prefix |
| get_sra | not ported | sra-tools | SRA-accession branch (fasterq-dump); not in the default path |
| fastp_se | not ported | fastp | single-end branch; not in the default path (paired-end fixtures) |
| rseqc_gtf2bed | `rseqc_gtf2bed` | gffutils 0.13 | gtf2bed.py ported to CLI args; `annotation.db` is `temp_output` (= upstream `temp()`) |
| rseqc_junction_annotation | `rseqc_junction_annotation` | RSeQC 5.0.4 | `junction_annotation.py -q 255 -i <bam> -r <bed> -o <prefix>` verbatim |
| rseqc_junction_saturation | `rseqc_junction_saturation` | RSeQC 5.0.4 | `junction_saturation.py -q 255 ...` verbatim |
| rseqc_stat | `rseqc_stat` | RSeQC 5.0.4 | `bam_stat.py -i <bam> > <out> 2> <log>` verbatim |
| rseqc_infer | `rseqc_infer` | RSeQC 5.0.4 | `infer_experiment.py -r <bed> -i <bam> > <out> 2> <log>` verbatim |
| rseqc_innerdis | `rseqc_innerdis` | RSeQC 5.0.4 | `inner_distance.py -r <bed> -i <bam> -o <prefix>` verbatim |
| rseqc_readdis | `rseqc_readdis` | RSeQC 5.0.4 | `read_distribution.py -r <bed> -i <bam>` verbatim |
| rseqc_readdup | `rseqc_readdup` | RSeQC 5.0.4 | `read_duplication.py -i <bam> -o <prefix>` verbatim |
| rseqc_readgc | `rseqc_readgc` | RSeQC 5.0.4 | `read_GC.py -i <bam> -o <prefix>` verbatim |
| multiqc | `multiqc` | MultiQC 1.29 | multiqc wrapper verbatim: parent dirs of all inputs (incl. the junction-annotation log dir), `--no-data-dir --outdir results/qc --filename multiqc_report`. (In upstream's `rule all` — ported, not excluded.) |
| count_matrix | `count_matrix` | pandas 2.3.2 | count-matrix.py logic identical (strandedness column pick 1/2/3, sample naming, `groupby(...).sum()` collapse of technical replicates); unit→(sample, strandedness) mapping read from `config/units.tsv` instead of snakemake params |
| gene_2_symbol | `gene_2_symbol_counts` / `gene_2_symbol_normcounts` / `gene_2_symbol_diffexp` | biomaRt 2.62.0, r-tidyverse 2.0.0 | upstream is one wildcard-generic rule over `{prefix}`; the port makes the three call sites explicit (oxo-flow has no arbitrary `{prefix}` wildcard). `{contrast}` variant scatters per contrast |
| deseq2_init | `deseq2_init` | DESeq2 1.46.0 | deseq2-init.R logic identical (relevel base levels, batch-effect factors, default interaction model, `rowSums>1` filter, normalized counts); config values passed as CLI args |
| pca | `pca` | DESeq2 1.46.0 | plot-pca.R verbatim (`rlog(blind=FALSE)`, `plotPCA(intgroup=variable)`); one instance per `pca_variables` entry; gated by `pca_activate` |
| deseq2 | `deseq2` | DESeq2 1.46.0, r-ashr 2.2_63 | deseq2.R logic identical (list-form contrast = vof + level + base_level, ashr `lfcShrink`, `order(padj)`, MA plot); complex string-form contrasts not ported (not in default config); one instance per `contrasts` entry |
| bwa_index | not ported | bwa | no consumer in the default path (upstream `rule all` never requests it) |
| genome_faidx | not ported | samtools | no consumer in the default path |
| report/ (`report/*.rst`) | not ported | — | Snakemake report artifacts; no oxo-flow equivalent |
| trimming.activate = False rewiring | not ported | — | upstream then feeds raw reads to star_align; the port keeps the default trimmed-read chain |

**Port-level conventions** (config-shape deviations, commands unchanged):
upstream wildcards are `(sample, unit)`; the port fans out over one composite
`{sample}` = `<sample>-<unit>` (e.g. `A-lane1`), so output paths are
byte-identical to upstream (`results/trimmed/A-lane1/A-lane1_R1.fastq.gz`,
`results/star/A-lane1/...`). Nested upstream config (`diffexp.*`, `ref.*`,
`trimming.activate`, `pca.*`) is flattened into flat `[config]` keys with the
same defaults (see `main.oxoflow` header). The upstream `config/samples.tsv`
demo sheet (samples A–E) and `config/units.tsv` (6 units) ship with the port;
raw reads live at `<raw_dir>/<unit-key>_R1.fastq.gz` / `_R2.fastq.gz`
(`[config] raw_dir` defaults to `test/fixtures/raw/`, which contains tiny real
reads so the dry-run resolves every input; point it at your data, e.g.
`raw_dir = "raw"`). Upstream demo-data FASTQ paths (`A.1.fq.gz` etc.) were
renamed to this convention — data-path substitution only.

## Quickstart

```bash
# 1. install oxo-flow (see Requirements)
# 2. prepare data: <raw_dir>/<unit-key>_R1.fastq.gz / _R2.fastq.gz per config/units.tsv
#    (raw_dir defaults to test/fixtures/raw — the tiny committed fixtures;
#    point it at your own data, e.g. raw_dir = "raw", in main.oxoflow)
# 3. preview the plan
oxo-flow dry-run main.oxoflow
# 4. run (downloads the Ensembl reference, then aligns and analyzes)
oxo-flow run main.oxoflow -j 8
# 5. run a subset (e.g. just the final DESeq2 contrast)
oxo-flow run main.oxoflow -t deseq2 --samples first:2
```

Outputs (identical to upstream): `resources/genome.fasta(.gtf)`,
`resources/star_genome/`, `results/trimmed/<unit-key>/…`,
`results/star/<unit-key>/…`, `results/qc/rseqc/…`,
`results/qc/multiqc_report.html`, `results/counts/all[.symbol].tsv`,
`results/deseq2/all.rds`, `results/deseq2/normcounts[.symbol].tsv`,
`results/diffexp/<contrast>.diffexp[.symbol].tsv` (+ MA plots),
`results/pca.<variable>.svg`.

## Requirements

- **oxo-flow ≥ 0.11.0** — install the prebuilt binary:

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/download/v0.11.0/oxo-flow-v0.11.0-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

- Conda users may alternatively `conda install -c bioconda oxo-flow-cli`
  (note: the bioconda package currently lags the release binary at 0.10.2 —
  some 0.11.0 format features may not validate).
- Conda at runtime for the tool environments declared in `main.oxoflow`
  (`envs/*.yaml`, all exact-pinned). Network access is required for the
  reference download (Ensembl), `gene_2_symbol` (Ensembl biomaRt) and conda
  env creation.

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. Upstream attribution in
[NOTICE.md](NOTICE.md).

## Community

https://oxo-flow-community.github.io/
