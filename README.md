# oxo-flow-rnaseq-star-deseq2 — RNA-seq: STAR alignment, DESeq2 differential expression and QC

[![CI](https://github.com/oxo-flow-community/oxo-flow-rnaseq-star-deseq2/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-rnaseq-star-deseq2/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

> ★ Verified · ⇄ Official port of [`snakemake-workflows/rna-seq-star-deseq2`](https://github.com/snakemake-workflows/rna-seq-star-deseq2) @ `v3.1.1` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

This workflow takes paired-end RNA-seq reads through the complete
differential-expression analysis: it downloads the Ensembl reference genome
and annotation (GRCh38, release 115), builds a STAR index, trims reads with
fastp, aligns with STAR to produce sorted BAMs and per-gene counts, runs
RSeQC QC plus a MultiQC report, builds the gene count matrix (collapsing
technical replicates), annotates gene symbols via Ensembl biomaRt, and runs
DESeq2 to produce normalized counts, PCA plots and per-contrast
differential-expression results with ashr shrinkage and MA plots. Every tool
is pinned to an exact version in conda environments, so results are
reproducible run after run.

## Installation

### 1. Install oxo-flow

Requires **oxo-flow >= 0.12.0**. Release binary (recommended):

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

Alternatively via conda: `conda install -c bioconda oxo-flow-cli` (note the
bioconda package may lag the release binary; other platform binaries are
available on the [releases page](https://github.com/Traitome/oxo-flow/releases)).

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-rnaseq-star-deseq2.git
cd oxo-flow-rnaseq-star-deseq2
```

### 3. Requirements

- **Input data** — paired-end FASTQ reads, one `_R1.fastq.gz` / `_R2.fastq.gz`
  pair per unit key at `<raw_dir>/<unit-key>_R1.fastq.gz` (e.g.
  `raw/A-lane1_R1.fastq.gz`), declared in `config/units.tsv`; sample
  conditions in `config/samples.tsv`. The repo ships tiny demo fixtures at
  `test/fixtures/raw/` so the dry-run resolves every input; point `raw_dir`
  at your own data (e.g. `raw_dir = "raw"` in `main.oxoflow`).
- **Reference data** — none to download manually: the Ensembl reference
  genome FASTA and annotation GTF (GRCh38, release 115, configurable via
  `ref_species`/`ref_release`/`ref_build`) are fetched automatically at run
  time, and `gene_2_symbol` queries Ensembl biomaRt — **network access
  required** for both.
- **Tools** — conda environments, all exact-pinned (`envs/*.yaml`): STAR
  2.7.11b, fastp 1.0.1, RSeQC 5.0.4, gffutils 0.13, pandas 2.3.2, MultiQC
  1.29, DESeq2 1.46.0 (r-stringr 1.5.1, r-ashr 2.2_63), biomaRt 2.62.0
  (r-tidyverse 2.0.0, r-dbplyr 2.5.0). Requires conda/mamba to create them;
  reference download uses system `curl`.
- **Compute** — up to 24 CPUs per rule (`star_align`); 8 for `fastp_pe`, 4
  for `star_index`, 1 elsewhere. No per-rule memory limits are configured in
  the workflow; budget RAM for STAR indexing and alignment of human data.
- **Disk** — several tens of GB: the GRCh38 STAR index
  (`resources/star_genome`) is ~30 GB, plus BAMs, trimmed reads and conda
  envs.

## Usage

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

Configuration lives in the `[config]` section of `main.oxoflow`:
`samples_file` / `units_file` point at the sample and unit sheets, `raw_dir`
at your FASTQ directory, and `contrasts` / `contrast_variables` /
`contrast_levels` (paired with `diffexp_variables`, `diffexp_base_levels`,
`diffexp_batch_effects`) drive the DESeq2 comparisons (`contrast_exprs` adds
upstream-style string-form contrasts); `pca_variables` sets the PCA groupings
and `pca_activate` / `trimming_activate` toggle the PCA and trimming rules.
The non-default branches are gated flags, all off by default: `single_end`
(single-end units, which provide only `<unit-key>_R1.fastq.gz`),
`sra_accessions` (comma-joined SRA accessions for the `get_sra` download),
`bwa_index_activate` and `genome_faidx_activate` (index rules with no
consumer in the default path). See the fidelity notes below for how these map
to the upstream config.

## Source

Upstream: **[snakemake-workflows/rna-seq-star-deseq2](https://github.com/snakemake-workflows/rna-seq-star-deseq2)**
@ `v3.1.1` (sha
`aa6b17edf3396230165c18709d04cd982bdaaa4c`), MIT license, Copyright (c) 2017
Johannes Köster. Created 2026-08-15; this workflow may lag behind upstream
releases. Upstream attribution and license details in
[NOTICE.md](NOTICE.md).

## Fidelity

Scope: the **default-parameters main execution path** (upstream `rule all`).
Rows cover every upstream rule; "not ported" rows carry a reason. Upstream
rules use snakemake wrappers v7.2.0 (`bio/fastp`, `bio/star/*`,
`bio/multiqc`, `bio/reference/ensembl-*`, `bio/samtools/faidx`,
`bio/bwa/index`, `bio/sra-tools/fasterq-dump`) whose conda pins were carried
over verbatim. Rules that upstream declares but never runs in the default
path (get_sra, fastp_se, bwa_index, genome_faidx) and the alternate
trimming/SE wiring are ported as config-gated rules (see the `when`/flag
notes per row); they appear as `skip` in the default dry-run plan.

| Upstream rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| get_genome | `get_genome` | curl (system) + Ensembl FTP/HTTPS | ensembl-sequence wrapper: primary_assembly URL with toplevel fallback; probe/fallback restructured into shell, HTTPS branch only (upstream also probes FTP) |
| get_annotation | `get_annotation` | curl (system) + Ensembl FTP/HTTPS | ensembl-annotation wrapper, identical URL + `gzip -d` logic |
| star_index | `star_index` | STAR 2.7.11b | star/index wrapper verbatim; tmpdir moved to `.oxo-flow/tmp/star_index` |
| fastp_pe | `fastp_pe` | fastp 1.0.1 | fastp wrapper verbatim (extra + adapters + reads + trimmed + json + html ordering); upstream per-unit `fastp_adapters`/`fastp_extra` columns → global `[config] fastp_adapters`/`fastp_extra` (defaults equal upstream defaults) |
| star_align | `star_align` (+ `star_align_raw`, `star_align_se`, `star_align_se_raw`) | STAR 2.7.11b | star/align wrapper verbatim: `--outSAMtype BAM SortedByCoordinate --quantMode GeneCounts --sjdbGTFfile "<gtf>"` in the upstream extra-string order, `--readFilesCommand gunzip -c`, `--outStd BAM_SortedByCoordinate` to the BAM, `cat` of ReadsPerGene/SJ/Logs out of the tmp prefix. Upstream's one rule takes trimmed or raw, one or two reads per sample (get_fq + units.tsv); engine rules have fixed input patterns, so the port makes the 2×2 matrix explicit: PE-trimmed (default), PE-raw (`trimming_activate = false`), SE-trimmed, SE-raw — each gated so exactly one variant is active |
| get_sra | `get_sra` | sra-tools 3.2.1 | fasterq-dump wrapper verbatim (`-x`, tmpdir via `mktemp -d`, outputs `sra/{accession}_1.fastq` / `_2.fastq`, log `logs/get-sra/{accession}.log`); gated on `sra_accessions` (comma-joined, default empty → skip). Upstream triggers it per unit from the units.tsv `sra` column and auto-feeds the reads into trimming/alignment (get_units_fastqs); the auto-feed needs per-unit input binding the engine cannot express — downloaded reads must be placed per the `raw_dir` naming convention to be consumed (see exclusions) |
| fastp_se | `fastp_se` | fastp 1.0.1 | fastp wrapper verbatim (single-end arg set: `--in1 --out1 --failed_out --json --html`); gated on `single_end && trimming_activate` (default off). Deviation: upstream writes the shared `{sample}-{unit}.json` for both SE and PE; the port names the SE report `{sample}_single.json` because two rules must not share an output file here |
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
| deseq2 | `deseq2` | DESeq2 1.46.0, r-ashr 2.2_63 | deseq2.R logic identical (list-form contrast = vof + level + base_level, ashr `lfcShrink`, `order(padj)`, MA plot); string-form contrasts ported via `contrast_exprs` (semicolon-joined R expressions parallel to `contrasts`, e.g. `list(c('a_vs_b', ...))`, evaluated `eval(parse(text = ...))` verbatim like upstream; entries must use single-quoted R strings and no semicolons); one instance per `contrasts` entry |
| bwa_index | `bwa_index` | bwa 0.7.19 | bwa/index wrapper verbatim: `-b <size/10 MB, clamped to [10, 51200]>M -p resources/genome.fasta` (the wrapper's block-size formula, replicated with `wc -c` + shell arithmetic), outputs `resources/genome.fasta.{amb,ann,bwt,pac,sa}`; gated on `bwa_index_activate` (default off — upstream `rule all` never requests it; snakemake lazy evaluation vs oxo-flow runs every rule) |
| genome_faidx | `genome_faidx` | samtools 1.22 | `samtools faidx` wrapper verbatim → `resources/genome.fasta.fai`; gated on `genome_faidx_activate` (same reasoning as bwa_index) |
| report/ (`report/*.rst`) | not ported | — | Snakemake report artifacts: the `.rst` captions are jinja templates rendered by the sphinx-based `snakemake --report` machinery (`report:` directive + `report()` output annotations); no oxo-flow equivalent |
| trimming.activate = False rewiring | `star_align_raw` / `star_align_se_raw` | STAR 2.7.11b | upstream then feeds raw reads to star_align; ported as explicit variants gated on `!trimming_activate` |
| per-unit fastp_adapters/fastp_extra | not ported | — | upstream `lookup()` reads per-unit TSV columns per wildcard; the engine has no per-wildcard data lookup — pipeline-level `fastp_adapters` / `fastp_adapters_se` / `fastp_extra` config keys instead (defaults equal upstream defaults) |
| SRA auto-feed | not ported | — | upstream `get_units_fastqs()` binds per-unit inputs (fq1/fq2 vs sra accession); the engine cannot bind per-unit inputs, so `get_sra` (ported) is standalone — its reads must be placed per the `raw_dir` naming convention to enter the pipeline |
| edger / kallisto / trimgalore | n/a | — | not present in upstream v3.1.1 (fastp is the trimmer, DESeq2 the DE tool) |

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

## Test

```bash
bash test/run.sh
```

Runs `validate`, `lint`, a dry-run with the default config and a debug check
that expanded commands contain no literal wildcards. Must exit 0 — the CI
badge above runs the same script.

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. This workflow is derived
from snakemake-workflows/rna-seq-star-deseq2, which is distributed under the
MIT license (see [LICENSE.upstream](LICENSE.upstream)); full attribution in
[NOTICE.md](NOTICE.md).
