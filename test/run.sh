#!/usr/bin/env bash
# Acceptance test for the oxo-flow-rnaseq-star-deseq2 port.
# Usage: ./test/run.sh            (uses ./main.oxoflow)
# Local: OXO=/Users/wsx/Documents/GitHub/oxo-community/bin/oxo-flow ./test/run.sh
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with default config"
"$OXO" dry-run main.oxoflow --samples first:1 > /tmp/oxo-dryrun-$$.txt 2>&1
grep -q "would execute" /tmp/oxo-dryrun-$$.txt

echo "==> debug: expanded commands contain no literal {wildcards}"
"$OXO" debug main.oxoflow | grep -q '{sample}' && { echo "unexpanded wildcards in debug output"; exit 1; } || true
"$OXO" debug main.oxoflow | grep -q '{config\.' && { echo "unexpanded config placeholders in debug output"; exit 1; } || true
"$OXO" debug main.oxoflow | grep -q '{contrast}' && { echo "unexpanded contrast scatter in debug output"; exit 1; } || true

echo "==> DRAFT-mode dry-run: SE + SRA + PCA labels + expression contrast activated"
python3 - "$PWD" <<'EOF' > /tmp/oxo-draft-$$.oxoflow
import sys
src = open(sys.argv[1] + "/main.oxoflow").read()
src = src.replace('name = "se_units"\nsamples = []', 'name = "se_units"\nsamples = ["X-lane1"]')
src = src.replace('name = "sra_units"\nsamples = []', 'name = "sra_units"\nsamples = ["Z-SRR999"]')
src = src.replace('pca_labels = ""', 'pca_labels = "batch"')
src = src.replace("contrast_expressions = \"\"",
                  "contrast_expressions = \"c('treatment_1','treated','untreated')\"")
sys.stdout.write(src)
EOF
"$OXO" dry-run /tmp/oxo-draft-$$.oxoflow > /tmp/oxo-draft-dry-$$.txt 2>&1
grep -q "fastp_se_X-lane1" /tmp/oxo-draft-dry-$$.txt || { echo "fastp_se instance missing"; exit 1; }
grep -q "star_align_se_X-lane1" /tmp/oxo-draft-dry-$$.txt || { echo "star_align_se instance missing"; exit 1; }
grep -q "get_sra_Z-SRR999" /tmp/oxo-draft-dry-$$.txt || { echo "get_sra instance missing"; exit 1; }
grep -q "pca_label_batch" /tmp/oxo-draft-dry-$$.txt || { echo "pca_label instance missing"; exit 1; }
grep -q "X-lane1_single.fastq.gz" /tmp/oxo-draft-dry-$$.txt || { echo "single-end trimmed naming missing"; exit 1; }
grep -q "acc=\"\${u##\*-\}\"" /tmp/oxo-draft-dry-$$.txt || { echo "get_sra accession derivation missing"; exit 1; }
grep -q 'Rscript scripts/deseq2.R .*treatment_1.diffexp.tsv' /tmp/oxo-draft-dry-$$.txt || { echo "deseq2 Rscript call missing"; exit 1; }
grep -q "logs/deseq2/treatment_1.diffexp.log" /tmp/oxo-draft-dry-$$.txt || { echo "per-contrast deseq2 log path missing"; exit 1; }

echo "PASS"
