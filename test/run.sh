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

echo "PASS"
