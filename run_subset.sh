#!/bin/bash
# ViCon — Viral Consensus Pipeline (subset: 3 samples, 20 cores)
# Builds reference sequences automatically via auto_select_refs.
# Execute with: ./vicon/run_subset.sh

set -euo pipefail

cd "$(dirname "$0")"

pixi run snakemake \
  --cores 20 \
  --use-singularity \
  --config \
    outdir=../results_subset \
    'samples=[B1,F1,K1]' \
    auto_select_refs=true \
  --printshellcmds \
  "$@"
