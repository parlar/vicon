#!/bin/bash
# ViCon — Viral Consensus Pipeline (subset samples, 20 cores)
# Execute with: ./vicon/run_subset.sh

set -euo pipefail

cd "$(dirname "$0")"

sg docker -c "pixi run snakemake \
  --cores 20 \
  --config fastq_dir=../subset_fastq outdir=../results_subset \
  --printshellcmds \
  $*"
