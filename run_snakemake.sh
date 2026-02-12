#!/bin/bash
# ViCon — Viral Consensus Pipeline (full run)
# Execute with: ./vicon/run_snakemake.sh

set -euo pipefail

cd "$(dirname "$0")"

# Override config on the command line as needed.
# Example: run only subset samples with 20 cores:
#   --config fastq_dir=../subset_fastq samples="[B1,F1]"

pixi run snakemake \
  --cores 72 \
  --use-singularity \
  --printshellcmds \
  "$@"
