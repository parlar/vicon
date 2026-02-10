#!/bin/bash
# Shared mapping helper: index reference, align reads, add read group, index BAM.
# Usage: map_reads.sh <ref> <r1> <r2> <threads> <sample> <out_bam> <aligner>
set -euo pipefail

REF="$1"; R1="$2"; R2="$3"; THREADS="$4"; SAMPLE="$5"; OUT_BAM="$6"; ALIGNER="$7"

# Index reference if needed
if [ "$ALIGNER" = "bwa-mem2" ] && [ ! -f "${REF}.bwt.2bit.64" ]; then
    bwa-mem2 index "$REF"
fi
[ -f "${REF}.fai" ] || samtools faidx "$REF"

# Align
TMP="$(dirname "$OUT_BAM")/mapped.unsorted.bam"
if [ "$ALIGNER" = "bwa-mem2" ]; then
    bwa-mem2 mem -t "$THREADS" "$REF" "$R1" "$R2" \
        | samtools sort -o "$TMP"
else
    strobealign -t "$THREADS" "$REF" "$R1" "$R2" \
        | samtools sort -o "$TMP"
fi

# Add read group and index
samtools addreplacerg \
    -r "ID:$SAMPLE" -r "SM:$SAMPLE" -r "PL:ILLUMINA" \
    -o "$OUT_BAM" "$TMP"
rm "$TMP"
samtools index "$OUT_BAM"
