#!/usr/bin/env python3
"""Summarize iVar minor variant output into per-segment statistics.

Reads the TSV produced by ``ivar variants`` and writes a per-segment summary
with variant counts, allele-frequency statistics, and a mixed-infection flag.

Usage::

    python summarize_minor_variants.py <ivar_variants.tsv> <summary.tsv> <seg1,seg2,...>
"""

import csv
import sys
from collections import defaultdict


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: summarize_minor_variants.py <variants.tsv> <summary.tsv> <segments>")

    variants_path = sys.argv[1]
    summary_path = sys.argv[2]
    segments = sys.argv[3].split(",")

    # ── Read iVar variants TSV ────────────────────────────────
    seg_variants = defaultdict(list)  # segment -> list of alt_freq

    with open(variants_path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if row.get("PASS") != "TRUE":
                continue
            region = row.get("REGION", "")
            alt_freq = float(row.get("ALT_FREQ", 0))
            # Match region to segment
            matched_seg = None
            for seg in segments:
                if seg in region:
                    matched_seg = seg
                    break
            if matched_seg is None:
                continue
            seg_variants[matched_seg].append(alt_freq)

    # ── Build per-segment summary ─────────────────────────────
    header = [
        "segment", "variant_count", "mean_alt_freq", "max_alt_freq",
        "sites_40_60_pct", "mixed_infection_flag",
    ]

    rows = []
    for seg in segments:
        freqs = seg_variants.get(seg, [])
        n = len(freqs)
        mean_af = sum(freqs) / n if n else 0.0
        max_af = max(freqs) if n else 0.0
        # Count sites with near-equal allele frequencies (40-60%) — mixed infection signal
        near_equal = sum(1 for f in freqs if 0.4 <= f <= 0.6)
        mixed_flag = "YES" if near_equal >= 10 else "NO"

        rows.append({
            "segment": seg,
            "variant_count": n,
            "mean_alt_freq": f"{mean_af:.4f}",
            "max_alt_freq": f"{max_af:.4f}",
            "sites_40_60_pct": near_equal,
            "mixed_infection_flag": mixed_flag,
        })

    # ── Write summary ─────────────────────────────────────────
    with open(summary_path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in header) + "\n")

    total = sum(r["variant_count"] for r in rows)
    print(f"[INFO] Minor variant summary: {total} variants across {len(segments)} segments → {summary_path}")


if __name__ == "__main__":
    main()
