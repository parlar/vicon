#!/usr/bin/env python3
"""
Validate CDS integrity of sequences in a FASTA file.

Checks each sequence for:
  1. Starts with ATG (start codon)
  2. Length is divisible by 3 (no frameshifts)
  3. No internal stop codons (TAA, TAG, TGA at codon positions, excluding final codon)

Writes passing sequences to output FASTA and a TSV report of all sequences.

Usage: validate_cds.py <input.fa> <output.fa> <report.tsv>
"""
import sys

STOP_CODONS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    """Parse FASTA -> list of (name, full_header_line, sequence)."""
    records = []
    name = None
    header = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    records.append((name, header, "".join(parts).upper()))
                name = line[1:].split()[0]
                header = line
                parts = []
            else:
                parts.append(line)
    if name:
        records.append((name, header, "".join(parts).upper()))
    return records


def validate_cds(seq):
    """
    Validate a single sequence as a proper CDS.

    Returns (is_valid, reasons) where reasons is a list of failure codes.
    """
    reasons = []

    if len(seq) < 3 or seq[:3] != "ATG":
        reasons.append("no_start_codon")

    if len(seq) % 3 != 0:
        reasons.append("length_not_divisible_by_3")

    n_codons = len(seq) // 3
    if n_codons > 1:
        for i in range(0, (n_codons - 1) * 3, 3):
            codon = seq[i : i + 3]
            if codon in STOP_CODONS:
                reasons.append(f"internal_stop_at_{i}")

    return (len(reasons) == 0, reasons)


def main():
    if len(sys.argv) != 4:
        print(
            "Usage: validate_cds.py <input.fa> <output.fa> <report.tsv>",
            file=sys.stderr,
        )
        sys.exit(1)

    input_fa = sys.argv[1]
    output_fa = sys.argv[2]
    report_tsv = sys.argv[3]

    records = read_fasta(input_fa)
    total = len(records)
    passed = 0
    failed = 0

    print(
        f"[validate_cds] Input: {total} sequences from {input_fa}",
        file=sys.stderr,
    )

    with open(output_fa, "w") as fout, open(report_tsv, "w") as rpt:
        rpt.write("sequence\tlength\tstatus\treasons\n")

        for name, header, seq in records:
            is_valid, reasons = validate_cds(seq)

            if is_valid:
                fout.write(header + "\n")
                for i in range(0, len(seq), 80):
                    fout.write(seq[i : i + 80] + "\n")
                passed += 1
                rpt.write(f"{name}\t{len(seq)}\tPASS\t\n")
            else:
                failed += 1
                reason_str = ";".join(reasons)
                rpt.write(f"{name}\t{len(seq)}\tFAIL\t{reason_str}\n")
                print(
                    f"[validate_cds] FILTERED: {name} "
                    f"(len={len(seq)}, reasons={reason_str})",
                    file=sys.stderr,
                )

    print(
        f"[validate_cds] Result: {passed} passed, {failed} filtered "
        f"out of {total} total",
        file=sys.stderr,
    )

    if failed > 0:
        pct = 100.0 * failed / total
        print(
            f"[validate_cds] WARNING: {pct:.1f}% of sequences were "
            f"filtered ({failed}/{total})",
            file=sys.stderr,
        )

    if passed == 0:
        print(
            "[validate_cds] ERROR: No sequences passed validation! "
            "Check your GenBank reference FASTA.",
            file=sys.stderr,
        )
        sys.exit(1)


if __name__ == "__main__":
    main()
