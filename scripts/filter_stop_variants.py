#!/usr/bin/env python3
"""
Filter VCF variants that would introduce internal stop codons.
Assumes reading frame starts at position 1 of each contig (ATG start).

Usage: filter_stop_variants.py <ref.fasta> <in.vcf.gz> <out.vcf.gz>
"""
import sys
import gzip

STOP_CODONS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    seqs = {}
    name = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
    if name:
        seqs[name] = "".join(parts).upper()
    return seqs


def internal_stops(seq):
    """Return set of codon start positions with internal stop codons."""
    stops = set()
    last_codon_start = len(seq) - (len(seq) % 3) - 3
    for i in range(0, last_codon_start, 3):
        if seq[i:i+3] in STOP_CODONS:
            stops.add(i)
    return stops


def variant_creates_stop(ref_seq, pos_1based, ref_allele, alt_allele):
    """Check if applying this variant introduces a new internal stop codon."""
    pos_0 = pos_1based - 1
    new_seq = ref_seq[:pos_0] + alt_allele.upper() + ref_seq[pos_0 + len(ref_allele):]

    ref_stops = internal_stops(ref_seq)
    new_stops = internal_stops(new_seq)

    return len(new_stops - ref_stops) > 0


def main():
    ref_path = sys.argv[1]
    vcf_in = sys.argv[2]
    vcf_out = sys.argv[3]

    refs = read_fasta(ref_path)

    filtered = 0
    kept = 0

    opener = gzip.open if vcf_in.endswith(".gz") else open

    with opener(vcf_in, "rt") as fin, open(vcf_out, "w") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue

            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref_allele = fields[3]
            alt_allele = fields[4]

            if chrom not in refs:
                fout.write(line)
                kept += 1
                continue

            if variant_creates_stop(refs[chrom], pos, ref_allele, alt_allele):
                filtered += 1
                print(
                    f"[filter_stops] Filtered {chrom}:{pos} "
                    f"{ref_allele}>{alt_allele} (introduces stop codon)",
                    file=sys.stderr,
                )
            else:
                fout.write(line)
                kept += 1

    print(f"[filter_stops] Kept {kept}, filtered {filtered} variants", file=sys.stderr)


if __name__ == "__main__":
    main()
