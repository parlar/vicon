#!/usr/bin/env python3
"""
Extend consensus sequences using soft-clipped reads at contig boundaries.

Parses a BAM file for reads with soft clips near the 5' and 3' ends of each
contig, builds majority-rule consensus extensions, and writes an extended FASTA.

Usage: extend_consensus.py <ref.fasta> <mapped.bam> <out.fasta> [min_clip=20] [min_support=3]
"""
import sys
import subprocess
import re
from collections import Counter


def read_fasta(path):
    """Parse FASTA → (dict of name→seq, ordered list of names)."""
    seqs = {}
    order = []
    name = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts).upper()
                name = line[1:].split()[0]
                order.append(name)
                parts = []
            else:
                parts.append(line)
    if name:
        seqs[name] = "".join(parts).upper()
    return seqs, order


def parse_cigar(cigar_str):
    """Parse CIGAR string into list of (length, operation) tuples."""
    return [(int(length), op) for length, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar_str)]


def get_soft_clips(bam_path, contig, ref_len, boundary_window=50, min_clip_len=20):
    """Extract soft-clipped sequences at contig boundaries.

    Returns (left_clips, right_clips) where each is a list of sequences.
    - left_clips: soft-clipped bases from reads mapping near position 1,
      reversed so index 0 = closest to contig start.
    - right_clips: soft-clipped bases from reads ending near the contig end,
      index 0 = closest to contig end.
    """
    left_clips = []
    right_clips = []

    result = subprocess.run(
        ["samtools", "view", bam_path, contig],
        capture_output=True, text=True,
    )

    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        flag = int(fields[1])
        if flag & 4:  # unmapped
            continue

        pos = int(fields[3])  # 1-based leftmost mapping position
        cigar_str = fields[5]
        seq = fields[9]

        if cigar_str == "*":
            continue

        cigar = parse_cigar(cigar_str)

        # Check for left soft clip (read maps near position 1)
        if cigar[0][1] == "S" and pos <= boundary_window + 1:
            clip_len = cigar[0][0]
            if clip_len >= min_clip_len:
                clipped_seq = seq[:clip_len]
                # Reverse so index 0 = base closest to contig start
                left_clips.append(clipped_seq[::-1])

        # Check for right soft clip (read ends near contig end)
        if cigar[-1][1] == "S":
            clip_len = cigar[-1][0]
            # Calculate aligned end position on reference
            ref_consumed = sum(length for length, op in cigar if op in "MDN=X")
            aln_end = pos - 1 + ref_consumed  # 0-based end
            if aln_end >= ref_len - boundary_window and clip_len >= min_clip_len:
                clipped_seq = seq[len(seq) - clip_len:]
                right_clips.append(clipped_seq)

    return left_clips, right_clips


def build_extension(clips, min_support=3, min_freq=0.6, max_extend=500):
    """Build consensus extension from soft-clipped sequences.

    clips: list of strings, all oriented so index 0 is closest to the
    contig boundary.  Returns consensus extension string.
    """
    if len(clips) < min_support:
        return ""

    extension = []
    max_len = min(max(len(c) for c in clips), max_extend)

    for i in range(max_len):
        bases = [c[i] for c in clips if i < len(c) and c[i] in "ACGT"]
        if len(bases) < min_support:
            break
        counts = Counter(bases)
        best_base, best_count = counts.most_common(1)[0]
        if best_count / len(bases) >= min_freq:
            extension.append(best_base)
        else:
            break

    return "".join(extension)


def main():
    ref_path = sys.argv[1]
    bam_path = sys.argv[2]
    out_path = sys.argv[3]
    min_clip_len = int(sys.argv[4]) if len(sys.argv) > 4 else 20
    min_support = int(sys.argv[5]) if len(sys.argv) > 5 else 3

    seqs, order = read_fasta(ref_path)
    total_extended = 0

    with open(out_path, "w") as fout:
        for name in order:
            seq = seqs[name]
            ref_len = len(seq)

            left_clips, right_clips = get_soft_clips(
                bam_path, name, ref_len, min_clip_len=min_clip_len,
            )

            left_ext = build_extension(left_clips, min_support=min_support)
            right_ext = build_extension(right_clips, min_support=min_support)

            # Left extension is reversed (index 0 = closest to contig start)
            new_seq = left_ext[::-1] + seq + right_ext

            if left_ext or right_ext:
                print(
                    f"[extend] {name}: +{len(left_ext)}bp 5' / +{len(right_ext)}bp 3' "
                    f"({ref_len} -> {len(new_seq)})",
                    file=sys.stderr,
                )
                total_extended += 1

            fout.write(f">{name}\n")
            for i in range(0, len(new_seq), 80):
                fout.write(new_seq[i : i + 80] + "\n")

    print(
        f"[extend] Extended {total_extended}/{len(order)} contigs",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
