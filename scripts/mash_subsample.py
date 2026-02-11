#!/usr/bin/env python3
"""
Subsample a FASTA file to N diverse representatives using Mash distances
and a greedy furthest-point-first (maximin) algorithm.

Usage: mash_subsample.py <input.fa> <output.fa> <max_seqs> [sketch_size] [kmer_size]
"""
import sys
import os
import subprocess
import tempfile


def read_fasta(path):
    """Parse FASTA -> list of (name, full_header, sequence)."""
    records = []
    name = None
    header = None
    parts = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    records.append((name, header, "".join(parts)))
                name = line[1:].split()[0]
                header = line
                parts = []
            else:
                parts.append(line)
    if name:
        records.append((name, header, "".join(parts)))
    return records


def mash_distances(fasta_path, tmpdir, sketch_size, kmer_size):
    """Sketch individual sequences and compute all-vs-all Mash distances."""
    sketch_prefix = os.path.join(tmpdir, "sketch")
    sketch_file = sketch_prefix + ".msh"

    subprocess.run(
        ["mash", "sketch", "-i", "-s", str(sketch_size),
         "-k", str(kmer_size), "-o", sketch_prefix, fasta_path],
        capture_output=True, text=True, check=True,
    )

    result = subprocess.run(
        ["mash", "dist", sketch_file, sketch_file],
        capture_output=True, text=True, check=True,
    )

    distances = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        fields = line.split("\t")
        distances[(fields[0], fields[1])] = float(fields[2])

    return distances


def greedy_furthest_first(names, distances, n):
    """
    Greedy furthest-point-first (maximin) selection.

    Seed: the sequence with the largest total distance to all others.
    Then iteratively add the sequence whose minimum distance to the
    selected set is largest.
    """

    def dist(a, b):
        return distances.get((a, b), distances.get((b, a), 0.0))

    # Seed: most distant from everything
    best_seed = max(names, key=lambda x: sum(dist(x, o) for o in names if o != x))
    selected = [best_seed]
    remaining = set(names) - {best_seed}

    # Track min distance from each remaining seq to the selected set
    min_dist = {name: dist(name, best_seed) for name in remaining}

    while len(selected) < n and remaining:
        best = max(remaining, key=lambda x: min_dist[x])
        selected.append(best)
        remaining.remove(best)
        # Update min distances against newly added sequence
        for name in remaining:
            d = dist(name, best)
            if d < min_dist[name]:
                min_dist[name] = d

    return selected


def main():
    if len(sys.argv) < 4:
        print("Usage: mash_subsample.py <input.fa> <output.fa> <max_seqs> "
              "[sketch_size] [kmer_size]", file=sys.stderr)
        sys.exit(1)

    input_fa = sys.argv[1]
    output_fa = sys.argv[2]
    max_seqs = int(sys.argv[3])
    sketch_size = int(sys.argv[4]) if len(sys.argv) > 4 else 1000
    kmer_size = int(sys.argv[5]) if len(sys.argv) > 5 else 21

    records = read_fasta(input_fa)
    num_seqs = len(records)
    print(f"[mash_subsample] Input: {num_seqs} sequences from {input_fa}",
          file=sys.stderr)

    # Passthrough if already at or below target, or subsampling disabled
    if num_seqs <= max_seqs or max_seqs <= 0:
        print(f"[mash_subsample] No subsampling needed ({num_seqs} <= {max_seqs})",
              file=sys.stderr)
        with open(output_fa, "w") as fout:
            for name, header, seq in records:
                fout.write(header + "\n")
                for i in range(0, len(seq), 80):
                    fout.write(seq[i:i + 80] + "\n")
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        print(f"[mash_subsample] Sketching (k={kmer_size}, s={sketch_size})...",
              file=sys.stderr)
        distances = mash_distances(input_fa, tmpdir, sketch_size, kmer_size)

    # Mash -i uses sequence names from the FASTA headers as identifiers.
    # Map record names to the names mash actually used in its output.
    dist_names = set()
    for a, b in distances:
        dist_names.add(a)
        dist_names.add(b)

    record_names = [r[0] for r in records]
    name_map = {}  # record_name -> mash_name
    for rec_name in record_names:
        for dn in dist_names:
            if rec_name == dn or dn.endswith(rec_name) or rec_name in dn:
                name_map[rec_name] = dn
                break

    # Remap distances to use record names
    reverse_map = {v: k for k, v in name_map.items()}
    remapped = {}
    for (a, b), d in distances.items():
        ra = reverse_map.get(a, a)
        rb = reverse_map.get(b, b)
        remapped[(ra, rb)] = d

    print(f"[mash_subsample] Selecting {max_seqs} diverse representatives...",
          file=sys.stderr)
    selected = greedy_furthest_first(record_names, remapped, max_seqs)
    selected_set = set(selected)

    with open(output_fa, "w") as fout:
        for name, header, seq in records:
            if name in selected_set:
                fout.write(header + "\n")
                for i in range(0, len(seq), 80):
                    fout.write(seq[i:i + 80] + "\n")

    print(f"[mash_subsample] Selected {len(selected)}/{num_seqs} sequences "
          f"-> {output_fa}", file=sys.stderr)


if __name__ == "__main__":
    main()
