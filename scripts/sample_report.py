#!/usr/bin/env python3
"""Generate a comprehensive HTML sample report for the ViCon pipeline.

Consolidates all per-sample results (coverage, depth, consensus quality,
homoplasy analysis, phylogenetic placement) into a single HTML file with
embedded coverage/depth plots.
"""

import argparse
import base64
import csv
import io
import os
import subprocess
import sys
from collections import defaultdict
from datetime import datetime

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ── Utilities ─────────────────────────────────────────────────


def read_tsv(path):
    """Read a TSV file, returning list of dicts. Returns [] if file is missing/empty/comment."""
    rows = []
    if not os.path.exists(path):
        return rows
    with open(path) as f:
        first = f.readline()
        if not first or first.startswith("#"):
            return rows
        reader = csv.DictReader(io.StringIO(first + f.read()), delimiter="\t")
        for row in reader:
            rows.append(row)
    return rows


def fig_to_base64(fig):
    """Convert a matplotlib figure to a base64-encoded PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def html_table(rows, headers=None, highlight_col=None, highlight_val=None):
    """Render a list of dicts as an HTML table."""
    if not rows:
        return "<p><em>No data available.</em></p>"
    if headers is None:
        headers = list(rows[0].keys())
    lines = ['<table>', '<thead><tr>']
    for h in headers:
        lines.append(f'<th>{h}</th>')
    lines.append('</tr></thead><tbody>')
    for row in rows:
        cls = ""
        if highlight_col and highlight_val:
            if str(row.get(highlight_col, "")) == str(highlight_val):
                cls = ' class="highlight"'
        lines.append(f'<tr{cls}>')
        for h in headers:
            lines.append(f'<td>{row.get(h, "")}</td>')
        lines.append('</tr>')
    lines.append('</tbody></table>')
    return "\n".join(lines)


def get_depth(bam_path):
    """Run samtools depth -a and return {contig: [(pos, depth), ...]}."""
    result = subprocess.run(
        ["samtools", "depth", "-a", "-J", bam_path],
        capture_output=True, text=True,
    )
    depth = defaultdict(list)
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        contig, pos, d = parts[0], int(parts[1]), int(parts[2])
        depth[contig].append((pos, d))
    return dict(depth)


# ── Plot functions ────────────────────────────────────────────


def plot_convergence(conv_rows):
    """Plot iterative mapping convergence: mapped reads and completeness per round."""
    if not conv_rows:
        return None

    rounds = [int(r["round"]) for r in conv_rows]
    mapped = [int(r["mapped_reads"]) for r in conv_rows]
    completeness = [float(r["completeness"]) for r in conv_rows]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    ax1.plot(rounds, mapped, "o-", color="#2563eb", linewidth=2, markersize=6)
    ax1.set_xlabel("Iteration")
    ax1.set_ylabel("Mapped reads")
    ax1.set_title("Mapped reads per iteration")
    ax1.grid(True, alpha=0.3)

    ax2.plot(rounds, completeness, "s-", color="#dc2626", linewidth=2, markersize=6)
    ax2.set_xlabel("Iteration")
    ax2.set_ylabel("Completeness (%)")
    ax2.set_title("Genome completeness per iteration")
    ax2.set_ylim(0, 105)
    ax2.grid(True, alpha=0.3)

    fig.tight_layout()
    return fig_to_base64(fig)


_PALETTE = [
    "#2563eb", "#16a34a", "#dc2626", "#9333ea", "#ea580c",
    "#0891b2", "#be185d", "#65a30d", "#6366f1", "#d97706",
]


def plot_depth(depth_data, segments=("L", "M", "S")):
    """Plot per-position depth for each segment as stacked subplots."""
    # Match segments to contigs (contigs may be named differently)
    contig_map = {}
    for contig in depth_data:
        for seg in segments:
            if seg in contig:
                contig_map[seg] = contig
                break

    available = [s for s in segments if s in contig_map]
    if not available:
        return None

    fig, axes = plt.subplots(len(available), 1, figsize=(14, 3.5 * len(available)),
                             squeeze=False)

    colors = {seg: _PALETTE[i % len(_PALETTE)] for i, seg in enumerate(segments)}

    for i, seg in enumerate(available):
        ax = axes[i][0]
        contig = contig_map[seg]
        positions = [p for p, _ in depth_data[contig]]
        depths = [d for _, d in depth_data[contig]]

        ax.fill_between(positions, depths, alpha=0.3, color=colors.get(seg, "#666"))
        ax.plot(positions, depths, linewidth=0.5, color=colors.get(seg, "#666"))

        mean_d = sum(depths) / len(depths) if depths else 0
        ax.axhline(y=mean_d, color="#f59e0b", linestyle="--", linewidth=1,
                    label=f"Mean: {mean_d:.1f}x")

        ax.set_xlabel("Position (bp)")
        ax.set_ylabel("Depth")
        ax.set_title(f"Segment {seg} ({contig}) — {len(positions):,} bp")
        ax.legend(loc="upper right")
        ax.grid(True, alpha=0.2)
        ax.set_xlim(positions[0], positions[-1])

    fig.tight_layout()
    return fig_to_base64(fig)


# ── Report builder ────────────────────────────────────────────


CSS = """
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
       max-width: 1200px; margin: 0 auto; padding: 20px; color: #1f2937;
       background: #f9fafb; }
h1 { color: #111827; border-bottom: 3px solid #2563eb; padding-bottom: 10px; }
h2 { color: #1f2937; margin-top: 40px; border-bottom: 1px solid #d1d5db; padding-bottom: 6px; }
h3 { color: #374151; margin-top: 25px; }
table { border-collapse: collapse; width: 100%; margin: 15px 0; font-size: 13px; }
th { background: #f3f4f6; padding: 8px 12px; text-align: left; border: 1px solid #d1d5db;
     font-weight: 600; white-space: nowrap; }
td { padding: 6px 12px; border: 1px solid #e5e7eb; }
tr:nth-child(even) { background: #f9fafb; }
tr.highlight { background: #dbeafe !important; font-weight: 600; }
img.plot { max-width: 100%; border: 1px solid #d1d5db; border-radius: 4px;
           margin: 10px 0; }
.meta { color: #6b7280; font-size: 14px; margin-bottom: 30px; }
.badge { display: inline-block; padding: 2px 8px; border-radius: 4px;
         font-size: 12px; font-weight: 600; margin: 0 2px; }
.badge-best { background: #dcfce7; color: #166534; }
.badge-warn { background: #fef3c7; color: #92400e; }
.badge-bad  { background: #fee2e2; color: #991b1b; }
.section { background: white; padding: 20px 25px; margin: 20px 0;
           border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
"""


def build_report(args):
    """Build the full HTML report."""
    sections = []

    # ── Section 1: Header ──
    sections.append(f"""
    <h1>ViCon Sample Report: {args.sample}</h1>
    <div class="meta">
        Analysis line: <strong>{args.line}</strong> &nbsp;|&nbsp;
        Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}
    </div>
    """)

    # ── Section 2: Iterative Mapping Convergence ──
    conv_rows = read_tsv(args.convergence)
    sections.append('<div class="section">')
    sections.append("<h2>Iterative Mapping Convergence</h2>")
    if conv_rows:
        img = plot_convergence(conv_rows)
        if img:
            sections.append(f'<img class="plot" src="data:image/png;base64,{img}" />')
        sections.append("<h3>Convergence statistics</h3>")
        sections.append(html_table(conv_rows))
    else:
        sections.append("<p><em>No convergence data available.</em></p>")
    sections.append("</div>")

    # ── Section 3: Coverage & Depth ──
    sections.append('<div class="section">')
    sections.append("<h2>Coverage &amp; Depth</h2>")

    if args.final_bam and os.path.exists(args.final_bam):
        depth_data = get_depth(args.final_bam)
        if depth_data:
            img = plot_depth(depth_data, segments=args.segments)
            if img:
                sections.append("<h3>Per-position depth (final mapping)</h3>")
                sections.append(f'<img class="plot" src="data:image/png;base64,{img}" />')
    else:
        sections.append("<p><em>No BAM file available for depth plotting.</em></p>")

    # Coverage summary table from summary.tsv
    summary_rows = read_tsv(args.summary)
    if summary_rows:
        sections.append("<h3>Coverage summary across pipeline stages</h3>")
        # Extract coverage columns dynamically
        if summary_rows:
            all_cols = list(summary_rows[0].keys())
            # Show segment info + all coverage columns (mapped_reads, mean_depth, completeness)
            cov_cols = ["segment", "segment_length"]
            cov_cols += [c for c in all_cols if any(
                c.endswith(s) for s in ("_mapped_reads", "_mean_depth", "_completeness")
            )]
            cov_display = [{c: r.get(c, "") for c in cov_cols} for r in summary_rows]
            sections.append(html_table(cov_display, cov_cols))
    sections.append("</div>")

    # ── Section 4: Consensus Quality ──
    sections.append('<div class="section">')
    sections.append("<h2>Consensus Quality</h2>")
    if summary_rows:
        all_cols = list(summary_rows[0].keys())
        # N-count columns
        n_cols = ["segment"]
        n_cols += [c for c in all_cols if c.endswith("_N") or c.endswith("_N_frac")]
        n_cols.append("internal_stop")
        n_display = [{c: r.get(c, "") for c in n_cols} for r in summary_rows]
        sections.append(html_table(n_display, n_cols))

        # Flag internal stops
        for row in summary_rows:
            if row.get("internal_stop") == "1":
                sections.append(
                    f'<p><span class="badge badge-bad">WARNING</span> '
                    f'Internal stop codon detected in iVar consensus for segment '
                    f'{row.get("segment", "?")}</p>'
                )
    else:
        sections.append("<p><em>No summary data available.</em></p>")
    sections.append("</div>")

    # ── Section 5: Consensus Ranking & Best Selection ──
    ranking_rows = read_tsv(args.ranking)
    sections.append('<div class="section">')
    sections.append("<h2>Consensus Ranking &amp; Best Selection</h2>")
    if ranking_rows:
        # Identify best methods
        best_methods = {}
        for r in ranking_rows:
            if r.get("rank") == "1":
                seg = r.get("segment", "")
                method = r.get("method", "")
                best_methods[seg] = method
                sections.append(
                    f'<p>Segment <strong>{seg}</strong>: best method = '
                    f'<span class="badge badge-best">{method}</span></p>'
                )

        sections.append("<h3>Full ranking</h3>")
        sections.append(html_table(ranking_rows, highlight_col="rank", highlight_val="1"))
    else:
        sections.append("<p><em>No ranking data available.</em></p>")
    sections.append("</div>")

    # ── Section 6: Homoplasy Analysis ──
    sections.append('<div class="section">')
    sections.append("<h2>Homoplasy Analysis</h2>")

    # Site-level homoplasies
    site_rows = read_tsv(args.site_homoplasy)
    if site_rows:
        # Summary counts per method per segment
        method_seg_counts = defaultdict(lambda: {"convergent": 0, "autapomorphic": 0})
        for r in site_rows:
            key = (r.get("segment", ""), r.get("method", ""))
            method_seg_counts[key][r.get("classification", "")] += 1

        summary_site = []
        for (seg, method), counts in sorted(method_seg_counts.items()):
            summary_site.append({
                "segment": seg, "method": method,
                "convergent_sites": counts["convergent"],
                "autapomorphic_sites": counts["autapomorphic"],
                "total_flagged": counts["convergent"] + counts["autapomorphic"],
            })
        sections.append("<h3>Homoplastic sites per consensus method</h3>")
        sections.append(html_table(summary_site))

        sections.append("<h3>Site-level details</h3>")
        if len(site_rows) > 100:
            sections.append(f"<p><em>Showing first 100 of {len(site_rows)} flagged sites.</em></p>")
            sections.append(html_table(site_rows[:100]))
        else:
            sections.append(html_table(site_rows))
    else:
        sections.append("<p><em>No site-level homoplasy data.</em></p>")

    # Delta score / pairwise distances
    homo_rows = read_tsv(args.homoplasy)
    if homo_rows:
        sections.append("<h3>Delta score (treelikeness)</h3>")
        sections.append(html_table(homo_rows))

    pw_rows = read_tsv(args.pairwise)
    if pw_rows:
        sections.append("<h3>Pairwise distances between consensus methods</h3>")
        if len(pw_rows) > 50:
            sections.append(f"<p><em>Showing first 50 of {len(pw_rows)} pairs.</em></p>")
            sections.append(html_table(pw_rows[:50]))
        else:
            sections.append(html_table(pw_rows))
    sections.append("</div>")

    # ── Section 7: Phylogenetic Placement Summary ──
    sections.append('<div class="section">')
    sections.append("<h2>Phylogenetic Placement Summary</h2>")

    # Read individual placement data from consensus_placements.tsv files
    placement_rows = []
    for tsv_path in args.placements:
        rows = read_tsv(tsv_path)
        for r in rows:
            # Parse the IQ-TREE report for basic metrics
            iq_path = r.get("iqtree_path", "")
            if iq_path and os.path.exists(iq_path):
                iq_info = parse_iqtree_report(iq_path)
                r.update(iq_info)
            placement_rows.append(r)

    if placement_rows:
        display_cols = [
            "segment", "method", "taxon_name",
            "best_model", "log_likelihood", "tree_length",
        ]
        display = [{c: r.get(c, "") for c in display_cols} for r in placement_rows]
        sections.append(html_table(display, display_cols))
    else:
        sections.append("<p><em>No placement data available.</em></p>")
    sections.append("</div>")

    # ── Section 8: Consensus Validation ──
    val_rows = read_tsv(args.validation)
    sections.append('<div class="section">')
    sections.append("<h2>Consensus Validation (Best &amp; Hybrid)</h2>")
    sections.append(
        "<p>Best and hybrid consensus sequences placed back on the backbone tree. "
        "Lower homoplastic site counts indicate better phylogenetic consistency.</p>"
    )
    if val_rows:
        # Group by segment for a clearer display
        seg_groups = defaultdict(list)
        for r in val_rows:
            seg_groups[r.get("segment", "")].append(r)

        for seg in sorted(seg_groups):
            sections.append(f"<h3>Segment {seg}</h3>")
            display_cols = [
                "consensus_type", "taxon", "homoplastic_sites",
                "autapomorphic_sites", "total_changes", "total_variable_sites",
            ]
            display = [{c: r.get(c, "") for c in display_cols} for r in seg_groups[seg]]
            sections.append(html_table(display, display_cols,
                                       highlight_col="consensus_type",
                                       highlight_val="rank1_original"))

            # Add comparison note
            rank1 = [r for r in seg_groups[seg] if r.get("consensus_type") == "rank1_original"]
            best = [r for r in seg_groups[seg] if r.get("consensus_type") == "best"]
            hybrid = [r for r in seg_groups[seg] if r.get("consensus_type") == "hybrid"]

            if rank1 and best:
                try:
                    orig_h = int(rank1[0].get("homoplastic_sites", 0))
                    best_h = int(best[0].get("homoplastic_sites", 0))
                    hybrid_h = int(hybrid[0].get("homoplastic_sites", 0)) if hybrid else None
                    if best_h < orig_h:
                        sections.append(
                            f'<p><span class="badge badge-best">IMPROVED</span> '
                            f"Best consensus reduced homoplastic sites from {orig_h} to {best_h}</p>"
                        )
                    elif best_h == orig_h:
                        sections.append(
                            f'<p><span class="badge badge-warn">SAME</span> '
                            f"Best consensus has same homoplastic sites ({best_h}) as rank-1 method</p>"
                        )
                    else:
                        sections.append(
                            f'<p><span class="badge badge-bad">WORSE</span> '
                            f"Best consensus has more homoplastic sites ({best_h}) than rank-1 ({orig_h})</p>"
                        )
                    if hybrid_h is not None:
                        if hybrid_h < orig_h:
                            sections.append(
                                f'<p><span class="badge badge-best">IMPROVED</span> '
                                f"Hybrid consensus reduced homoplastic sites from {orig_h} to {hybrid_h}</p>"
                            )
                        elif hybrid_h == orig_h:
                            sections.append(
                                f'<p><span class="badge badge-warn">SAME</span> '
                                f"Hybrid consensus has same homoplastic sites ({hybrid_h}) as rank-1 method</p>"
                            )
                except (ValueError, TypeError):
                    pass
    else:
        sections.append("<p><em>No validation data available.</em></p>")
    sections.append("</div>")

    # ── Section 9: Read Support at Disputed Sites ──
    rs_rows = read_tsv(args.read_support)
    sections.append('<div class="section">')
    sections.append("<h2>Read Support at Disputed Sites</h2>")
    sections.append(
        "<p>Allele frequencies from the final BAM at positions flagged as homoplastic. "
        "Verdicts: <strong>supported</strong> (&gt;80% reads agree), "
        "<strong>weak</strong> (50-80%), <strong>contradicted</strong> (&lt;50%).</p>"
    )
    if rs_rows:
        # Summary counts
        verdict_counts = defaultdict(int)
        for r in rs_rows:
            verdict_counts[r.get("read_support_verdict", "unknown")] += 1

        summary_items = []
        for verdict in ("supported", "weak", "contradicted"):
            count = verdict_counts.get(verdict, 0)
            badge_cls = {"supported": "badge-best", "weak": "badge-warn",
                         "contradicted": "badge-bad"}.get(verdict, "")
            summary_items.append(
                f'<span class="badge {badge_cls}">{verdict}: {count}</span>'
            )
        sections.append(f'<p>{"  ".join(summary_items)} '
                        f"(total: {len(rs_rows)} disputed sites)</p>")

        # Per-segment tables
        seg_groups = defaultdict(list)
        for r in rs_rows:
            seg_groups[r.get("segment", "")].append(r)

        display_cols = [
            "segment", "method", "aln_position", "contig_position",
            "sample_state", "ancestral_state", "classification",
            "depth", "A_count", "C_count", "G_count", "T_count",
            "consensus_allele_freq", "read_support_verdict",
        ]
        for seg in sorted(seg_groups):
            sections.append(f"<h3>Segment {seg}</h3>")
            seg_rows = seg_groups[seg]
            if len(seg_rows) > 50:
                sections.append(f"<p><em>Showing first 50 of {len(seg_rows)} sites.</em></p>")
                seg_rows = seg_rows[:50]
            display = [{c: r.get(c, "") for c in display_cols} for r in seg_rows]
            sections.append(html_table(display, display_cols))
    else:
        sections.append("<p><em>No read support data available.</em></p>")
    sections.append("</div>")

    # ── Assemble HTML ──
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ViCon Report: {args.sample} ({args.line})</title>
<style>{CSS}</style>
</head>
<body>
{"".join(sections)}
</body>
</html>"""

    return html


def parse_iqtree_report(path):
    """Extract key metrics from an IQ-TREE .iqtree report file."""
    import re
    info = {}
    try:
        with open(path) as f:
            content = f.read()
        if content.startswith("#"):
            return info
        m = re.search(r"Best-fit model according to BIC:\s+(.+)", content)
        if m:
            info["best_model"] = m.group(1).strip()
        m = re.search(r"Log-likelihood of the tree:\s+([-\d.]+)", content)
        if m:
            info["log_likelihood"] = m.group(1)
        m = re.search(r"Total tree length.*?:\s+([\d.]+)", content)
        if m:
            info["tree_length"] = m.group(1)
    except Exception:
        pass
    return info


# ── CLI ───────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(description="Generate ViCon sample report")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--line", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--convergence", required=True)
    parser.add_argument("--final-bam", required=True)
    parser.add_argument("--homoplasy", required=True)
    parser.add_argument("--pairwise", required=True)
    parser.add_argument("--ranking", required=True)
    parser.add_argument("--site-homoplasy", required=True)
    parser.add_argument("--best-fasta", required=True)
    parser.add_argument("--validation", required=True)
    parser.add_argument("--read-support", required=True)
    parser.add_argument("--placements", nargs="+", required=True)
    parser.add_argument("--segments", nargs="+", default=["L", "M", "S"],
                        help="Segment names (default: L M S)")
    parser.add_argument("--output", required=True)

    args = parser.parse_args()
    html = build_report(args)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as f:
        f.write(html)

    print(f"Report written to {args.output}")


if __name__ == "__main__":
    main()
