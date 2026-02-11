#!/usr/bin/env python3
"""Generate a cross-sample HTML overview report for the ViCon pipeline.

Compares results across all samples: coverage heatmaps, method frequency,
homoplasy distributions, and validation summaries.
"""

import argparse
import base64
import csv
import io
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ── Utilities ─────────────────────────────────────────────────


def read_tsv(path):
    """Read a TSV file, returning list of dicts."""
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
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode("ascii")


def html_table(rows, headers=None, highlight_col=None, highlight_val=None):
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


CSS = """
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
       max-width: 1400px; margin: 0 auto; padding: 20px; color: #1f2937;
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
.section { background: white; padding: 20px 25px; margin: 20px 0;
           border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
"""


# ── Plot functions ────────────────────────────────────────────


def plot_heatmap(data, row_labels, col_labels, title, cmap="YlOrRd", vmin=None, vmax=None,
                 fmt=".1f", label=""):
    """Generic heatmap plot."""
    fig, ax = plt.subplots(figsize=(max(6, len(col_labels) * 1.5), max(3, len(row_labels) * 0.5)))
    arr = np.array(data, dtype=float)

    im = ax.imshow(arr, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=10)

    # Annotate cells
    for i in range(len(row_labels)):
        for j in range(len(col_labels)):
            val = arr[i, j]
            if not np.isnan(val):
                text = f"{val:{fmt}}" if fmt else str(int(val))
                color = "white" if val > (arr[~np.isnan(arr)].max() * 0.6) else "black"
                ax.text(j, i, text, ha="center", va="center", fontsize=9, color=color)

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    if label:
        cbar.set_label(label)
    ax.set_title(title, fontsize=13, fontweight="bold")
    fig.tight_layout()
    return fig_to_base64(fig)


def plot_method_frequency(method_counts, title="Best method frequency"):
    """Bar chart of method frequency."""
    if not method_counts:
        return None
    methods = sorted(method_counts.keys(), key=lambda m: -method_counts[m])
    counts = [method_counts[m] for m in methods]

    fig, ax = plt.subplots(figsize=(max(6, len(methods) * 0.8), 4))
    colors = plt.cm.Set2(np.linspace(0, 1, len(methods)))
    ax.bar(methods, counts, color=colors)
    ax.set_ylabel("Times selected as best")
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.set_xticklabels(methods, rotation=45, ha="right")
    for i, c in enumerate(counts):
        ax.text(i, c + 0.1, str(c), ha="center", fontsize=10, fontweight="bold")
    fig.tight_layout()
    return fig_to_base64(fig)


# ── Report builder ────────────────────────────────────────────


def build_report(args):
    sections = []

    # Parse all data
    all_summaries = []
    for path in args.summaries:
        all_summaries.extend(read_tsv(path))

    all_rankings = []
    for path in args.rankings:
        all_rankings.extend(read_tsv(path))

    all_validations = []
    for path in args.validations:
        all_validations.extend(read_tsv(path))

    # Determine unique samples, lines, segments
    samples = sorted(set(r.get("sample", "") for r in all_summaries))
    lines = sorted(set(r.get("line", "") for r in all_summaries))
    segments = sorted(set(r.get("segment", "") for r in all_summaries))

    n_samples = len(samples)
    n_lines = len(lines)

    # ── Header ──
    sections.append(f"""
    <h1>ViCon Cross-Sample Overview</h1>
    <div class="meta">
        {n_samples} samples &times; {n_lines} analysis lines &times; {len(segments)} segments
        &nbsp;|&nbsp; Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}
    </div>
    """)

    # Process per (sample, line) pairs
    sample_line_pairs = sorted(set(
        (r.get("sample", ""), r.get("line", "")) for r in all_summaries
    ))

    # ── Section 1: Coverage heatmap ──
    sections.append('<div class="section">')
    sections.append("<h2>Coverage Overview</h2>")

    for line_name in lines:
        # Build heatmap: samples (rows) × segments (cols) → final completeness
        row_labels = [s for s, l in sample_line_pairs if l == line_name]
        if not row_labels:
            continue

        data = []
        for s in row_labels:
            row_data = []
            for seg in segments:
                match = [r for r in all_summaries
                         if r.get("sample") == s and r.get("line") == line_name
                         and r.get("segment") == seg]
                if match:
                    # Find final completeness column
                    r = match[0]
                    # Try final_completeness first, then ivar_completeness
                    val = r.get("final_completeness", r.get("ivar_completeness", "0"))
                    try:
                        row_data.append(float(val))
                    except (ValueError, TypeError):
                        row_data.append(float("nan"))
                else:
                    row_data.append(float("nan"))
            data.append(row_data)

        if data:
            img = plot_heatmap(data, row_labels, segments,
                               f"Final completeness (%) — {line_name}",
                               cmap="RdYlGn", vmin=0, vmax=100, label="%")
            sections.append(f'<img class="plot" src="data:image/png;base64,{img}" />')
    sections.append("</div>")

    # ── Section 2: Mean depth heatmap ──
    sections.append('<div class="section">')
    sections.append("<h2>Mean Depth Overview</h2>")

    for line_name in lines:
        row_labels = [s for s, l in sample_line_pairs if l == line_name]
        if not row_labels:
            continue

        data = []
        for s in row_labels:
            row_data = []
            for seg in segments:
                match = [r for r in all_summaries
                         if r.get("sample") == s and r.get("line") == line_name
                         and r.get("segment") == seg]
                if match:
                    r = match[0]
                    val = r.get("final_mean_depth", r.get("ivar_mean_depth", "0"))
                    try:
                        row_data.append(float(val))
                    except (ValueError, TypeError):
                        row_data.append(float("nan"))
                else:
                    row_data.append(float("nan"))
            data.append(row_data)

        if data:
            img = plot_heatmap(data, row_labels, segments,
                               f"Final mean depth — {line_name}",
                               cmap="YlOrRd", fmt=".0f", label="Depth")
            sections.append(f'<img class="plot" src="data:image/png;base64,{img}" />')
    sections.append("</div>")

    # ── Section 3: Best method frequency ──
    sections.append('<div class="section">')
    sections.append("<h2>Best Method Frequency</h2>")

    rank1 = [r for r in all_rankings if r.get("rank") == "1"]
    method_counts = Counter(r.get("method", "") for r in rank1)
    if method_counts:
        img = plot_method_frequency(dict(method_counts))
        if img:
            sections.append(f'<img class="plot" src="data:image/png;base64,{img}" />')

    # Per-line breakdown
    for line_name in lines:
        line_rank1 = [r for r in rank1 if r.get("line") == line_name]
        if not line_rank1:
            continue
        sections.append(f"<h3>{line_name} line</h3>")
        overview = []
        for r in sorted(line_rank1, key=lambda x: (x.get("sample", ""), x.get("segment", ""))):
            overview.append({
                "sample": r.get("sample", ""),
                "segment": r.get("segment", ""),
                "best_method": r.get("method", ""),
                "homoplastic_sites": r.get("homoplastic_sites", ""),
                "autapomorphic_sites": r.get("autapomorphic_sites", ""),
                "total_changes": r.get("total_changes", ""),
            })
        sections.append(html_table(overview))
    sections.append("</div>")

    # ── Section 4: Homoplasy heatmap ──
    sections.append('<div class="section">')
    sections.append("<h2>Homoplasy Overview (Best Method)</h2>")

    for line_name in lines:
        row_labels = [s for s, l in sample_line_pairs if l == line_name]
        if not row_labels:
            continue

        data = []
        for s in row_labels:
            row_data = []
            for seg in segments:
                match = [r for r in rank1
                         if r.get("sample") == s and r.get("line") == line_name
                         and r.get("segment") == seg]
                if match:
                    try:
                        row_data.append(int(match[0].get("homoplastic_sites", 0)))
                    except (ValueError, TypeError):
                        row_data.append(float("nan"))
                else:
                    row_data.append(float("nan"))
            data.append(row_data)

        if data:
            img = plot_heatmap(data, row_labels, segments,
                               f"Homoplastic sites (best method) — {line_name}",
                               cmap="YlOrRd", fmt="d", label="Sites")
            sections.append(f'<img class="plot" src="data:image/png;base64,{img}" />')
    sections.append("</div>")

    # ── Section 5: Validation summary ──
    sections.append('<div class="section">')
    sections.append("<h2>Consensus Validation Summary</h2>")
    sections.append("<p>Comparing best/hybrid consensus to original best method.</p>")

    if all_validations:
        val_table = []
        for r in sorted(all_validations,
                        key=lambda x: (x.get("sample", ""), x.get("line", ""),
                                       x.get("segment", ""), x.get("consensus_type", ""))):
            val_table.append({
                "sample": r.get("sample", ""),
                "line": r.get("line", ""),
                "segment": r.get("segment", ""),
                "type": r.get("consensus_type", ""),
                "homoplastic": r.get("homoplastic_sites", ""),
                "autapomorphic": r.get("autapomorphic_sites", ""),
                "total_changes": r.get("total_changes", ""),
            })
        sections.append(html_table(val_table))
    else:
        sections.append("<p><em>No validation data available.</em></p>")
    sections.append("</div>")

    # ── Assemble HTML ──
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ViCon Cross-Sample Overview</title>
<style>{CSS}</style>
</head>
<body>
{"".join(sections)}
</body>
</html>"""

    return html


def main():
    parser = argparse.ArgumentParser(description="Generate ViCon cross-sample overview report")
    parser.add_argument("--summaries", nargs="+", required=True)
    parser.add_argument("--rankings", nargs="+", required=True)
    parser.add_argument("--validations", nargs="+", required=True)
    parser.add_argument("--output", required=True)

    args = parser.parse_args()
    html = build_report(args)

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w") as f:
        f.write(html)

    print(f"Cross-sample overview written to {args.output}")


if __name__ == "__main__":
    main()
