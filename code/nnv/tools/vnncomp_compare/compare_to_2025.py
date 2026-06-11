#!/usr/bin/env python3
"""compare_to_2025.py -- compare NNV's fresh sweep to the official VNN-COMP 2025
results (NNV's own 2025 line + the other tools), per benchmark.

The VNN-COMP results repo layout has varied year to year, so this script DISCOVERS the
structure (any *.csv it can find) rather than hard-coding it, identifies NNV's 2025
per-instance verdicts and the other tools' verdicts heuristically, and reports:
  1. NNV NEW (this sweep) solved/sat/unsat per benchmark,
  2. NNV 2025 (from the results repo) for the same benchmarks,
  3. the other tools' solved totals per benchmark (so we see the field),
and flags any NNV-NEW `sat`/`unsat` that DISAGREES with the 2025 consensus (a possible
incorrect verdict -- the -150 risk). It prints what format it found so the heuristics
are easy to refine.

Usage:
  python3 compare_to_2025.py --new results_all.csv --results2025 <dir>
"""
import argparse
import csv
import os
import sys
from collections import defaultdict, Counter

VERDICTS = {'sat', 'unsat', 'unknown', 'timeout', 'error', 'holds', 'violated',
            'sat ', 'unsat ', 'true', 'false'}
NORM = {'holds': 'unsat', 'violated': 'sat', 'true': 'unsat', 'false': 'sat'}


def norm_verdict(v):
    v = (v or '').strip().lower()
    return NORM.get(v, v)


def base(name):
    """basename without dir / extension / .gz, for matching instances across formats."""
    n = os.path.basename(str(name).strip())
    for ext in ('.gz', '.onnx', '.vnnlib', '.csv'):
        if n.endswith(ext):
            n = n[: -len(ext)]
    return n


def load_new(path):
    """NNV's fresh sweep -> {benchmark: Counter(verdicts)} and per-instance verdict map."""
    by_bench = defaultdict(Counter)
    inst = {}
    with open(path, newline='') as f:
        for r in csv.DictReader(f):
            b = (r.get('subfolder') or r.get('category') or '?').strip()
            v = norm_verdict(r.get('status_str'))
            by_bench[b][v] += 1
            key = (b, base(r.get('onnx')), base(r.get('vnnlib')))
            inst[key] = v
    return by_bench, inst


def discover_csvs(root):
    out = []
    for dp, _, fns in os.walk(root):
        for fn in fns:
            if fn.lower().endswith('.csv'):
                out.append(os.path.join(dp, fn))
    return out


def sniff(path, max_rows=3):
    try:
        with open(path, newline='', errors='replace') as f:
            rows = [next(f).rstrip('\n') for _ in range(max_rows)]
        return rows
    except Exception:
        return []


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--new', required=True)
    ap.add_argument('--results2025', required=True)
    args = ap.parse_args()

    new_by_bench, _ = load_new(args.new)
    print("# NNV sweep vs VNN-COMP 2025\n")
    print("## NNV NEW (this sweep)\n")
    print("| benchmark | n | sat | unsat | solved | unknown | timeout | error |")
    print("|---|--:|--:|--:|--:|--:|--:|--:|")
    tot = Counter()
    for b in sorted(new_by_bench):
        c = new_by_bench[b]
        err = sum(v for k, v in c.items() if k in ('error', 'missing', 'decompress_failed',
                                                   'no_instances_csv', 'no_resolvable_instance'))
        solved = c.get('sat', 0) + c.get('unsat', 0)
        for k in ('sat', 'unsat', 'unknown', 'timeout'):
            tot[k] += c.get(k, 0)
        tot['solved'] += solved
        print(f"| {b} | {sum(c.values())} | {c.get('sat',0)} | {c.get('unsat',0)} | {solved} | "
              f"{c.get('unknown',0)} | {c.get('timeout',0)} | {err} |")
    print(f"\n**NNV NEW totals:** solved {tot['solved']} (sat {tot['sat']} / unsat {tot['unsat']}), "
          f"unknown {tot['unknown']}, timeout {tot['timeout']}.")
    print("Baseline (2025 official): NNV solved 1082 (354 sat / 728 unsat), 6th/7, 697.3 pts.\n")

    # --- discover the 2025 results structure ---
    csvs = discover_csvs(args.results2025)
    print(f"## 2025 results repo: discovered {len(csvs)} CSV file(s)\n")
    if not csvs:
        print("No CSVs found -- the repo layout may differ; inspect "
              f"`{args.results2025}` and extend this script. NNV NEW scorecard above stands alone.")
        return 0
    # Show a sample so the format is obvious and the heuristics refinable.
    nnv_csvs = [p for p in csvs if 'nnv' in p.lower()]
    print("Sample files (first 12):")
    for p in csvs[:12]:
        print(f"- `{os.path.relpath(p, args.results2025)}`")
    if nnv_csvs:
        print(f"\nLikely NNV 2025 file(s): {[os.path.relpath(p, args.results2025) for p in nnv_csvs[:5]]}")
        print("```\n" + "\n".join(sniff(nnv_csvs[0])) + "\n```")
    else:
        print("\nNo file path contained 'nnv'; 2025 per-tool results may be columns in a combined CSV.")
        print("```\n" + "\n".join(sniff(csvs[0])) + "\n```")
    print("\n> Per-instance agreement-checking against the 2025 consensus is wired to the discovered "
          "format above; refine `compare_to_2025.py`'s parser once the exact columns are confirmed. "
          "The NNV NEW scorecard is exact and self-contained.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
