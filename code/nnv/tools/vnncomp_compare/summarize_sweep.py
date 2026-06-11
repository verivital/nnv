#!/usr/bin/env python3
"""summarize_sweep.py -- turn a run_all_benchmarks results CSV into a VNN-COMP-style
scorecard (per-category and overall solved counts + a rough score estimate).

Input CSV columns (from run_all_benchmarks.m):
    subfolder,category,onnx,vnnlib,status,status_str,time_s,error_message

VNN-COMP scoring (per the rules): a correct sat/unsat = +10, an INCORRECT verdict =
-150, and timeout/unknown/error = 0. We do NOT have ground truth here, so this reports
SOLVED counts (sat+unsat) and an OPTIMISTIC score (assumes every emitted sat/unsat is
correct -- the real score needs the official checker / 2025 ground truth, which
compare_to_2025.py cross-checks). Use this for a quick directional read.

Usage:  python3 summarize_sweep.py results_all.csv [> summary.md]
"""
import csv
import sys
from collections import defaultdict, Counter


def main():
    if len(sys.argv) < 2:
        print("usage: summarize_sweep.py <results.csv>", file=sys.stderr)
        return 2
    path = sys.argv[1]
    rows = []
    with open(path, newline='') as f:
        for r in csv.DictReader(f):
            rows.append(r)

    by_cat = defaultdict(Counter)
    overall = Counter()
    times = []
    for r in rows:
        st = (r.get('status_str') or '').strip()
        cat = (r.get('subfolder') or r.get('category') or '?').strip()
        by_cat[cat][st] += 1
        overall[st] += 1
        try:
            times.append(float(r.get('time_s') or 0))
        except ValueError:
            pass

    n = len(rows)
    sat = overall.get('sat', 0)
    unsat = overall.get('unsat', 0)
    solved = sat + unsat
    unknown = overall.get('unknown', 0)
    timeout = overall.get('timeout', 0)
    error = sum(v for k, v in overall.items() if k in ('error', 'missing', 'decompress_failed',
                                                       'no_instances_csv', 'no_resolvable_instance'))

    print("# NNV VNN-COMP sweep scorecard\n")
    print(f"- **Instances:** {n}")
    print(f"- **Solved:** {solved}  (sat {sat} / unsat {unsat})")
    print(f"- **Unknown:** {unknown}  |  **Timeout:** {timeout}  |  **Error/missing:** {error}")
    if times:
        print(f"- **Time/instance:** mean {sum(times)/len(times):.1f}s, max {max(times):.1f}s, total {sum(times)/60:.1f} min")
    print(f"- **Optimistic score** (all emitted verdicts assumed correct; +10 each): **{10*solved}** "
          f"(real score needs the official checker -- see comparison)\n")

    print("## Per-benchmark\n")
    print("| benchmark | n | sat | unsat | unknown | timeout | error |")
    print("|---|--:|--:|--:|--:|--:|--:|")
    for cat in sorted(by_cat):
        c = by_cat[cat]
        err = sum(v for k, v in c.items() if k in ('error', 'missing', 'decompress_failed',
                                                   'no_instances_csv', 'no_resolvable_instance'))
        print(f"| {cat} | {sum(c.values())} | {c.get('sat',0)} | {c.get('unsat',0)} | "
              f"{c.get('unknown',0)} | {c.get('timeout',0)} | {err} |")
    return 0


if __name__ == '__main__':
    sys.exit(main())
