#!/usr/bin/env python3
"""compare_to_2025.py -- compare NNV's fresh sweep to the official VNN-COMP 2025
results, per benchmark AND per instance, and FLAG disagreements (the soundness check).

The 2025 results repo (github.com/VNN-COMP/vnncomp2025_results) stores, per tool, one
CSV per benchmark: <tool>/2025_<benchmark>/results.csv with rows
    category, onnx_path, vnnlib_path, prepare_time, result, run_time
where result is sat / unsat / timeout / error_*. Tools: alpha_beta_crown, neuralsat,
nnenum, pyrat, cora, rover, sobolbox, nnv.

This script:
  1. Reports NNV NEW (this sweep) vs NNV 2025 solved/sat/unsat per benchmark.
  2. Builds a per-instance REFERENCE verdict: alpha-beta-CROWN's verdict when definitive
     (it won 2025 with 0 incorrect; sound + complete), else the majority definitive
     verdict among the other tools.
  3. FLAGS every NNV-NEW sat/unsat that DISAGREES with the reference -- a NNV `sat` where
     the reference is `unsat` is a probable FALSE counterexample (-150 in competition);
     a NNV `unsat` where the reference is `sat` is a probable unsound proof. This is the
     check that surfaced the load_vnnlib (or ...)-conjunction bug.

Usage:
  python3 compare_to_2025.py --new results_all.csv --results2025 <dir>
"""
import argparse
import csv
import os
import sys
from collections import defaultdict, Counter

REF_TOOL = 'alpha_beta_crown'
TOOLS = ['alpha_beta_crown', 'neuralsat', 'nnenum', 'pyrat', 'cora', 'rover', 'sobolbox']
DEFINITIVE = {'sat', 'unsat'}


def base(name):
    n = os.path.basename(str(name).strip())
    for ext in ('.gz', '.onnx', '.vnnlib'):
        if n.endswith(ext):
            n = n[: -len(ext)]
    return n


def norm(v):
    v = (v or '').strip().lower()
    if v.startswith('error') or v in ('exception', 'crash'):
        return 'error'
    if v in ('holds',):
        return 'unsat'
    if v in ('violated',):
        return 'sat'
    return v


def load_new(path):
    by_bench = defaultdict(Counter)
    inst = {}
    with open(path, newline='') as f:
        for r in csv.DictReader(f):
            b = (r.get('subfolder') or r.get('category') or '?').strip()
            v = norm(r.get('status_str'))
            by_bench[b][v] += 1
            inst[(b, base(r.get('onnx')), base(r.get('vnnlib')))] = v
    return by_bench, inst


def load_tool_bench(results_dir, tool, bench):
    """{(onnx_base, vnnlib_base): verdict} from <tool>/2025_<bench>/results.csv (or {})."""
    p = os.path.join(results_dir, tool, '2025_' + bench, 'results.csv')
    out = {}
    if not os.path.isfile(p):
        return out
    with open(p, newline='', errors='replace') as f:
        for row in csv.reader(f):
            if len(row) >= 5:
                out[(base(row[1]), base(row[2]))] = norm(row[4])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--new', required=True)
    ap.add_argument('--results2025', required=True)
    args = ap.parse_args()

    new_by_bench, new_inst = load_new(args.new)
    rd = args.results2025

    print("# NNV sweep vs VNN-COMP 2025\n")

    # --- per-benchmark NNV-new vs NNV-2025 ---
    print("## Per-benchmark: NNV NEW vs NNV 2025\n")
    print("| benchmark | NNV-new solved (sat/unsat) | NNV-2025 solved (sat/unsat) | NNV-new unk/to/err |")
    print("|---|---|---|---|")
    tot_new = Counter()
    tot_2025 = Counter()
    for b in sorted(new_by_bench):
        c = new_by_bench[b]
        nsat, nunsat = c.get('sat', 0), c.get('unsat', 0)
        nsolved = nsat + nunsat
        nerr = sum(v for k, v in c.items() if k not in ('sat', 'unsat', 'unknown', 'timeout'))
        nnv2025 = load_tool_bench(rd, 'nnv', b)
        s25 = sum(1 for v in nnv2025.values() if v == 'sat')
        u25 = sum(1 for v in nnv2025.values() if v == 'unsat')
        tot_new['sat'] += nsat; tot_new['unsat'] += nunsat; tot_new['solved'] += nsolved
        tot_2025['sat'] += s25; tot_2025['unsat'] += u25; tot_2025['solved'] += s25 + u25
        print(f"| {b} | {nsolved} ({nsat}/{nunsat}) | {s25+u25} ({s25}/{u25}) | "
              f"{c.get('unknown',0)}/{c.get('timeout',0)}/{nerr} |")
    print(f"\n**Totals -- NNV NEW:** solved {tot_new['solved']} (sat {tot_new['sat']} / unsat {tot_new['unsat']})")
    print(f"**Totals -- NNV 2025 (this results repo):** solved {tot_2025['solved']} "
          f"(sat {tot_2025['sat']} / unsat {tot_2025['unsat']}). Official 2025: 1082 solved, 697.3 pts, 6th/7.\n")

    # --- per-instance disagreement check (the soundness gate) ---
    print("## Soundness check: NNV-NEW verdicts vs the cross-tool reference\n")
    benches = sorted(new_by_bench)
    ref_cache = {b: {t: load_tool_bench(rd, t, b) for t in TOOLS} for b in benches}
    false_sat, false_unsat, agree, no_ref = [], [], 0, 0
    for (b, onnx, vl), v in new_inst.items():
        if v not in DEFINITIVE:
            continue
        per_tool = ref_cache.get(b, {})
        ref = per_tool.get(REF_TOOL, {}).get((onnx, vl))
        if ref not in DEFINITIVE:
            votes = Counter(per_tool[t].get((onnx, vl)) for t in TOOLS
                            if per_tool[t].get((onnx, vl)) in DEFINITIVE)
            ref = votes.most_common(1)[0][0] if votes else None
        if ref is None:
            no_ref += 1
        elif ref == v:
            agree += 1
        elif v == 'sat':
            false_sat.append((b, onnx, vl, ref))
        else:
            false_unsat.append((b, onnx, vl, ref))

    print(f"- NNV-NEW definitive verdicts: {agree + len(false_sat) + len(false_unsat) + no_ref}")
    print(f"- **AGREE with reference:** {agree}")
    print(f"- **NO reference** (no tool gave a definitive verdict): {no_ref}")
    print(f"- !! **NNV sat but reference unsat (probable FALSE SAT, -150 risk):** {len(false_sat)}")
    print(f"- !! **NNV unsat but reference sat (probable unsound proof):** {len(false_unsat)}\n")
    for tag, lst in (('FALSE-SAT', false_sat), ('FALSE-UNSAT', false_unsat)):
        if lst:
            print(f"### {tag} (NNV vs reference)\n")
            for b, onnx, vl, ref in lst[:60]:
                print(f"- `{b}` {onnx} / {vl}: NNV={'sat' if tag=='FALSE-SAT' else 'unsat'} vs reference={ref}")
            if len(lst) > 60:
                print(f"- ... and {len(lst)-60} more")
            print()
    if not false_sat and not false_unsat:
        print("OK: No NNV-NEW verdict disagreed with the cross-tool reference -- sound on the checked instances.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
