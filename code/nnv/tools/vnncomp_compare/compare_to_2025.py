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
            if len(row) < 5 or row[0].strip().lower() == 'category':
                continue                       # skip blank rows and the header line
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
    # Reference verdict = MAJORITY VOTE of the other tools' definitive (sat/unsat) results
    # (per the soundness strategy: when tools disagree, trust the majority). We also report
    # where the field itself splits sat-vs-unsat (a likely soundness issue in SOME tool).
    print("## Soundness check: NNV-NEW verdicts vs the cross-tool MAJORITY\n")
    benches = sorted(new_by_bench)
    ref_cache = {b: {t: load_tool_bench(rd, t, b) for t in TOOLS} for b in benches}
    false_sat, false_unsat, agree, no_ref = [], [], 0, 0
    contested, ties, gold_disagree = [], [], []
    for (b, onnx, vl), v in new_inst.items():
        if v not in DEFINITIVE:
            continue
        per_tool = ref_cache.get(b, {})
        votes = Counter(per_tool[t].get((onnx, vl)) for t in TOOLS
                        if per_tool[t].get((onnx, vl)) in DEFINITIVE)
        n_sat, n_unsat = votes.get('sat', 0), votes.get('unsat', 0)
        # alpha-beta-CROWN was sound+complete with 0 incorrect verdicts in 2025: a
        # disagreement with it is high-signal even when there's no clear majority.
        gold = per_tool.get(REF_TOOL, {}).get((onnx, vl))
        if gold in DEFINITIVE and gold != v:
            gold_disagree.append((b, onnx, vl, v, gold, dict(votes)))
        if n_sat and n_unsat:                      # the field itself disagrees sat-vs-unsat
            contested.append((b, onnx, vl, dict(votes), v))
        if n_sat == 0 and n_unsat == 0:
            no_ref += 1
            continue
        if n_sat == n_unsat:                       # a tie is NOT a majority -- don't hard-flag
            ties.append((b, onnx, vl, dict(votes), v))
            continue
        ref = 'sat' if n_sat > n_unsat else 'unsat'   # strict majority
        if ref == v:
            agree += 1
        elif v == 'sat':
            false_sat.append((b, onnx, vl, dict(votes)))
        else:
            false_unsat.append((b, onnx, vl, dict(votes)))

    print(f"- NNV-NEW definitive verdicts checked: {agree + len(false_sat) + len(false_unsat) + len(ties) + no_ref}")
    print(f"- **AGREE with majority:** {agree}")
    print(f"- **NO reference** (no other tool gave a definitive verdict): {no_ref}")
    print(f"- **NO majority** (tools tied sat-vs-unsat): {len(ties)} -- see contested list, not hard-flagged")
    print(f"- !! **NNV sat but majority unsat (probable FALSE SAT, -150 risk):** {len(false_sat)}")
    print(f"- !! **NNV unsat but majority sat (probable unsound proof):** {len(false_unsat)}")
    print(f"- !! **NNV disagrees with alpha-beta-CROWN (gold standard):** {len(gold_disagree)}\n")
    for tag, lst in (('FALSE-SAT', false_sat), ('FALSE-UNSAT', false_unsat)):
        if lst:
            print(f"### {tag}: NNV={'sat' if tag=='FALSE-SAT' else 'unsat'} vs the majority\n")
            for b, onnx, vl, votes in lst[:80]:
                print(f"- `{b}` {onnx} / {vl}: tool votes {dict(votes)}")
            if len(lst) > 80:
                print(f"- ... and {len(lst)-80} more")
            print()
    if gold_disagree:
        print("### NNV disagrees with alpha-beta-CROWN (review first -- gold standard)\n")
        for b, onnx, vl, nv, g, votes in gold_disagree[:80]:
            print(f"- `{b}` {onnx} / {vl}: NNV={nv} vs a-b-CROWN={g} (all votes {dict(votes)})")
        if len(gold_disagree) > 80:
            print(f"- ... and {len(gold_disagree)-80} more")
        print()
    if not false_sat and not false_unsat and not gold_disagree:
        print("OK: No NNV-NEW verdict disagreed with the majority or with alpha-beta-CROWN -- "
              "sound on the checked instances.\n")

    # --- coverage / opportunity: where the FIELD solves but NNV does not ---
    print("## Coverage gap: instances the FIELD solved but NNV did not (minimize these)\n")
    print("| benchmark | field-solved & NNV-unknown/to/err | NNV-solved |")
    print("|---|--:|--:|")
    for b in benches:
        per_tool = ref_cache.get(b, {})
        keys = set().union(*[set(per_tool[t].keys()) for t in TOOLS]) if per_tool else set()
        gap = 0
        nnv_solved = sum(1 for (bb, o, vv), x in new_inst.items()
                         if bb == b and x in DEFINITIVE)
        for (o, vl) in keys:
            votes = Counter(per_tool[t].get((o, vl)) for t in TOOLS
                            if per_tool[t].get((o, vl)) in DEFINITIVE)
            if votes and new_inst.get((b, o, vl)) not in DEFINITIVE:
                gap += 1
        if gap or nnv_solved:
            print(f"| {b} | {gap} | {nnv_solved} |")
    return 0


if __name__ == '__main__':
    sys.exit(main())
