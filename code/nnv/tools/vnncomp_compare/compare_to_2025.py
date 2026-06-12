#!/usr/bin/env python3
"""compare_to_2025.py -- soundness gate + scorecard for an NNV sweep against the official
VNN-COMP 2025 results.

Reads `<results2025>/<tool>/2025_<bench>/results.csv` for each reference tool and, for every
NNV definitive (sat/unsat) verdict, compares against the cross-tool result:
  1. Per-benchmark NNV-NEW vs NNV-2025 solved (sat/unsat) + the NNV-new unknown/timeout/error
     split, against the published 2025 NNV baseline (1082 solved / 697.3 pts).
  2. Soundness check (classify_instances): reference = STRICT majority of the other tools'
     definitive verdicts. NNV sat vs majority unsat -> FALSE-SAT (-150 risk); NNV unsat vs
     majority sat -> FALSE-UNSAT. A sat-vs-unsat TIE is reported as "no majority" (NOT
     hard-flagged). Plus a separate "disagrees with alpha-beta-CROWN (gold standard)" list
     (a-b-CROWN was 0-incorrect in 2025), high-signal even on a tie.
  3. Coverage gap: instances the field solved but NNV left unknown/timeout/error.

The per-instance classification lives in the pure function `classify_instances` (unit-tested
in test_compare.py); main() only formats its output.

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


def norm(v):
    """Alias of norm_verdict (verdict-string normalizer)."""
    return norm_verdict(v)


def load_tool_bench(root, tool, bench):
    """Load one tool's per-benchmark results.

    Expects:  <root>/<tool>/2025_<bench>/results.csv
    with the VNN-COMP per-tool column layout:
        category,onnx_path,vnnlib_path,prepare_time,result,run_time
    Returns {(onnx_base, vnnlib_base): normalized_verdict}. The HEADER row is
    skipped (so the literal ('onnx_path','vnnlib_path') key never appears); only
    data rows are keyed. Missing file -> {}.
    """
    path = os.path.join(root, tool, '2025_%s' % bench, 'results.csv')
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path, newline='') as f:
        for r in csv.reader(f):
            if len(r) < 6:
                continue
            onnx_col, vnnlib_col, verdict_col = r[1], r[2], r[4]
            # skip the header row (literal column names), not a real instance
            if onnx_col.strip() == 'onnx_path' and vnnlib_col.strip() == 'vnnlib_path':
                continue
            out[(base(onnx_col), base(vnnlib_col))] = norm_verdict(verdict_col)
    return out


def _majority(verdicts):
    """Given an iterable of solved (sat/unsat) verdicts, return ('sat'|'unsat', strict)
    where strict is True iff there is a UNIQUE majority. A tie (equal sat/unsat) ->
    (None, False): there is no majority, so NNV must NOT be hard-flagged against it.
    Verdicts that are neither sat nor unsat are ignored.
    """
    c = Counter(v for v in verdicts if v in ('sat', 'unsat'))
    sat = c.get('sat', 0)
    unsat = c.get('unsat', 0)
    if sat == 0 and unsat == 0:
        return None, False
    if sat == unsat:
        return None, False           # tie -> no majority
    return ('sat' if sat > unsat else 'unsat'), True


def classify_instances(new_inst, ref_cache, tools, ref_tool='alpha_beta_crown'):
    """Pure soundness classifier for an NNV sweep against reference-tool verdicts.

    Args:
      new_inst   : {(bench, onnx_base, vnnlib_base): nnv_verdict}  (from load_new)
      ref_cache  : {(tool, bench): {(onnx_base, vnnlib_base): verdict}}  (load_tool_bench)
      tools      : list of reference tool names to form the majority over.
      ref_tool   : the gold-standard tool (default alpha_beta_crown). NNV verdicts that
                   DISAGREE with this tool specifically are collected separately, even
                   on a majority tie.

    Returns a dict of lists of keys:
      agree        : NNV solved verdict == majority.
      false_sat    : NNV says sat, strict majority says unsat (HARD false flag, -150 risk).
      false_unsat  : NNV says unsat, strict majority says sat (HARD false flag).
      ties         : NNV solved but reference is a sat/unsat TIE -> NO majority, NOT a
                     hard false-flag (counted separately).
      no_ref       : no reference solved verdict available for the instance.
      contested    : the reference tools THEMSELVES disagree (both sat and unsat appear
                     among them) -> the instance is contested in the field, independent of
                     what NNV said. (A superset of `ties`: a tie is the special case of
                     an even split; an odd split is contested with a strict majority.)
      gold_disagree: NNV's solved verdict disagrees specifically with ref_tool's solved
                     verdict (the gold standard), regardless of majority/tie.
      coverage_gap : the field (majority) SOLVED the instance but NNV did NOT (unknown/
                     timeout/error/missing) -> a coverage gap, not a soundness error.
    """
    res = {
        'agree': [], 'false_sat': [], 'false_unsat': [], 'ties': [],
        'no_ref': [], 'contested': [], 'gold_disagree': [], 'coverage_gap': [],
    }
    for key, nnv_v in new_inst.items():
        bench, onnx_b, vnnlib_b = key
        ikey = (onnx_b, vnnlib_b)

        # gather reference verdicts for this instance across the requested tools
        ref_verdicts = []
        for t in tools:
            d = ref_cache.get((t, bench))
            if d and ikey in d:
                ref_verdicts.append(d[ikey])

        # gold-standard (ref_tool) verdict, if present
        gold = None
        gd = ref_cache.get((ref_tool, bench))
        if gd and ikey in gd:
            gold = gd[ikey]

        nnv_solved = nnv_v in ('sat', 'unsat')
        maj, strict = _majority(ref_verdicts)
        solved_ref = [v for v in ref_verdicts if v in ('sat', 'unsat')]

        # reference tools themselves split (both polarities present) -> contested instance
        if 'sat' in solved_ref and 'unsat' in solved_ref:
            res['contested'].append(key)

        # gold-standard disagreement: both sides solved and differ (independent of majority)
        if nnv_solved and gold in ('sat', 'unsat') and gold != nnv_v:
            res['gold_disagree'].append(key)

        if not solved_ref:
            res['no_ref'].append(key)
            continue

        # coverage gap: field solved it, NNV did not
        if strict and not nnv_solved:
            res['coverage_gap'].append(key)
            continue

        if not nnv_solved:
            # NNV unknown and no strict field majority either -> nothing to flag
            continue

        if not strict:
            # reference is a sat/unsat tie -> no majority -> NOT a hard false-flag
            res['ties'].append(key)
            continue

        if maj == nnv_v:
            res['agree'].append(key)
        elif nnv_v == 'sat':           # NNV sat vs majority unsat
            res['false_sat'].append(key)
        else:                          # NNV unsat vs majority sat
            res['false_unsat'].append(key)

    return res


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


# The reference tools whose verdicts form the cross-tool majority (alpha-beta-CROWN first
# -- it was sound+complete with 0 incorrect verdicts in 2025, so it is also the gold standard).
TOOLS = ['alpha_beta_crown', 'neuralsat', 'nnenum', 'pyrat', 'cora', 'rover', 'sobolbox']
REF_TOOL = 'alpha_beta_crown'


def _solved_counts(inst_map):
    """(#sat, #unsat) over a {key: verdict} map."""
    sat = sum(1 for v in inst_map.values() if v == 'sat')
    unsat = sum(1 for v in inst_map.values() if v == 'unsat')
    return sat, unsat


def _votes(ref_cache, bench, onnx, vl, tools):
    """{tool: verdict} for the tools that gave a definitive verdict on this instance."""
    out = {}
    for t in tools:
        v = ref_cache.get((t, bench), {}).get((onnx, vl))
        if v in ('sat', 'unsat'):
            out[t] = v
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--new', required=True)
    ap.add_argument('--results2025', required=True)
    args = ap.parse_args()

    new_by_bench, new_inst = load_new(args.new)
    rd = args.results2025
    benches = sorted(new_by_bench)

    print("# NNV sweep vs VNN-COMP 2025\n")

    # --- per-benchmark: NNV NEW vs NNV 2025 ---
    print("## Per-benchmark: NNV NEW vs NNV 2025\n")
    print("| benchmark | NNV-new solved (sat/unsat) | NNV-2025 solved (sat/unsat) | NNV-new unk/to/err |")
    print("|---|---|---|---|")
    tot_new = Counter()
    for b in benches:
        c = new_by_bench[b]
        nsat, nunsat = c.get('sat', 0), c.get('unsat', 0)
        err = sum(v for k, v in c.items() if k in ('error', 'missing', 'decompress_failed',
                                                   'no_instances_csv', 'no_resolvable_instance'))
        osat, ounsat = _solved_counts(load_tool_bench(rd, 'nnv', b))
        tot_new['sat'] += nsat
        tot_new['unsat'] += nunsat
        print(f"| {b} | {nsat+nunsat} ({nsat}/{nunsat}) | {osat+ounsat} ({osat}/{ounsat}) | "
              f"{c.get('unknown',0)}/{c.get('timeout',0)}/{err} |")
    print(f"\n**Totals -- NNV NEW:** solved {tot_new['sat']+tot_new['unsat']} "
          f"(sat {tot_new['sat']} / unsat {tot_new['unsat']}).")
    print("Baseline (2025 official): NNV solved 1082 (354 sat / 728 unsat), 6th/7, 697.3 pts.\n")

    # --- soundness check: NNV verdicts vs the cross-tool MAJORITY (and gold standard) ---
    ref_cache = {(t, b): load_tool_bench(rd, t, b) for t in TOOLS for b in benches}
    res = classify_instances(new_inst, ref_cache, TOOLS, REF_TOOL)
    nchecked = (len(res['agree']) + len(res['false_sat']) + len(res['false_unsat'])
                + len(res['ties']) + len(res['no_ref']))
    print("## Soundness check: NNV-NEW verdicts vs the cross-tool MAJORITY\n")
    print(f"- NNV-NEW definitive verdicts checked: {nchecked}")
    print(f"- **AGREE with majority:** {len(res['agree'])}")
    print(f"- **NO reference** (no other tool gave a definitive verdict): {len(res['no_ref'])}")
    print(f"- **NO majority** (tools tied sat-vs-unsat): {len(res['ties'])} -- contested, not hard-flagged")
    print(f"- !! **NNV sat but majority unsat (probable FALSE SAT, -150 risk):** {len(res['false_sat'])}")
    print(f"- !! **NNV unsat but majority sat (probable unsound proof):** {len(res['false_unsat'])}")
    print(f"- !! **NNV disagrees with alpha-beta-CROWN (gold standard):** {len(res['gold_disagree'])}\n")

    def show(title, keys):
        if not keys:
            return
        print(f"### {title}\n")
        for (b, onnx, vl) in keys[:80]:
            print(f"- `{b}` {onnx} / {vl}: NNV={new_inst[(b, onnx, vl)]}, tool votes {_votes(ref_cache, b, onnx, vl, TOOLS)}")
        if len(keys) > 80:
            print(f"- ... and {len(keys)-80} more")
        print()
    show("FALSE-SAT: NNV=sat vs the majority", res['false_sat'])
    show("FALSE-UNSAT: NNV=unsat vs the majority", res['false_unsat'])
    show("NNV disagrees with alpha-beta-CROWN (review first -- gold standard)", res['gold_disagree'])
    if not res['false_sat'] and not res['false_unsat'] and not res['gold_disagree']:
        print("OK: No NNV-NEW verdict disagreed with the majority or with alpha-beta-CROWN -- "
              "sound on the checked instances.\n")

    # --- coverage gap: the field solved it but NNV did not (minimize these) ---
    print("## Coverage gap: instances the FIELD solved but NNV did not (minimize these)\n")
    gap_by_bench = Counter(b for (b, o, vl) in res['coverage_gap'])
    solved_by_bench = Counter(b for (b, o, vl), v in new_inst.items() if v in ('sat', 'unsat'))
    print("| benchmark | field-solved & NNV-unknown/to/err | NNV-solved |")
    print("|---|--:|--:|")
    for b in benches:
        if gap_by_bench[b] or solved_by_bench[b]:
            print(f"| {b} | {gap_by_bench[b]} | {solved_by_bench[b]} |")
    return 0


if __name__ == '__main__':
    sys.exit(main())
