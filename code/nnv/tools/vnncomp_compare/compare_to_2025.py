#!/usr/bin/env python3
"""compare_to_2025.py -- put NNV's fresh sweep next to the official VNN-COMP 2025
results for context.

CURRENT behaviour (the VNN-COMP results repo layout has varied year to year, so this is
deliberately conservative):
  1. Prints NNV NEW (this sweep) solved/sat/unsat per benchmark + totals, against the
     published 2025 NNV baseline (1082 solved / 697.3 pts).
  2. DISCOVERS the 2025 results repo structure (walks it for *.csv, identifies the likely
     NNV file, and prints a sample) so the exact column format is visible.

NOT yet implemented (the format above must be confirmed first): per-instance parsing of
the 2025 verdicts, per-tool totals, and flagging NNV-NEW sat/unsat that disagree with the
2025 consensus (the -150 risk). The discovery output makes wiring those a small, targeted
follow-up; the NNV NEW scorecard is exact and self-contained.

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


def discover_csvs(root):
    out = []
    for dp, _, fns in os.walk(root):
        for fn in fns:
            if fn.lower().endswith('.csv'):
                out.append(os.path.join(dp, fn))
    return out


def sniff(path, max_rows=3):
    # Read up to max_rows lines; a short file (fewer lines) must still return what it has
    # -- don't let StopIteration (or any read error) swallow the whole sample.
    rows = []
    try:
        with open(path, newline='', errors='replace') as f:
            for line in f:
                rows.append(line.rstrip('\n'))
                if len(rows) >= max_rows:
                    break
    except OSError:
        pass
    return rows


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
    print("\n> NEXT: per-instance parsing of the 2025 verdicts and agreement-checking against the "
          "consensus are NOT yet implemented -- add them to `compare_to_2025.py` using the format "
          "shown above (a small, targeted follow-up). The NNV NEW scorecard is exact and self-contained.")
    return 0


if __name__ == '__main__':
    sys.exit(main())
