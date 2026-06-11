#!/usr/bin/env python3
"""Aggregate NNV matrix-CI shard results into one report.

Reads downloaded shard artifacts (JUnit XML + progress-*.txt) from a directory, then:
  - totals tests passed/failed/errored,
  - lists failures WITH their messages (preserves the "what happened" logs),
  - detects crashed shards (last 'START <t>' with no 'DONE <t>' in a progress file =
    the test that killed MATLAB; GitHub discards a crashed job's stdout),
  - detects MISSING shards (a shard that produced no artifact at all = a hang/crash whose
    upload was skipped by job cancellation),
  - classifies failures against ci_allowed_failures.txt (known WIP = visible/non-blocking)
    vs NEW regressions,
  - writes a Markdown summary to stdout and $GITHUB_STEP_SUMMARY,
  - exits non-zero on any NEW failure, crash, or missing shard.

Usage: python3 ci_report.py <artifacts-dir> [--allowed ci_allowed_failures.txt]
"""
import sys, os, re, glob, argparse, xml.etree.ElementTree as ET


def load_allowed(path):
    allowed = set()
    if path and os.path.isfile(path):
        with open(path, encoding='utf-8') as fh:
            for line in fh:
                s = line.strip()
                if s and not s.startswith('#'):
                    allowed.add(s)
    return allowed


def parse_junit(d):
    results = []  # (name, status, message)
    for f in glob.glob(os.path.join(d, '**', '*.xml'), recursive=True):
        try:
            root = ET.parse(f).getroot()
        except Exception:
            continue
        for tc in root.iter('testcase'):
            name = (tc.get('classname', '') or '') + '/' + (tc.get('name', '') or '')
            fe = tc.find('failure')
            er = tc.find('error')
            sk = tc.find('skipped')
            if fe is not None:
                results.append((name, 'fail', (fe.get('message') or fe.text or '').strip()[:500]))
            elif er is not None:
                results.append((name, 'error', (er.get('message') or er.text or '').strip()[:500]))
            elif sk is not None:
                # A <skipped> testcase (MATLAB assumption failure / filtered) is
                # NOT a pass. Counting it as pass let a reachability test that
                # assumeFail'd masquerade as green (a soundness-gate hole).
                results.append((name, 'skip', (sk.get('message') or sk.text or '').strip()[:500]))
            else:
                results.append((name, 'pass', ''))
    return results


def parse_crashes(d):
    crashes = []  # (progress-file, crashing-test)
    for pf in glob.glob(os.path.join(d, '**', 'progress-*.txt'), recursive=True):
        starts, done = [], set()
        with open(pf, encoding='utf-8', errors='replace') as fh:
            for line in fh:
                if line.startswith('START '):
                    starts.append(line[6:].strip())
                elif line.startswith('DONE  '):
                    done.add(line[6:].strip())
        if starts and starts[-1] not in done:
            crashes.append((os.path.basename(pf), starts[-1]))
    return crashes


def detect_missing_shards(d):
    """Determine expected shard count from artifact filenames (results-*-shardKK-ofNN)
    and return the list of shard indices that produced NO artifact (their job hung/crashed
    and the upload was skipped)."""
    present, total = set(), None
    for f in (glob.glob(os.path.join(d, '**', '*.xml'), recursive=True)
              + glob.glob(os.path.join(d, '**', 'progress-*.txt'), recursive=True)):
        m = re.search(r'shard(\d+)-of(\d+)', os.path.basename(f))
        if m:
            present.add(int(m.group(1)))
            total = int(m.group(2))
    if total is None:
        return []
    return [k for k in range(1, total + 1) if k not in present]


def main():
    try:
        sys.stdout.reconfigure(encoding='utf-8')   # emojis print on Windows consoles too
    except Exception:
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument('artifacts')
    ap.add_argument('--allowed', default='')
    a = ap.parse_args()

    allowed = load_allowed(a.allowed)
    res = parse_junit(a.artifacts)
    crashes = parse_crashes(a.artifacts)
    missing = detect_missing_shards(a.artifacts)

    total = len(res)
    npass = sum(1 for _, s, _ in res if s == 'pass')
    fails = [(n, s, m) for n, s, m in res if s in ('fail', 'error')]
    skips = [(n, s, m) for n, s, m in res if s == 'skip']
    new = sorted(x for x in fails if x[0] not in allowed)
    known = sorted(x for x in fails if x[0] in allowed)
    # A skipped soundness/reach test is a gate hole (assumeFail masking a real
    # error). Surface these loudly; they should be verifyError/containment, not
    # assumptions. Not auto-red (legitimate GPU/toolbox skips exist), but visible.
    susp_skips = sorted(n for n, _, _ in skips
                        if any(k in n.lower() for k in ('sound', 'reach', 'containment', 'attention')))

    out = ["## NNV matrix-CI report", ""]
    out.append(f"- **{npass}/{total}** passed, **{len(skips)}** skipped across recorded shards")
    out.append(f"- **{len(fails)}** failed/errored — {len(new)} NEW, {len(known)} known-allowed")
    out.append(f"- **{len(crashes)}** shard crash(es), **{len(missing)}** missing shard(s)")
    out.append("")
    if total == 0:
        out.append("### ⛔ NO test results recorded — the gate is meaningless; failing the build.")
        out.append("")
    if susp_skips:
        out.append("### ⚠️ SKIPPED soundness/reach tests (assumeFail masking? should be verifyError/containment)")
        for n in susp_skips:
            out.append(f"- {n}")
        out.append("")
    if missing:
        out.append("### ⛔ MISSING shards — produced no results (hang/crash whose upload was skipped)")
        out.append(f"- shard(s): {', '.join(map(str, missing))}. Their tests did not run/report; "
                   f"check the job's live log for the last `START` (the hanger).")
        out.append("")
    if crashes:
        out.append("### 💥 Crashes — test that killed the shard (exit 255 / OOM; add to confirmedCrashPatterns)")
        for pf, t in crashes:
            out.append(f"- `{t}`  _({pf})_")
        out.append("")
    if new:
        out.append("### ❌ NEW failures (not in ci_allowed_failures.txt) — these red the build")
        for n, s, m in new:
            out.append(f"- **{n}** ({s}): {m or '(no message)'}")
        out.append("")
    if known:
        out.append("### 🟡 Known/allowed failures (WIP — visible but non-blocking)")
        for n, s, _ in known:
            out.append(f"- {n} ({s})")
        out.append("")
    report = "\n".join(out)
    print(report)
    gs = os.environ.get('GITHUB_STEP_SUMMARY')
    if gs:
        with open(gs, 'a', encoding='utf-8') as fh:
            fh.write(report + "\n")

    # total==0 means no shard reported anything (all wiped out / all skipped) ->
    # the gate would otherwise pass "0/0" green. Treat as failure.
    sys.exit(1 if (new or crashes or missing or total == 0) else 0)


if __name__ == '__main__':
    main()
