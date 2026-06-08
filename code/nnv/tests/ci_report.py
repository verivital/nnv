#!/usr/bin/env python3
"""Aggregate NNV matrix-CI shard results into one report.

Reads downloaded shard artifacts (JUnit XML + progress-*.txt) from a directory, then:
  - totals tests passed/failed/errored,
  - lists failures WITH their messages (preserves the "what happened" logs),
  - detects crashed shards (last 'START <t>' with no 'DONE <t>' in a progress file =
    the test that killed MATLAB; GitHub discards a crashed job's stdout),
  - classifies failures against ci_allowed_failures.txt (known WIP, kept visible but
    non-blocking) vs NEW regressions,
  - writes a Markdown summary to stdout and to $GITHUB_STEP_SUMMARY,
  - exits non-zero on any NEW failure or any crash (so CI is red only for real problems,
    while known-WIP/transformer failures stay green and visible).

Usage: python3 ci_report.py <artifacts-dir> [--allowed ci_allowed_failures.txt]
"""
import sys, os, glob, argparse, xml.etree.ElementTree as ET


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
            if fe is not None:
                results.append((name, 'fail', (fe.get('message') or fe.text or '').strip()[:500]))
            elif er is not None:
                results.append((name, 'error', (er.get('message') or er.text or '').strip()[:500]))
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

    total = len(res)
    npass = sum(1 for _, s, _ in res if s == 'pass')
    fails = [(n, s, m) for n, s, m in res if s in ('fail', 'error')]
    new = sorted(x for x in fails if x[0] not in allowed)
    known = sorted(x for x in fails if x[0] in allowed)

    out = ["## NNV matrix-CI report", ""]
    out.append(f"- **{npass}/{total}** passed across recorded shards")
    out.append(f"- **{len(fails)}** failed/errored — {len(new)} NEW, {len(known)} known-allowed")
    out.append(f"- **{len(crashes)}** shard crash(es)")
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

    sys.exit(1 if (new or crashes) else 0)


if __name__ == '__main__':
    main()
