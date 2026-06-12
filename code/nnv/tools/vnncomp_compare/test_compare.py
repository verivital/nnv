#!/usr/bin/env python3
"""Unit tests for the VNN-COMP comparison tooling (stdlib unittest only).

Run from this directory:
    python -m unittest test_compare

Covers:
  - compare_to_2025.base / norm_verdict / norm
  - compare_to_2025.load_tool_bench  (header skipping + keying)
  - compare_to_2025.classify_instances  (per-instance soundness logic)
  - summarize_sweep tallies + time-stat exclusion of 0-time rows
"""
import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stdout

# import the modules under test (this file lives alongside them)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_to_2025 as cmp  # noqa: E402
import summarize_sweep as summ  # noqa: E402


class TestBase(unittest.TestCase):
    def test_strip_extensions(self):
        self.assertEqual(cmp.base('foo.onnx'), 'foo')
        self.assertEqual(cmp.base('foo.vnnlib'), 'foo')
        self.assertEqual(cmp.base('foo.csv'), 'foo')
        self.assertEqual(cmp.base('foo.gz'), 'foo')

    def test_strips_directory(self):
        self.assertEqual(cmp.base('a/b/c/model.onnx'), 'model')
        # os.path.basename uses the platform separator; backslash on Windows.
        self.assertEqual(cmp.base(os.path.join('x', 'y', 'spec.vnnlib')), 'spec')

    def test_ext_strip_is_single_pass_outer_first(self):
        # base() strips ONE matching extension per ext in order (.gz then .onnx ...).
        # 'model.onnx.gz' -> '.gz' stripped -> 'model.onnx' -> '.onnx' stripped -> 'model'
        self.assertEqual(cmp.base('model.onnx.gz'), 'model')
        # '.vnnlib.gz' -> 'spec.vnnlib' -> 'spec'
        self.assertEqual(cmp.base('spec.vnnlib.gz'), 'spec')

    def test_no_known_extension_kept(self):
        self.assertEqual(cmp.base('weird.txt'), 'weird.txt')
        self.assertEqual(cmp.base('plain'), 'plain')

    def test_whitespace_trimmed(self):
        self.assertEqual(cmp.base('  spaced.onnx  '), 'spaced')


class TestNormVerdict(unittest.TestCase):
    def test_holds_to_unsat(self):
        self.assertEqual(cmp.norm_verdict('holds'), 'unsat')

    def test_violated_to_sat(self):
        self.assertEqual(cmp.norm_verdict('violated'), 'sat')

    def test_true_to_unsat(self):
        self.assertEqual(cmp.norm_verdict('true'), 'unsat')

    def test_false_to_sat(self):
        self.assertEqual(cmp.norm_verdict('false'), 'sat')

    def test_trim_and_lowercase(self):
        self.assertEqual(cmp.norm_verdict('  SAT '), 'sat')
        self.assertEqual(cmp.norm_verdict('UNSAT'), 'unsat')
        self.assertEqual(cmp.norm_verdict(' Holds '), 'unsat')

    def test_passthrough_unknown(self):
        self.assertEqual(cmp.norm_verdict('unknown'), 'unknown')
        self.assertEqual(cmp.norm_verdict('timeout'), 'timeout')
        self.assertEqual(cmp.norm_verdict('error'), 'error')

    def test_none_and_empty(self):
        self.assertEqual(cmp.norm_verdict(None), '')
        self.assertEqual(cmp.norm_verdict(''), '')

    def test_norm_is_alias(self):
        for v in ('holds', 'violated', ' TRUE ', 'false', 'sat', None):
            self.assertEqual(cmp.norm(v), cmp.norm_verdict(v))


class TestLoadToolBench(unittest.TestCase):
    HEADER = 'category,onnx_path,vnnlib_path,prepare_time,result,run_time\n'

    def _write(self, root, tool, bench, body):
        d = os.path.join(root, tool, '2025_%s' % bench)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, 'results.csv'), 'w', newline='') as f:
            f.write(self.HEADER)
            f.write(body)

    def test_header_skipped_and_keyed_by_bases(self):
        with tempfile.TemporaryDirectory() as root:
            self._write(root, 'nnv', 'acasxu',
                        'acasxu,model_a.onnx,prop_1.vnnlib,0.1,holds,2.3\n'
                        'acasxu,model_b.onnx.gz,prop_2.vnnlib.gz,0.2,violated,1.1\n')
            d = cmp.load_tool_bench(root, 'nnv', 'acasxu')

        # the header row's literal column names must NOT appear as a key
        self.assertNotIn(('onnx_path', 'vnnlib_path'), d)
        # data rows keyed by (onnx_base, vnnlib_base) with normalized verdicts
        self.assertEqual(d[('model_a', 'prop_1')], 'unsat')   # holds -> unsat
        self.assertEqual(d[('model_b', 'prop_2')], 'sat')     # violated -> sat ; .gz/.onnx stripped
        self.assertEqual(len(d), 2)

    def test_missing_file_returns_empty(self):
        with tempfile.TemporaryDirectory() as root:
            self.assertEqual(cmp.load_tool_bench(root, 'nnv', 'nope'), {})

    def test_short_rows_ignored(self):
        with tempfile.TemporaryDirectory() as root:
            self._write(root, 't', 'b',
                        'b,m.onnx,s.vnnlib,0.0,sat,0.5\n'
                        'truncated,row\n')              # < 6 cols -> skipped
            d = cmp.load_tool_bench(root, 't', 'b')
        self.assertEqual(d, {('m', 's'): 'sat'})


class TestClassifyInstances(unittest.TestCase):
    REF_TOOL = 'alpha_beta_crown'

    def _cache(self, mapping):
        # mapping: {(tool, bench): {(onnx, vnnlib): verdict}}
        return dict(mapping)

    def test_nnv_sat_vs_majority_unsat_is_false_sat(self):
        new = {('b', 'm', 's'): 'sat'}
        ref = self._cache({
            ('t1', 'b'): {('m', 's'): 'unsat'},
            ('t2', 'b'): {('m', 's'): 'unsat'},
            ('t3', 'b'): {('m', 's'): 'unsat'},
        })
        out = cmp.classify_instances(new, ref, ['t1', 't2', 't3'], ref_tool=self.REF_TOOL)
        self.assertIn(('b', 'm', 's'), out['false_sat'])
        self.assertEqual(out['false_unsat'], [])
        self.assertEqual(out['agree'], [])

    def test_nnv_unsat_vs_majority_sat_is_false_unsat(self):
        new = {('b', 'm', 's'): 'unsat'}
        ref = self._cache({
            ('t1', 'b'): {('m', 's'): 'sat'},
            ('t2', 'b'): {('m', 's'): 'sat'},
        })
        out = cmp.classify_instances(new, ref, ['t1', 't2'], ref_tool=self.REF_TOOL)
        self.assertIn(('b', 'm', 's'), out['false_unsat'])
        self.assertEqual(out['false_sat'], [])

    def test_agreement(self):
        new = {('b', 'm', 's'): 'unsat'}
        ref = self._cache({
            ('t1', 'b'): {('m', 's'): 'unsat'},
            ('t2', 'b'): {('m', 's'): 'unsat'},
        })
        out = cmp.classify_instances(new, ref, ['t1', 't2'], ref_tool=self.REF_TOOL)
        self.assertIn(('b', 'm', 's'), out['agree'])
        self.assertEqual(out['false_sat'], [])
        self.assertEqual(out['false_unsat'], [])

    def test_tie_is_not_a_hard_false_flag(self):
        # 1 sat vs 1 unsat -> NO majority. NNV sat must be in 'ties', not 'false_sat'.
        new = {('b', 'm', 's'): 'sat'}
        ref = self._cache({
            ('t1', 'b'): {('m', 's'): 'sat'},
            ('t2', 'b'): {('m', 's'): 'unsat'},
        })
        out = cmp.classify_instances(new, ref, ['t1', 't2'], ref_tool=self.REF_TOOL)
        self.assertIn(('b', 'm', 's'), out['ties'])
        self.assertNotIn(('b', 'm', 's'), out['false_sat'])
        self.assertNotIn(('b', 'm', 's'), out['false_unsat'])
        # a tie with both polarities is also a 'contested' instance
        self.assertIn(('b', 'm', 's'), out['contested'])

    def test_gold_disagree_listed_even_on_a_tie(self):
        # NNV sat; alpha_beta_crown says unsat; the field is a tie overall.
        new = {('b', 'm', 's'): 'sat'}
        ref = self._cache({
            (self.REF_TOOL, 'b'): {('m', 's'): 'unsat'},
            ('other', 'b'): {('m', 's'): 'sat'},
        })
        out = cmp.classify_instances(new, ref, [self.REF_TOOL, 'other'], ref_tool=self.REF_TOOL)
        # tie -> not a hard flag
        self.assertIn(('b', 'm', 's'), out['ties'])
        self.assertNotIn(('b', 'm', 's'), out['false_sat'])
        # but it DOES disagree with the gold standard specifically
        self.assertIn(('b', 'm', 's'), out['gold_disagree'])

    def test_gold_agree_not_flagged(self):
        new = {('b', 'm', 's'): 'sat'}
        ref = self._cache({
            (self.REF_TOOL, 'b'): {('m', 's'): 'sat'},
            ('other', 'b'): {('m', 's'): 'unsat'},
        })
        out = cmp.classify_instances(new, ref, [self.REF_TOOL, 'other'], ref_tool=self.REF_TOOL)
        self.assertNotIn(('b', 'm', 's'), out['gold_disagree'])

    def test_no_ref_when_reference_unknown(self):
        new = {('b', 'm', 's'): 'sat'}
        ref = self._cache({
            ('t1', 'b'): {('m', 's'): 'unknown'},
            ('t2', 'b'): {('m', 's'): 'timeout'},
        })
        out = cmp.classify_instances(new, ref, ['t1', 't2'], ref_tool=self.REF_TOOL)
        self.assertIn(('b', 'm', 's'), out['no_ref'])
        self.assertEqual(out['false_sat'], [])

    def test_no_ref_when_instance_absent(self):
        new = {('b', 'm', 's'): 'sat'}
        ref = self._cache({('t1', 'b'): {('other', 'spec'): 'unsat'}})
        out = cmp.classify_instances(new, ref, ['t1'], ref_tool=self.REF_TOOL)
        self.assertIn(('b', 'm', 's'), out['no_ref'])

    def test_coverage_gap_field_solved_nnv_unknown(self):
        # field has a strict majority (solved); NNV says unknown -> coverage gap, not error.
        new = {('b', 'm', 's'): 'unknown'}
        ref = self._cache({
            ('t1', 'b'): {('m', 's'): 'unsat'},
            ('t2', 'b'): {('m', 's'): 'unsat'},
        })
        out = cmp.classify_instances(new, ref, ['t1', 't2'], ref_tool=self.REF_TOOL)
        self.assertIn(('b', 'm', 's'), out['coverage_gap'])
        self.assertEqual(out['false_sat'], [])
        self.assertEqual(out['false_unsat'], [])

    def test_nnv_unknown_and_field_tie_not_coverage_gap(self):
        # field is a tie (no strict majority) and NNV unknown -> nothing flagged,
        # specifically NOT a coverage gap (the field didn't "solve" it either).
        new = {('b', 'm', 's'): 'unknown'}
        ref = self._cache({
            ('t1', 'b'): {('m', 's'): 'sat'},
            ('t2', 'b'): {('m', 's'): 'unsat'},
        })
        out = cmp.classify_instances(new, ref, ['t1', 't2'], ref_tool=self.REF_TOOL)
        self.assertNotIn(('b', 'm', 's'), out['coverage_gap'])
        self.assertIn(('b', 'm', 's'), out['contested'])

    def test_contested_odd_split_with_majority(self):
        # 2 unsat vs 1 sat: strict majority unsat, but the field is contested.
        new = {('b', 'm', 's'): 'unsat'}
        ref = self._cache({
            ('t1', 'b'): {('m', 's'): 'unsat'},
            ('t2', 'b'): {('m', 's'): 'unsat'},
            ('t3', 'b'): {('m', 's'): 'sat'},
        })
        out = cmp.classify_instances(new, ref, ['t1', 't2', 't3'], ref_tool=self.REF_TOOL)
        self.assertIn(('b', 'm', 's'), out['agree'])       # NNV matches the majority
        self.assertIn(('b', 'm', 's'), out['contested'])   # but tools split


class TestSummarizeSweep(unittest.TestCase):
    def _run(self, csv_text):
        with tempfile.NamedTemporaryFile('w', suffix='.csv', delete=False, newline='') as f:
            f.write(csv_text)
            path = f.name
        try:
            buf = io.StringIO()
            old = sys.argv
            sys.argv = ['summarize_sweep.py', path]
            try:
                with redirect_stdout(buf):
                    rc = summ.main()
            finally:
                sys.argv = old
            return rc, buf.getvalue()
        finally:
            os.unlink(path)

    CSV = (
        'subfolder,category,onnx,vnnlib,status,status_str,time_s,error_message\n'
        'acasxu,acasxu,a.onnx,p1.vnnlib,0,sat,1.5,\n'
        'acasxu,acasxu,a.onnx,p2.vnnlib,1,unsat,2.5,\n'
        'acasxu,acasxu,a.onnx,p3.vnnlib,2,unknown,3.0,\n'
        'cersyve,cersyve,b.onnx,q1.vnnlib,-1,timeout,0,exceeded\n'
        'cersyve,cersyve,b.onnx,q2.vnnlib,-1,error,0,boom\n'
        'cersyve,cersyve,b.onnx,q3.vnnlib,-2,missing,0,\n'
    )

    def test_tallies(self):
        rc, out = self._run(self.CSV)
        self.assertEqual(rc, 0)
        self.assertIn('**Instances:** 6', out)
        self.assertIn('**Solved:** 2  (sat 1 / unsat 1)', out)
        # unknown 1, timeout 1, error/missing = error + missing = 2
        self.assertIn('**Unknown:** 1', out)
        self.assertIn('**Timeout:** 1', out)
        self.assertIn('**Error/missing:** 2', out)

    def test_time_stats_exclude_zero_time_rows(self):
        # only the 3 acasxu rows have time>0 (1.5, 2.5, 3.0) -> mean 2.3, max 3.0.
        # the 0-time timeout/error/missing rows must NOT drag the mean down.
        rc, out = self._run(self.CSV)
        self.assertIn('mean 2.3s', out)
        self.assertIn('max 3.0s', out)

    def test_no_time_line_when_all_zero(self):
        csv_text = (
            'subfolder,category,onnx,vnnlib,status,status_str,time_s,error_message\n'
            'x,x,a.onnx,p.vnnlib,-2,missing,0,\n'
        )
        rc, out = self._run(csv_text)
        self.assertEqual(rc, 0)
        self.assertNotIn('Time/instance', out)   # no time>0 rows -> stat suppressed

    def test_optimistic_score(self):
        rc, out = self._run(self.CSV)
        # 2 solved * 10 = 20
        self.assertIn('**20**', out)

    def test_per_benchmark_rows(self):
        rc, out = self._run(self.CSV)
        self.assertIn('| acasxu |', out)
        self.assertIn('| cersyve |', out)


if __name__ == '__main__':
    unittest.main()
