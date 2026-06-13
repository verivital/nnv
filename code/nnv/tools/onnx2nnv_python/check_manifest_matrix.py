"""check_manifest_matrix.py — end-to-end layout-model regression matrix for
the ONNX -> NNV manifest importer.

For each benchmark model:
  1. preprocess the ONNX exactly like onnx2nnv.py (opset upgrade, onnxsim,
     onnxoptimizer, dead-node elimination),
  2. emit the .nnv.mat manifest (onnx2nnv.walk + manifest_to_mat),
  3. run onnxruntime on the SAME preprocessed graph with every intermediate
     exposed as an output,
  4. execute the manifest with manifest_sim.ManifestSim (a numpy replica of
     the MATLAB load_nnv_from_mat + layer evaluate semantics),
  5. compare every mapped tensor ORIENTATION-STRICTLY (per-layer MlLayout
     convention; no permutation guessing).

Gate: max|diff| < 1e-4 * max(1, ||ref||_inf) on every clean mapped tensor,
< 1e-3 on the final output. The scale factor exists ONLY because the ground
truth itself is float32 (onnxruntime CPU has no float64 Conv kernel): an
IBP-trained ViT carries ~5e2-magnitude residual activations whose float32
ULP is ~6e-5, so a few e-4 of ABSOLUTE noise is the reference's own
rounding (~2e-7 relative), demonstrably not an orientation error (layout
bugs produce O(||x||) errors). Shape/orientation checks remain STRICT —
any mismatch is reported as inf.

Layers wrapping ops MATLAB refuses to evaluate (UnsupportedOp:*, Resize)
are oracle-patched from onnxruntime and tracked separately, as are tensors
transitively downstream of a patch ("tainted") — those mirror what the
MATLAB xval also cannot check (collins).

Usage:
  python3 check_manifest_matrix.py                # full matrix
  python3 check_manifest_matrix.py vit_pgd cgan   # subset
  VNNCOMP_BENCH=<path> overrides the benchmark root (default: the
  vnncomp2026_benchmarks clone next to the nnv repo).
"""

import os
import sys
import tempfile

import numpy as np
import onnx
import onnxruntime as ort

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import onnx2nnv
from manifest_sim import ManifestSim, ort_to_matlab, strict_diff

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_BENCH = os.path.normpath(os.path.join(
    _HERE, '..', '..', '..', '..', '..', 'vnncomp2026_benchmarks', 'benchmarks'))
BENCH = os.environ.get('VNNCOMP_BENCH', _DEFAULT_BENCH)

MODELS = {
    'vit_pgd': 'vit_2023/1.0/onnx/pgd_2_3_16.onnx',
    'vit_ibp': 'vit_2023/1.0/onnx/ibp_3_3_8.onnx',
    'soundnessbench': 'soundnessbench/1.0/onnx/model.onnx',
    'cgan': 'cgan_2023/1.0/onnx/cGAN_imgSz32_nCh_1.onnx',
    'lsnc': 'lsnc_relu/1.0/onnx/relu_quadrotor2d_state.onnx',
    'traffic': ('traffic_signs_recognition_2023/1.0/onnx/'
                '3_30_30_QConv_16_3_QConv_32_2_Dense_43_ep_30.onnx'),
    'collins': 'collins_aerospace_benchmark/1.0/onnx/yolov5nano_LRelu_640.onnx',
}

TOL_TENSOR = 1e-4
TOL_FINAL = 1e-3
SEED = 20260612


def _maybe_gunzip(path):
    if os.path.exists(path):
        return path
    gz = path + '.gz'
    if os.path.exists(gz):
        import gzip, shutil
        with gzip.open(gz, 'rb') as f, open(path, 'wb') as g:
            shutil.copyfileobj(f, g)
        return path
    raise FileNotFoundError(path)


def _ort_session_with_intermediates(model, wanted):
    """Build an ORT session on `model` with `wanted` tensor names added as
    graph outputs (only those that have value_info / are already outputs)."""
    m = onnx.ModelProto()
    m.CopyFrom(model)
    have = {o.name for o in m.graph.output}
    vi_map = {vi.name: vi for vi in m.graph.value_info}
    for name in wanted:
        if name in have or name not in vi_map:
            continue
        m.graph.output.append(vi_map[name])
        have.add(name)
    sess = ort.InferenceSession(m.SerializeToString(),
                                providers=['CPUExecutionProvider'])
    return sess, [o.name for o in sess.get_outputs()]


def check_model(key, onnx_rel, verbose=False):
    path = _maybe_gunzip(os.path.join(BENCH, onnx_rel))
    model = onnx.load(path)
    model = onnx2nnv.preprocess_onnx(model)

    layers, weights, inp_name, inp_shape, out_name = onnx2nnv.walk(model)
    tmpdir = tempfile.mkdtemp(prefix=f'nnvman_{key}_')
    mat_path = os.path.join(tmpdir, key + '.nnv.mat')
    onnx2nnv.manifest_to_mat(layers, weights, inp_name, inp_shape, out_name,
                             mat_path, 0, path)

    sim = ManifestSim(mat_path)

    # tensors we want from ORT: every real (graph) tensor a manifest layer
    # produces. Synthetic decomposition outputs ('*_out', '*_flat') are
    # internal and not present in the ONNX graph.
    graph_tensors = set()
    for n in model.graph.node:
        graph_tensors.update(n.output)
    graph_tensors.update(o.name for o in model.graph.output)
    layer_out = {}
    for L in sim.layers:
        if L.outputs and L.outputs[0] in graph_tensors:
            layer_out[L.name] = L.outputs[0]

    sess, out_names = _ort_session_with_intermediates(model, set(layer_out.values()))
    rng = np.random.default_rng(SEED)
    x = rng.uniform(-1.0, 1.0, size=[d if d > 0 else 1 for d in inp_shape])
    x = x.astype(np.float32)
    ort_vals = dict(zip(out_names, sess.run(None, {inp_name: x})))

    def oracle(L):
        t = L.outputs[0] if L.outputs else None
        if t is None or t not in ort_vals:
            return None
        return ort_to_matlab(ort_vals[t], L.layout)

    values, patched, tainted = sim.run(x, oracle=oracle)

    rows = []
    clean_max, clean_scaled_max, taint_max = 0.0, 0.0, 0.0
    n_clean = n_taint = 0
    final_diff = None
    for L in sim.layers:
        nm = L.name
        if nm not in layer_out or nm in patched:
            continue
        t = layer_out[nm]
        if t not in ort_vals:
            continue
        d = strict_diff(values[nm], ort_vals[t], L.layout)
        scale = max(1.0, float(np.max(np.abs(ort_vals[t]))) if np.asarray(ort_vals[t]).size else 1.0)
        d_scaled = d / scale
        is_taint = nm in tainted
        rows.append((nm, t, L.layout, d, d_scaled, is_taint))
        if t == out_name:
            final_diff = d
        if is_taint:
            taint_max = max(taint_max, d); n_taint += 1
        else:
            clean_max = max(clean_max, d)
            clean_scaled_max = max(clean_scaled_max, d_scaled)
            n_clean += 1
        if verbose:
            flag = ' TAINTED' if is_taint else ''
            mark = '  <<<<' if (d_scaled > TOL_TENSOR and not is_taint) else ''
            print(f'  {nm[:34]:34s} {L.layout:4s} {str(t)[:42]:42s} {d:9.2e}{flag}{mark}')

    ok = clean_scaled_max < TOL_TENSOR and (final_diff is None or final_diff < TOL_FINAL
                                            or out_name in {layer_out.get(n) for n in tainted})
    return {
        'key': key, 'ok': ok, 'clean_max': clean_max,
        'clean_scaled_max': clean_scaled_max, 'n_clean': n_clean,
        'taint_max': taint_max, 'n_taint': n_taint, 'n_patched': len(patched),
        'final': final_diff, 'n_layers': len(sim.layers), 'mat': mat_path,
        'rows': rows,
    }


def main(argv):
    verbose = '-v' in argv
    keys = [a for a in argv if not a.startswith('-')] or list(MODELS)
    results = []
    for k in keys:
        print(f'=== {k} ({MODELS[k]})', flush=True)
        try:
            r = check_model(k, MODELS[k], verbose=verbose)
        except Exception as e:
            print(f'  ERROR: {type(e).__name__}: {e}')
            results.append({'key': k, 'ok': False, 'error': str(e)})
            continue
        f = 'n/a' if r['final'] is None else f"{r['final']:.2e}"
        print(f"  layers={r['n_layers']} compared={r['n_clean']} clean_max={r['clean_max']:.2e} "
              f"(scaled {r['clean_scaled_max']:.2e}) "
              f"final={f} patched={r['n_patched']} tainted_cmp={r['n_taint']}"
              + ("" if r['n_taint'] == 0 else " taint_max={:.2e}".format(r['taint_max']))
              + "  -> " + ("PASS" if r['ok'] else "FAIL"), flush=True)
        results.append(r)
    n_bad = sum(1 for r in results if not r.get('ok'))
    print(f"\n{len(results) - n_bad}/{len(results)} models PASS")
    return 1 if n_bad else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
