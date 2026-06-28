#!/usr/bin/env python3
"""Sound batched-IBP decider for nn4sys lindex / lindex_deep (VNN-COMP 2026 early route).

The lindex spec is one big OR over MANY ultra-thin X-boxes (width ~8e-8), each disjunct
  (and (>= X[0,0] lo) (<= X[0,0] hi) (<= Y[0,0] a))   OR   (... (>= Y[0,0] b))
i.e. for X in [lo,hi] the UNSAFE region is Y<=a or Y>=b (safe band = (a,b)). The whole instance is
SAT iff ANY box has an X whose Y is unsafe; UNSAT iff EVERY box keeps Y strictly inside its band.

The net is a tiny FC+ReLU stack (1->128->128->128->1; deep = 5 layers). Over a box this thin, sound
Interval-Bound-Propagation is essentially tight. NNV's general star/LP reach instead walks the 200-4000
OR clauses one at a time and TIMES OUT (lindex_200 ran 444s). This does batched IBP (a sound
over-approximation) in milliseconds, plus a concrete onnxruntime witness search for any non-safe box.

SOUND-OR-UNKNOWN: prints exactly one of
  unsat   -- EVERY box's IBP output interval is provably inside its safe band  (sound proof)
  sat     -- a concrete onnxruntime witness lands in an unsafe region          (real counterexample)
  unknown -- some box's IBP interval touches unsafe but no witness found (IBP too loose) / parse fail
Never emits unsat/sat without the proof/witness above. A self-check asserts IBP contains the true ORT
forward at sampled points (a wrong-direction interval would abort to unknown).
"""
import sys, re, os, gzip, tempfile, atexit
import numpy as np

# CRASH-SAFE exit codes: NEVER reuse 1 (an uncaught Python exception exits 1) for a verdict, so an
# unexpected crash maps to fail-open unknown in the MATLAB caller, never a wrong unsat (-150).
EXIT_SAT, EXIT_UNSAT, EXIT_UNKNOWN = 10, 11, 12

# .gz inputs are extracted to temp files; remove them on exit so an eval box running many instances
# does not leak /tmp.
_TMPFILES = []


@atexit.register
def _cleanup_tmpfiles():
    for p in _TMPFILES:
        try:
            os.unlink(p)
        except OSError:
            pass


def _gunzip(p):
    if os.path.exists(p):
        return p
    if os.path.exists(p + ".gz"):
        with gzip.open(p + ".gz", "rb") as f:
            data = f.read()
        t = tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(p)[1])
        t.write(data); t.close()
        _TMPFILES.append(t.name)
        return t.name
    return p


def cant(msg):
    print("unknown")
    print("CANT " + msg, file=sys.stderr)
    sys.exit(EXIT_UNKNOWN)


def load_fc_layers(onnx_path):
    """Return an ordered list of ('gemm', W[out,in], b[out]) / ('relu',) for a pure Gemm/Relu net."""
    import onnx
    from onnx import numpy_helper
    m = onnx.load(onnx_path)
    init = {i.name: numpy_helper.to_array(i) for i in m.graph.initializer}
    layers = []
    for n in m.graph.node:
        if n.op_type == "Gemm":
            A = init[n.input[1]].astype(np.float64)
            b = (init[n.input[2]].astype(np.float64) if len(n.input) > 2
                 else np.zeros(A.shape[0]))
            attr = {a.name: a for a in n.attribute}
            transB = bool(attr["transB"].i) if "transB" in attr else False
            # NNV/torch export uses y = x @ W^T + b (transB=1), so W=[out,in] is A as-is; else transpose.
            W = A if transB else A.T
            alpha = attr["alpha"].f if "alpha" in attr else 1.0
            beta = attr["beta"].f if "beta" in attr else 1.0
            layers.append(("gemm", alpha * W, beta * b))
        elif n.op_type == "Relu":
            layers.append(("relu",))
        elif n.op_type in ("Flatten", "Identity", "Reshape"):
            continue  # shape-only on a [N,1] feature vector -- no-op for IBP
        else:
            cant("unsupported op for IBP: " + n.op_type)
    if not any(L[0] == "gemm" for L in layers):
        cant("no Gemm layers found")
    return layers


def ibp(layers, lb, ub):
    """Sound interval propagation. lb,ub: [N,d_in]; returns [N,d_out] interval bounds."""
    for L in layers:
        if L[0] == "gemm":
            W, b = L[1], L[2]                 # y = x @ W^T + b, W=[out,in]
            Wp = np.maximum(W, 0.0); Wn = np.minimum(W, 0.0)
            nlb = lb @ Wp.T + ub @ Wn.T + b
            nub = ub @ Wp.T + lb @ Wn.T + b
            lb, ub = nlb, nub
        else:                                  # relu (monotone -> interval is [relu(lb),relu(ub)])
            lb, ub = np.maximum(lb, 0.0), np.maximum(ub, 0.0)
    return lb, ub


def parse_boxes(vnn_path):
    """Return boxes [(lo,hi)], a[N], b[N] (unsafe = Y<=a or Y>=b; defaults open)."""
    with open(vnn_path, encoding="utf-8") as f:
        txt = f.read()
    # Accept BOTH vnnlib index conventions: the 2.0 `X[0,0]`/`Y[0,0]` form AND the 1.0 `X_0`/`Y_0`
    # form. nn4sys is a CARRY-OVER cat whose competition vnnlib is served from the 1.0/ dir (X_0),
    # but this parser was written only for X[0,0] -> it silently parsed 0 sub-queries -> CANT ->
    # every lindex instance fell to unknown (the +24 lindex recovery, PR #426, was lost). Widening
    # the index pattern is purely a PARSE change; the sound IBP+witness decision logic is unchanged.
    pat = (r"\(and\s*\(>=\s*X(?:\[0,\s*0\]|_0)\s*([-\d.eE]+)\)\s*\(<=\s*X(?:\[0,\s*0\]|_0)\s*([-\d.eE]+)\)\s*"
           r"\((<=|>=)\s*Y(?:\[0,\s*0\]|_0)\s*([-\d.eE]+)\)")
    from collections import OrderedDict
    g = OrderedDict()
    for lo, hi, sense, bnd in re.findall(pat, txt):
        k = (float(lo), float(hi))
        g.setdefault(k, {})
        g[k]["a" if sense == "<=" else "b"] = float(bnd)
    if not g:
        cant("no lindex sub-queries parsed")
    boxes = list(g.keys())
    a = np.array([g[k].get("a", -np.inf) for k in boxes])
    b = np.array([g[k].get("b", np.inf) for k in boxes])
    return boxes, a, b


def main():
    if len(sys.argv) < 3:
        cant("usage: nn4sys_lindex_decide.py <onnx> <vnnlib>")
    onnx_p = _gunzip(sys.argv[1]); vnn_p = _gunzip(sys.argv[2])
    try:
        layers = load_fc_layers(onnx_p)
        boxes, a, b = parse_boxes(vnn_p)
    except SystemExit:
        raise
    except Exception as e:
        cant("setup failed: " + str(e))

    los = np.array([lo for lo, _ in boxes], dtype=np.float64)
    his = np.array([hi for _, hi in boxes], dtype=np.float64)

    def subdiv_bounds(idx, K):
        """Sound output interval per box (worst case over K equal sub-boxes). IBP is exact at a point,
        so subdividing the thin box shrinks IBP's correlation-looseness ~linearly -> tight + sound."""
        lo = los[idx]; hi = his[idx]; n = len(idx)
        Ylb = np.full(n, np.inf); Yub = np.full(n, -np.inf)
        e = np.linspace(0.0, 1.0, K + 1)
        for j in range(K):
            sl = (lo + (hi - lo) * e[j]).reshape(-1, 1)
            su = (lo + (hi - lo) * e[j + 1]).reshape(-1, 1)
            ylb, yub = ibp(layers, sl, su)
            Ylb = np.minimum(Ylb, ylb.ravel()); Yub = np.maximum(Yub, yub.ravel())
        return Ylb, Yub

    # ADAPTIVE REFINEMENT: certify all boxes at the smallest K; re-subdivide ONLY the still-unsafe boxes
    # at higher K (keeps memory bounded -- only the few band-boundary boxes need depth). All sound.
    pending = np.arange(len(boxes))
    YlbAll = np.full(len(boxes), np.nan); YubAll = np.full(len(boxes), np.nan)
    for K in (1, 16, 64, 256, 1024, 4096):
        if pending.size == 0:
            break
        ylb, yub = subdiv_bounds(pending, K)
        YlbAll[pending] = ylb; YubAll[pending] = yub
        safe = (ylb > a[pending]) & (yub < b[pending])
        pending = pending[~safe]

    # SOUNDNESS SELF-CHECK: the sound interval MUST contain the true onnxruntime forward at sampled points.
    try:
        import onnxruntime as ort
        # pin CPU (like the other VNN_COMP2026 ORT scripts): the session backs the soundness self-check +
        # witness validation, so behaviour must not depend on the local ORT GPU build, and it must not
        # contend for the GPU during a sweep.
        s = ort.InferenceSession(onnx_p, providers=["CPUExecutionProvider"])
        iname = s.get_inputs()[0].name; ish = s.get_inputs()[0].shape
        shp = tuple(1 if (d is None or isinstance(d, str)) else d for d in ish)

        def fwd(x):
            xv = np.array([[x]], dtype=np.float32).reshape(shp)
            return float(np.array(s.run(None, {iname: xv})[0]).ravel()[0])
        chk = np.random.default_rng(0).choice(len(boxes), size=min(50, len(boxes)), replace=False)
        for i in chk:
            for x in (boxes[i][0], boxes[i][1], 0.5 * (boxes[i][0] + boxes[i][1])):
                y = fwd(x)
                if y < YlbAll[i] - 1e-4 or y > YubAll[i] + 1e-4:
                    cant("soundness self-check FAILED (interval does not contain forward) -> refuse verdict")
    except SystemExit:
        raise
    except Exception as e:
        cant("onnxruntime self-check unavailable: " + str(e))

    if pending.size == 0:
        print("unsat"); sys.exit(EXIT_UNSAT)

    # boxes still un-certified at K_max -> hunt a concrete witness (densely) before claiming anything
    for i in pending:
        lo, hi = boxes[i]
        for x in np.linspace(lo, hi, 64):
            y = fwd(float(x))
            if y <= a[i] or y >= b[i]:
                print("sat"); print("(X_0 %.17g)" % float(x)); print("(Y_0 %.17g)" % y); sys.exit(EXIT_SAT)
    print("unknown"); sys.exit(EXIT_UNKNOWN)


if __name__ == "__main__":
    main()
