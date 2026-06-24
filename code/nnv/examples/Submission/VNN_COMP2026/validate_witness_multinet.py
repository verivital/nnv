#!/usr/bin/env python3
"""Authoritative SAT-witness gate for VNN-LIB 2.0 TWO-NETWORK specs (the `equal-to` / monotonicity
family: isomorphic_acasxu, monotonic_acasxu) -- the multi-network sibling of
validate_witness_authoritative.py.

WHY: verify_multinet emits a STACKED [x_f; x_g] witness for a two-network spec. The single-
InferenceSession gate (validate_witness_authoritative.py) opens ONE network and cannot replay two
coupled inputs/outputs, so run_vnncomp_instance used to BLANKET-downgrade every multinet sat to
unknown. That is sound but OVER-conservative: it discards REAL counterexamples (verified 2026-06-24:
a monotonic_acasxu witness replayed through onnxruntime matches NNV's forward to 1e-9 and genuinely
violates Y_f[3] < Y_g[3]). This gate forwards BOTH halves through the network(s) via onnxruntime and
checks the LITERAL two-network vnnlib on the witness, so a real sat is TRUSTED and only a
genuinely-spurious one is downgraded.

Usage: validate_witness_multinet.py <onnx_f[.gz]> <vnnlib[.gz]> <result_file> [--onnx-g G] [--tol T]
Exit:  0 = VALID     (every assertion holds on the ORT replay -> emit sat)
       1 = SPURIOUS  (some assertion fails -> downgrade to unknown; a would-be -150 -> 0 points)
       2 = CANNOT-CHECK (missing deps / parse fail / shape mismatch / non-sat / not 2 nets) -> FAIL OPEN

SOUND DIRECTION: the gate can only ever turn sat->unknown (lose +10), never manufacture a verdict.
"""
import sys, os, re, gzip, argparse
import numpy as np

EXIT_VALID, EXIT_SPURIOUS, EXIT_CANT = 0, 1, 2


def cant(msg):
    print("CANNOT-CHECK " + msg, file=sys.stderr)
    sys.exit(EXIT_CANT)


def _maybe_gunzip_path(p):
    if os.path.exists(p):
        return p
    if os.path.exists(p + ".gz"):
        import tempfile
        with gzip.open(p + ".gz", "rb") as f:
            data = f.read()
        t = tempfile.NamedTemporaryFile(delete=False, suffix=os.path.splitext(p)[1])
        t.write(data); t.close()
        return t.name
    return p


def read_stacked_witness(resfile):
    """Parse the runner's sat result file: 'sat' then (X_i val) lines, the STACKED [x_f; x_g] input
    written by write_counterexample. Returns a contiguous X array, None (not sat), or 'noncontig'."""
    txt = open(resfile).read()
    if not txt.lstrip().lower().startswith("sat"):
        return None
    fnum = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
    xs = {}
    for m in re.finditer(r"\(\s*X_(\d+)\s+(" + fnum + r")\s*\)", txt):
        xs[int(m.group(1))] = float(m.group(2))
    if not xs:
        return None
    n = max(xs) + 1
    if len(xs) != n:
        return "noncontig"
    return np.array([xs[i] for i in range(n)], dtype=np.float64)


def _flat(e):
    idx = list(e.indices)
    try:
        shp = list(e.shape)
        if len(shp) == len(idx) and int(np.prod(shp)) > 1:
            return int(np.ravel_multi_index(idx, shp))
    except Exception:
        pass
    return idx[-1] if idx else 0


def ev_arith(e, env):
    """env maps (network_name, 'Input'|'Output', flat_index) -> float."""
    t = type(e).__name__
    if t == "Var":
        kind = "Input" if str(e.kind).endswith("Input") else "Output"
        key = (e.network_name, kind, _flat(e))
        if key not in env:
            raise ValueError("unknown var " + str(key))
        return float(env[key])
    if t in ("Float", "Int", "IntExpr"):
        return float(e.value)
    if t == "Literal":
        return float(e.lexeme)
    if t == "Negate":
        return -ev_arith(e.expr, env)
    if t == "Plus":
        return sum(ev_arith(a, env) for a in e.args)
    if t == "Minus":
        return ev_arith(e.head, env) - sum(ev_arith(a, env) for a in e.rest)
    if t == "Multiply":
        r = 1.0
        for a in e.args:
            r *= ev_arith(a, env)
        return r
    raise ValueError("arith " + t)


def ev_cmp(c, env, tol):
    """A comparison holds within the LENIENT tolerance so a real boundary witness is not falsely
    downgraded (the spec stores <=-relaxations of possibly-strict asserts; tol >= official CE tol)."""
    t = type(c).__name__
    lo = ev_arith(c.lhs, env); ro = ev_arith(c.rhs, env)
    if t == "GreaterEqual": return lo >= ro - tol
    if t == "GreaterThan":  return lo >  ro - tol
    if t == "LessEqual":    return lo <= ro + tol
    if t == "LessThan":     return lo <  ro + tol
    if t == "Equal":        return abs(lo - ro) <= tol
    if t == "NotEqual":     return abs(lo - ro) >  0.0
    raise ValueError("cmp " + t)


def holds(dnf, env, tol):
    return any(all(ev_cmp(c, env, tol) for c in clause) for clause in dnf)


def net_io_size(net, which):
    defs = net.inputs if which == "in" else net.outputs
    sz = 0
    for d in defs:
        shp = getattr(d, "shape", None)
        sz += int(np.prod(shp)) if shp else 1
    return sz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("onnx"); ap.add_argument("vnnlib"); ap.add_argument("resultfile")
    ap.add_argument("--onnx-g", default=None,
                    help="second-network onnx (for isomorphic / different-weight g; default = shared <onnx>)")
    ap.add_argument("--tol", type=float, default=1e-4)   # >= official VNN-COMP CE-checker tol
    a = ap.parse_args()
    try:
        import vnnlib
        import onnxruntime as ort
    except Exception as e:
        cant("missing vnnlib/onnxruntime: " + str(e))

    X = read_stacked_witness(a.resultfile)
    if X is None:
        cant("result file is not a parseable sat witness")
    if isinstance(X, str):
        cant("non-contiguous witness indices")

    vnnlib_p = _maybe_gunzip_path(a.vnnlib)
    try:
        q = vnnlib.parse_query_file(vnnlib_p)
        nets = list(q.networks)
        dnf = [assn.expr.to_dnf() for assn in q.assertions]
    except Exception as e:
        cant("vnnlib parse failed: " + str(e))
    if len(nets) != 2:
        cant("expected exactly 2 networks, got %d" % len(nets))

    # split the stacked witness by each network's input size (declaration order [f, g], which is the
    # order write_counterexample stacks counterEx{1}=x_f then counterEx{2}=x_g).
    insz = [net_io_size(n, "in") for n in nets]
    if sum(insz) != X.size:
        cant("stacked witness size %d != sum of network input sizes %s" % (X.size, insz))

    onnx_f = _maybe_gunzip_path(a.onnx)
    onnx_g = _maybe_gunzip_path(a.onnx_g) if a.onnx_g else onnx_f

    def net_onnx(net, is_first):
        # `equal-to` (g references f) shares the SAME onnx as f; a genuine different-weight g needs
        # --onnx-g. verify_multinet only emits sat for equal-to today, so the default (shared) is right.
        if is_first or getattr(net, "equal_to", ""):
            return onnx_f
        return onnx_g

    try:
        env = {}
        off = 0
        for ni, net in enumerate(nets):
            xi = X[off:off + insz[ni]]; off += insz[ni]
            s = ort.InferenceSession(net_onnx(net, ni == 0))
            iname = s.get_inputs()[0].name
            ish = s.get_inputs()[0].shape
            shp = tuple(1 if (d is None or isinstance(d, str)) else d for d in ish)
            if int(np.prod(shp)) != xi.size:
                shp = (1, xi.size)
            yi = s.run(None, {iname: xi.reshape(shp).astype(np.float32)})[0].ravel()
            for k in range(xi.size):
                env[(net.name, "Input", k)] = float(xi[k])
            for k in range(yi.size):
                env[(net.name, "Output", k)] = float(yi[k])
    except Exception as e:
        cant("onnx forward failed: " + str(e))

    try:
        ok = all(holds(d, env, a.tol) for d in dnf)
    except Exception as e:
        cant("assertion eval failed: " + str(e))

    if ok:
        print("VALID")
        sys.exit(EXIT_VALID)
    print("SPURIOUS")
    sys.exit(EXIT_SPURIOUS)


if __name__ == "__main__":
    main()
