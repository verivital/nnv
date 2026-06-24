#!/usr/bin/env python3
"""Authoritative SAT-witness gate: re-validate an emitted `sat` against the AUTHORITATIVE vnnlib
parser (`vnnlib` pypi) + a real-onnx onnxruntime forward, before the verdict is trusted.

This is the inline, per-witness form of status-repo tools/vnnlib_audit.py (whose flagging arm was
validated on T1 out-of-box + T_PRIMARY in-box/output-safe). It reads the EXACT `(X_i val)` result-file
format the runner's write_counterexample emits (the same format the auditor was validated on), so the
flat<->ravel order matching is the validated path.

Usage:  validate_witness_authoritative.py <onnx[.gz]> <vnnlib[.gz]> <result_file> [--in-tol T] [--out-tol T]
Exit:   0 = VALID (every assertion holds; emit sat)         <- the witness is a real counterexample
        1 = SPURIOUS (some assertion fails; downgrade->unknown)  <- a would-be -150, converted to 0 points
        2 = CANNOT-CHECK (missing deps / parse fail / shape mismatch / non-sat file) -> FAIL OPEN
                          (caller keeps today's verdict; the gate never makes things worse)

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


def read_witness(resfile):
    """Parse the runner's sat result file: 'sat' then (X_i val) lines. Returns X array or None."""
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


def ev_arith(e, X, Y):
    t = type(e).__name__
    if t == "Var":
        flat = _flat(e)
        arr = X if str(e.kind).endswith("Input") else Y
        return float(arr[flat])
    if t in ("Float", "Int", "IntExpr"):
        return float(e.value)
    if t == "Literal":
        return float(e.lexeme)
    if t == "Negate":
        return -ev_arith(e.expr, X, Y)
    if t == "Plus":
        return sum(ev_arith(a, X, Y) for a in e.args)
    if t == "Minus":
        return ev_arith(e.head, X, Y) - sum(ev_arith(a, X, Y) for a in e.rest)
    if t == "Multiply":
        r = 1.0
        for a in e.args:
            r *= ev_arith(a, X, Y)
        return r
    raise ValueError("arith " + t)


def ev_cmp(c, X, Y, in_tol, out_tol):
    """A comparison holds if it is satisfied within the LENIENT tolerance (so a real boundary
    witness is not falsely downgraded). out_tol is used (>= the official checker tol)."""
    t = type(c).__name__
    lo = ev_arith(c.lhs, X, Y); ro = ev_arith(c.rhs, X, Y)
    tol = out_tol
    if t == "GreaterEqual": return lo >= ro - tol
    if t == "GreaterThan":  return lo >  ro - tol
    if t == "LessEqual":    return lo <= ro + tol
    if t == "LessThan":     return lo <  ro + tol
    if t == "Equal":        return abs(lo - ro) <= tol
    if t == "NotEqual":     return abs(lo - ro) >  0.0   # lenient: any non-equality holds
    raise ValueError("cmp " + t)


def holds(dnf, X, Y, in_tol, out_tol):
    return any(all(ev_cmp(c, X, Y, in_tol, out_tol) for c in clause) for clause in dnf)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("onnx"); ap.add_argument("vnnlib"); ap.add_argument("resultfile")
    ap.add_argument("--in-tol", type=float, default=1e-4)
    ap.add_argument("--out-tol", type=float, default=1e-4)   # >= official VNN-COMP CE-checker tol
    a = ap.parse_args()
    try:
        import vnnlib
        import onnxruntime as ort
    except Exception as e:
        cant("missing vnnlib/onnxruntime: " + str(e))
    X = read_witness(a.resultfile)
    if X is None:
        cant("result file is not a parseable sat witness")
    if isinstance(X, str):
        cant("non-contiguous witness indices")
    onnx_p = _maybe_gunzip_path(a.onnx)
    vnnlib_p = _maybe_gunzip_path(a.vnnlib)
    try:
        q = vnnlib.parse_query_file(vnnlib_p)
        dnf = [assn.expr.to_dnf() for assn in q.assertions]
    except Exception as e:
        cant("vnnlib parse failed: " + str(e))
    try:
        s = ort.InferenceSession(onnx_p)
        iname = s.get_inputs()[0].name
        ish = s.get_inputs()[0].shape
        shp = tuple(1 if (d is None or isinstance(d, str)) else d for d in ish)
        if int(np.prod(shp)) != X.size:
            shp = (1, X.size)
        Y = s.run(None, {iname: X.reshape(shp).astype(np.float32)})[0].ravel()
    except Exception as e:
        cant("onnx forward failed: " + str(e))
    try:
        ok = all(holds(d, X, Y, a.in_tol, a.out_tol) for d in dnf)
    except Exception as e:
        cant("assertion eval failed: " + str(e))
    if ok:
        print("VALID")
        sys.exit(EXIT_VALID)
    print("SPURIOUS")
    sys.exit(EXIT_SPURIOUS)


if __name__ == "__main__":
    main()
