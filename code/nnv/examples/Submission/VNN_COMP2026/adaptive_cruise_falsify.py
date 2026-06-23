#!/usr/bin/env python3
"""Sound nonlinear-property falsifier for adaptive_cruise_control_non_linear_2026.

NNV gates nonlinear vnnlib (X*Y, X*X) as unsupported, so this shells-out path finds a SAT
witness by sampling the input box + forwarding the REAL onnx via onnxruntime, then accepting
a witness ONLY if it satisfies EVERY assertion under the AUTHORITATIVE vnnlib parser (the
`vnnlib` pypi package: parse_query_file + to_dnf for boolean structure, exact arith eval on
the parsed AST) with a strict margin. SAT-or-unknown: never emits unsat (sampling can't prove
it), never emits a witness the authoritative parser doesn't confirm -> the only -150 surface
(the property semantics) is delegated to the official parser.

Exit codes (mirrors cctsdb_enumerate.py protocol): 10=SAT (+ witness csv), 12=unknown, 1=error.
Usage: adaptive_cruise_falsify.py <onnx> <vnnlib> <witness_out_csv> [--n N] [--seed S] [--tol T]
"""
import sys, argparse
import numpy as np

def eprint(*a): print(*a, file=sys.stderr)

def ev_arith(e, X, Y):
    t = type(e).__name__
    if t == 'Var':
        idx = list(e.indices); flat = idx[-1] if idx else 0
        return float(X[flat]) if str(e.kind).endswith('Input') else float(Y[flat])
    if t in ('Float', 'Int', 'IntExpr'): return float(e.value)
    if t == 'Literal': return float(e.lexeme)
    if t == 'Negate': return -ev_arith(e.expr, X, Y)
    if t == 'Plus': return sum(ev_arith(a, X, Y) for a in e.args)
    if t == 'Minus': return ev_arith(e.head, X, Y) - sum(ev_arith(a, X, Y) for a in e.rest)
    if t == 'Multiply':
        r = 1.0
        for a in e.args: r *= ev_arith(a, X, Y)
        return r
    raise ValueError("arith node " + t)

def ev_cmp(c, X, Y, tol=0.0):
    """Evaluate a single Comparison with a strict margin `tol` on inequalities."""
    t = type(c).__name__
    lo = ev_arith(c.lhs, X, Y); ro = ev_arith(c.rhs, X, Y)
    if t == 'GreaterEqual': return lo >= ro + tol
    if t == 'GreaterThan':  return lo >  ro + tol
    if t == 'LessEqual':    return lo <= ro - tol
    if t == 'LessThan':     return lo <  ro - tol
    if t == 'Equal':        return lo == ro
    if t == 'NotEqual':     return lo != ro
    raise ValueError("cmp node " + t)

def holds_dnf(dnf, X, Y, tol=0.0):
    """dnf = List[List[Comparison]] (precomputed). Holds iff some clause has ALL comparisons true."""
    return any(all(ev_cmp(c, X, Y, tol) for c in clause) for clause in dnf)

def split_assertions(query):
    """Partition assertions into input-only (filters/box) vs output-involving."""
    import vnnlib
    def uses_output(e):
        if type(e).__name__ == 'Var': return str(e.kind).endswith('Output')
        for ch in e.children():
            if uses_output(ch): return True
        return False
    inp, out = [], []
    for a in query.assertions:
        (out if uses_output(a.expr) else inp).append(a.expr)
    return inp, out

def input_box(query, nin):
    """Tight box bounds per input feature from simple Var>=c / Var<=c comparisons (recursively)."""
    lb = np.full(nin, -np.inf); ub = np.full(nin, np.inf)
    def walk(e):
        t = type(e).__name__
        if t in ('And', 'Or'):
            for a in e.args: walk(a)
        elif t in ('GreaterEqual', 'GreaterThan', 'LessEqual', 'LessThan'):
            L, R = e.lhs, e.rhs
            # ONLY simple bounds on INPUT vars (an output Var can share a feature index with an input)
            if type(L).__name__ == 'Var' and str(L.kind).endswith('Input') \
                    and type(R).__name__ in ('Float', 'Int', 'Literal'):
                i = list(L.indices)[-1]; v = ev_arith(R, [0]*nin, [0]*nin)
                if t.startswith('Greater'): lb[i] = max(lb[i], v)
                else: ub[i] = min(ub[i], v)
    for a in query.assertions: walk(a.expr)
    return lb, ub

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('onnx'); ap.add_argument('vnnlib'); ap.add_argument('witness_csv')
    ap.add_argument('--n', type=int, default=1_000_000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--tol', type=float, default=1e-2)     # ROBUST margin: excludes thin saturation
    # artifacts (the net saturates ~100.0011, only 1.2e-4 above the 100.001 threshold -> a fp-noise-grade
    # "violation" that could flip under the official checker). 1e-2 is ~100x a typical checker tolerance, so
    # an accepted witness robustly clears the boundary and stays a real counterexample under any runtime.
    args = ap.parse_args()
    try:
        import vnnlib, onnxruntime as ort
    except Exception as e:
        eprint("missing vnnlib/onnxruntime ->", e); sys.exit(12)   # cannot check -> unknown (sound)
    try:
        q = vnnlib.parse_query_file(args.vnnlib)
    except Exception as e:
        eprint("vnnlib parse failed ->", e); sys.exit(12)
    # input dim from the network declaration
    net = q.networks[0]
    nin = int(np.prod(net.inputs[0].shape))
    nout = int(np.prod(net.outputs[0].shape))
    inp_asserts, out_asserts = split_assertions(q)
    lb, ub = input_box(q, nin)
    if not (np.all(np.isfinite(lb)) and np.all(np.isfinite(ub)) and np.all(ub >= lb)):
        eprint("no finite input box -> unknown"); sys.exit(12)
    # Prefer a dynamic-batch session (one-time onnx surgery: batch dim -> symbolic) so we can forward
    # millions of samples in one call; safe for a plain Gemm/Relu FFNN. Fall back to per-sample if that
    # fails (correctness preserved either way -- same real ONNX).
    sess = ort.InferenceSession(args.onnx)
    iname = sess.get_inputs()[0].name; ishape = sess.get_inputs()[0].shape
    shp = tuple(1 if (x is None or isinstance(x, str)) else x for x in ishape)
    batched_ok = True
    try:
        import onnx as _onnx
        _m = _onnx.load(args.onnx)
        _m.graph.input[0].type.tensor_type.shape.dim[0].dim_param = 'N'
        if _m.graph.output:
            _m.graph.output[0].type.tensor_type.shape.dim[0].dim_param = 'N'
        sess = ort.InferenceSession(_m.SerializeToString())
    except Exception as e:
        eprint("dynamic-batch surgery failed, per-sample fallback:", e); batched_ok = False
    rng = np.random.default_rng(args.seed)

    def forward(Xb):
        Xb = Xb.astype(np.float32)
        if batched_ok:
            try:
                return sess.run(None, {iname: Xb.reshape((-1,) + shp[1:])})[0].reshape(Xb.shape[0], -1)
            except Exception:
                pass
        return np.vstack([sess.run(None, {iname: Xb[r].reshape(shp)})[0].reshape(1, -1)
                          for r in range(Xb.shape[0])])

    inp_dnf = [e.to_dnf() for e in inp_asserts]      # precompute once (to_dnf is expensive)
    out_dnf = [e.to_dnf() for e in out_asserts]
    Yz = [0.0] * nout

    def all_hold(X, Y, tol):
        return all(holds_dnf(d, X, Y, tol) for d in inp_dnf) and \
               all(holds_dnf(d, X, Y, tol) for d in out_dnf)

    found = None
    BS = 200_000
    done = 0
    while done < args.n and found is None:
        m = min(BS, args.n - done); done += m
        X = lb + (ub - lb) * rng.random((m, nin))
        # keep only samples passing the (possibly nonlinear) input-only assertions
        keep = np.array([all(holds_dnf(d, X[r], Yz, 0.0) for d in inp_dnf) for r in range(m)])
        if not keep.any(): continue
        Xk = X[keep]; Yk = forward(Xk)
        for r in range(Xk.shape[0]):
            if all(holds_dnf(d, Xk[r], Yk[r], args.tol) for d in out_dnf):
                found = (Xk[r].copy(), Yk[r].copy()); break

    if found is None:
        eprint("no strict witness in %d samples -> unknown" % args.n); sys.exit(12)

    Xw, _ = found
    # AUTHORITATIVE re-validation: fresh ORT forward + every assertion must hold (strict margin)
    Yw = forward(Xw.reshape(1, -1))[0]
    if not all_hold(Xw, Yw, args.tol):
        eprint("witness failed authoritative re-validation -> unknown"); sys.exit(12)
    # also require it to hold at tol=0 trivially (sanity) and report
    np.savetxt(args.witness_csv, Xw.reshape(-1, 1), fmt='%.10g')
    print("SAT " + " ".join("%.10g" % v for v in Yw))
    sys.exit(10)

if __name__ == '__main__':
    main()
