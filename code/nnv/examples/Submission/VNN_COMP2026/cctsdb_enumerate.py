#!/usr/bin/env python3
"""Complete-enumeration verifier for the VNN-COMP cctsdb_yolo_2023 benchmark.

Usage:
    python cctsdb_enumerate.py <model.onnx[.gz]> <spec.vnnlib[.gz]> <out_witness.csv>
                               [--timeout S] [--margin M]

Why enumeration is SOUND AND COMPLETE for this benchmark (and ONLY under the
structural guards below, all re-verified per instance):

  * Every cctsdb_yolo_2023 spec fixes ALL inputs (lb == ub) except the two patch
    position coordinates X_12288 and X_12289, each ranging over [0, 62].
  * In the ONNX graph the ONLY consumers of those two inputs are Gather(idx 12288)
    and Gather(idx 12289), each feeding directly into Cast(to=int64) -- i.e. the
    network sees only trunc(x). The output is therefore PIECEWISE CONSTANT on the
    unit cells [p, p+1) x [q, q+1), and evaluating one representative point per
    cell (the integer points trunc(lb)..trunc(ub) per axis, clamped into the box)
    covers the ENTIRE input set exactly. 63 x 63 = 3969 onnxruntime evaluations.
  * The output spec is a single half-space on Y_0 (e.g. "(assert (<= Y_0 0.5))",
    the UNSAFE region). Any cell with Y_0 in the unsafe region => SAT with that
    cell's representative point as witness; every cell outside it => UNSAT, a
    complete proof.

Any deviation from this structure (different free variables, non-Cast consumers,
complex asserts, subgraphs, borderline margins, timeout, any error at all) prints
UNKNOWN and exits 2 -- never a guessed verdict. VNN-COMP scores a wrong verdict
at -150; this script is sound-or-unknown by construction.

Output protocol (stdout last line / exit code):
    SAT p q y       exit 10   witness CSV written: the FULL input vector, one
                              value per line, in FLAT vnnlib order X_0..X_{N-1},
                              with (X_12288, X_12289) = the violating cell point;
                              y is the onnxruntime output Y_0 at that point.
    UNSAT ymin ymax exit 11   complete proof; ymin/ymax = grid output range.
    UNKNOWN <why>   exit 2    guard violation / borderline / timeout / error.

Margin: verdicts within --margin (default 1e-4) of the threshold are reported
UNKNOWN. onnxruntime results can differ across versions/CPUs in the last ulps;
the margin keeps a borderline instance from flipping verdict on another machine.
"""

import argparse
import gzip
import math
import os
import re
import sys
import time

import numpy as np

EXIT_SAT = 10
EXIT_UNSAT = 11
EXIT_UNKNOWN = 2

# The only free-input structure this verifier understands (guard (a)).
EXPECTED_FREE = frozenset({12288, 12289})
MAX_GRID = 100_000  # refuse absurd enumeration sizes (future-proofing guard)


def unknown(why):
    print(f"UNKNOWN {why}")
    sys.exit(EXIT_UNKNOWN)


def read_maybe_gz(path):
    """Return file bytes, transparently gunzipping *.gz.

    If a non-.gz path is missing but a sibling ``<path>.gz`` exists, read that:
    some cctsdb_yolo 2.0 specs ship gz-only while instances.csv names the plain
    path. The fallback only triggers when the literal path is absent, so it never
    changes behavior for an existing file (sound).
    """
    p = str(path)
    if p.endswith(".gz"):
        with gzip.open(p, "rb") as f:
            return f.read()
    if not os.path.exists(p) and os.path.exists(p + ".gz"):
        with gzip.open(p + ".gz", "rb") as f:
            return f.read()
    with open(p, "rb") as f:
        return f.read()


# ---------------------------------------------------------------------------
# vnnlib parsing (simple regex parse; ANY unrecognized assert => UNKNOWN)
# ---------------------------------------------------------------------------

_NUM = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
# Accept BOTH the 1.0 token `X_0`/`Y_0` and the 2.0 bracket token `X[0]`/`Y[0]`.
# The vnnlib bumped to 2.0 (declare-network/declare-input header + bracket indices) but the
# ONNX is byte-identical and the spec numerically identical to 1.0, so this is a cosmetic
# syntax change: widening the variable token recovers the category WITHOUT touching any of the
# four soundness guards below (the count guard, complete-bound-set, single-Y0, free-set check),
# which is what keeps the enumeration sound. Group numbering is unchanged (op, var, idx, val) so
# parse_vnnlib needs no other edit, and the 1.0 underscore branch still matches (no regression).
# Verified on real specs: 1.0 24593/24593 (unchanged), 2.0 0 -> 24593, parsed box byte-identical
# 1.0<->2.0, all 39 referenced 2.0 specs pass every guard. (Deliberately NOT a swap to the vnnlib
# pkg: that would accept the full grammar and DELETE the count guard -> a -150 surface.)
_ASSERT_RE = re.compile(
    r"\(\s*assert\s+\(\s*(<=|>=)\s+([XY])(?:_|\[)(\d+)\]?\s+(" + _NUM + r")\s*\)\s*\)"
)


def parse_vnnlib(text):
    """Parse the simple cctsdb vnnlib form.

    Returns (lb, ub, out_op, out_thr) where lb/ub are dicts {index: float} and
    the single output constraint is "out_op" in {'<=', '>='} applied to Y_0
    (the asserted = UNSAFE region). UNKNOWN on any non-simple construct.
    """
    text = re.sub(r";[^\n]*", "", text)  # strip comments
    n_asserts = text.count("(assert")
    matches = list(_ASSERT_RE.finditer(text))
    if len(matches) != n_asserts:
        unknown(
            f"vnnlib has {n_asserts} asserts but only {len(matches)} match the "
            "simple (assert (<=|>= V_i const)) form"
        )
    lb, ub, ycons = {}, {}, []
    for m in matches:
        op, var, idx, val = m.group(1), m.group(2), int(m.group(3)), float(m.group(4))
        if var == "X":
            if op == ">=":  # X_i >= v  -> lower bound
                lb[idx] = max(lb.get(idx, -math.inf), val)
            else:           # X_i <= v  -> upper bound
                ub[idx] = min(ub.get(idx, math.inf), val)
        else:
            ycons.append((op, idx, val))

    # guard (c): EXACTLY one output constraint, on Y_0, <= or >= a constant
    if len(ycons) != 1 or ycons[0][1] != 0:
        unknown(f"output spec is not exactly one <=/>= threshold on Y_0: {ycons}")
    out_op, _, out_thr = ycons[0]

    if not lb:
        unknown("no input bounds parsed")
    n = max(max(lb), max(ub)) + 1
    if set(lb) != set(range(n)) or set(ub) != set(range(n)):
        unknown("input box is not a complete X_0..X_{N-1} bound set")
    for i in range(n):
        if not (math.isfinite(lb[i]) and math.isfinite(ub[i]) and lb[i] <= ub[i]):
            unknown(f"invalid bounds on X_{i}: [{lb[i]}, {ub[i]}]")
    return lb, ub, out_op, out_thr


# ---------------------------------------------------------------------------
# ONNX structural guard (guard (b)): the free inputs may ONLY flow through
# Gather -> Cast(int64) (at most one intermediate Gather, i.e. 2 hops).
# ---------------------------------------------------------------------------

def _attr_i(node, name, default=None):
    for a in node.attribute:
        if a.name == name:
            return a.i
    return default


def _build_const_map(graph):
    from onnx import numpy_helper

    cmap = {t.name: numpy_helper.to_array(t) for t in graph.initializer}
    for nd in graph.node:
        if nd.op_type == "Constant":
            for a in nd.attribute:
                if a.name == "value":
                    cmap[nd.output[0]] = numpy_helper.to_array(a.t)
    return cmap


def _slice_read_set(node, cmap, dim):
    """Indices (axis-0) a 1-D Slice reads, or None if not statically resolvable.

    ONNX Slice clamping/negative-index semantics match Python slicing, so we
    evaluate via range(dim)[start:end:step]. Handles both the input form
    (opset >= 10) and the attribute form (opset < 10).
    """
    if len(node.input) >= 3:  # input form: data, starts, ends[, axes[, steps]]
        vals = []
        for k in range(1, len(node.input)):
            name = node.input[k]
            if name == "":
                vals.append(None)
                continue
            if name not in cmap:
                return None
            vals.append(cmap[name].flatten().astype(np.int64))
        while len(vals) < 4:
            vals.append(None)
        starts, ends, axes, steps = vals[0], vals[1], vals[2], vals[3]
    else:  # attribute form
        starts = ends = axes = steps = None
        for a in node.attribute:
            arr = np.array(a.ints, dtype=np.int64)
            if a.name == "starts":
                starts = arr
            elif a.name == "ends":
                ends = arr
            elif a.name == "axes":
                axes = arr
    if starts is None or ends is None or starts.size != 1 or ends.size != 1:
        return None
    if axes is not None and (axes.size != 1 or int(axes[0]) not in (0, -1)):
        return None  # data is 1-D; only axis 0 is meaningful
    step = 1
    if steps is not None:
        if steps.size != 1 or int(steps[0]) == 0:
            return None
        step = int(steps[0])
    return set(range(dim)[int(starts[0]):int(ends[0]):step])


def _paths_hit_int_cast(tname, consumers, out_names, hops):
    """True iff EVERY consumer path of tensor `tname` reaches Cast(to=int64)
    within `hops` hops, with only Gather permitted as an intermediate node."""
    from onnx import TensorProto

    if tname in out_names:
        return False  # free value would flow to a graph output un-truncated
    for c in consumers.get(tname, []):
        if c.op_type == "Cast" and _attr_i(c, "to") == TensorProto.INT64:
            continue  # this path is truncated: OK
        if (
            c.op_type == "Gather"
            and hops > 1
            and c.input[0] == tname
            and tname not in list(c.input)[1:]
        ):
            if all(
                _paths_hit_int_cast(o, consumers, out_names, hops - 1)
                for o in c.output
            ):
                continue
        return False
    return True  # all paths truncated (no consumers = dead value: also safe)


def check_onnx_structure(model, n_inputs, free_set):
    """UNKNOWN unless the graph provably consumes the free inputs only through
    Cast(int64) truncation. Returns (input_name, input_shape)."""
    import onnx
    from onnx import TensorProto

    g = model.graph

    # no control-flow subgraphs: outer-scope capture could smuggle the input past
    # this walk, so refuse them outright
    for nd in g.node:
        for a in nd.attribute:
            if a.type in (onnx.AttributeProto.GRAPH, onnx.AttributeProto.GRAPHS):
                unknown(f"graph contains a subgraph node ({nd.op_type})")

    init_names = {t.name for t in g.initializer}
    real_inputs = [i for i in g.input if i.name not in init_names]
    if len(real_inputs) != 1:
        unknown(f"expected exactly 1 graph input, found {len(real_inputs)}")
    inp = real_inputs[0]
    if inp.type.tensor_type.elem_type != TensorProto.FLOAT:
        unknown(f"input is not float32 (elem_type={inp.type.tensor_type.elem_type})")
    dims = [d.dim_value for d in inp.type.tensor_type.shape.dim]
    if any(d <= 0 for d in dims) or int(np.prod(dims)) != n_inputs:
        unknown(f"input shape {dims} does not match the {n_inputs} vnnlib inputs")

    cmap = _build_const_map(g)
    consumers = {}
    for nd in g.node:
        for i in nd.input:
            consumers.setdefault(i, []).append(nd)
    out_names = {o.name for o in g.output}

    for nd in g.node:
        if inp.name not in nd.input:
            continue
        if (
            nd.op_type == "Gather"
            and nd.input[0] == inp.name
            and inp.name not in list(nd.input)[1:]
        ):
            if _attr_i(nd, "axis", 0) != 0:
                unknown(f"Gather {nd.name} on input has axis != 0")
            idx = cmap.get(nd.input[1])
            if idx is None:
                unknown(f"Gather {nd.name} on input has non-constant indices")
            idxs = {int(v) + (n_inputs if int(v) < 0 else 0) for v in idx.flatten()}
            if idxs & free_set:
                # 2 hops: Gather -> Cast(int64), or Gather -> Gather -> Cast(int64)
                if not _paths_hit_int_cast(nd.output[0], consumers, out_names, 2):
                    unknown(
                        f"free input flows through Gather {nd.name} without "
                        "hitting Cast(int64) within 2 hops"
                    )
        elif (
            nd.op_type == "Slice"
            and nd.input[0] == inp.name
            and inp.name not in list(nd.input)[1:]
        ):
            read = _slice_read_set(nd, cmap, n_inputs)
            if read is None:
                unknown(f"Slice {nd.name} on input is not statically resolvable")
            if read & free_set:
                unknown(f"Slice {nd.name} reads the free inputs {sorted(read & free_set)}")
        else:
            unknown(f"unexpected consumer of the input tensor: {nd.op_type} {nd.name}")

    return inp.name, dims


# ---------------------------------------------------------------------------
# enumeration
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("onnx_path")
    ap.add_argument("vnnlib_path")
    ap.add_argument("out_witness_csv")
    ap.add_argument("--timeout", type=float, default=300.0,
                    help="wall-clock budget in seconds (default 300)")
    ap.add_argument("--margin", type=float, default=1e-4,
                    help="threshold margin below which the verdict is UNKNOWN")
    args = ap.parse_args()
    t0 = time.monotonic()
    deadline = t0 + args.timeout

    import onnx
    import onnxruntime as ort

    # ---- parse spec --------------------------------------------------------
    try:
        spec_text = read_maybe_gz(args.vnnlib_path).decode("utf-8", errors="replace")
    except OSError as e:
        unknown(f"cannot read vnnlib: {e}")
    lb, ub, out_op, out_thr = parse_vnnlib(spec_text)
    n = len(lb)

    # guard (a): the free variables are EXACTLY {X_12288, X_12289}
    free = {i for i in range(n) if ub[i] > lb[i]}
    if free != set(EXPECTED_FREE):
        unknown(f"free variables {sorted(free)[:10]}... != {sorted(EXPECTED_FREE)}")
    i0, i1 = sorted(free)
    # the unit-cell argument below uses trunc-toward-zero == floor, which needs
    # the free range to be nonnegative
    if lb[i0] < 0 or lb[i1] < 0:
        unknown("free-input lower bound is negative; truncation cells differ")

    # ---- guard (b): structural check on the graph --------------------------
    try:
        model_bytes = read_maybe_gz(args.onnx_path)
    except OSError as e:
        unknown(f"cannot read onnx: {e}")
    model = onnx.load_model_from_string(model_bytes)
    in_name, in_shape = check_onnx_structure(model, n, free)

    # ---- enumerate the integer grid (one cell representative each) ---------
    p_lo, p_hi = int(math.trunc(lb[i0])), int(math.trunc(ub[i0]))
    q_lo, q_hi = int(math.trunc(lb[i1])), int(math.trunc(ub[i1]))
    n_evals = (p_hi - p_lo + 1) * (q_hi - q_lo + 1)
    if n_evals > MAX_GRID:
        unknown(f"grid of {n_evals} cells exceeds the {MAX_GRID} safety cap")
    print(f"enumerating {n_evals} cells: X_{i0} in [{p_lo},{p_hi}], "
          f"X_{i1} in [{q_lo},{q_hi}]", file=sys.stderr)

    so = ort.SessionOptions()
    so.log_severity_level = 3
    sess = ort.InferenceSession(model_bytes, sess_options=so,
                                providers=["CPUExecutionProvider"])

    base = np.array([lb[i] for i in range(n)], dtype=np.float64)

    def eval_point(session, xp, xq):
        x = base.copy()
        x[i0], x[i1] = xp, xq
        out = session.run(None, {in_name: x.astype(np.float32).reshape(in_shape)})[0]
        arr = np.asarray(out).flatten()
        if arr.size < 1:
            unknown("model produced an empty output")
        y = float(arr[0])  # Y_0 = first element in flat order
        if not math.isfinite(y):
            unknown(f"non-finite output {y} at ({xp},{xq})")
        return y

    ymin = math.inf
    ymax = -math.inf
    argmin = argmax = None  # (p, q, xp, xq)
    for p in range(p_lo, p_hi + 1):
        # cell [p, p+1) representative INSIDE the box: max(lb, p); for p > trunc(lb)
        # this is p itself, and for p == trunc(lb) it is lb (which truncates to p)
        xp = max(lb[i0], float(p))
        for q in range(q_lo, q_hi + 1):
            if time.monotonic() > deadline:
                unknown(f"timeout after {time.monotonic() - t0:.1f}s "
                        f"({(p - p_lo) * (q_hi - q_lo + 1) + (q - q_lo)}/{n_evals} cells)")
            xq = max(lb[i1], float(q))
            y = eval_point(sess, xp, xq)
            if y < ymin:
                ymin, argmin = y, (p, q, xp, xq)
            if y > ymax:
                ymax, argmax = y, (p, q, xp, xq)

    print(f"grid done in {time.monotonic() - t0:.2f}s: "
          f"Y_0 in [{ymin:.17g}, {ymax:.17g}], unsafe region Y_0 {out_op} {out_thr:g}",
          file=sys.stderr)

    # ---- verdict (asserted region = UNSAFE; entering it = SAT) -------------
    m = args.margin
    if out_op == "<=":
        viol = argmin
        clearly_sat = ymin <= out_thr - m
        clearly_unsat = ymin > out_thr + m
    else:  # '>='
        viol = argmax
        clearly_sat = ymax >= out_thr + m
        clearly_unsat = ymax < out_thr - m

    if clearly_sat:
        p, q, xp, xq = viol
        # double-check the single witness point through a FRESH session
        sess2 = ort.InferenceSession(model_bytes, sess_options=so,
                                     providers=["CPUExecutionProvider"])
        y2 = eval_point(sess2, xp, xq)
        recheck = (y2 <= out_thr - m) if out_op == "<=" else (y2 >= out_thr + m)
        if not recheck:
            unknown(f"witness re-evaluation disagreed (y={y2:.17g} at ({p},{q}))")
        x = base.copy()
        x[i0], x[i1] = xp, xq
        with open(args.out_witness_csv, "w") as f:
            for v in x:
                f.write(f"{v:.17g}\n")
        print(f"SAT {p} {q} {y2:.17g}")
        sys.exit(EXIT_SAT)
    elif clearly_unsat:
        print(f"UNSAT {ymin:.17g} {ymax:.17g}")
        sys.exit(EXIT_UNSAT)
    else:
        unknown(f"verdict within margin {m:g} of threshold {out_thr:g} "
                f"(Y_0 range [{ymin:.17g}, {ymax:.17g}])")


if __name__ == "__main__":
    try:
        main()
    except SystemExit:
        raise
    except BaseException as e:  # any unexpected failure must stay sound
        unknown(f"{type(e).__name__}: {e}")
