#!/usr/bin/env python3
"""collins_falsify.py -- sound falsification-only path for collins_aerospace_benchmark.

Usage:
    python collins_falsify.py <onnx> <vnnlib> <out_witness_csv> <budget_seconds> [--check-center-only]

Exit codes:
    10  SAT     -- a numerically re-validated counterexample was written to <out_witness_csv>
     2  UNKNOWN -- no (validated) counterexample found within budget, or any inconsistency
                   (NEVER emits unsat; set-based reach on this 640x640x3 YOLOv5 is out of scope)

Benchmark facts this script relies on (verified 2026-06 against the official 1.0 specs):
  * Model: yolov5nano_LRelu_640.onnx, input [1,3,640,640] (CHW), output [1,25200,11].
  * The vnnlib X_i order is **HWC flat** (NOT CHW):
        x_chw[c, h, w] = flat[h*640*3 + w*3 + c]
    Reading the spec as CHW makes every instance look trivially 'violated' at the box
    center, which the official checker (which replays HWC) rejects -> guaranteed -150.
    So before doing ANYTHING we validate the mapping with a center-consistency check:
    at the box center (the nominal image) the output property must NOT hold (the
    objectness sits strictly inside its band). If the center looks violated, the
    mapping is wrong -> print UNKNOWN and never proceed.
  * Spec shape: per-pixel bounds for all 1,228,800 inputs (most fixed lb==ub,
    216-3040 free), and ONE output assert of the form
        (assert (or (and atom) ... ))
    over ~7 consecutive Y indices of a single anchor row r: an objectness band
    (Y[r,4] outside [lower,upper]), per-class [0,1] escapes, and class-argmax swaps.
    SAT semantics (VNN-LIB): the property holds (counterexample exists) iff some
    disjunct's atoms ALL hold.

Witness CSV format (consumed by run_vnncomp_instance.m's collins branch):
    line 1: nx        (number of inputs, flat vnnlib order = HWC)
    line 2: ny        (number of outputs, flat row-major = vnnlib Y order)
    then nx x-values, then ny y-values, one per line, full round-trip precision.
"""

import gzip
import os
import re
import sys
import tempfile
import time

import numpy as np
import onnxruntime as ort

H = W = 640
C = 3
MARGIN = 1e-6          # required violation margin for accepting a counterexample
FD_STEP = 1e-3         # finite-difference step (pixels live in [0,1])


def log(msg):
    print(msg, flush=True)


def open_maybe_gz(path, mode="rt"):
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, mode)
    return open(path, mode)


def materialize_onnx(path):
    """onnxruntime needs a real .onnx file; transparently gunzip if needed."""
    with open(path, "rb") as f:
        magic = f.read(2)
    if magic != b"\x1f\x8b":
        return path, None
    fd, tmp = tempfile.mkstemp(suffix=".onnx")
    with os.fdopen(fd, "wb") as out, gzip.open(path, "rb") as g:
        out.write(g.read())
    return tmp, tmp


# ---------------------------------------------------------------------------
# vnnlib parsing (the standard 1.0 grammar subset these 6 specs actually use)
# ---------------------------------------------------------------------------

_X_BOUND = re.compile(r"^\(assert \((<=|>=) X_(\d+) ([-+0-9.eE]+)\)\)\s*$")


def _strides(shape):
    """row-major (C) strides for a shape: [1,3,640,640] -> [1228800,409600,640,1]."""
    st = [1] * len(shape)
    for i in range(len(shape) - 2, -1, -1):
        st[i] = st[i + 1] * shape[i + 1]
    return st


def _strip_top_sexpr(text, head):
    """Remove every balanced `(head ...)` s-expression from text (head e.g. 'declare-network')."""
    out, i, pat = [], 0, "(" + head
    while i < len(text):
        j = text.find(pat, i)
        if j == -1:
            out.append(text[i:]); break
        out.append(text[i:j])
        depth, k = 0, j
        while k < len(text):
            if text[k] == "(":
                depth += 1
            elif text[k] == ")":
                depth -= 1
                if depth == 0:
                    k += 1; break
            k += 1
        i = k
    return "".join(out)


def _normalize_vnnlib_2x(text):
    """Rewrite VNN-LIB 2.0 tensor notation into the flat 1.0 form parse_vnnlib already handles, WITHOUT
    touching any proven downstream logic. 2.0 declares e.g.
        (declare-input  X float32 [1,3,640,640])   (declare-output Y float32 [1,25200,11])
    and indexes tensors X[0,c,h,w] / Y[0,r,k]; map each to the SAME row-major flat index the 1.0 path used
    (X[0,c,h,w] -> X_(c*409600+h*640+w); Y[0,r,k] -> Y_(r*11+k)) and return n_x from the input shape.
    Returns (rewritten_text, n_x). SOUNDNESS: this is a falsify-only path -- any mis-rewrite yields a witness
    that the INDEPENDENT authoritative_witness_gate (vnnlib-pypi parse + onnxruntime replay, a different
    parser) rejects -> 'unknown', never an unsound UNSAT/false-robust.
    """
    def shape_of(kind):
        m = re.search(r"\(declare-" + kind + r"\s+\w+\s+\w+\s*\[\s*([0-9,\s]+?)\s*\]\s*\)", text)
        return [int(x) for x in m.group(1).split(",")] if m else None
    in_shape = shape_of("input")
    if in_shape is None:
        raise ValueError("vnnlib-2.0: no (declare-input ...) found")
    out_shape = shape_of("output")
    n_x = 1
    for d in in_shape:
        n_x *= d
    st_in = _strides(in_shape)
    st_out = _strides(out_shape) if out_shape else None

    def _sub(var, strides):
        def f(m):
            idx = [int(x) for x in m.group(1).split(",")]
            if len(idx) != len(strides):
                raise ValueError("vnnlib-2.0: index rank mismatch vs declared shape")
            return f"{var}_{sum(i * s for i, s in zip(idx, strides))}"
        return f
    t = re.sub(r"X\[\s*([0-9,\s]+?)\s*\]", _sub("X", st_in), text)
    if st_out is not None:
        t = re.sub(r"Y\[\s*([0-9,\s]+?)\s*\]", _sub("Y", st_out), t)
    # strip the 2.0 preamble so only (assert ...) remain (declare-network wraps declare-input/output when nested)
    for head in ("vnnlib-version", "declare-network", "declare-input", "declare-output", "declare-hidden"):
        t = _strip_top_sexpr(t, head)
    return t, n_x


def parse_vnnlib(path):
    """Returns (lb, ub, groups).

    lb/ub: float64 arrays over X (flat vnnlib order).
    groups: list of OR-groups; each group is a list of disjuncts; each disjunct is a
    list of atoms (op, lhs, rhs) with op in {'<=','>='} and lhs/rhs one of
    ('Y', idx) or ('const', val). The property (= counterexample condition) holds
    iff EVERY group has at least one disjunct whose atoms ALL hold.
    Raises ValueError on anything outside the supported subset (caller -> UNKNOWN).
    """
    le_idx, le_val, ge_idx, ge_val = [], [], [], []
    n_x = 0
    rest_lines = []
    with open_maybe_gz(path) as f:
        raw = f.read()
    n_x_override = None
    if ("(declare-input" in raw) or (re.search(r"\b[XY]\[", raw) is not None):
        raw, n_x_override = _normalize_vnnlib_2x(raw)   # VNN-LIB 2.0 tensor notation -> flat 1.0 form
    for line in raw.splitlines(keepends=True):
            if line.startswith("(assert ("):
                m = _X_BOUND.match(line)
                if m is not None:
                    op, idx, val = m.group(1), int(m.group(2)), float(m.group(3))
                    if op == "<=":
                        le_idx.append(idx)
                        le_val.append(val)
                    else:
                        ge_idx.append(idx)
                        ge_val.append(val)
                    continue
                rest_lines.append(line)
            elif line.startswith("(declare-const X_"):
                n_x += 1
            elif line.startswith("(declare-const"):
                continue
            else:
                s = line.strip()
                if not s or s.startswith(";"):
                    continue
                rest_lines.append(line)

    if n_x_override is not None:
        n_x = n_x_override
    if n_x == 0:
        raise ValueError("no X declarations found")
    lb = np.full(n_x, np.nan)
    ub = np.full(n_x, np.nan)
    ub[np.asarray(le_idx, dtype=np.int64)] = np.asarray(le_val)
    lb[np.asarray(ge_idx, dtype=np.int64)] = np.asarray(ge_val)
    if np.isnan(lb).any() or np.isnan(ub).any():
        raise ValueError("input box is not fully bounded (missing lb/ub for some X_i)")
    if (lb > ub).any():
        raise ValueError("inconsistent input box (lb > ub)")

    groups = parse_output_asserts("".join(rest_lines))
    if not groups:
        raise ValueError("no output constraints found")
    return lb, ub, groups


def tokenize(text):
    return re.findall(r"\(|\)|[^\s()]+", text)


def parse_sexprs(tokens):
    """Token list -> list of nested-list s-expressions."""
    out, stack = [], []
    for tok in tokens:
        if tok == "(":
            stack.append([])
        elif tok == ")":
            if not stack:
                raise ValueError("unbalanced ')'")
            done = stack.pop()
            if stack:
                stack[-1].append(done)
            else:
                out.append(done)
        else:
            if stack:
                stack[-1].append(tok)
            else:
                out.append(tok)
    if stack:
        raise ValueError("unbalanced '('")
    return out


def parse_term(tok):
    if isinstance(tok, str):
        if tok.startswith("Y_"):
            return ("Y", int(tok[2:]))
        if tok.startswith("X_"):
            raise ValueError("X variable inside output property -- unsupported (mixed in/out)")
        return ("const", float(tok))
    raise ValueError("nested term in atom -- unsupported (arithmetic)")


def parse_atom(sx):
    if not (isinstance(sx, list) and len(sx) == 3 and sx[0] in ("<=", ">=")):
        raise ValueError(f"unsupported atom: {sx!r}")
    return (sx[0], parse_term(sx[1]), parse_term(sx[2]))


def parse_output_asserts(text):
    groups = []
    for sx in parse_sexprs(tokenize(text)):
        if not (isinstance(sx, list) and sx and sx[0] == "assert" and len(sx) == 2):
            raise ValueError(f"unsupported top-level form: {sx!r}")
        body = sx[1]
        if isinstance(body, list) and body and body[0] == "or":
            disjuncts = []
            for cl in body[1:]:
                if isinstance(cl, list) and cl and cl[0] == "and":
                    disjuncts.append([parse_atom(a) for a in cl[1:]])
                else:
                    disjuncts.append([parse_atom(cl)])
            groups.append(disjuncts)
        elif isinstance(body, list) and body and body[0] == "and":
            for a in body[1:]:
                groups.append([[parse_atom(a)]])
        else:
            groups.append([[parse_atom(body)]])
    return groups


# ---------------------------------------------------------------------------
# property evaluation
# ---------------------------------------------------------------------------

def term_val(term, y):
    kind, v = term
    return y[v] if kind == "Y" else v


def atom_margin(atom, y):
    """<= 0 means the atom HOLDS; value = how far from holding."""
    op, a, b = atom
    va, vb = term_val(a, y), term_val(b, y)
    return (va - vb) if op == "<=" else (vb - va)


def clause_margin(atoms, y):
    """<= 0 means ALL atoms of the (and ...) hold."""
    return max(atom_margin(a, y) for a in atoms)


def property_margin(groups, y):
    """<= 0 means the WHOLE output property holds (y is in the counterexample region).
    sum over groups of min over disjuncts (each group must have a satisfied disjunct)."""
    total = 0.0
    for g in groups:
        m = min(clause_margin(cl, y) for cl in g)
        total += max(m, 0.0) if len(groups) > 1 else m
    return total


def property_holds_with_margin(groups, y, margin):
    """True iff every group has a disjunct with ALL atoms satisfied by > margin."""
    return all(any(clause_margin(cl, y) < -margin for cl in g) for g in groups)


# ---------------------------------------------------------------------------
# model evaluation (HWC flat -> CHW tensor)
# ---------------------------------------------------------------------------

class Model:
    def __init__(self, onnx_path):
        # ORT's default thread heuristic beat an explicit ncpu-1 override on this
        # nano model in measurement; keep defaults.
        self.sess = ort.InferenceSession(onnx_path, providers=["CPUExecutionProvider"])
        self.in_name = self.sess.get_inputs()[0].name
        self.n_fwd = 0

    def forward_flat_hwc(self, x_flat):
        """x_flat: float64/32 (1228800,) in HWC order -> y_flat (277200,) row-major."""
        x = np.asarray(x_flat, dtype=np.float32).reshape(H, W, C)
        x_chw = np.ascontiguousarray(np.transpose(x, (2, 0, 1)))[None]  # [1,3,640,640]
        y = self.sess.run(None, {self.in_name: x_chw})[0]
        self.n_fwd += 1
        return y.reshape(-1).astype(np.float64)


# ---------------------------------------------------------------------------
# falsification
# ---------------------------------------------------------------------------

def fd_gradient(model, x, free, loss_fn, f0, deadline):
    """One-sided finite-difference gradient of loss_fn(forward(x)) over indices `free`.
    Steps stay inside the box implicitly: caller clips after the update anyway.
    Returns (grad over free, completed_flag)."""
    g = np.zeros(free.size)
    for k, i in enumerate(free):
        if time.time() > deadline:
            return g, False
        xp = x.copy()
        xp[i] += FD_STEP
        g[k] = (loss_fn(model.forward_flat_hwc(xp)) - f0) / FD_STEP
    return g, True


def try_candidate(model, x, lb, ub, groups, best):
    """Evaluate x (clipped to the box) and track the best property margin.
    Returns (sat, x, y): sat means the property HOLDS with > MARGIN -- the explicit
    check, not a margin-arithmetic proxy, so multi-group (conjunctive) specs work."""
    x = np.clip(x, lb, ub)
    y = model.forward_flat_hwc(x)
    m = property_margin(groups, y)
    sat = property_holds_with_margin(groups, y, MARGIN)
    if sat or m < best[0]:
        best[0], best[1], best[2] = (-np.inf if sat else m), x, y
    return sat, x, y


def falsify(model, lb, ub, groups, deadline):
    """Search the box for a property-satisfying point. Returns (x, y) or None."""
    free = np.flatnonzero(ub - lb > 0)
    log(f"free pixels: {free.size} / {lb.size}")
    center = 0.5 * (lb + ub)
    best = [np.inf, None, None]  # [margin, x, y]

    # 0) center (already known consistent by the caller's gate) + box corners
    for probe in (center, lb.copy(), ub.copy()):
        sat, _, _ = try_candidate(model, probe, lb, ub, groups, best)
        if sat:
            return best[1], best[2]
        if time.time() > deadline:
            return None

    if free.size == 0:
        return None

    # 1) rank single-clause targets by achievability at the center, attack each with
    #    FD-gradient FGSM steps (one gradient + a line search along the box corner),
    #    then a few PGD refinement rounds while budget remains.
    y_center = model.forward_flat_hwc(center)
    order = []
    for gi, g in enumerate(groups):
        for ci, cl in enumerate(g):
            cm = clause_margin(cl, y_center)
            if cm >= 0.5:
                continue  # far targets are a last resort; random corners cover them
            # Search-order heuristic ONLY (never affects verdicts): the YOLO head is
            # sigmoid-activated, so the (<= Y 0.0)/(>= Y 1.0) escapes can show a
            # near-zero center margin (a class score ~1e-4) while being
            # asymptotically unreachable. Attack band / argmax-swap clauses first,
            # [0,1]-boundary escapes last.
            boundary = all(
                atom[1][0] == "Y" and atom[2][0] == "const" and atom[2][1] in (0.0, 1.0)
                for atom in cl
            )
            order.append((1 if boundary else 0, cm, gi, ci))
    order.sort()
    if len(groups) > 1:
        # multi-group conjunction: optimize the global property margin directly
        order = [(0, property_margin(groups, y_center), -1, -1)]

    # NB: iterate ALL plausible targets (deadline-guarded), not a fixed top-k; the
    # per-target stall detection below abandons dead ends and moves on.
    for _, cm0, gi, ci in order:
        if time.time() > deadline:
            break
        if gi < 0:
            loss_fn = lambda y: property_margin(groups, y)
        else:
            cl = groups[gi][ci]
            loss_fn = lambda y, cl=cl: clause_margin(cl, y)
        log(f"target clause group={gi} disjunct={ci} center-margin={cm0:.6g}")

        x_cur = center.copy()
        f_cur = loss_fn(model.forward_flat_hwc(x_cur))
        for it in range(8):  # PGD-style: FD gradient -> corner line search -> repeat
            if time.time() > deadline:
                break
            grad, complete = fd_gradient(model, x_cur, free, loss_fn, f_cur, deadline)
            if not complete:
                log("budget exhausted during gradient")
                break
            corner = x_cur.copy()
            dn = grad < 0
            corner[free[dn]] = ub[free[dn]]   # loss decreases as x increases -> go up
            corner[free[~dn]] = lb[free[~dn]]
            improved = False
            for alpha in (1.0, 0.5, 0.25, 0.125, 0.0625):
                x_try = x_cur + alpha * (corner - x_cur)
                sat, x_try, y_try = try_candidate(model, x_try, lb, ub, groups, best)
                if sat:
                    return best[1], best[2]
                f_try = loss_fn(y_try)
                if f_try < f_cur - 1e-9:
                    x_cur, f_cur = x_try, f_try
                    improved = True
                    break
                if time.time() > deadline:
                    break
            log(f"  iter {it}: clause-loss {f_cur:.6g} (property margin best {best[0]:.6g})")
            if not improved:
                break

    # 2) leftover budget: random corner sampling on the free pixels
    rng = np.random.default_rng(0)
    n_rand = 0
    while time.time() < deadline:
        x_try = center.copy()
        pick = rng.random(free.size) < 0.5
        x_try[free[pick]] = ub[free[pick]]
        x_try[free[~pick]] = lb[free[~pick]]
        sat, _, _ = try_candidate(model, x_try, lb, ub, groups, best)
        n_rand += 1
        if sat:
            return best[1], best[2]
    log(f"random corner sampling: {n_rand} samples, best property margin {best[0]:.6g}")
    return None


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def describe_band(groups, y):
    """Diagnostic: find the objectness-band disjuncts (a Y index with both a >=c and
    a <=c single-atom disjunct, 0<c<1) and report nominal vs band."""
    if len(groups) != 1:
        return
    lo = {}
    hi = {}
    for cl in groups[0]:
        if len(cl) != 1:
            continue
        op, a, b = cl[0]
        if a[0] == "Y" and b[0] == "const" and 0.0 < b[1] < 1.0:
            (hi if op == ">=" else lo)[a[1]] = b[1]
    for idx in sorted(set(lo) & set(hi)):
        log(f"band check Y_{idx}: lower={lo[idx]!r} <= value={float(y[idx])!r} <= upper={hi[idx]!r} "
            f"-> {'INSIDE (no violation here)' if lo[idx] < y[idx] < hi[idx] else 'OUTSIDE'}")


def main(argv):
    if len(argv) < 5:
        log(__doc__)
        return 2
    onnx_path, vnnlib_path, out_csv = argv[1], argv[2], argv[3]
    budget = float(argv[4])
    center_only = "--check-center-only" in argv[5:]
    t0 = time.time()
    deadline = t0 + max(10.0, budget)

    try:
        lb, ub, groups = parse_vnnlib(vnnlib_path)
    except Exception as e:  # malformed/unsupported spec -> never guess
        log(f"vnnlib parse failed/unsupported: {e}")
        log("UNKNOWN")
        return 2
    log(f"parsed vnnlib in {time.time()-t0:.1f}s: {lb.size} inputs, "
        f"{sum(len(g) for g in groups)} output disjunct(s) in {len(groups)} group(s)")
    if lb.size != H * W * C:
        log(f"unexpected input count {lb.size} (expected {H*W*C})")
        log("UNKNOWN")
        return 2

    tmp = None
    try:
        model_path, tmp = materialize_onnx(onnx_path)
        model = Model(model_path)

        # --- center-consistency gate: validates the HWC->CHW mapping --------------
        center = 0.5 * (lb + ub)
        y_center = model.forward_flat_hwc(center)
        if y_center.size != 25200 * 11:
            log(f"unexpected output count {y_center.size}")
            log("UNKNOWN")
            return 2
        describe_band(groups, y_center)
        if property_holds_with_margin(groups, y_center, 0.0):  # holds even at zero margin
            # Under the verified HWC mapping the nominal center NEVER satisfies the
            # property (bands are built around the nominal output). A 'violated'
            # center means the input mapping is WRONG -> emitting it would be the
            # exact -150 failure mode. Refuse.
            log("center-consistency FAILED: property already 'holds' at the box center "
                "-> input-order mapping suspect; refusing to falsify")
            log("UNKNOWN")
            return 2
        log("center-consistency OK: property does not hold at the nominal center")
        if center_only:
            log("CENTER_CONSISTENT")
            log("UNKNOWN")  # mode is diagnostic only; never claims sat
            return 2

        # --- falsify ---------------------------------------------------------------
        res = falsify(model, lb, ub, groups, deadline - 5.0)
        if res is None:
            log(f"no counterexample found ({model.n_fwd} forward passes, "
                f"{time.time()-t0:.1f}s)")
            log("UNKNOWN")
            return 2

        x_cand = np.clip(res[0], lb, ub)  # exact final witness, float64, in-box

        # --- independent re-validation of the EXACT witness ------------------------
        y_cand = model.forward_flat_hwc(x_cand)  # fresh forward on the exact values
        if not ((x_cand >= lb).all() and (x_cand <= ub).all()):
            log("witness left the input box after clipping?! refusing")
            log("UNKNOWN")
            return 2
        if not property_holds_with_margin(groups, y_cand, MARGIN):
            log("re-validation failed: candidate does not violate by > "
                f"{MARGIN:g} on a fresh forward pass; refusing")
            log("UNKNOWN")
            return 2
        for gi, g in enumerate(groups):
            bi = int(np.argmin([clause_margin(cl, y_cand) for cl in g]))
            log(f"violated disjunct group={gi} idx={bi} margin={-clause_margin(g[bi], y_cand):.6g}: "
                f"{g[bi]}")
        describe_band(groups, y_cand)

        with open(out_csv, "w") as f:
            f.write(f"{x_cand.size}\n{y_cand.size}\n")
            # NB: float(v) BEFORE repr -- numpy 2.x repr of a np.float64 scalar is
            # 'np.float64(0.447...)', which MATLAB's readmatrix cannot parse.
            # repr(float) is shortest-round-trip, so the witness replays exactly.
            for v in x_cand:
                f.write(f"{float(v)!r}\n")
            for v in y_cand:
                f.write(f"{float(v)!r}\n")
        log(f"witness written: {out_csv} ({model.n_fwd} forward passes, "
            f"{time.time()-t0:.1f}s)")
        log("SAT")
        return 10
    except Exception as e:
        log(f"error: {type(e).__name__}: {e}")
        log("UNKNOWN")
        return 2
    finally:
        if tmp is not None:
            try:
                os.remove(tmp)
            except OSError:
                pass


if __name__ == "__main__":
    sys.exit(main(sys.argv))
