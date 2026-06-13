"""manifest_sim.py — a pure-numpy interpreter for the .nnv.mat manifests
produced by onnx2nnv.py, replicating the MATLAB-side semantics of
engine/utils/load_nnv_from_mat.m + each NNV layer's evaluate() EXACTLY
(column-major reshapes, trailing-singleton trimming, MATLAB implicit
expansion, the FullyConnectedLayer token fast path, ...).

Purpose: iterate on the importer's tensor-LAYOUT model without a MATLAB
session. If this simulator reproduces onnxruntime on manifests whose MATLAB
evaluate is KNOWN-GOOD (soundnessbench, cgan), it faithfully models the
MATLAB semantics; the same simulator then validates layout fixes (vit_2023)
ORIENTATION-STRICTLY before the human re-checks in MATLAB.

The companion matrix runner lives in check_manifest_matrix.py.

MATLAB-array model
------------------
Values are numpy float64 arrays kept in "MATLAB normal form": ndim >= 2 and
no trailing singleton dimensions beyond ndim 2 (MATLAB's size() semantics).
All reshapes are column-major (order='F'); broadcasting pads TRAILING
singleton dims (MATLAB implicit expansion), unlike numpy's leading-dim rule.
"""

import numpy as np
from scipy.io import loadmat


class SimError(Exception):
    """A MATLAB-side evaluate() would have errored the same way."""


class OracleNeeded(Exception):
    """Layer wraps an op MATLAB refuses to evaluate (UnsupportedOp:* /
    Resize); the harness substitutes onnxruntime ground truth."""


# ---------------------------------------------------------------------------
# MATLAB array helpers
# ---------------------------------------------------------------------------

def mn(a):
    """MATLAB normal form: float64, ndim >= 2, trailing singletons trimmed."""
    a = np.asarray(a, dtype=np.float64)
    if a.ndim == 0:
        return a.reshape(1, 1)
    if a.ndim == 1:
        return a.reshape(1, -1)   # scipy/savemat convention for 1-D: row
    shp = list(a.shape)
    while len(shp) > 2 and shp[-1] == 1:
        shp.pop()
    return a.reshape(shp)


def m_size(a):
    return list(a.shape)


def m_ndims(a):
    return a.ndim


def m_numel(a):
    return a.size


def m_isvector(a):
    return a.ndim == 2 and (a.shape[0] == 1 or a.shape[1] == 1)


def m_isscalar(a):
    return a.size == 1


def m_reshape(a, tgt):
    """MATLAB reshape: column-major; numel must match; allows trailing 1s."""
    tgt = [int(t) for t in tgt]
    if -1 in tgt:
        known = int(np.prod([t for t in tgt if t != -1]))
        tgt[tgt.index(-1)] = a.size // max(known, 1)
    if int(np.prod(tgt)) != a.size:
        raise SimError(f"reshape: numel mismatch {a.size} -> {tgt}")
    return mn(np.reshape(a, tgt, order='F'))


def m_permute(a, perm1):
    """MATLAB permute with a 1-indexed perm; perm may exceed ndims (pads the
    array with trailing singletons) and is padded with identity if shorter."""
    perm1 = [int(p) for p in perm1]
    k = max(len(perm1), a.ndim)
    perm1 = perm1 + list(range(len(perm1) + 1, k + 1))
    a2 = a.reshape(a.shape + (1,) * (k - a.ndim))
    return mn(np.transpose(a2, [p - 1 for p in perm1]))


def m_broadcast(a, b, f):
    """MATLAB implicit expansion: align LEADING dims (pad TRAILING 1s)."""
    a = mn(a); b = mn(b)
    k = max(a.ndim, b.ndim)
    a2 = a.reshape(a.shape + (1,) * (k - a.ndim))
    b2 = b.reshape(b.shape + (1,) * (k - b.ndim))
    for da, db in zip(a2.shape, b2.shape):
        if da != db and da != 1 and db != 1:
            raise SimError(f"implicit expansion size mismatch {a.shape} vs {b.shape}")
    return mn(f(a2, b2))


def m_col(a):
    """MATLAB a(:) — column-major flatten to a column vector."""
    return mn(np.reshape(a, (-1, 1), order='F'))


# ---------------------------------------------------------------------------
# .mat decoding helpers (scipy loadmat, struct_as_record=False)
# ---------------------------------------------------------------------------

def _to_str(x):
    a = np.asarray(x)
    if a.size == 0:
        return ''
    return str(a.flatten()[0])


def _to_str_list(x):
    if x is None:
        return []
    a = np.asarray(x)
    if a.size == 0:
        return []
    if a.dtype == object:
        return [_to_str(e) for e in a.flatten()]
    if a.dtype.kind in 'US':
        return [str(e) for e in a.flatten()]
    return [_to_str(a)]


def _attr(attrs, name, default=None):
    if attrs is None or not hasattr(attrs, name):
        return default
    return getattr(attrs, name)


def _num(x, default=None):
    if x is None:
        return default
    a = np.asarray(x)
    if a.size == 0:
        return default
    return float(np.asarray(a).flatten()[0])


def _ivec(x):
    return [int(v) for v in np.asarray(x).flatten()]


# ---------------------------------------------------------------------------
# ElementwiseAffineLayer helpers (mirror ElementwiseAffineLayer.m)
# ---------------------------------------------------------------------------

def _ea_broadcastable(sz_p, sz_x):
    if len(sz_p) != len(sz_x):
        return False
    return all(p == x or p == 1 for p, x in zip(sz_p, sz_x))


def _ea_trailing_align(p, x):
    # onnx_trailing_align: pad LEADING singletons up to ndims(x)
    if p.ndim < x.ndim:
        p = mn(p.reshape((1,) * (x.ndim - p.ndim) + p.shape))
    return p


def _ea_padshape_left(arr, target_ndims):
    sz = m_size(arr)
    while len(sz) > target_ndims and sz[0] == 1:
        sz = sz[1:]
    need = target_ndims - len(sz)
    if need > 0:
        sz = [1] * need + sz
    if len(sz) < 2:
        sz = sz + [1] * (2 - len(sz))
    return m_reshape(arr, sz)


def _ea_align_to_input(arr, x_shape):
    sz_a = m_size(arr)
    ns_a = [d for d in sz_a if d > 1]
    if not ns_a:
        return arr
    if len(ns_a) == 1:
        target = [1] * max(2, len(x_shape))
        placed = False
        val = ns_a[0]
        for k, d in enumerate(x_shape):
            if d == val and not placed:
                target[k] = val
                placed = True
        if placed:
            return m_reshape(arr, target)
    target = [1] * max(2, len(x_shape))
    remaining = list(ns_a)
    for k, d in enumerate(x_shape):
        if remaining and d == remaining[0]:
            target[k] = remaining[0]
            remaining.pop(0)
    if not remaining:
        return m_reshape(arr, target)
    return _ea_padshape_left(arr, len(x_shape))


# ---------------------------------------------------------------------------
# Layer evaluates
# ---------------------------------------------------------------------------

def _fc_eval(W, b_col, x):
    """FullyConnectedLayer.evaluate — including the transformer batched-FC
    fast path (rank>=3 input whose LAST dim equals in_features)."""
    n = m_size(x)
    in_feats = W.shape[1]
    if len(n) >= 3 and n[-1] == in_feats and int(np.prod(n[:-1])) > 1:
        lead = n[:-1]
        x_flat = m_reshape(x, [int(np.prod(lead)), in_feats])
        y_flat = (W @ x_flat.T).T + b_col.T   # [N, out] + [1, out]
        return m_reshape(y_flat, list(lead) + [W.shape[0]])
    if len(n) > 4:
        raise SimError('FC: invalid input rank')
    x_col = m_col(x)
    if x_col.shape[0] != in_feats:
        raise SimError(f'FC: inconsistent input vector ({x_col.shape[0]} vs {in_feats})')
    return mn(W @ x_col + b_col)


def _conv2d_eval(x, Wm, b_col, pads_tblr, strides, dils, groups=1):
    """Conv2DLayer.evaluate (dlconv 'SSCB' on HWC; correlation, no flip)."""
    if groups != 1:
        raise SimError('Conv2D: grouped conv not supported by NNV manifest path')
    x = mn(x)
    if x.ndim == 2:
        x = x.reshape(x.shape + (1,))
    t, bo, l, r = [int(v) for v in pads_tblr]
    xp = np.pad(x, ((t, bo), (l, r), (0, 0)))
    kH, kW, C, F = Wm.shape
    if xp.shape[2] != C:
        raise SimError(f'Conv2D: channel mismatch {xp.shape[2]} vs {C}')
    sh, sw = [int(v) for v in strides]
    dh, dw = [int(v) for v in dils]
    eKH = (kH - 1) * dh + 1
    eKW = (kW - 1) * dw + 1
    Ho = (xp.shape[0] - eKH) // sh + 1
    Wo = (xp.shape[1] - eKW) // sw + 1
    if Ho <= 0 or Wo <= 0:
        raise SimError('Conv2D: empty output')
    y = np.tile(b_col.reshape(1, 1, F), (Ho, Wo, 1)).astype(np.float64)
    for i in range(kH):
        for j in range(kW):
            xs = xp[i * dh: i * dh + (Ho - 1) * sh + 1: sh,
                    j * dw: j * dw + (Wo - 1) * sw + 1: sw, :]
            y += xs @ Wm[i, j]          # [Ho,Wo,C] @ [C,F]
    return mn(y)


def _transp_conv2d_eval(x, Wm, b_col, crops_tblr, strides):
    """TransposedConv2DLayer.evaluate (dltranspconv 'SSCU')."""
    x = mn(x)
    if x.ndim == 2:
        x = x.reshape(x.shape + (1,))
    kH, kW, F, C = Wm.shape
    H, W_, Cin = x.shape
    if Cin != C:
        raise SimError(f'TranspConv2D: channel mismatch {Cin} vs {C}')
    sh, sw = [int(v) for v in strides]
    t, bo, l, r = [int(v) for v in crops_tblr]
    Hf = (H - 1) * sh + kH
    Wf = (W_ - 1) * sw + kW
    y = np.zeros((Hf, Wf, F), dtype=np.float64)
    for i in range(kH):
        for j in range(kW):
            # scatter: y[i + h*sh, j + w*sw, f] += sum_c x[h,w,c] W[i,j,f,c]
            y[i: i + (H - 1) * sh + 1: sh,
              j: j + (W_ - 1) * sw + 1: sw, :] += x @ Wm[i, j].T
    y = y[t: Hf - bo, l: Wf - r, :]
    return mn(y + b_col.reshape(1, 1, F))


def _maxpool_eval(x, pool, stride, pads_tblr):
    x = mn(x)
    if x.ndim == 2:
        x = x.reshape(x.shape + (1,))
    t, bo, l, r = [int(v) for v in pads_tblr]
    xp = np.pad(x, ((t, bo), (l, r), (0, 0)), constant_values=-np.inf)
    kH, kW = [int(v) for v in pool]
    sh, sw = [int(v) for v in stride]
    Ho = (xp.shape[0] - kH) // sh + 1
    Wo = (xp.shape[1] - kW) // sw + 1
    y = np.full((Ho, Wo, x.shape[2]), -np.inf)
    for i in range(kH):
        for j in range(kW):
            y = np.maximum(y, xp[i: i + (Ho - 1) * sh + 1: sh,
                                 j: j + (Wo - 1) * sw + 1: sw, :])
    return mn(y)


def _avgpool_eval(x, pool, stride, pads_tblr):
    x = mn(x)
    if x.ndim == 2:
        x = x.reshape(x.shape + (1,))
    t, bo, l, r = [int(v) for v in pads_tblr]
    xp = np.pad(x, ((t, bo), (l, r), (0, 0)))
    kH, kW = [int(v) for v in pool]
    sh, sw = [int(v) for v in stride]
    Ho = (xp.shape[0] - kH) // sh + 1
    Wo = (xp.shape[1] - kW) // sw + 1
    y = np.zeros((Ho, Wo, x.shape[2]))
    for i in range(kH):
        for j in range(kW):
            y += xp[i: i + (Ho - 1) * sh + 1: sh,
                    j: j + (Wo - 1) * sw + 1: sw, :]
    return mn(y / (kH * kW))


def _softmax_layer_eval(x):
    """SoftmaxLayer.evaluate: isvector -> column-wise softmax (NN toolbox
    softmax()); rank-3 -> 'SSC' channel softmax over MATLAB dim 3; else error."""
    x = mn(x)
    if m_isvector(x):
        e = np.exp(x - np.max(x, axis=0, keepdims=True))
        return mn(e / np.sum(e, axis=0, keepdims=True))
    if x.ndim == 3:
        e = np.exp(x - np.max(x, axis=2, keepdims=True))
        return mn(e / np.sum(e, axis=2, keepdims=True))
    raise SimError('SoftmaxLayer: input is not a vector or multi-channel image')


def _flatten_into2d_eval(x):
    """FlattenLayer.evaluate, Type='nnet.onnx.layer.FlattenInto2dLayer'."""
    x = mn(x)
    n = m_size(x)
    if len(n) == 2:
        xi = m_permute(x, [2, 1])
        return m_reshape(xi, [1, 1, n[0] * n[1]])
    if len(n) == 3:
        xi = m_permute(x, [2, 1, 3])
        return m_reshape(xi, [1, 1, n[0] * n[1] * n[2]])
    if len(n) == 4:
        xi = m_permute(x, [3, 2, 1, 4])
        return m_reshape(xi, [1, 1, int(np.prod(n))])
    if len(n) == 5:
        xi = m_permute(x, [1, 5, 4, 3, 2])
        return m_reshape(xi, [1, 1, int(np.prod(n))])
    raise SimError('FlattenLayer: invalid input rank')


def _reshape_layer_eval(x, tgt, onnx_bchw):
    x = mn(x)
    tgt = [int(t) for t in tgt]
    neg = [k for k, t in enumerate(tgt) if t < 0]
    if len(neg) > 1:
        raise SimError('ReshapeLayer: more than one -1 sentinel')
    if len(neg) == 1:
        known = list(tgt)
        known[neg[0]] = 1
        rest = int(np.prod(known))
        tgt[neg[0]] = x.size // rest if rest > 0 else 1
    if onnx_bchw and len(tgt) == 3:
        H, W_, C = tgt
        tmp = m_reshape(x, [W_, H, C])
        return m_permute(tmp, [2, 1, 3])
    return m_reshape(x, tgt)


def _dynmatmul_eval(A, B):
    A = mn(A); B = mn(B)
    if A.ndim <= 2 and B.ndim <= 2:
        if A.shape[1] != B.shape[0]:
            raise SimError(f'DynamicMatmul: inner dims mismatch ({A.shape[1]} vs {B.shape[0]})')
        return mn(A @ B)
    sa = m_size(A); sb = m_size(B)
    k_a, k_a_in = sa[-2], sa[-1]
    k_b, k_b_in = sb[-2], sb[-1]
    if k_a_in != k_b:
        raise SimError(f'DynamicMatmul: inner dims mismatch ({k_a_in} vs {k_b})')
    lead_a, lead_b = sa[:-2], sb[:-2]
    if lead_a != lead_b:
        raise SimError(f'DynamicMatmul: leading dims differ ({lead_a} vs {lead_b})')
    n_lead = max(int(np.prod(lead_a)), 1) if lead_a else 1
    A_flat = m_reshape(A, [n_lead, k_a, k_a_in])
    B_flat = m_reshape(B, [n_lead, k_b, k_b_in])
    # guard: m_reshape may have trimmed trailing singletons
    A_flat = np.reshape(A_flat, (n_lead, k_a, k_a_in), order='F')
    B_flat = np.reshape(B_flat, (n_lead, k_b, k_b_in), order='F')
    out_flat = np.zeros((n_lead, k_a, k_b_in))
    for i in range(n_lead):
        out_flat[i] = A_flat[i] @ B_flat[i]
    return m_reshape(out_flat, list(lead_a) + [k_a, k_b_in])


# ---------------------------------------------------------------------------
# Layer objects (mirror load_nnv_from_mat.build_layer + each evaluate)
# ---------------------------------------------------------------------------

class SimLayer:
    def __init__(self, type_str, name, attrs, wkeys, weights, inputs, outputs):
        self.type = type_str
        self.name = name
        self.attrs = attrs
        self.inputs = inputs
        self.outputs = outputs
        self.layout = _to_str(_attr(attrs, 'MlLayout', 'FLAT')) or 'FLAT'
        w = [mn(np.asarray(getattr(weights, k))) for k in wkeys] if wkeys else []

        t = type_str
        if t == 'FullyConnectedLayer':
            self.W = w[0]
            self.b = m_col(w[1])
        elif t == 'Conv2DLayer':
            # MATLAB trims trailing singleton dims of the stored [F,C,kH,kW]
            # tensor (1x1 kernels); permute() pads them back — replicate.
            W4 = np.asarray(w[0], dtype=np.float64)
            W4 = W4.reshape(W4.shape + (1,) * (4 - W4.ndim))
            self.Wm = np.transpose(W4, (2, 3, 1, 0))   # [kH,kW,C,F]
            self.b = m_col(w[1])
            self.pads = self._pads(_attr(attrs, 'Pads'))
            self.strides = _ivec(_attr(attrs, 'Strides', [1, 1]))
            self.dils = _ivec(_attr(attrs, 'Dilations', [1, 1]))
            self.groups = int(_num(_attr(attrs, 'Groups'), 1))
        elif t == 'TransposedConv2DLayer':
            W4 = np.asarray(w[0], dtype=np.float64)
            W4 = W4.reshape(W4.shape + (1,) * (4 - W4.ndim))
            self.Wm = np.transpose(W4, (2, 3, 1, 0))   # [kH,kW,F,C]
            self.b = m_col(w[1])
            self.pads = self._pads(_attr(attrs, 'Pads'))
            self.strides = _ivec(_attr(attrs, 'Strides', [1, 1]))
        elif t == 'BatchNormalizationLayer':
            nC = w[0].size
            self.nC = nC
            self.scale = np.reshape(w[0], (1, 1, nC), order='F')
            self.bias = np.reshape(w[1], (1, 1, nC), order='F')
            self.mean = np.reshape(w[2], (1, 1, nC), order='F')
            self.var = np.reshape(w[3], (1, 1, nC), order='F')
            self.eps = _num(_attr(attrs, 'Epsilon'), 1e-5)
        elif t == 'ElementwiseAffineLayer':
            s, o = w[0], w[1]
            if m_isvector(s) or m_isscalar(s):
                s = m_col(s)
            if m_isvector(o) or m_isscalar(o):
                o = m_col(o)
            self.scale, self.offset = s, o
            self.do_scale = bool(_num(_attr(attrs, 'DoScale'), 0))
            self.do_offset = bool(_num(_attr(attrs, 'DoOffset'), 0))
        elif t == 'LeakyReluLayer':
            self.gamma = _num(_attr(attrs, 'Scale'), 0.01)
        elif t == 'ReshapeLayer':
            self.tgt = _ivec(_attr(attrs, 'TargetShape'))
            self.onnx_bchw = bool(_num(_attr(attrs, 'OnnxBCHW'), 0))
        elif t == 'TransposeLayer':
            perm = _ivec(_attr(attrs, 'Perm', []))
            mperm = [p + 1 for p in perm]
            self.perm = mperm if (mperm and mperm != list(range(1, len(mperm) + 1))) else None
        elif t == 'ConcatenationLayer':
            self.dim = self._concat_dim(attrs)
        elif t in ('MaxPooling2DLayer', 'AveragePooling2DLayer'):
            self.pool = _ivec(_attr(attrs, 'KernelSize'))
            self.strides = _ivec(_attr(attrs, 'Strides', [1, 1]))
            self.pads = self._pads(_attr(attrs, 'Pads'))
        elif t == 'SignLayer':
            self.mode = _to_str(_attr(attrs, 'Mode', 'polar_zero_to_pos_one'))
        elif t == 'SoftmaxLayer':
            pass
        elif t == 'PlaceholderLayer':
            # mirror the loader's tag logic
            op = _to_str(_attr(attrs, 'OriginalOp', ''))
            if op in ('Sign', 'Abs', 'Constant', 'Floor', 'Ceil', 'Round',
                      'Sin', 'Cos', 'Tan', 'Exp', 'Log', 'Sqrt'):
                tag = op
            elif op in ('', 'Identity', 'Cast', 'Dropout', 'extra_graph_input',
                        'Transpose_2D_no_op_vector', 'Transpose_BCHW_to_BHWC',
                        'Transpose_BHWC_to_BCHW'):
                tag = 'Identity'
            else:
                tag = 'UnsupportedOp:' + op
            self.tag = tag
            self.constant = None
            if tag == 'Constant':
                if w:
                    self.constant = w[0]
                else:
                    self.tag = 'UnsupportedOp:Constant'
        self.in_size = _attr(attrs, 'InputSize')

    @staticmethod
    def _pads(p):
        pf = _ivec(p) if p is not None else []
        if len(pf) == 4:
            return [pf[0], pf[2], pf[1], pf[3]]   # ONNX [t,l,b,r] -> [t,b,l,r]
        return [0, 0, 0, 0]

    @staticmethod
    def _concat_dim(attrs):
        axis = int(_num(_attr(attrs, 'Axis'), 0))
        in_rank = int(_num(_attr(attrs, 'InRank'), 0))
        if axis < 0:
            onnx_pos = axis + in_rank if in_rank > 0 else max(0, axis + 4)
        else:
            onnx_pos = axis
        if in_rank == 2:
            d = (onnx_pos == 1) * 1 + (onnx_pos == 0) * 2
            return d if d != 0 else 1
        if in_rank == 3:
            return onnx_pos + 1
        if in_rank == 4:
            map4 = [3, 3, 1, 2]
            return map4[min(max(onnx_pos + 1, 1), 4) - 1]
        d = axis
        if d < 0:
            d = max(1, d + 4)
        return max(d, 1)

    # ---- evaluate dispatch -------------------------------------------------

    def evaluate(self, vals):
        t = self.type
        if t in ('FeatureInputLayer', 'ImageInputLayer'):
            return vals[0]
        if t == 'FullyConnectedLayer':
            return _fc_eval(self.W, self.b, vals[0])
        if t == 'Conv2DLayer':
            return _conv2d_eval(vals[0], self.Wm, self.b, self.pads,
                                self.strides, self.dils, self.groups)
        if t == 'TransposedConv2DLayer':
            return _transp_conv2d_eval(vals[0], self.Wm, self.b, self.pads,
                                       self.strides)
        if t == 'BatchNormalizationLayer':
            return self._bn_eval(vals[0])
        if t == 'ReluLayer':
            return mn(np.maximum(vals[0], 0))
        if t == 'LeakyReluLayer':
            x = mn(vals[0])
            return mn(np.where(x >= 0, x, self.gamma * x))
        if t == 'SigmoidLayer':
            return mn(1.0 / (1.0 + np.exp(-mn(vals[0]))))
        if t == 'TanhLayer':
            return mn(np.tanh(vals[0]))
        if t == 'SoftmaxLayer':
            return _softmax_layer_eval(vals[0])
        if t == 'SignLayer':
            y = np.sign(mn(vals[0]))
            if self.mode == 'polar_zero_to_pos_one':
                y = y + (y == 0)
            elif self.mode == 'nonnegative_zero_to_pos_one':
                y = y + (y == 0)
                y = y + (y == -1)
            return mn(y)
        if t == 'ElementwiseAffineLayer':
            return self._ea_eval(vals[0])
        if t == 'AdditionLayer':
            out = mn(vals[0])
            for v in vals[1:]:
                out = m_broadcast(out, v, np.add)
            return out
        if t == 'ElementwiseProductLayer':
            out = mn(vals[0])
            for v in vals[1:]:
                out = m_broadcast(out, v, np.multiply)
            return out
        if t == 'ElementwiseDivisionLayer':
            return m_broadcast(vals[0], vals[1], np.divide)
        if t == 'DynamicMatmulLayer':
            return _dynmatmul_eval(vals[0], vals[1])
        if t == 'ConcatenationLayer':
            out = mn(vals[0])
            for v in vals[1:]:
                v = mn(v)
                k = max(out.ndim, v.ndim, self.dim)
                o2 = out.reshape(out.shape + (1,) * (k - out.ndim))
                v2 = v.reshape(v.shape + (1,) * (k - v.ndim))
                out = mn(np.concatenate([o2, v2], axis=self.dim - 1))
            return out
        if t == 'FlattenLayer':
            return _flatten_into2d_eval(vals[0])
        if t == 'ReshapeLayer':
            return _reshape_layer_eval(vals[0], self.tgt, self.onnx_bchw)
        if t == 'TransposeLayer':
            # loader maps TransposeLayer -> PlaceholderLayer with Perm
            if self.perm is None:
                return mn(vals[0])
            return m_permute(vals[0], self.perm)
        if t == 'MaxPooling2DLayer':
            return _maxpool_eval(vals[0], self.pool, self.strides, self.pads)
        if t == 'AveragePooling2DLayer':
            return _avgpool_eval(vals[0], self.pool, self.strides, self.pads)
        if t == 'GlobalAveragePooling2DLayer':
            x = mn(vals[0])
            return mn(np.mean(x, axis=(0, 1), keepdims=True))
        if t == 'Resize2DLayer':
            raise OracleNeeded(self.name)
        if t == 'PlaceholderLayer':
            return self._placeholder_eval(vals)
        raise SimError(f'unknown layer type {t}')

    def _placeholder_eval(self, vals):
        if self.tag.startswith('UnsupportedOp:'):
            raise OracleNeeded(self.name)
        if self.constant is not None:
            return mn(self.constant)
        x = mn(vals[0]) if vals else None
        f = {'Sign': np.sign, 'Abs': np.abs, 'Floor': np.floor,
             'Ceil': np.ceil, 'Round': np.round, 'Sin': np.sin,
             'Cos': np.cos, 'Tan': np.tan, 'Exp': np.exp, 'Log': np.log,
             'Sqrt': np.sqrt}.get(self.tag)
        if f is not None:
            return mn(f(x))
        return x

    def _bn_eval(self, x):
        x = mn(x)
        scale, bias, mean, var, eps = self.scale, self.bias, self.mean, self.var, self.eps
        nC = self.nC
        if m_isvector(x):
            m = m_col(mean); v = m_col(var); s = m_col(scale); o = m_col(bias)
            y = (m_col(x) - m) / np.sqrt(v + eps)
            y = s * y + o
            if x.shape[0] == 1:   # row in -> row out
                y = y.T
            return mn(y)
        # axis-aware fast path: exactly one input dim equals NumChannels
        if nC > 1:
            sz = m_size(x)
            ch_dims = [k for k, d in enumerate(sz) if d == nC]
            if len(ch_dims) == 1:
                ax = ch_dims[0]
                tgt = [1] * max(2, len(sz))
                tgt[ax] = nC
                m = np.reshape(mean, tgt, order='F')
                v = np.reshape(var, tgt, order='F')
                s = np.reshape(scale, tgt, order='F')
                o = np.reshape(bias, tgt, order='F')
                return mn(((x - m) / np.sqrt(v + eps)) * s + o)
        # generic image path: [1,1,C] params broadcast over MATLAB dim 3
        if x.ndim >= 3 and x.shape[2] == nC:
            return mn(((x - mean) / np.sqrt(var + eps)) * scale + bias)
        if nC == 1:
            return mn(((x - mean.flatten()[0]) / np.sqrt(var.flatten()[0] + eps))
                      * scale.flatten()[0] + bias.flatten()[0])
        raise SimError(f'BatchNorm: ambiguous/unsupported input shape {m_size(x)} for C={nC}')

    def _ea_eval(self, x):
        x = mn(x)
        y = x
        if self.do_scale:
            s = self.scale
            if not m_isscalar(s) and not _ea_broadcastable(m_size(s), m_size(x)):
                s = _ea_trailing_align(s, x)
                if not _ea_broadcastable(m_size(s), m_size(x)):
                    s = _ea_align_to_input(s, m_size(x))
            y = m_broadcast(y, s, np.multiply)
        if self.do_offset:
            o = self.offset
            if not m_isscalar(o) and not _ea_broadcastable(m_size(o), m_size(x)):
                o = _ea_trailing_align(o, x)
                if not _ea_broadcastable(m_size(o), m_size(x)):
                    o = _ea_align_to_input(o, m_size(x))
            y = m_broadcast(y, o, np.add)
        return y


# ---------------------------------------------------------------------------
# Manifest network
# ---------------------------------------------------------------------------

class ManifestSim:
    """Loads a .nnv.mat manifest (scipy) and evaluates it layer-by-layer the
    way NNV's NN.evaluate would after load_nnv_from_mat."""

    def __init__(self, mat_path):
        s = loadmat(mat_path, struct_as_record=False, squeeze_me=False)
        weights = s['weights'][0, 0] if isinstance(s['weights'], np.ndarray) else s['weights']
        raw_layers = np.asarray(s['layers']).flatten()
        self.layers = []
        for L in raw_layers:
            if isinstance(L, np.ndarray):
                L = L.flatten()[0]
            attrs = getattr(L, 'attrs', None)
            if isinstance(attrs, np.ndarray):
                attrs = attrs.flatten()[0] if attrs.size else None
            self.layers.append(SimLayer(
                _to_str(L.type), _to_str(L.name), attrs,
                _to_str_list(getattr(L, 'weight_keys', None)), weights,
                _to_str_list(getattr(L, 'inputs', None)),
                _to_str_list(getattr(L, 'outputs', None))))
        self.input_shape = _ivec(s['input_shape'])
        self.input_name = _to_str(s['input_name'])
        self.output_name = _to_str(s['output_name'])

    def prepare_input(self, x_onnx):
        """Convert the ONNX-layout input array to the MATLAB-side form the
        manifest's input layer expects."""
        L0 = self.layers[0]
        x = np.asarray(x_onnx, dtype=np.float64)
        if L0.type == 'ImageInputLayer':
            want = _ivec(L0.in_size)
            if x.ndim == 4:
                x = x[0]
            if list(x.shape) == want:           # BHWC source: already HWC
                return mn(x)
            if x.ndim == 3 and [x.shape[1], x.shape[2], x.shape[0]] == want:
                return mn(np.transpose(x, (1, 2, 0)))   # CHW -> HWC
            raise SimError(f'input shape {x.shape} does not map to InputSize {want}')
        return mn(x.reshape(-1, 1))             # FLAT column (ONNX C-order)

    def run(self, x_onnx, oracle=None):
        """Execute all layers. Returns (values, patched, tainted):
        values: layer_name -> MATLAB-form numpy array
        patched: set of layer names whose output came from the oracle
        tainted: set of layer names transitively downstream of a patch."""
        values, patched, tainted = {}, set(), set()
        name2layer = {L.name: L for L in self.layers}
        for L in self.layers:
            if L.type in ('FeatureInputLayer', 'ImageInputLayer'):
                values[L.name] = self.prepare_input(x_onnx)
                continue
            vals = []
            missing = None
            for p in L.inputs:
                if p in values:
                    vals.append(values[p])
                else:
                    missing = p
            is_tainted = any(p in patched or p in tainted for p in L.inputs)
            try:
                if missing is not None:
                    raise OracleNeeded(f'{L.name}: missing input {missing}')
                out = L.evaluate(vals)
            except Exception as e:
                # An evaluation error on a CLEAN path is a real importer/sim
                # bug — never mask it. On a TAINTED path (downstream of an
                # unsupported-op placeholder, e.g. the collins detection head)
                # garbage shapes are expected; MATLAB xval cannot check there
                # either. OracleNeeded is patchable on any path.
                if not isinstance(e, OracleNeeded) and not is_tainted:
                    raise
                if oracle is None:
                    raise
                out = oracle(L)
                if out is None:
                    raise SimError(f'{L.name}: oracle has no value for {L.outputs} '
                                   f'(after {type(e).__name__}: {e})')
                patched.add(L.name)
                is_tainted = False   # ground truth resets the taint at this node
            if is_tainted:
                tainted.add(L.name)
            values[L.name] = out
        return values, patched, tainted


# ---------------------------------------------------------------------------
# ORT <-> MATLAB layout conversion + strict comparison
# ---------------------------------------------------------------------------

def ort_to_matlab(arr, lay):
    """Convert an onnxruntime tensor to the MATLAB-side form the manifest
    DEFINES for that layout. No permutation guessing."""
    a = np.asarray(arr, dtype=np.float64)
    if lay == 'HWC':
        if a.ndim == 4 and a.shape[0] == 1:
            a = a[0]
        if a.ndim == 3:
            return mn(np.transpose(a, (1, 2, 0)))
        return mn(a)        # non-conforming rank: treat as RAW
    if lay == 'HWCB':
        # ONNX tensor is [B,H,W,C]; the MATLAB value is the same [H,W,C]
        if a.ndim == 4 and a.shape[0] == 1:
            a = a[0]
        return mn(a)
    if lay == 'FLAT':
        return mn(a.reshape(-1, 1))     # ONNX C-order sequence as a column
    return mn(a)            # RAW: dims literally equal (trailing 1s trimmed)


def strict_diff(sim_val, ort_val, lay):
    """Orientation-STRICT max|diff|: shapes must agree under the declared
    layout convention; no permutation search. Returns np.inf on mismatch."""
    ref = ort_to_matlab(ort_val, lay)
    a = mn(np.asarray(sim_val, dtype=np.float64))
    if lay == 'FLAT':
        # contract: <=1 non-singleton dim; the SEQUENCE is the comparison
        if sum(1 for d in a.shape if d > 1) > 1 or a.size != ref.size:
            return np.inf
        return float(np.max(np.abs(a.flatten() - ref.flatten()))) if a.size else 0.0
    if a.shape != ref.shape:
        return np.inf
    if a.size == 0:
        return 0.0
    return float(np.max(np.abs(a - ref)))
