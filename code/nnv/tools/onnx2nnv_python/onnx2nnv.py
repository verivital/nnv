"""onnx2nnv.py — preprocess an ONNX file into a .mat manifest that NNV can
load directly, bypassing MATLAB's importNetworkFromONNX (which collapses
many subgraphs into opaque custom layers).

Pipeline:
  1. Load ONNX
  2. Set static input shape (from VNN-LIB if provided; else use ONNX-declared)
  3. Convert opset >= 11 if needed (some passes only work on newer opsets)
  4. shape_inference + onnxsim simplify  (constant-fold dynamic-shape noise)
  5. onnxoptimizer passes: fuse_matmul_add_bias_into_gemm, fuse_bn_into_conv,
     fuse_consecutive_transposes, eliminate_nop_*, fuse_add_bias_into_conv
  6. Walk simplified graph, emit a manifest of layer records + weights
  7. scipy.io.savemat to disk

Output structure (in MATLAB after loadmat):
  manifest.layers   — struct array, fields: type, name, attrs (struct), inputs, outputs, weight_keys
  manifest.weights  — struct, each value is a numeric array (float32)
  manifest.input_name, manifest.input_shape, manifest.input_dtype
  manifest.output_name, manifest.output_shape

Usage:
  python onnx2nnv.py path/to/model.onnx [path/to/out.mat] [--vnnlib path/to.vnnlib]

Layer record types we emit (matching NNV layer classes):
  FeatureInputLayer, ImageInputLayer (auto)
  FullyConnectedLayer
  Conv2DLayer
  BatchNormalizationLayer
  ReluLayer, LeakyReluLayer, SigmoidLayer, TanhLayer
  SoftmaxLayer (FLAGGED unsound for mid-network use)
  ElementwiseAffineLayer  (Mul/Sub scaling, ScalingLayer)
  AdditionLayer  (residual connection)
  FlattenLayer
  ReshapeLayer  (constant-shape only)
  ConcatenationLayer
  MaxPooling2DLayer, AveragePooling2DLayer, GlobalAveragePooling2DLayer
  TransposedConv2DLayer
  PlaceholderLayer (Dropout, Identity)
"""

import argparse, os, sys, hashlib, re
from dataclasses import dataclass, field
from typing import Any
import numpy as np
import onnx
from onnx import numpy_helper, shape_inference, version_converter
import onnxsim
import onnxoptimizer
from scipy.io import savemat


# ---------- tensor-layout model ----------
#
# Every tensor we emit lives in exactly one MATLAB-side layout:
#
#   'HWC'  — image tensor. ONNX [B,C,H,W] is stored MATLAB-side as an
#            [H,W,C] array (NNV's image convention; Conv/Pool/ImageInput
#            layers all assume it).
#   'HWCB' — image tensor whose ONNX form is [B,H,W,C] (keras/larq BHWC
#            graphs around the no-op layout transposes). MATLAB-side it is
#            the SAME [H,W,C] array as 'HWC'; the tag only records the ONNX
#            dim order for strict ground-truth comparison.
#   'FLAT' — feature vector. ONNX [B,F] (or [F]) is stored as a column
#            [F,1] (or row/[1,1,F] from FlattenLayer); the element SEQUENCE
#            equals the ONNX C-order flattening.
#   'RAW'  — everything else (attention/token tensors: [B,T,C], [B,h,T,d],
#            ...). Stored MATLAB-side with dims IDENTICAL to the ONNX dims,
#            element-by-element (no permutation). MATLAB trims trailing
#            singleton dims, which is harmless here (batch leads).
#
# The old importer pushed the HWC/flat machinery through transformer
# pipelines too, silently rotating token tensors relative to what the
# Transpose perms / DynamicMatmul shapes expect (vit_2023: pos-emb added in
# the wrong orientation, then 'inner dims mismatch' at attn x V). The walk()
# below tracks the layout of every tensor and only applies HWC conventions
# on actual image flows; attention/token tensors stay RAW with explicit
# C-order reshapes.
#
# Key consequences for RAW tensors:
#   * ONNX Reshape/Flatten must use C-order semantics. MATLAB reshape is
#     column-major, so we emit the standard decomposition
#         permute(reverse-dims) -> reshape(flip(target)) -> permute(reverse)
#     using the existing TransposeLayer/ReshapeLayer (exact, no new MATLAB
#     machinery).
#   * ONNX Transpose perms apply literally (TransposeLayer; never the
#     BHWC/BCHW "no-op" placeholders, which are an HWC-only convention).
#   * FullyConnectedLayer's batched fast path (last dim == in_features)
#     evaluates per-token FCs on [1,T,C] without reordering.
#   * Selector matrices (Split/Slice/Gather/Reduce*) must index the
#     COLUMN-MAJOR flattening of the RAW array (MATLAB reshape(x,[],1)),
#     not the ONNX C-order one used for FLAT flows.
#   * Per-channel biases are emitted trailing-aligned ([1,...,1,C]) so ONNX
#     broadcasting semantics hold under MATLAB implicit expansion.
#   * Softmax over the last axis of a rank-3 RAW tensor maps to NNV's
#     SoftmaxLayer 'SSC' channel softmax; other ranks are wrapped in
#     reshape-[1,L,N] / softmax / reshape-back (exact: softmax is
#     per-leading-slice over the last dim).

IMAGE_CONSUMER_OPS = {'Conv', 'ConvTranspose', 'MaxPool', 'AveragePool',
                      'GlobalAveragePool', 'BatchNormalization', 'Resize',
                      'Upsample'}

# Ops whose NNV layer always produces an HWC image output.
IMAGE_PRODUCER_OPS = {'Conv', 'ConvTranspose', 'MaxPool', 'AveragePool',
                      'GlobalAveragePool', 'Resize', 'Upsample'}


# ---------- helpers ----------

def to_array(t):
    """Convert ONNX TensorProto to numpy array."""
    return np.asarray(numpy_helper.to_array(t)).astype(np.float32)

def get_attr(node, name, default=None):
    for a in node.attribute:
        if a.name == name:
            if a.type == onnx.AttributeProto.INT:    return a.i
            if a.type == onnx.AttributeProto.FLOAT:  return a.f
            if a.type == onnx.AttributeProto.STRING: return a.s.decode('utf-8')
            if a.type == onnx.AttributeProto.INTS:   return list(a.ints)
            if a.type == onnx.AttributeProto.FLOATS: return list(a.floats)
            if a.type == onnx.AttributeProto.TENSOR: return to_array(a.t)
    return default

def _trailing_align_param(arr, x_shape):
    """Pad an ONNX broadcast parameter with LEADING singleton dims so its
    dims line up with the trailing dims of an x of rank len(x_shape) — the
    ONNX/NumPy broadcasting rule, made explicit for MATLAB implicit
    expansion (which aligns LEADING dims). Used for RAW-layout consumers."""
    a = np.asarray(arr)
    if x_shape is None:
        return a
    r = len(x_shape)
    if a.ndim < r:
        a = a.reshape((1,) * (r - a.ndim) + a.shape)
    return a


_NAME_HASH_TABLE = {}
def safe_name(s):
    """sanitize for MATLAB struct field names. MATLAB field names are limited
    to 63 characters (savemat allows up to 31 in some old code paths)."""
    s = re.sub(r'[^A-Za-z0-9_]', '_', s)
    if not s or not s[0].isalpha():
        s = 'L_' + s
    if len(s) > 28:
        # truncate + append a short hash for uniqueness
        import hashlib as _h
        digest = _h.md5(s.encode()).hexdigest()[:6]
        s = s[:21] + '_' + digest
    return s

# ---------- preprocessing ----------

def _try_emit_reduce(node, nm, value_shapes, initializers, producer,
                     add_weight, emit, op_name, raw_order=False):
    """Try to emit a ReduceSum / ReduceMean as a sparse FullyConnectedLayer.

    For an input of known shape S = [d_0, ..., d_{n-1}] with axes A and
    keepdims=K, the reduction is a linear map. Build the FC matrix:
       W[r, c] = (1 / |reduce_volume|)  if input flat-index c contributes to
                                        output flat-index r else 0
    where the "reduce_volume" is the product of dims being reduced.
    Returns True on successful emit, False on fallback.

    op_name: 'sum' or 'mean'. For sum, no division.

    raw_order: when the producer tensor is in 'RAW' layout, the MATLAB-side
    FullyConnectedLayer flattens it COLUMN-MAJOR (F-order over the ONNX
    dims), not C-order. Build the selector with F-order strides on both
    sides, and follow the FC with a ReshapeLayer restoring the RAW output
    shape (the plain MATLAB column-major reshape inverts the F-order
    flatten exactly).
    """
    in_name = node.input[0]
    in_sh = value_shapes.get(in_name)
    if in_sh is None:
        return False
    eff_sh = [int(d) if d > 0 else 1 for d in in_sh]
    n = len(eff_sh)

    # Axes from attribute (older opset) or input[1] initializer (opset>=13).
    axes_attr = None
    for a in node.attribute:
        if a.name == 'axes':
            axes_attr = list(a.ints); break
    if axes_attr is None and len(node.input) > 1 and node.input[1] in initializers:
        axes_attr = [int(x) for x in initializers[node.input[1]].flatten()]
    if axes_attr is None:
        return False

    keepdims = 1
    for a in node.attribute:
        if a.name == 'keepdims': keepdims = int(a.i); break
    noop_with_empty_axes = 0
    for a in node.attribute:
        if a.name == 'noop_with_empty_axes': noop_with_empty_axes = int(a.i); break
    if not axes_attr and not noop_with_empty_axes:
        # reduce all dims
        axes_attr = list(range(n))
    elif not axes_attr:
        # noop
        return False

    axes = sorted(((a + n) if a < 0 else a) for a in axes_attr)
    if any(ax < 0 or ax >= n for ax in axes):
        return False

    # Output shape:
    if keepdims:
        out_sh = list(eff_sh)
        for ax in axes: out_sh[ax] = 1
    else:
        out_sh = [d for i, d in enumerate(eff_sh) if i not in set(axes)]
        if not out_sh: out_sh = [1]

    flat_in = 1
    for d in eff_sh: flat_in *= d
    flat_out = 1
    for d in out_sh: flat_out *= d

    if flat_in <= 0 or flat_out <= 0 or flat_in > 100000 or flat_out > 100000:
        return False  # too big to materialize as dense FC

    # Build selector matrix.
    import numpy as _np
    W = _np.zeros((flat_out, flat_in), dtype=_np.float32)

    if raw_order:
        # F-order (column-major) strides over the RAW dims — matches how
        # MATLAB's FullyConnectedLayer flattens a RAW [d0,...,dn-1] array.
        in_strides = [1] * n
        for k in range(1, n):
            in_strides[k] = in_strides[k-1] * eff_sh[k-1]
        out_strides = [1] * len(out_sh)
        for k in range(1, len(out_sh)):
            out_strides[k] = out_strides[k-1] * out_sh[k-1]
    else:
        # Input strides (row-major / C-order, matching how we treat NNV's flat layout).
        in_strides = [1] * n
        for k in range(n - 2, -1, -1):
            in_strides[k] = in_strides[k+1] * eff_sh[k+1]
        out_strides = [1] * len(out_sh)
        for k in range(len(out_sh) - 2, -1, -1):
            out_strides[k] = out_strides[k+1] * out_sh[k+1]

    reduce_vol = 1
    for ax in axes: reduce_vol *= eff_sh[ax]
    coeff = 1.0 / reduce_vol if op_name == 'mean' else 1.0

    # For each input flat index c, find the corresponding output position.
    from itertools import product
    ranges = [range(d) for d in eff_sh]
    for in_coord in product(*ranges):
        in_flat = sum(c * s for c, s in zip(in_coord, in_strides))
        # Output coord: drop or zero the reduced axes
        if keepdims:
            out_coord = tuple(0 if i in set(axes) else in_coord[i] for i in range(n))
        else:
            out_coord = tuple(in_coord[i] for i in range(n) if i not in set(axes))
            if not out_coord: out_coord = (0,)
        if keepdims:
            out_flat = sum(c * s for c, s in zip(out_coord, out_strides))
        else:
            out_flat = sum(c * s for c, s in zip(out_coord, out_strides))
        W[out_flat, in_flat] += coeff

    b = _np.zeros(flat_out, dtype=_np.float32)
    if raw_order:
        # FC emits an F-order flat column; restore the RAW output shape with
        # a plain (column-major) ReshapeLayer so downstream RAW consumers
        # (BatchNorm, FC, ...) see the ONNX dims.
        fc_nm = safe_name(nm + '_red')
        wkey = add_weight(fc_nm, 'W', W)
        bkey = add_weight(fc_nm, 'b', b)
        emit(LayerSpec('FullyConnectedLayer', fc_nm,
                       attrs={'OutputSize': int(flat_out)},
                       inputs=[producer[in_name]],
                       outputs=[fc_nm + '_out'],
                       weight_keys=[wkey, bkey]))
        spec = LayerSpec('ReshapeLayer', nm,
                         attrs={'TargetShape': [int(d) for d in out_sh]},
                         inputs=[fc_nm],
                         outputs=[node.output[0]])
        emit(spec)
        return True
    wkey = add_weight(nm, 'W', W)
    bkey = add_weight(nm, 'b', b)
    spec = LayerSpec('FullyConnectedLayer', nm,
                     attrs={'OutputSize': int(flat_out)},
                     inputs=[producer[in_name]],
                     outputs=[node.output[0]],
                     weight_keys=[wkey, bkey])
    emit(spec)
    return True


def find_data_input(model):
    """Return the graph.input entry that is the actual data input, not an
    initializer-also-listed-as-input (e.g. weights). Some PyTorch exports
    list every initializer in graph.input too; the real input is the one
    that is NOT also an initializer."""
    init_names = {t.name for t in model.graph.initializer}
    for inp in model.graph.input:
        if inp.name not in init_names:
            return inp
    # fallback: first one
    return model.graph.input[0]

def static_input_shape(model, batch_size=1):
    """Set batch dim to a concrete value if unknown. Return new shape list."""
    inp = find_data_input(model)
    dims = inp.type.tensor_type.shape.dim
    shape = []
    for i, d in enumerate(dims):
        if d.dim_value > 0:
            shape.append(d.dim_value)
        else:
            shape.append(batch_size if i == 0 else 1)
            d.dim_value = shape[-1]
            d.dim_param = ''
    return shape

def preprocess_onnx(model, batch_size=1, target_opset=17):
    """Static shape, version convert, simplify, optimize."""
    # Static input shape. Use the ACTUAL data input (find_data_input) for BOTH the
    # shape and the name: some PyTorch exports list every initializer in graph.input,
    # so graph.input[0] can be a weight, not the data input. Pairing in_shape (from
    # find_data_input) with a graph.input[0] name would apply the shape to the wrong
    # tensor in onnxsim.simplify(overwrite_input_shapes=...). (Copilot #290 finding.)
    in_shape = static_input_shape(model, batch_size)
    in_name  = find_data_input(model).name

    # Version convert if needed
    cur_opset = model.opset_import[0].version if model.opset_import else 9
    if cur_opset < target_opset and cur_opset >= 7:
        try:
            model = version_converter.convert_version(model, target_version=target_opset)
        except Exception as e:
            print(f"  [info] opset upgrade {cur_opset} -> {target_opset} skipped: {e}", file=sys.stderr)

    # Shape inference
    try:
        model = shape_inference.infer_shapes(model)
    except Exception:
        pass

    # Simplify
    try:
        model_sim, ok = onnxsim.simplify(model, perform_optimization=True,
            overwrite_input_shapes={in_name: in_shape})
        if ok:
            model = model_sim
    except Exception as e:
        print(f"  [info] onnxsim skipped: {e}", file=sys.stderr)

    # Optimizer passes (fuses MatMul+Add->Gemm, BN+Conv, etc.)
    passes = [
        'eliminate_identity', 'eliminate_nop_dropout', 'eliminate_nop_flatten',
        'eliminate_nop_transpose', 'eliminate_nop_pad', 'eliminate_nop_concat',
        'fuse_matmul_add_bias_into_gemm',
        'fuse_add_bias_into_conv',
        'fuse_pad_into_conv',
        'fuse_consecutive_transposes',
        'eliminate_deadend',
    ]
    try:
        model = onnxoptimizer.optimize(model, passes)
    except Exception as e:
        print(f"  [info] onnxoptimizer skipped: {e}", file=sys.stderr)

    # Manual dead-node elimination: some optimizer passes leave producer
    # nodes whose outputs no other node references and which aren't graph
    # outputs (e.g. fuse_matmul_add_bias_into_gemm leaves the original
    # MatMul behind). Walk backwards from graph outputs and keep only
    # reachable nodes.
    used = set(o.name for o in model.graph.output)
    keep_node = [False] * len(model.graph.node)
    # First pass: mark direct producers of graph outputs
    output_names = {o.name for o in model.graph.output}
    for i, n in enumerate(model.graph.node):
        if any(o in output_names for o in n.output):
            keep_node[i] = True
            for inp in n.input:
                used.add(inp)
    # Iterate backward
    changed = True
    while changed:
        changed = False
        for i, n in enumerate(model.graph.node):
            if keep_node[i]: continue
            if any(o in used for o in n.output):
                keep_node[i] = True
                for inp in n.input:
                    if inp not in used:
                        used.add(inp); changed = True
    new_nodes = [n for i, n in enumerate(model.graph.node) if keep_node[i]]
    if len(new_nodes) < len(model.graph.node):
        del model.graph.node[:]
        model.graph.node.extend(new_nodes)

    return model


# ---------- graph walking ----------

@dataclass
class LayerSpec:
    type: str
    name: str
    attrs: dict = field(default_factory=dict)
    inputs: list = field(default_factory=list)   # list of producer_layer_names (output_index always 0 in our IR)
    outputs: list = field(default_factory=list)  # list of output tensor names produced by this layer
    weight_keys: list = field(default_factory=list)

def walk(model):
    """Walk simplified ONNX graph and emit list of LayerSpecs + dict of weights."""
    initializers = {t.name: to_array(t) for t in model.graph.initializer}

    # Hoist any Constant nodes' values into the initializers dict — onnxsim
    # sometimes leaves Constants behind (especially for opset < 13 models).
    for node in model.graph.node:
        if node.op_type == 'Constant':
            v = get_attr(node, 'value', None)
            if v is None:
                v = get_attr(node, 'value_float', None)
                if v is not None:
                    v = np.array([v], dtype=np.float32)
            if v is not None:
                for out in node.output:
                    initializers[out] = v
    value_shapes = {}
    for vi in list(model.graph.input) + list(model.graph.value_info) + list(model.graph.output):
        if vi.type.tensor_type.shape.dim:
            value_shapes[vi.name] = [d.dim_value if d.dim_value > 0 else -1
                                      for d in vi.type.tensor_type.shape.dim]

    data_inp = find_data_input(model)
    inp_name  = data_inp.name
    inp_shape = value_shapes.get(inp_name,
        [d.dim_value if d.dim_value > 0 else 1 for d in data_inp.type.tensor_type.shape.dim])

    # Track which tensor name is produced by which layer
    producer = {}    # tensor_name -> layer_name

    layers = []
    weights = {}    # weight_key -> ndarray

    # Input layer: pick FeatureInputLayer if input is flat (or the model has
    # no spatial/conv ops), else ImageInputLayer.
    chw = [d for d in inp_shape if d > 0]
    has_spatial_op = any(n.op_type in ('Conv', 'ConvTranspose', 'MaxPool', 'AveragePool',
                                       'GlobalAveragePool', 'Resize')
                         for n in model.graph.node)
    # Detect BHWC layout: rank-4 input + first non-Constant node is a Transpose
    # with perm [0, 3, 1, 2] (BHWC -> BCHW). In NNV's HWC convention, that
    # leading transpose is a no-op, so drop it and keep input shape as HWC.
    bhwc_input = False
    drop_first_transpose = False
    if len(inp_shape) == 4 and has_spatial_op:
        for n in model.graph.node:
            if n.op_type == 'Constant':
                continue
            if n.op_type == 'Transpose' and inp_name in n.input:
                perm_attr = get_attr(n, 'perm', [])
                if list(perm_attr) == [0, 3, 1, 2]:
                    bhwc_input = True
                    drop_first_transpose = n.output[0]
            break
    if len([d for d in inp_shape[1:] if d > 1]) <= 1 or not has_spatial_op:
        # flat feature input (single non-batch dim, or no spatial ops in graph)
        # For rank-1 ONNX inputs (no batch dim), use the full shape; for
        # higher-rank inputs, drop the leading batch dim.
        if len(inp_shape) <= 1:
            in_size = int(np.prod([d for d in inp_shape if d > 0]))
        else:
            in_size = int(np.prod([d for d in inp_shape[1:] if d > 0]))
        if in_size <= 0:
            in_size = 1
        spec = LayerSpec('FeatureInputLayer', 'input',
                         attrs={'InputSize': in_size}, outputs=[inp_name])
    else:
        # image input
        if len(inp_shape) == 4:
            if bhwc_input:
                # [B, H, W, C] -> NNV [H, W, C]
                sz = [int(inp_shape[1]), int(inp_shape[2]), int(inp_shape[3])]
            else:
                # [B, C, H, W] -> NNV [H, W, C]
                sz = [int(inp_shape[2]), int(inp_shape[3]), int(inp_shape[1])]
        else:
            sz = [int(d) for d in inp_shape[1:]]
        spec = LayerSpec('ImageInputLayer', 'input',
                         attrs={'InputSize': sz}, outputs=[inp_name])
    layers.append(spec)
    producer[inp_name] = 'input'

    # Per-tensor MATLAB-side layout tracking ('HWC' | 'FLAT' | 'RAW') — see
    # the layout-model comment at the top of this file.
    layout = {}
    if spec.type == 'ImageInputLayer':
        layout[inp_name] = 'HWCB' if bhwc_input else 'HWC'
    else:
        layout[inp_name] = 'FLAT'

    def tensor_layout(tname):
        return layout.get(tname, 'FLAT')

    def node_in_layout(node):
        """Layout of the first input we know about (dynamic tensors only)."""
        for x in node.input:
            if x in layout:
                return layout[x]
        return 'FLAT'

    # If we detected a BHWC input layout, the leading BHWC->BCHW Transpose
    # is a no-op in NNV's HWC convention. Skip it by mapping its output
    # name to 'input' directly so downstream Conv layers wire correctly.
    if drop_first_transpose:
        producer[drop_first_transpose] = 'input'
        layout[drop_first_transpose] = 'HWC'   # post-transpose ONNX form is BCHW
        # value_shapes for the transpose's output should match input HWC.
    # Also register any other graph inputs (some ONNX exports have multiple
    # inputs, e.g. multi-input attention models). They map to PlaceholderLayer
    # since NNV networks have a single data input; downstream ops referencing
    # them get a sensible passthrough.
    init_names = {t.name for t in model.graph.initializer}
    for extra_inp in model.graph.input:
        if extra_inp.name in init_names: continue
        if extra_inp.name == inp_name: continue
        ph_name = safe_name('extra_input_' + extra_inp.name)
        ph_spec = LayerSpec('PlaceholderLayer', ph_name,
                            attrs={'OriginalOp': 'extra_graph_input'},
                            outputs=[extra_inp.name])
        layers.append(ph_spec)
        producer[extra_inp.name] = ph_name
        layout[extra_inp.name] = 'FLAT'

    # Helper: register a layer + its outputs
    def emit(spec):
        layers.append(spec)
        for out in spec.outputs:
            producer[out] = spec.name
        return spec

    # Helper: ensure a producer exists for a tensor name. If the tensor is
    # only known as an initializer (e.g. CLS token fed as data input to
    # Concat), emit a synthetic PlaceholderLayer that returns the constant.
    # Returns the producer layer name (or None if the tensor truly cannot be
    # produced — caller can choose to drop it).
    def ensure_producer(tensor_name):
        if tensor_name in producer:
            return producer[tensor_name]
        if tensor_name in initializers:
            arr = initializers[tensor_name]
            const_name = safe_name('const_' + tensor_name)
            wkey = add_weight(const_name, 'value', np.asarray(arr, dtype=np.float32))
            spec = LayerSpec('PlaceholderLayer', const_name,
                             attrs={'OriginalOp': 'Constant'},
                             outputs=[tensor_name],
                             weight_keys=[wkey])
            emit(spec)
            return const_name
        return None

    _wkey_counter = [0]
    def add_weight(layer_name, role, arr):
        key = safe_name(f"{layer_name}__{role}")
        if key in weights:
            # collision — append counter
            _wkey_counter[0] += 1
            key = f"{key[:25]}_{_wkey_counter[0]}"
        weights[key] = np.ascontiguousarray(arr.astype(np.float32))
        return key

    def emit_conorder_reshape(nm, src_tensor, tgt_static, out_tensor, pre_perm0):
        """Emit an EXACT ONNX (C-order) reshape for a RAW-destined tensor as a
        composition of existing NNV layers:

            [pre-permute] -> reshape(flip(tgt)) -> permute(reverse(tgt-rank))

        MATLAB's reshape is column-major; flattening the REVERSE-permuted
        array column-major equals the ONNX C-order flattening, and refilling
        flip(tgt) then reverse-permuting realizes the C-order refill.

        pre_perm0: 0-indexed permutation applied to the MATLAB-side input
        array so its column-major flatten equals the ONNX C-order flatten of
        the source tensor: reversed(range(rank)) for RAW inputs, [1,0,2] for
        HWC image inputs ([H,W,C] -> [W,H,C]), None for FLAT columns (their
        element sequence is already the C-order one).

        The final layer carries `nm` and `out_tensor` so downstream wiring
        and debug-name mapping stay 1:1 with the ONNX node.
        """
        cur = producer[src_tensor]
        if pre_perm0 is not None and list(pre_perm0) != list(range(len(pre_perm0))):
            pre_nm = safe_name(nm + '_cpre')
            emit(LayerSpec('TransposeLayer', pre_nm,
                           attrs={'Perm': [int(p) for p in pre_perm0]},
                           inputs=[cur], outputs=[pre_nm + '_out']))
            cur = pre_nm
        rs_nm = safe_name(nm + '_crs')
        emit(LayerSpec('ReshapeLayer', rs_nm,
                       attrs={'TargetShape': [int(d) for d in reversed(tgt_static)]},
                       inputs=[cur], outputs=[rs_nm + '_out']))
        cur = rs_nm
        post = list(range(len(tgt_static) - 1, -1, -1))
        spec = LayerSpec('TransposeLayer', nm, attrs={'Perm': post},
                         inputs=[cur], outputs=[out_tensor])
        emit(spec)
        return spec

    def static_reshape_target(node, tgt):
        """Resolve an ONNX Reshape target ([-1]/0 sentinels) to a fully
        static shape, preferring the inferred output shape."""
        out_sh = value_shapes.get(node.output[0])
        if out_sh and len(out_sh) == len(tgt) and all(d > 0 for d in out_sh):
            return [int(d) for d in out_sh]
        in_sh = value_shapes.get(node.input[0])
        res = [int(t) for t in tgt]
        if in_sh and all(d > 0 for d in in_sh):
            for k, t in enumerate(res):
                if t == 0:
                    res[k] = int(in_sh[k]) if k < len(in_sh) else 1
            if res.count(-1) == 1:
                known = 1
                for t in res:
                    if t != -1:
                        known *= t
                tot = int(np.prod(in_sh))
                if known > 0 and tot % known == 0:
                    res[res.index(-1)] = tot // known
            if all(t > 0 for t in res):
                return res
        return None

    def to_hwc_if_image_bias(arr):
        """Convert ONNX-shape (broadcasting) bias arrays to HWC layout if they
        look like image-shaped per-pixel constants.

        Recognizes:
          - rank-4 [1, C, H, W]  -> [H, W, C] when C in (1, 3, 4) AND H==W
            (avoids mis-firing on transformer per-token biases like [1, 1, T, C]).
        Leaves rank-3, 1D, and scalar arrays alone — rank-3 transformer biases
        like [1, T, C] (positional embeddings) and [1, 1, C] (CLS tokens) need
        to broadcast naturally without permutation.
        """
        a = np.asarray(arr)
        if (a.ndim == 4 and a.shape[0] == 1 and a.shape[1] in (1, 3, 4)
                and a.shape[2] == a.shape[3]):
            # [1, C, H, W] -> [H, W, C], only when H==W (genuine spatial bias)
            return np.transpose(a[0], (1, 2, 0))
        return a

    # Pre-pass: detect Larq STE-Sign patterns
    #   Sign(x) -> Add(_, +small_const) -> Sign(_)
    # Larq emits this for binarized weights/activations; the inner Sign
    # gives -1/0/+1, then Add(>0) shifts so the second Sign yields +1
    # for x >= 0 (polar/zero-to-pos-one mode). Fold the trio into a single
    # SignLayer in polar mode acting on the ORIGINAL x.
    ste_skip = set()        # node names to skip (inner Sign and Add)
    ste_polar = {}          # outer Sign name -> producer of original x
    nodes_by_name = {n.name or f"_n{i}": n for i, n in enumerate(model.graph.node)}
    out_to_node = {}
    for n in model.graph.node:
        for o in n.output:
            out_to_node[o] = n
    for n in model.graph.node:
        if n.op_type != 'Sign':
            continue
        # Look at producer of n.input[0]; should be Add
        prev = out_to_node.get(n.input[0])
        if prev is None or prev.op_type != 'Add':
            continue
        # Add must have one initializer-side and one Sign-side
        a0, a1 = prev.input[0], prev.input[1]
        if a0 in initializers and a1 in initializers:
            continue
        if a0 in initializers:
            sign_in_name, const_name = a1, a0
        elif a1 in initializers:
            sign_in_name, const_name = a0, a1
        else:
            continue
        const_val = initializers[const_name]
        if const_val.size != 1:
            continue
        if float(const_val.flatten()[0]) <= 0:
            continue
        inner = out_to_node.get(sign_in_name)
        if inner is None or inner.op_type != 'Sign':
            continue
        # Pattern matched. Mark inner Sign and Add to skip; remember the
        # outer Sign should rewrite its input to the inner Sign's input.
        ste_skip.add(inner.name or f"_n_inner_{id(inner)}")
        ste_skip.add(prev.name or f"_n_add_{id(prev)}")
        ste_polar[n.name or f"_n_outer_{id(n)}"] = inner.input[0]

    # Pre-pass: detect ml4acopf-style Sigmoid PWL pattern:
    #   Unsqueeze(input, axis=last) -> Sub(_, thresholds[N]) -> Relu(_)
    #   -> MatMul(_, delta_m[N]) -> Add(_, b)
    # which approximates Sigmoid as a sum of N relu kinks. ONNX broadcasting
    # does outer subtraction [in_size,1] - [N] = [in_size, N], which NNV's
    # ElementwiseAffineLayer can't do. Detect the trio and replace with an
    # explicit expansion + per-element-PWL chain.
    pwl_skip = set()
    pwl_emit = {}
    for n in model.graph.node:
        if n.op_type != 'Add':
            continue
        # Add(_, b_scalar) preceded by MatMul(...)
        b_init = None
        x_in = None
        for inp in n.input:
            if inp in initializers and initializers[inp].size == 1:
                b_init = initializers[inp]
            else:
                x_in = inp
        if b_init is None or x_in is None:
            continue
        mm = out_to_node.get(x_in)
        if mm is None or mm.op_type != 'MatMul':
            continue
        # delta_m: one MatMul input must be initializer with shape (N,) or (N,1)
        dm_init = None
        relu_in = None
        for inp in mm.input:
            if inp in initializers and initializers[inp].ndim <= 2 and initializers[inp].size > 1:
                dm_init = initializers[inp]
                dm_name = inp
            else:
                relu_in = inp
        if dm_init is None or relu_in is None:
            continue
        relu = out_to_node.get(relu_in)
        if relu is None or relu.op_type != 'Relu':
            continue
        sub = out_to_node.get(relu.input[0])
        if sub is None or sub.op_type != 'Sub':
            continue
        # thresholds (init) and unsqueeze input (data)
        th_init = None
        unsq_in = None
        for inp in sub.input:
            if inp in initializers:
                th_init = initializers[inp]
                th_name = inp
            else:
                unsq_in = inp
        if th_init is None or unsq_in is None:
            continue
        if th_init.size != dm_init.size:
            continue
        unsq = out_to_node.get(unsq_in)
        if unsq is None or unsq.op_type != 'Unsqueeze':
            continue
        # Pattern matched. Skip Unsqueeze, Sub, Relu, MatMul, Add inner;
        # emit explicit PWL chain when we hit the Add (last node).
        pwl_skip.add(unsq.name or f"_n_unsq_{id(unsq)}")
        pwl_skip.add(sub.name or f"_n_sub_{id(sub)}")
        pwl_skip.add(relu.name or f"_n_relu_{id(relu)}")
        pwl_skip.add(mm.name or f"_n_mm_{id(mm)}")
        # The outer Add will trigger the PWL emit
        pwl_emit[n.name or f"_n_add_{id(n)}"] = {
            'data_in': unsq.input[0],
            'thresholds': th_init.flatten(),
            'delta_m': dm_init.flatten(),
            'b': float(b_init.flatten()[0]),
            'output': n.output[0],
        }

    # Build the layer chain
    for i, node in enumerate(model.graph.node):
        op = node.op_type
        nm = safe_name(node.name) if node.name else f"{op}_{i}"
        node_key = node.name or f"_n_{op}_{i}"
        # Default layout propagation: image producers emit HWC; everything
        # else inherits the first known input layout. Handlers that change
        # representation (Reshape boundary, FC, Flatten, ...) overwrite
        # these entries below. Assign BEFORE any skip/rewrite so fused
        # patterns (STE Sign, PWL) still propagate to their consumers.
        _in_lay = node_in_layout(node)
        for _o in node.output:
            layout[_o] = 'HWC' if op in IMAGE_PRODUCER_OPS else _in_lay
        # Larq STE pattern: skip the inner Sign + Add entirely; the outer
        # Sign takes the original x and emits a polar SignLayer.
        if node_key in ste_skip:
            continue
        # ml4acopf Sigmoid-via-PWL pattern: skip inner ops, emit chain at Add
        if node_key in pwl_skip:
            continue
        if node_key in pwl_emit:
            info = pwl_emit[node_key]
            data_in = info['data_in']
            thresholds = info['thresholds']
            delta_m = info['delta_m']
            b_val = info['b']
            out_name = info['output']
            # Determine input size from value_shapes / initializer
            in_sh = value_shapes.get(data_in)
            if in_sh is None and data_in in initializers:
                in_sh = list(initializers[data_in].shape)
            if in_sh is None:
                # Can't determine; bail to placeholder
                spec = LayerSpec('PlaceholderLayer', nm,
                                 inputs=[producer[data_in]],
                                 outputs=[out_name])
                emit(spec)
                print(f"  [warn] PWL pattern at {nm} could not infer input size; placeholder", file=sys.stderr)
                continue
            in_size = int(np.prod([d for d in in_sh if d > 0]))
            N = int(thresholds.size)
            # Step 1: Expand input via FC. W_expand[i*N + j, i] = 1 for all j
            # so output[i*N + j] = x[i].
            W_exp = np.zeros((in_size * N, in_size), dtype=np.float32)
            for ii in range(in_size):
                for jj in range(N):
                    W_exp[ii * N + jj, ii] = 1.0
            b_exp = np.zeros(in_size * N, dtype=np.float32)
            exp_nm = safe_name(f"{nm}_pwl_expand")
            wkey = add_weight(exp_nm, 'W', W_exp)
            bkey = add_weight(exp_nm, 'b', b_exp)
            spec = LayerSpec('FullyConnectedLayer', exp_nm,
                             attrs={'OutputSize': int(in_size * N)},
                             inputs=[producer[data_in]],
                             outputs=[exp_nm + '_out'],
                             weight_keys=[wkey, bkey])
            emit(spec)
            # Step 2: subtract tiled thresholds
            tiled_th = np.tile(thresholds.astype(np.float32), in_size)
            sub_nm = safe_name(f"{nm}_pwl_sub")
            sk = add_weight(sub_nm, 'scale', np.ones(1, dtype=np.float32))
            bk = add_weight(sub_nm, 'bias', (-tiled_th).astype(np.float32))
            spec = LayerSpec('ElementwiseAffineLayer', sub_nm,
                             attrs={'DoScale': False, 'DoOffset': True},
                             inputs=[exp_nm], outputs=[sub_nm + '_out'],
                             weight_keys=[sk, bk])
            emit(spec)
            # Step 3: ReLU
            relu_nm = safe_name(f"{nm}_pwl_relu")
            spec = LayerSpec('ReluLayer', relu_nm,
                             inputs=[sub_nm], outputs=[relu_nm + '_out'])
            emit(spec)
            # Step 4: sum-with-weights via FC. W_sum[i, i*N + j] = delta_m[j]
            W_sum = np.zeros((in_size, in_size * N), dtype=np.float32)
            for ii in range(in_size):
                for jj in range(N):
                    W_sum[ii, ii * N + jj] = float(delta_m[jj])
            b_sum = np.full(in_size, b_val, dtype=np.float32)
            wkey = add_weight(nm, 'W', W_sum)
            bkey = add_weight(nm, 'b', b_sum)
            spec = LayerSpec('FullyConnectedLayer', nm,
                             attrs={'OutputSize': int(in_size)},
                             inputs=[relu_nm],
                             outputs=[out_name],
                             weight_keys=[wkey, bkey])
            emit(spec)
            continue

        try:
            if op == 'Gemm':
                # y = alpha * A @ B + beta * C   (with optional transA/transB)
                alpha = get_attr(node, 'alpha', 1.0)
                beta  = get_attr(node, 'beta', 1.0)
                transA = get_attr(node, 'transA', 0)
                transB = get_attr(node, 'transB', 0)
                inputs = node.input
                A_name, B_name = inputs[0], inputs[1]
                C_name = inputs[2] if len(inputs) > 2 else None
                B = initializers.get(B_name)
                C = initializers.get(C_name) if C_name else None
                if B is None:
                    raise RuntimeError(f"Gemm {nm}: B is not an initializer; need MatMul-style handling")
                # NNV FullyConnectedLayer expects W of shape [out, in]. ONNX Gemm with
                # transB=1 means B is [out, in] already; transB=0 means B is [in, out] -> transpose.
                W = (alpha * B).astype(np.float32)
                if transB == 0:
                    W = W.T
                b = (beta * C).astype(np.float32) if C is not None else np.zeros(W.shape[0], np.float32)
                wkey = add_weight(nm, 'W', W)
                bkey = add_weight(nm, 'b', b)
                spec = LayerSpec('FullyConnectedLayer', nm,
                                 attrs={'OutputSize': int(W.shape[0])},
                                 inputs=[producer[A_name]],
                                 outputs=[node.output[0]],
                                 weight_keys=[wkey, bkey])
                emit(spec)
                layout[node.output[0]] = 'FLAT'
                continue

            if op == 'MatMul':
                # If either input is initializer -> FullyConnectedLayer
                A_name, B_name = node.input[0], node.input[1]
                A = initializers.get(A_name)
                B = initializers.get(B_name)
                if B is None and A is None:
                    # Both dynamic: real DynamicMatmulLayer (attention Q*K^T,
                    # attn@V). Use ensure_producer so any initializer that's
                    # somehow not in `initializers` (e.g. a Constant op output)
                    # gets a synthetic PlaceholderLayer(Constant) producer.
                    ins = []
                    for x in node.input:
                        p = ensure_producer(x)
                        if p is not None:
                            ins.append(p)
                    spec = LayerSpec('DynamicMatmulLayer', nm,
                                     attrs={'NumInputs': len(ins)},
                                     inputs=ins,
                                     outputs=[node.output[0]])
                    emit(spec); continue
                if B is not None:
                    # ONNX MatMul: C = A @ B; A is [..., K]
                    #   * B rank-2 [K, M]  -> NNV W = B.T (shape [M, K])
                    #   * B rank-1 [K]     -> W shape [1, K]
                    if B.ndim == 1:
                        W = B.reshape(1, -1).astype(np.float32)
                    else:
                        W = B.T.astype(np.float32)
                    in_producer = producer[A_name]
                else:
                    # A is constant, B is dynamic (e.g. ml4acopf
                    # Constant_47[14,11] @ Transpose_output_0[11,1]).
                    # Output: [14, 1] = A applied to B as a column vector.
                    # NNV FC: W = A, applied to B (treated as input).
                    if A.ndim == 1:
                        W = A.reshape(1, -1).astype(np.float32)
                    else:
                        W = A.astype(np.float32)
                    in_producer = producer[B_name]
                b = np.zeros(W.shape[0], np.float32)
                wkey = add_weight(nm, 'W', W)
                bkey = add_weight(nm, 'b', b)
                spec = LayerSpec('FullyConnectedLayer', nm,
                                 attrs={'OutputSize': int(W.shape[0])},
                                 inputs=[in_producer],
                                 outputs=[node.output[0]],
                                 weight_keys=[wkey, bkey])
                emit(spec)
                # Token-batched FC on a RAW [B,T,C] input keeps RAW layout
                # (MATLAB FC's last-dim fast path); everything else flattens
                # to a feature column.
                out_sh_mm = value_shapes.get(node.output[0])
                if (B is not None and tensor_layout(A_name) == 'RAW'
                        and out_sh_mm is not None and len(out_sh_mm) >= 3):
                    layout[node.output[0]] = 'RAW'
                else:
                    layout[node.output[0]] = 'FLAT'
                continue

            if op == 'Conv':
                W = initializers[node.input[1]]
                b = initializers[node.input[2]] if len(node.input) > 2 else np.zeros(W.shape[0], np.float32)
                kernel_shape = get_attr(node, 'kernel_shape')
                strides      = get_attr(node, 'strides', [1]*len(kernel_shape))
                pads         = get_attr(node, 'pads', [0]*(2*len(kernel_shape)))
                dilations    = get_attr(node, 'dilations', [1]*len(kernel_shape))
                groups       = get_attr(node, 'group', 1)
                wkey = add_weight(nm, 'W', W)
                bkey = add_weight(nm, 'b', b)
                spec = LayerSpec('Conv2DLayer', nm,
                                 attrs={'KernelSize': list(map(int, kernel_shape)),
                                        'Strides': list(map(int, strides)),
                                        'Pads': list(map(int, pads)),
                                        'Dilations': list(map(int, dilations)),
                                        'Groups': int(groups),
                                        'NumFilters': int(W.shape[0]),
                                        'NumChannels': int(W.shape[1] * groups)},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]],
                                 weight_keys=[wkey, bkey])
                emit(spec)
                continue

            if op == 'ConvTranspose':
                W = initializers[node.input[1]]
                b = initializers[node.input[2]] if len(node.input) > 2 else np.zeros(W.shape[1], np.float32)
                kernel_shape = get_attr(node, 'kernel_shape')
                strides = get_attr(node, 'strides', [1]*len(kernel_shape))
                pads = get_attr(node, 'pads', [0]*(2*len(kernel_shape)))
                wkey = add_weight(nm, 'W', W); bkey = add_weight(nm, 'b', b)
                spec = LayerSpec('TransposedConv2DLayer', nm,
                                 attrs={'KernelSize': list(map(int, kernel_shape)),
                                        'Strides': list(map(int, strides)),
                                        'Pads': list(map(int, pads))},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]],
                                 weight_keys=[wkey, bkey])
                emit(spec); continue

            if op == 'BatchNormalization':
                scale = initializers[node.input[1]]
                B     = initializers[node.input[2]]
                mean  = initializers[node.input[3]]
                var   = initializers[node.input[4]]
                eps_v = get_attr(node, 'epsilon', 1e-5)
                k_s = add_weight(nm, 'scale', scale); k_b = add_weight(nm, 'bias', B)
                k_m = add_weight(nm, 'mean', mean);   k_v = add_weight(nm, 'var', var)
                spec = LayerSpec('BatchNormalizationLayer', nm,
                                 attrs={'Epsilon': float(eps_v)},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]],
                                 weight_keys=[k_s, k_b, k_m, k_v])
                emit(spec); continue

            if op == 'Relu':
                spec = LayerSpec('ReluLayer', nm,
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'LeakyRelu':
                alpha = get_attr(node, 'alpha', 0.01)
                spec = LayerSpec('LeakyReluLayer', nm, attrs={'Scale': float(alpha)},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Sigmoid':
                spec = LayerSpec('SigmoidLayer', nm,
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Tanh':
                spec = LayerSpec('TanhLayer', nm,
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Softmax':
                axis = get_attr(node, 'axis', -1)
                in_sh_sm = value_shapes.get(node.input[0])
                if tensor_layout(node.input[0]) == 'RAW' and in_sh_sm \
                        and all(d > 0 for d in in_sh_sm):
                    # RAW tensor: NNV's SoftmaxLayer computes a channel ('SSC')
                    # softmax over MATLAB dim 3 for rank-3 inputs and errors on
                    # other ranks. ONNX last-axis softmax on a rank-3 RAW
                    # [1,L,N] maps to it directly; for other ranks wrap in
                    # reshape-[1,L,N] / softmax / reshape-back (exact: softmax
                    # is per-leading-slice over the last dim, and the plain
                    # column-major reshape repacks leading dims bijectively
                    # while preserving the last dim).
                    r_sm = len(in_sh_sm)
                    ax_sm = int(axis) % r_sm
                    if ax_sm == r_sm - 1:
                        if r_sm == 3:
                            spec = LayerSpec('SoftmaxLayer', nm,
                                             attrs={'Axis': int(axis)},
                                             inputs=[producer[node.input[0]]],
                                             outputs=[node.output[0]])
                            emit(spec)
                        else:
                            Lp = int(np.prod(in_sh_sm[:-1]))
                            Nl = int(in_sh_sm[-1])
                            pre_nm = safe_name(nm + '_smpre')
                            emit(LayerSpec('ReshapeLayer', pre_nm,
                                           attrs={'TargetShape': [1, Lp, Nl]},
                                           inputs=[producer[node.input[0]]],
                                           outputs=[pre_nm + '_out']))
                            sm_nm = safe_name(nm + '_smax')
                            emit(LayerSpec('SoftmaxLayer', sm_nm,
                                           attrs={'Axis': int(axis)},
                                           inputs=[pre_nm],
                                           outputs=[sm_nm + '_out']))
                            spec = LayerSpec('ReshapeLayer', nm,
                                             attrs={'TargetShape': [int(d) for d in in_sh_sm]},
                                             inputs=[sm_nm],
                                             outputs=[node.output[0]])
                            emit(spec)
                        layout[node.output[0]] = 'RAW'
                        continue
                    print(f"  [warn] Softmax {nm}: RAW input with non-last axis "
                          f"{axis}; legacy emission (likely wrong dim)", file=sys.stderr)
                spec = LayerSpec('SoftmaxLayer', nm, attrs={'Axis': int(axis)},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Add':
                # Two cases: (a) bias add (one operand is initializer): fold as ElementwiseAffine
                # (b) residual: emit AdditionLayer
                a_init = node.input[0] in initializers
                b_init = node.input[1] in initializers
                if a_init or b_init:
                    raw_bias = initializers[node.input[1]] if b_init else initializers[node.input[0]]
                    other_in = node.input[0] if b_init else node.input[1]
                    if tensor_layout(other_in) == 'RAW':
                        # RAW consumer: store the parameter trailing-aligned
                        # ([1,...,1,*shape]) so ONNX trailing-dim broadcasting
                        # holds under MATLAB implicit expansion (e.g. a [T,C]
                        # pos-emb or a [C] per-channel bias against [1,T,C]).
                        bias_arr = _trailing_align_param(
                            raw_bias, value_shapes.get(other_in))
                    else:
                        # Preserve multi-dim shape for broadcasting (e.g. positional
                        # embedding [1, T, C], CLS bias [1, 1, C]); permute genuine
                        # image biases via to_hwc_if_image_bias.
                        bias_arr = to_hwc_if_image_bias(raw_bias)
                    bkey = add_weight(nm, 'bias', bias_arr)
                    scale_key = add_weight(nm, 'scale', np.ones(1, dtype=np.float32))
                    spec = LayerSpec('ElementwiseAffineLayer', nm,
                                     attrs={'DoScale': False, 'DoOffset': True},
                                     inputs=[producer[other_in]],
                                     outputs=[node.output[0]],
                                     weight_keys=[scale_key, bkey])
                else:
                    spec = LayerSpec('AdditionLayer', nm,
                                     inputs=[producer[node.input[0]],
                                             producer[node.input[1]]],
                                     outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Sub':
                # bias subtract: one operand is initializer -> ElementwiseAffine
                if node.input[1] in initializers:
                    if tensor_layout(node.input[0]) == 'RAW':
                        bias_arr = _trailing_align_param(
                            -initializers[node.input[1]],
                            value_shapes.get(node.input[0]))
                    else:
                        bias_arr = to_hwc_if_image_bias(-initializers[node.input[1]])
                    bkey = add_weight(nm, 'bias', bias_arr)
                    scale_key = add_weight(nm, 'scale', np.ones_like(bias_arr).flatten()[:1])
                    spec = LayerSpec('ElementwiseAffineLayer', nm,
                                     attrs={'DoScale': False, 'DoOffset': True},
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]],
                                     weight_keys=[scale_key, bkey])
                    emit(spec); continue
                if node.input[0] in initializers:
                    # const - x  =  -1*x + const
                    if tensor_layout(node.input[1]) == 'RAW':
                        bias_arr = _trailing_align_param(
                            initializers[node.input[0]],
                            value_shapes.get(node.input[1]))
                    else:
                        bias_arr = to_hwc_if_image_bias(initializers[node.input[0]])
                    bkey = add_weight(nm, 'bias', bias_arr)
                    scale_key = add_weight(nm, 'scale', -np.ones(1, dtype=np.float32))
                    spec = LayerSpec('ElementwiseAffineLayer', nm,
                                     attrs={'DoScale': True, 'DoOffset': True},
                                     inputs=[producer[node.input[1]]],
                                     outputs=[node.output[0]],
                                     weight_keys=[scale_key, bkey])
                    emit(spec); continue
                # dynamic - dynamic: emit AdditionLayer with first input passthrough
                # and second input scaled by -1 via an inserted ElementwiseAffine
                neg_name = f"{nm}__neg"
                # we don't know the size, so use a scalar -1 broadcast (OK only for
                # element-wise; reach() may be inexact for shape-changing ops)
                neg_scale = -np.ones(1, dtype=np.float32)
                neg_bias  = np.zeros(1, dtype=np.float32)
                sk = add_weight(neg_name, 'scale', neg_scale)
                bk = add_weight(neg_name, 'bias', neg_bias)
                neg_spec = LayerSpec('ElementwiseAffineLayer', neg_name,
                    attrs={'DoScale': True, 'DoOffset': False},
                    inputs=[producer[node.input[1]]], outputs=[neg_name + '_out'],
                    weight_keys=[sk, bk])
                emit(neg_spec)
                spec = LayerSpec('AdditionLayer', nm,
                    inputs=[producer[node.input[0]], neg_name],
                    outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Mul':
                # scaling: one operand is initializer -> ElementwiseAffine with offset=0
                a_init = node.input[0] in initializers
                b_init = node.input[1] in initializers
                if a_init or b_init:
                    scale_arr = initializers[node.input[1]] if b_init else initializers[node.input[0]]
                    other_in = node.input[0] if b_init else node.input[1]
                    if tensor_layout(other_in) == 'RAW' and scale_arr.size > 1:
                        scale_out = _trailing_align_param(
                            scale_arr, value_shapes.get(other_in))
                    else:
                        scale_out = scale_arr.flatten()
                    skey = add_weight(nm, 'scale', scale_out)
                    bkey = add_weight(nm, 'bias', np.zeros(1, dtype=np.float32))
                    spec = LayerSpec('ElementwiseAffineLayer', nm,
                                     attrs={'DoScale': True, 'DoOffset': False},
                                     inputs=[producer[other_in]],
                                     outputs=[node.output[0]],
                                     weight_keys=[skey, bkey])
                    emit(spec); continue
                # both operands dynamic — emit ElementwiseProductLayer for
                # element-wise (Hadamard) product. Used in Lyapunov bilinear
                # forms, attention scaling, etc.
                ins = [producer[x] for x in node.input if x in producer]
                spec = LayerSpec('ElementwiseProductLayer', nm,
                                 attrs={'NumInputs': len(ins)},
                                 inputs=ins,
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Flatten':
                axis = get_attr(node, 'axis', 1)
                if tensor_layout(node.input[0]) == 'RAW':
                    # ONNX Flatten(axis) == C-order reshape to
                    # [prod(d[:axis]), prod(d[axis:])]; the legacy
                    # FlattenLayer assumes image (HWC) input and flattens ALL
                    # dims, which is wrong for RAW tensors (e.g. the ViT
                    # [1,h,T,T] -> [h*T, T] pre-softmax flatten).
                    in_sh_f = value_shapes.get(node.input[0])
                    out_sh_f = value_shapes.get(node.output[0])
                    if in_sh_f and all(d > 0 for d in in_sh_f):
                        if not (out_sh_f and len(out_sh_f) == 2
                                and all(d > 0 for d in out_sh_f)):
                            ax_f = int(axis) % (len(in_sh_f) + 1)
                            out_sh_f = [int(np.prod(in_sh_f[:ax_f])) if ax_f > 0 else 1,
                                        int(np.prod(in_sh_f[ax_f:])) if ax_f < len(in_sh_f) else 1]
                        pre0 = list(range(len(in_sh_f) - 1, -1, -1))
                        emit_conorder_reshape(nm, node.input[0],
                                              [int(d) for d in out_sh_f],
                                              node.output[0], pre0)
                        layout[node.output[0]] = 'RAW'
                        continue
                    print(f"  [warn] Flatten {nm}: RAW input without static shape; "
                          f"legacy emission (likely wrong)", file=sys.stderr)
                spec = LayerSpec('FlattenLayer', nm, attrs={'Axis': int(axis)},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                layout[node.output[0]] = 'FLAT'
                emit(spec); continue

            if op == 'Reshape':
                # Need shape constant
                shape_arr = initializers.get(node.input[1])
                if shape_arr is None:
                    # Dynamic shape input — fall back to passthrough placeholder.
                    # Common in YOLO-style nets where shape is computed via Shape+Concat ops.
                    spec = LayerSpec('PlaceholderLayer', nm,
                                     attrs={'OriginalOp': 'Reshape_dynamic'},
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]])
                    emit(spec)
                    print(f"  [warn] Reshape {nm}: dynamic shape -> PlaceholderLayer", file=sys.stderr)
                    continue
                tgt = [int(x) for x in shape_arr]

                in_lay = tensor_layout(node.input[0])
                consumer_ops = {n2.op_type for n2 in model.graph.node
                                if node.output[0] in n2.input}
                feeds_image = bool(consumer_ops & IMAGE_CONSUMER_OPS)
                # Flatten-like rank-2 targets ([B,F] / [-1,F] / [F,1]) stay on
                # the proven legacy path even when leaving an image flow (the
                # FlattenLayer's C-style flatten IS the ONNX semantics there,
                # and a column target is the NNV feature convention already).
                flatten_like = len(tgt) == 2 and (tgt[0] in (-1, 1) or tgt[1] in (-1, 1))

                # ---- RAW (attention/token) path: exact ONNX C-order reshape.
                # Taken when the producer is already RAW, or when an HWC/FLAT
                # producer feeds NON-image consumers with a non-flatten target
                # (the image<->token boundary, e.g. ViT patch-embed
                # [1,C,H,W] -> [1,C,T]).
                want_raw = (in_lay == 'RAW') or \
                           (bool(consumer_ops) and not feeds_image and not flatten_like)
                if want_raw:
                    tgt_static = static_reshape_target(node, tgt)
                    in_sh = value_shapes.get(node.input[0])
                    pre0 = None
                    ok = tgt_static is not None
                    if ok:
                        if in_lay == 'HWC':
                            pre0 = [1, 0, 2]      # [H,W,C] -> [W,H,C]
                        elif in_lay == 'HWCB':
                            pre0 = [2, 1, 0]      # [H,W,C] -> [C,W,H] (ONNX BHWC C-order)
                        elif in_lay == 'RAW':
                            if in_sh is not None and all(d > 0 for d in in_sh):
                                pre0 = list(range(len(in_sh) - 1, -1, -1))
                            else:
                                ok = False
                        # FLAT: column-major flatten already equals the ONNX
                        # C-order sequence; no pre-permute needed.
                    if ok:
                        emit_conorder_reshape(nm, node.input[0], tgt_static,
                                              node.output[0], pre0)
                        layout[node.output[0]] = 'RAW'
                        continue
                    print(f"  [warn] Reshape {nm}: RAW path needs static shapes; "
                          f"falling back to legacy emission", file=sys.stderr)

                # ---- legacy image/feature-flow paths (behavior unchanged).
                # Common cases:
                #   * [-1, F]  or  [B, F]  -> FlattenLayer (rank-2, the typical
                #     Conv->FC adapter; numel preserved)
                #   * [B, C, H, W]   -> ReshapeLayer with NNV-native [H, W, C]
                #     because NNV ImageInputLayer convention is HWC.
                if len(tgt) == 2:
                    # Two flavors:
                    #   * input rank >= 3 (image-derived): genuine flatten with
                    #     permute → emit FlattenLayer (rank-3 [1,1,N] result is
                    #     fine because NNV's downstream FC treats trailing-only
                    #     non-1 dim as the feature length).
                    #   * input rank 2 (already feature flow [B, F]): identity
                    #     reshape; emit ReshapeLayer with target [F, 1] so the
                    #     column-vector representation is preserved.
                    in_sh_r = value_shapes.get(node.input[0])
                    # Feature flow if the input has at most one non-singleton
                    # dim (rank-2 [B, F] or rank-3 [1, 1, F] post-onnxsim).
                    if in_sh_r is not None:
                        non_singleton = sum(1 for d in in_sh_r if d > 1)
                        is_feature = non_singleton <= 1
                    else:
                        is_feature = False
                    if is_feature:
                        # Feature size: prefer the inferred output shape; fall
                        # back to the explicit target dim if available, else
                        # numel(input).
                        out_sh = value_shapes.get(node.output[0])
                        if out_sh is not None and len(out_sh) >= 1:
                            F = int(np.prod([d for d in out_sh if d > 0]))
                        elif in_sh_r is not None:
                            F = int(np.prod([d for d in in_sh_r if d > 0]))
                        else:
                            F = max(int(t) for t in tgt if int(t) > 0)
                        spec = LayerSpec('ReshapeLayer', nm,
                                         attrs={'TargetShape': [F, 1]},
                                         inputs=[producer[node.input[0]]],
                                         outputs=[node.output[0]])
                    else:
                        spec = LayerSpec('FlattenLayer', nm,
                                         inputs=[producer[node.input[0]]],
                                         outputs=[node.output[0]])
                    layout[node.output[0]] = 'FLAT'
                elif len(tgt) == 4:
                    # IMAGE restore ([B,C,H,W] feeding Conv/Pool): drop batch,
                    # swap from CHW -> HWC for NNV convention. OnnxBCHW=1
                    # tells ReshapeLayer to apply ONNX C-order semantics
                    # (reshape to [W,H,C] then permute to [H,W,C]). Non-image
                    # consumers were already diverted to the RAW path above.
                    _, c, h, w = tgt
                    nnv_shape = [int(h), int(w), int(c)]
                    spec = LayerSpec('ReshapeLayer', nm,
                                     attrs={'TargetShape': nnv_shape, 'OnnxBCHW': 1},
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]])
                    layout[node.output[0]] = 'HWC'
                else:
                    spec = LayerSpec('ReshapeLayer', nm, attrs={'TargetShape': tgt},
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Transpose':
                # Skip the leading BHWC->BCHW transpose if the input layer
                # is already HWC (we detected & dropped it above).
                if drop_first_transpose and node.output[0] == drop_first_transpose:
                    continue
                perm = get_attr(node, 'perm')
                perm_l = list(map(int, perm)) if perm else []
                # RAW (attention/token) tensors carry the ONNX dims literally,
                # so EVERY permutation is a real data movement — never apply
                # the HWC-convention "no-op" shortcuts below (treating e.g.
                # the K^T perm [0,2,3,1] as a no-op transposed the key matrix
                # and broke Q@K^T in vit_2023).
                if tensor_layout(node.input[0]) == 'RAW':
                    spec = LayerSpec('TransposeLayer', nm, attrs={'Perm': perm_l},
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]])
                    emit(spec)
                    layout[node.output[0]] = 'RAW'
                    continue
                # Common BHWC<->BCHW transposes (perm=[0,2,3,1] or [0,3,1,2])
                # are no-ops in NNV's HWC convention — *unless* the next op is a
                # Flatten/Reshape, in which case the ONNX flatten happens on
                # BHWC data (C fastest in C-order). NNV's FlattenLayer in
                # ONNX-style assumes BCHW input (W fastest), so we need to
                # explicitly permute HWC->CWH so MATLAB column-major flatten
                # yields the same C-fastest order ONNX would.
                if perm_l == [0, 2, 3, 1]:
                    # Look for a downstream Flatten/Reshape consumer of this
                    # transpose's output before any other "real" op.
                    out_name = node.output[0]
                    next_is_flatten = False
                    for later in model.graph.node[i+1:]:
                        if out_name in later.input:
                            if later.op_type in ('Flatten', 'Reshape', 'Squeeze',
                                                 'GlobalAveragePool', 'ReduceMean'):
                                next_is_flatten = True
                            break
                    if next_is_flatten:
                        # The downstream Flatten (ONNX C-style) does its own
                        # permute([2,1,3]) + reshape internally. To make the
                        # final flat order match ONNX BHWC C-order (c-fastest,
                        # w, h), pre-permute HWC with MATLAB perm [2, 3, 1] so:
                        #   HWC --perm[2,3,1]--> [W, C, H]
                        #   --Flatten perm[2,1,3]--> [C, W, H]
                        #   --col-major flat--> c-fastest, then w, then h. ✓
                        # The loader adds 1 to convert 0-indexed ONNX perm
                        # to MATLAB 1-indexed, so emit [1, 2, 0].
                        spec = LayerSpec('TransposeLayer', nm, attrs={'Perm': [1, 2, 0]},
                                         inputs=[producer[node.input[0]]],
                                         outputs=[node.output[0]])
                        emit(spec); continue
                # rank-2 transpose [1, 0] on a feature vector (one of the dims
                # is 1) is a row-vector ↔ column-vector swap. NNV stores rank-2
                # feature data as column vector [F, 1], so this is semantically
                # a no-op. Only fire for vector-shaped tensors where the
                # singleton dim makes the transpose feature-preserving.
                if perm_l == [1, 0]:
                    in_sh_t = value_shapes.get(node.input[0])
                    if in_sh_t is not None and len(in_sh_t) == 2:
                        a, b = in_sh_t
                        if a == 1 or b == 1:
                            spec = LayerSpec('PlaceholderLayer', nm,
                                             attrs={'OriginalOp': 'Transpose_2D_no_op_vector'},
                                             inputs=[producer[node.input[0]]],
                                             outputs=[node.output[0]])
                            emit(spec); continue
                if perm_l in ([0, 2, 3, 1], [0, 3, 1, 2]):
                    tag = 'Transpose_BCHW_to_BHWC' if perm_l == [0, 2, 3, 1] else 'Transpose_BHWC_to_BCHW'
                    spec = LayerSpec('PlaceholderLayer', nm,
                                     attrs={'OriginalOp': tag},
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]])
                    # The MATLAB value stays the same [H,W,C] array; only the
                    # ONNX-side dim order flips (for strict comparison).
                    layout[node.output[0]] = 'HWCB' if perm_l == [0, 2, 3, 1] else 'HWC'
                else:
                    spec = LayerSpec('TransposeLayer', nm, attrs={'Perm': perm_l},
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Concat':
                axis = get_attr(node, 'axis', 0)
                # Use ensure_producer so initializer-as-data inputs (e.g.
                # CLS tokens, positional embeddings) get synthetic constant
                # producers instead of being dropped.
                ins = []
                for x in node.input:
                    p = ensure_producer(x)
                    if p is not None:
                        ins.append(p)
                # Stash input rank so the loader can map ONNX axis to NNV's
                # MATLAB axis correctly. Try value_shapes first; if the
                # input is an initializer (e.g. a CLS token), fall back to
                # its shape. Try all inputs in order until one resolves.
                in_rank = 0
                for inp in node.input:
                    in_sh_c = value_shapes.get(inp)
                    if in_sh_c is None and inp in initializers:
                        in_sh_c = list(initializers[inp].shape)
                    if in_sh_c:
                        in_rank = len(in_sh_c)
                        break
                spec = LayerSpec('ConcatenationLayer', nm,
                                 attrs={'Axis': int(axis), 'InRank': int(in_rank)},
                                 inputs=ins,
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'MaxPool':
                kshape = get_attr(node, 'kernel_shape')
                strides = get_attr(node, 'strides', [1]*len(kshape))
                pads = get_attr(node, 'pads', [0]*(2*len(kshape)))
                spec = LayerSpec('MaxPooling2DLayer', nm,
                                 attrs={'KernelSize': list(map(int, kshape)),
                                        'Strides': list(map(int, strides)),
                                        'Pads': list(map(int, pads))},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'AveragePool':
                kshape = get_attr(node, 'kernel_shape')
                strides = get_attr(node, 'strides', [1]*len(kshape))
                pads = get_attr(node, 'pads', [0]*(2*len(kshape)))
                spec = LayerSpec('AveragePooling2DLayer', nm,
                                 attrs={'KernelSize': list(map(int, kshape)),
                                        'Strides': list(map(int, strides)),
                                        'Pads': list(map(int, pads))},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'GlobalAveragePool':
                spec = LayerSpec('GlobalAveragePooling2DLayer', nm,
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op in ('Dropout', 'Identity'):
                # passthrough at inference time
                spec = LayerSpec('PlaceholderLayer', nm,
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Constant':
                # Already hoisted into initializers above; just skip
                continue

            if op == 'Cast':
                # Inference-time no-op for our purposes (NNV uses doubles regardless)
                spec = LayerSpec('PlaceholderLayer', nm,
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Squeeze':
                # Constant-fold: emit Reshape to the de-squeezed shape if input shape known.
                # If input is in value_shapes, compute the squeezed shape; else, treat as
                # PlaceholderLayer (works iff downstream doesn't care about shape).
                axes = get_attr(node, 'axes', None)
                in_shape_v = value_shapes.get(node.input[0])
                if in_shape_v is not None:
                    if axes is None:
                        new_shape = [d for d in in_shape_v if d != 1]
                    else:
                        axset = set(int(a) % len(in_shape_v) for a in axes)
                        new_shape = [d for i, d in enumerate(in_shape_v) if i not in axset]
                    spec = LayerSpec('ReshapeLayer', nm,
                                     attrs={'TargetShape': list(map(int, new_shape))},
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]])
                else:
                    spec = LayerSpec('PlaceholderLayer', nm,
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Unsqueeze':
                # Symmetric: insert singleton dims. If input shape known, build
                # a static Reshape; else PlaceholderLayer.
                axes = get_attr(node, 'axes', None)
                in_shape_v = value_shapes.get(node.input[0])
                if in_shape_v is not None and axes is not None:
                    new_shape = list(in_shape_v)
                    for a in sorted(int(x) for x in axes):
                        new_shape.insert(a, 1)
                    spec = LayerSpec('ReshapeLayer', nm,
                                     attrs={'TargetShape': list(map(int, new_shape))},
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]])
                else:
                    spec = LayerSpec('PlaceholderLayer', nm,
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Slice':
                # Static slice: when starts/ends/axes/steps are all initializer
                # constants and the input shape is known, emit a sparse
                # FullyConnectedLayer that selects the indexed elements. This
                # works for the common feature-input slicing case (e.g. taking
                # the first k of N features).
                in_name = node.input[0]
                in_sh = value_shapes.get(in_name)
                emitted = False
                if in_sh is not None and len(node.input) >= 3:
                    starts = initializers.get(node.input[1])
                    ends   = initializers.get(node.input[2])
                    axes_v = initializers.get(node.input[3]) if len(node.input) > 3 else None
                    steps  = initializers.get(node.input[4]) if len(node.input) > 4 else None
                    if starts is not None and ends is not None:
                        # Compute the slice for each axis
                        try:
                            starts_i = [int(x) for x in starts.flatten()]
                            ends_i   = [int(x) for x in ends.flatten()]
                            axes_i   = [int(x) for x in axes_v.flatten()] if axes_v is not None else list(range(len(starts_i)))
                            steps_i  = [int(x) for x in steps.flatten()]  if steps is not None else [1]*len(starts_i)
                            # Build the output shape and the index sets
                            out_shape = list(in_sh)
                            indices_per_axis = {}
                            n = len(in_sh)
                            for k, ax in enumerate(axes_i):
                                if ax < 0: ax += n
                                dim = in_sh[ax] if in_sh[ax] > 0 else 1
                                start = starts_i[k] if starts_i[k] >= 0 else dim + starts_i[k]
                                end   = ends_i[k]   if ends_i[k]   >= 0 else dim + ends_i[k]
                                start = max(0, min(start, dim))
                                end   = max(0, min(end, dim))
                                step  = steps_i[k]
                                idxs = list(range(start, end, step))
                                indices_per_axis[ax] = idxs
                                out_shape[ax] = len(idxs)
                            # If the slice is on a flat feature input (rank=2 or rank=1)
                            # with axis=last, build an FC selector.
                            sliced_axes = list(indices_per_axis.keys())
                            # General case: build a [flat_out, flat_in] selector
                            # by enumerating output coordinates and mapping back
                            # to source flat indices. Treat zero-dim as 1.
                            eff_sh = [d if d > 0 else 1 for d in in_sh]
                            n = len(eff_sh)
                            flat_in = 1
                            for d in eff_sh:
                                flat_in *= d
                            # out_dims: replace each sliced axis with len(idxs)
                            out_dims = list(eff_sh)
                            for ax2, idxs2 in indices_per_axis.items():
                                out_dims[ax2] = len(idxs2)
                            flat_out = 1
                            for d in out_dims:
                                flat_out *= d
                            # Size guard: a [flat_out, flat_in] dense selector
                            # matrix would be O(flat_out * flat_in) bytes. For
                            # large feature maps (640x640 spatial inputs in
                            # YOLO etc.) this can blow past machine memory.
                            # Cap at ~200k * 200k => 4e10 bytes, but actually
                            # use 50k * 50k = 2.5e9 = 10 GB float32. Be more
                            # aggressive: refuse if flat_in or flat_out > 200k.
                            # Memory cap: float32 dense [flat_out, flat_in].
                            # Cap at 256M elements (~1 GB float32) per selector.
                            # Below that, even if savemat compresses, it's fine.
                            if (flat_out > 200000 or flat_in > 200000
                                or flat_out * flat_in > 256_000_000):
                                raise RuntimeError(
                                    f"slice selector too large "
                                    f"(flat_out={flat_out}, flat_in={flat_in}); "
                                    f"falling back to placeholder")
                            # RAW producers flatten column-major in MATLAB's
                            # FullyConnectedLayer; FLAT/HWC keep legacy C-order.
                            raw_in_sl = tensor_layout(in_name) == 'RAW'
                            if raw_in_sl:
                                strides = [1]*n
                                for k in range(1, n):
                                    strides[k] = strides[k-1] * eff_sh[k-1]
                                out_strides_sl = [1]*n
                                for k in range(1, n):
                                    out_strides_sl[k] = out_strides_sl[k-1] * out_dims[k-1]
                            else:
                                strides = [1]*n
                                for k in range(n-2, -1, -1):
                                    strides[k] = strides[k+1] * eff_sh[k+1]
                                out_strides_sl = None
                            from itertools import product
                            ranges = [range(d) for d in out_dims]
                            W_sel = np.zeros((flat_out, flat_in), dtype=np.float32)
                            for r, coord in enumerate(product(*ranges)):
                                src = list(coord)
                                for ax2, idxs2 in indices_per_axis.items():
                                    src[ax2] = idxs2[coord[ax2]]
                                src_flat = sum(s*st for s, st in zip(src, strides))
                                if out_strides_sl is not None:
                                    r = sum(c*st for c, st in zip(coord, out_strides_sl))
                                W_sel[r, src_flat] = 1.0
                            b_sel = np.zeros(flat_out, dtype=np.float32)
                            if raw_in_sl:
                                fc_nm = safe_name(nm + '_sel')
                                wkey = add_weight(fc_nm, 'W', W_sel)
                                bkey = add_weight(fc_nm, 'b', b_sel)
                                emit(LayerSpec('FullyConnectedLayer', fc_nm,
                                               attrs={'OutputSize': int(flat_out)},
                                               inputs=[producer[in_name]],
                                               outputs=[fc_nm + '_out'],
                                               weight_keys=[wkey, bkey]))
                                spec = LayerSpec('ReshapeLayer', nm,
                                                 attrs={'TargetShape': list(map(int, out_dims))},
                                                 inputs=[fc_nm],
                                                 outputs=[node.output[0]])
                                emit(spec)
                                layout[node.output[0]] = 'RAW'
                            else:
                                wkey = add_weight(nm, 'W', W_sel)
                                bkey = add_weight(nm, 'b', b_sel)
                                spec = LayerSpec('FullyConnectedLayer', nm,
                                                 attrs={'OutputSize': int(flat_out)},
                                                 inputs=[producer[in_name]],
                                                 outputs=[node.output[0]],
                                                 weight_keys=[wkey, bkey])
                                emit(spec)
                            emitted = True
                        except Exception as e:
                            print(f"  [warn] Slice {nm}: selector build failed ({e}), falling back to placeholder", file=sys.stderr)

                if not emitted:
                    spec = LayerSpec('PlaceholderLayer', nm,
                                     inputs=[producer[in_name]],
                                     outputs=[node.output[0]])
                    emit(spec)
                    print(f"  [warn] Slice op {nm} mapped to PlaceholderLayer (no-op); reach may be inexact",
                          file=sys.stderr)
                continue

            if op == 'Pad':
                # Constant-pad: skip; downstream Conv usually has compensating padding
                spec = LayerSpec('PlaceholderLayer', nm,
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec)
                print(f"  [warn] Pad op {nm} mapped to PlaceholderLayer; check downstream Conv has compensating padding",
                      file=sys.stderr)
                continue

            if op == 'Split':
                # Real multi-output split: emit one selector FC per output
                # that picks the right channel slice. Each selector layer
                # has its own output name so downstream consumers connect
                # to the correct slice.
                axis = get_attr(node, 'axis', 0)
                in_name = node.input[0]
                in_sh = value_shapes.get(in_name)
                if in_sh is None and in_name in initializers:
                    in_sh = list(initializers[in_name].shape)

                # Try real selector emission if shape is known.
                emitted_real = False
                if in_sh is not None:
                    eff_sh = [d if d > 0 else 1 for d in in_sh]
                    n_dims = len(eff_sh)
                    ax = axis if axis >= 0 else axis + n_dims
                    if 0 <= ax < n_dims:
                        # Determine split sizes: explicit `split` input or attr,
                        # else equal split.
                        split_sizes = None
                        if len(node.input) > 1 and node.input[1] in initializers:
                            split_sizes = [int(s) for s in initializers[node.input[1]].flatten()]
                        else:
                            split_attr = None
                            for a in node.attribute:
                                if a.name == 'split':
                                    split_attr = list(a.ints); break
                            if split_attr:
                                split_sizes = split_attr
                            else:
                                n_out = len(node.output)
                                if eff_sh[ax] % n_out == 0:
                                    split_sizes = [eff_sh[ax] // n_out] * n_out

                        if split_sizes and sum(split_sizes) == eff_sh[ax]:
                            try:
                                # ATOMIC size pre-check: validate EVERY part
                                # against the memory cap BEFORE emitting any
                                # layer. Emitting part0 and then bailing on
                                # part1 used to leave orphaned, mis-wired
                                # selector FCs alongside the fallback
                                # placeholder (collins /model.24/Split_2),
                                # which crash evaluate on both the MATLAB and
                                # simulator sides.
                                flat_in_all = int(np.prod(eff_sh))
                                for ssz in split_sizes:
                                    od = list(eff_sh); od[ax] = ssz
                                    fo = int(np.prod(od))
                                    if (fo > 200000 or flat_in_all > 200000
                                            or fo * flat_in_all > 64_000_000):
                                        raise RuntimeError("split selector too large")
                                # Flat-order of the selector must match how the
                                # MATLAB FullyConnectedLayer flattens its input:
                                # column-major (F-order over the ONNX dims) for
                                # RAW producers; legacy C-order otherwise.
                                raw_in = tensor_layout(in_name) == 'RAW'
                                offset = 0
                                for k, (out_name, ssz) in enumerate(zip(node.output, split_sizes)):
                                    sub_nm = safe_name(f"{nm}_part{k}")
                                    # Build selector matrix [flat_out, flat_in] like Slice
                                    # for a single-axis range select.
                                    out_dims = list(eff_sh)
                                    out_dims[ax] = ssz
                                    flat_in = int(np.prod(eff_sh))
                                    flat_out = int(np.prod(out_dims))
                                    # Memory cap: float32 dense [flat_out, flat_in]
                                    # is 4 * flat_out * flat_in bytes. Cap at 256 MB
                                    # (~64M float32) per selector.
                                    if (flat_out > 200000 or flat_in > 200000
                                        or flat_out * flat_in > 64_000_000):
                                        raise RuntimeError("split selector too large")
                                    if raw_in:
                                        strides = [1] * n_dims
                                        for j in range(1, n_dims):
                                            strides[j] = strides[j-1] * eff_sh[j-1]
                                        out_strides = [1] * n_dims
                                        for j in range(1, n_dims):
                                            out_strides[j] = out_strides[j-1] * out_dims[j-1]
                                    else:
                                        strides = [1] * n_dims
                                        for j in range(n_dims - 2, -1, -1):
                                            strides[j] = strides[j+1] * eff_sh[j+1]
                                        out_strides = None
                                    W_sel = np.zeros((flat_out, flat_in), dtype=np.float32)
                                    from itertools import product
                                    ranges = [range(d) for d in out_dims]
                                    for r, coord in enumerate(product(*ranges)):
                                        src = list(coord)
                                        src[ax] = coord[ax] + offset
                                        src_flat = sum(s*st for s, st in zip(src, strides))
                                        if out_strides is not None:
                                            r = sum(c*st for c, st in zip(coord, out_strides))
                                        W_sel[r, src_flat] = 1.0
                                    b_sel = np.zeros(flat_out, dtype=np.float32)
                                    wkey = add_weight(sub_nm, 'W', W_sel)
                                    bkey = add_weight(sub_nm, 'b', b_sel)
                                    sel_internal_out = f"{out_name}_flat"
                                    spec = LayerSpec('FullyConnectedLayer', sub_nm,
                                                     attrs={'OutputSize': int(flat_out)},
                                                     inputs=[producer[in_name]],
                                                     outputs=[sel_internal_out],
                                                     weight_keys=[wkey, bkey])
                                    emit(spec)
                                    # Reshape flat output back to the original
                                    # logical shape so downstream ops with
                                    # rank-3 broadcasting (Add, etc.) see [1, T, C].
                                    # (For RAW producers the F-order selector
                                    # output + column-major reshape compose
                                    # exactly back to the RAW dims.)
                                    rs_nm = safe_name(f"{sub_nm}_rs")
                                    rs_spec = LayerSpec('ReshapeLayer', rs_nm,
                                                        attrs={'TargetShape': list(map(int, out_dims))},
                                                        inputs=[sub_nm],
                                                        outputs=[out_name])
                                    emit(rs_spec)
                                    if raw_in:
                                        layout[out_name] = 'RAW'
                                    offset += ssz
                                emitted_real = True
                            except Exception as e:
                                print(f"  [warn] Split {nm}: selector emit failed ({e}); fallback to placeholder", file=sys.stderr)

                if not emitted_real:
                    spec = LayerSpec('PlaceholderLayer', nm,
                                     attrs={'Axis': int(axis)},
                                     inputs=[producer[in_name]],
                                     outputs=list(node.output))
                    emit(spec)
                    print(f"  [warn] Split op {nm} mapped to PlaceholderLayer; outputs share one tensor", file=sys.stderr)
                continue

            if op == 'Sign':
                # Larq STE-Sign pattern (Sign + Add(+const) + Sign): emit a
                # single SignLayer in polar mode that maps 0 -> +1.
                if node_key in ste_polar:
                    real_in_name = ste_polar[node_key]
                    real_in_producer = ensure_producer(real_in_name)
                    if real_in_producer is None:
                        real_in_producer = producer.get(real_in_name)
                    spec = LayerSpec('SignLayer', nm,
                                     attrs={'Mode': 'polar_zero_to_pos_one'},
                                     inputs=[real_in_producer],
                                     outputs=[node.output[0]])
                    emit(spec)
                    continue
                # Map to NNV's SignLayer. Sign.m returns plain sign(x) (ONNX
                # semantics: -1, 0, +1) when mode is anything not in the
                # known polar/nonnegative set. Use 'onnx' as a clear sentinel.
                spec = LayerSpec('SignLayer', nm,
                                 attrs={'Mode': 'onnx'},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec)
                continue

            if op == 'Resize':
                # Image resize: map to NNV Resize2DLayer if scale/size factors
                # are constants. For now, emit Resize2DLayer with mode='nearest'
                # by default; downstream evaluate may compensate if attrs match.
                mode = get_attr(node, 'mode', 'nearest')
                spec = LayerSpec('Resize2DLayer', nm, attrs={'Mode': mode},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec)
                continue

            if op == 'ReduceMean':
                _raw_red = tensor_layout(node.input[0]) == 'RAW'
                if _try_emit_reduce(node, nm, value_shapes, initializers, producer,
                                    add_weight, emit, op_name='mean',
                                    raw_order=_raw_red):
                    layout[node.output[0]] = 'RAW' if _raw_red else 'FLAT'
                    continue
                axes = get_attr(node, 'axes', None)
                spec = LayerSpec('PlaceholderLayer', nm, attrs={'OriginalOp':'ReduceMean',
                    'Axes': list(map(int, axes)) if axes else []},
                    inputs=[producer[node.input[0]]], outputs=[node.output[0]])
                emit(spec)
                print(f"  [warn] ReduceMean op {nm} fallback to PlaceholderLayer", file=sys.stderr)
                continue

            if op == 'Pow':
                # x^k for constant k -> PlaceholderLayer (xval will fail; mark)
                spec = LayerSpec('PlaceholderLayer', nm, attrs={'OriginalOp':'Pow'},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec)
                print(f"  [warn] Pow op {nm} mapped to PlaceholderLayer", file=sys.stderr)
                continue

            if op == 'Neg':
                # y = -x  =  -1 * x + 0
                # Encode as ElementwiseAffineLayer with scale=-1.
                shape_in = value_shapes.get(node.input[0])
                fan = max(1, int(np.prod([d for d in (shape_in or [1]) if d > 0])))
                scale_arr = -np.ones(1, dtype=np.float32)
                bias_arr  = np.zeros(1, dtype=np.float32)
                sk = add_weight(nm, 'scale', scale_arr)
                bk = add_weight(nm, 'bias', bias_arr)
                spec = LayerSpec('ElementwiseAffineLayer', nm,
                    attrs={'DoScale': True, 'DoOffset': False},
                    inputs=[producer[node.input[0]]],
                    outputs=[node.output[0]],
                    weight_keys=[sk, bk])
                emit(spec); continue

            if op == 'Gather':
                # Static-index gather on a tracked-shape data tensor: build
                # a sparse FullyConnectedLayer (selector matrix) that pulls
                # the gathered elements out of the flattened input. This
                # turns Gather into a sound, exact linear op for verification.
                emitted = False
                data_name = node.input[0]
                idx_name  = node.input[1] if len(node.input) > 1 else None
                in_sh = value_shapes.get(data_name)
                idx_arr = initializers.get(idx_name) if idx_name else None
                axis = get_attr(node, 'axis', 0)
                if in_sh is not None and idx_arr is not None:
                    try:
                        # Compute total flattened input length, treating
                        # zero/negative dims as 1 (single-batch).
                        flat_in = 1
                        eff_sh = []
                        for d in in_sh:
                            d2 = d if d > 0 else 1
                            flat_in *= d2
                            eff_sh.append(d2)
                        n = len(eff_sh)
                        ax = axis if axis >= 0 else axis + n
                        if 0 <= ax < n:
                            # ONNX Gather: output[I_0,...,I_{ax-1}, q_0,...,q_{r-1}, I_{ax+1},...]
                            # = data[I_0,..., indices[q_0,...,q_{r-1}], I_{ax+1},...]
                            idxs = np.asarray(idx_arr).flatten().astype(int).tolist()
                            dim_ax = eff_sh[ax]
                            # Negative indices wrap
                            idxs = [(i if i >= 0 else i + dim_ax) for i in idxs]
                            # Build flat row->col map. For each output position
                            # (lex over I_0..ax-1, q_0..q_{r-1}, I_{ax+1}..),
                            # compute the source flat index.
                            # Strides over eff_sh: C-order for legacy FLAT/HWC
                            # flows; F-order (MATLAB column-major) for RAW.
                            raw_in_g = tensor_layout(data_name) == 'RAW'
                            if raw_in_g:
                                strides = [1]*n
                                for k in range(1, n):
                                    strides[k] = strides[k-1] * eff_sh[k-1]
                            else:
                                strides = [1]*n
                                for k in range(n-2, -1, -1):
                                    strides[k] = strides[k+1] * eff_sh[k+1]
                            # Iterate output coordinates
                            outer_dims = eff_sh[:ax]
                            inner_dims = eff_sh[ax+1:]
                            out_dims_full = list(outer_dims) + [len(idxs)] + list(inner_dims)
                            flat_out = int(np.prod(out_dims_full)) if out_dims_full else 1
                            if raw_in_g:
                                out_strides_g = [1]*len(out_dims_full)
                                for k in range(1, len(out_dims_full)):
                                    out_strides_g[k] = out_strides_g[k-1] * out_dims_full[k-1]
                            else:
                                out_strides_g = None
                            W_sel = np.zeros((flat_out, flat_in), dtype=np.float32)
                            from itertools import product
                            ranges = [range(d) for d in outer_dims] + [range(len(idxs))] + [range(d) for d in inner_dims]
                            for r, coord in enumerate(product(*ranges)):
                                # Build src index using idxs[q] for the gather axis
                                src = list(coord)
                                src[ax] = idxs[coord[ax]]
                                src_flat = sum(s*st for s, st in zip(src, strides))
                                if out_strides_g is not None:
                                    r = sum(c*st for c, st in zip(coord, out_strides_g))
                                W_sel[r, src_flat] = 1.0
                            b_sel = np.zeros(flat_out, dtype=np.float32)
                            if raw_in_g:
                                fc_nm = safe_name(nm + '_sel')
                                wkey = add_weight(fc_nm, 'W', W_sel)
                                bkey = add_weight(fc_nm, 'b', b_sel)
                                emit(LayerSpec('FullyConnectedLayer', fc_nm,
                                               attrs={'OutputSize': int(flat_out)},
                                               inputs=[producer[data_name]],
                                               outputs=[fc_nm + '_out'],
                                               weight_keys=[wkey, bkey]))
                                spec = LayerSpec('ReshapeLayer', nm,
                                                 attrs={'TargetShape': list(map(int, out_dims_full))},
                                                 inputs=[fc_nm],
                                                 outputs=[node.output[0]])
                                emit(spec)
                                layout[node.output[0]] = 'RAW'
                            else:
                                wkey = add_weight(nm, 'W', W_sel)
                                bkey = add_weight(nm, 'b', b_sel)
                                spec = LayerSpec('FullyConnectedLayer', nm,
                                                 attrs={'OutputSize': int(flat_out)},
                                                 inputs=[producer[data_name]],
                                                 outputs=[node.output[0]],
                                                 weight_keys=[wkey, bkey])
                                emit(spec)
                            emitted = True
                    except Exception as e:
                        print(f"  [warn] Gather {nm}: selector build failed ({e}), falling back to placeholder", file=sys.stderr)
                if not emitted:
                    spec = LayerSpec('PlaceholderLayer', nm, attrs={'OriginalOp':'Gather'},
                                     inputs=[producer[x] for x in node.input if x in producer],
                                     outputs=[node.output[0]])
                    emit(spec)
                    print(f"  [warn] Gather op {nm} mapped to PlaceholderLayer", file=sys.stderr)
                continue

            if op == 'MatMul' and not (node.input[1] in initializers):
                # Dynamic MatMul (e.g. attention Q*K^T). Map to PlaceholderLayer.
                spec = LayerSpec('PlaceholderLayer', nm, attrs={'OriginalOp':'MatMul_dynamic'},
                                 inputs=[producer[x] for x in node.input if x in producer],
                                 outputs=[node.output[0]])
                emit(spec)
                print(f"  [warn] dynamic MatMul {nm} mapped to PlaceholderLayer", file=sys.stderr)
                continue

            if op == 'Mul' and not any(x in initializers for x in node.input):
                # Dynamic-by-dynamic mul (e.g. Lyapunov bilinear, attention scaling).
                # Emit ElementwiseProductLayer so evaluate computes the
                # Hadamard product. reach uses an interval over-approximation.
                ins = [producer[x] for x in node.input if x in producer]
                spec = LayerSpec('ElementwiseProductLayer', nm,
                                 attrs={'NumInputs': len(ins)},
                                 inputs=ins,
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op == 'Div':
                # If divisor is a constant, fold to ElementwiseAffineLayer
                # with scale = 1/divisor. Otherwise emit dynamic division.
                if len(node.input) >= 2 and node.input[1] in initializers:
                    div_arr = initializers[node.input[1]].astype(np.float32)
                    if np.any(div_arr == 0):
                        # zero divisor — leave as placeholder
                        spec = LayerSpec('PlaceholderLayer', nm,
                                         attrs={'OriginalOp':'Div_zero'},
                                         inputs=[producer[x] for x in node.input if x in producer],
                                         outputs=[node.output[0]])
                        emit(spec)
                        print(f"  [warn] Div {nm} has zero divisor, placeholder", file=sys.stderr)
                        continue
                    scale = (1.0 / div_arr).flatten()
                    bias  = np.zeros_like(scale)
                    skey = add_weight(nm, 'scale', scale)
                    bkey = add_weight(nm, 'bias',  bias)
                    spec = LayerSpec('ElementwiseAffineLayer', nm,
                                     attrs={'DoScale': True, 'DoOffset': False},
                                     inputs=[producer[node.input[0]]],
                                     outputs=[node.output[0]],
                                     weight_keys=[skey, bkey])
                    emit(spec); continue
                # dynamic / dynamic
                ins = [producer[x] for x in node.input if x in producer]
                spec = LayerSpec('ElementwiseDivisionLayer', nm,
                                 attrs={'NumInputs': len(ins)},
                                 inputs=ins,
                                 outputs=[node.output[0]])
                emit(spec); continue

            if op in ('Min', 'Max', 'Clip'):
                # Saturating activations; map to ReluLayer-like via PlaceholderLayer
                spec = LayerSpec('PlaceholderLayer', nm, attrs={'OriginalOp':op},
                                 inputs=[producer[x] for x in node.input if x in producer],
                                 outputs=[node.output[0]])
                emit(spec)
                continue

            if op == 'Floor':
                spec = LayerSpec('PlaceholderLayer', nm, attrs={'OriginalOp':'Floor'},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec)
                print(f"  [warn] Floor op {nm} mapped to PlaceholderLayer", file=sys.stderr)
                continue

            if op in ('Cos','Sin','Tan','Exp','Log','Sqrt','Abs'):
                spec = LayerSpec('PlaceholderLayer', nm, attrs={'OriginalOp':op},
                                 inputs=[producer[node.input[0]]],
                                 outputs=[node.output[0]])
                emit(spec)
                print(f"  [warn] {op} op {nm} mapped to PlaceholderLayer", file=sys.stderr)
                continue

            if op == 'ReduceSum':
                _raw_red = tensor_layout(node.input[0]) == 'RAW'
                if _try_emit_reduce(node, nm, value_shapes, initializers, producer,
                                    add_weight, emit, op_name='sum',
                                    raw_order=_raw_red):
                    layout[node.output[0]] = 'RAW' if _raw_red else 'FLAT'
                    continue
                axes = get_attr(node, 'axes', None)
                spec = LayerSpec('PlaceholderLayer', nm, attrs={'OriginalOp':'ReduceSum',
                    'Axes': list(map(int, axes)) if axes else []},
                    inputs=[producer[node.input[0]]], outputs=[node.output[0]])
                emit(spec)
                print(f"  [warn] ReduceSum op {nm} fallback to PlaceholderLayer", file=sys.stderr)
                continue

            # Best-effort fallbacks for ops we don't fully model. Each is
            # mapped to a PlaceholderLayer that passes its first input
            # through. xval will diverge for uses where the op actually
            # changes values, but the network at least imports.
            if op in ('Cast', 'Identity', 'Unsqueeze', 'Squeeze', 'Expand',
                      'Shape', 'Equal', 'Where', 'ScatterND', 'ArgMax',
                      'Min', 'Max', 'Clip', 'Pow', 'Resize'):
                # For multi-input ops (Where: 3 inputs, ScatterND: 3, etc.)
                # we still only forward the first known producer.
                in_first = None
                for x in node.input:
                    if x in producer:
                        in_first = x; break
                ins = [producer[in_first]] if in_first else []
                spec = LayerSpec('PlaceholderLayer', nm,
                                 attrs={'OriginalOp': op},
                                 inputs=ins,
                                 outputs=[node.output[0]])
                emit(spec)
                print(f"  [warn] {op} op {nm} -> PlaceholderLayer (passthrough)", file=sys.stderr)
                continue

            raise RuntimeError(f"unsupported op {op} (node {nm})")

        except KeyError as e:
            # Unknown producer — could be an unsupported topology op upstream
            # that didn't register itself. Emit a PlaceholderLayer with the
            # producers we DO know about, so the rest of the graph still
            # loads. xval will likely fail downstream but the net imports.
            ins_known = [producer[x] for x in node.input if x in producer]
            spec = LayerSpec('PlaceholderLayer', nm,
                attrs={'OriginalOp': op, 'MissingProducerFor': str(e)},
                inputs=ins_known,
                outputs=list(node.output))
            emit(spec)
            print(f"  [warn] node {nm} ({op}): missing producer for {e} — mapped to PlaceholderLayer", file=sys.stderr)

    # Output spec
    out_name = model.graph.output[0].name

    # Stamp each layer with the MATLAB-side layout of its (first) output —
    # metadata only (the MATLAB loader ignores unknown attrs); used by
    # manifest_sim.py for orientation-STRICT comparison against onnxruntime.
    for L in layers:
        if L.outputs:
            L.attrs = dict(L.attrs) if L.attrs else {}
            L.attrs.setdefault('MlLayout', layout.get(L.outputs[0], 'FLAT'))

    return layers, weights, inp_name, inp_shape, out_name


# ---------- mat writer ----------

def manifest_to_mat(layers, weights, inp_name, inp_shape, out_name, mat_path, opset, src_path):
    # Convert layers to a struct array via a list of dicts. Encode list-of-strings
    # fields as numpy object arrays of fixed-length strings so scipy.io.savemat
    # produces a cell-of-char on the MATLAB side rather than mangled char tensors.
    def _str_array(lst):
        if not lst: return np.array([], dtype=object)
        arr = np.empty(len(lst), dtype=object)
        for i, s in enumerate(lst):
            arr[i] = str(s)
        return arr

    layer_dicts = []
    for L in layers:
        d = {
            'type': L.type,
            'name': L.name,
            'attrs': L.attrs if L.attrs else {'_empty_': 0},
            'inputs':  _str_array(L.inputs),
            'outputs': _str_array(L.outputs),
            'weight_keys': _str_array(L.weight_keys),
        }
        layer_dicts.append(d)

    with open(src_path, 'rb') as f:
        checksum = hashlib.md5(f.read()).hexdigest()

    manifest = {
        'layers': np.array(layer_dicts, dtype=object),
        'weights': weights,
        'input_name': inp_name,
        'input_shape': np.array(inp_shape, dtype=np.int64),
        'output_name': out_name,
        'opset': int(opset),
        'src_path': src_path,
        'checksum': checksum,
    }

    # scipy savemat (v5 format) cannot write tensors > 2^31 bytes. For large
    # models (YOLO at 640x640, etc.) fall back to HDF5 (v7.3) via h5py.
    total_w = sum(int(np.asarray(v).nbytes) for v in weights.values()) if weights else 0
    use_hdf5 = total_w > 2_000_000_000
    if use_hdf5:
        try:
            import hdf5storage  # type: ignore
            hdf5storage.write(manifest, '.', mat_path, matlab_compatible=True,
                              store_python_metadata=False, truncate_existing=True)
        except ImportError:
            print(f"  [warn] manifest is large ({total_w/1e9:.1f} GB) and hdf5storage "
                  "is not installed; install via `pip install hdf5storage` for "
                  "automatic v7.3 .mat support.", file=sys.stderr)
            savemat(mat_path, manifest, do_compression=True)
    else:
        savemat(mat_path, manifest, do_compression=True)


# ---------- main ----------

def main():
    p = argparse.ArgumentParser()
    p.add_argument('onnx', help='input ONNX path')
    p.add_argument('mat', nargs='?', help='output .mat path (default: alongside ONNX)')
    p.add_argument('--vnnlib', help='optional VNN-LIB path (for input shape inference)')
    p.add_argument('--no-simplify', action='store_true')
    p.add_argument('--target-opset', type=int, default=17)
    args = p.parse_args()

    if args.mat is None:
        args.mat = os.path.splitext(args.onnx)[0] + '.nnv.mat'

    print(f"Loading: {args.onnx}", file=sys.stderr)
    model = onnx.load(args.onnx)
    opset = model.opset_import[0].version if model.opset_import else 9
    print(f"  opset={opset}, nodes={len(model.graph.node)}, initializers={len(model.graph.initializer)}", file=sys.stderr)

    if not args.no_simplify:
        model = preprocess_onnx(model, target_opset=args.target_opset)
        print(f"  after preprocess: nodes={len(model.graph.node)}", file=sys.stderr)

    layers, weights, inp_name, inp_shape, out_name = walk(model)
    print(f"  emitted {len(layers)} layers, {len(weights)} weight tensors", file=sys.stderr)
    op_counts = {}
    for L in layers:
        op_counts[L.type] = op_counts.get(L.type, 0) + 1
    print(f"  layer kinds: {op_counts}", file=sys.stderr)

    manifest_to_mat(layers, weights, inp_name, inp_shape, out_name,
                    args.mat, opset, args.onnx)
    print(f"Wrote: {args.mat}", file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
