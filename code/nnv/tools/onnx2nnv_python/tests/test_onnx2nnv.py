"""Regression tests for the V04/V05 onnx2nnv.py importer.

Run:  python -m pytest tools/onnx2nnv_python/tests/

(or just `python tools/onnx2nnv_python/tests/test_onnx2nnv.py` to run as a script.)
"""

import os, sys, tempfile, hashlib
import numpy as np
import onnx
from onnx import helper, TensorProto

# Allow running from repo root
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import onnx2nnv as O


def _make_onnx(nodes, inputs, outputs, initializers, opset=14, name='m'):
    """Build a minimal ONNX ModelProto for testing."""
    g = helper.make_graph(nodes, name, inputs, outputs, initializer=initializers)
    m = helper.make_model(g, opset_imports=[helper.make_opsetid("", opset)])
    m.ir_version = 7
    return m


def _const(name, arr):
    return helper.make_tensor(name, TensorProto.FLOAT, list(arr.shape), arr.flatten().tolist())


def test_safe_name_truncates_and_hashes():
    short = O.safe_name("a_short_name")
    assert short == "a_short_name"
    long_n = O.safe_name("/some/very/long/onnx/style/node/name/that/is_way_too_long_for_matlab")
    assert len(long_n) <= 28
    # deterministic
    assert long_n == O.safe_name("/some/very/long/onnx/style/node/name/that/is_way_too_long_for_matlab")
    print(f"  safe_name short={short!r}, long={long_n!r}")


def test_walk_emits_fc_for_gemm():
    W = np.random.randn(3, 5).astype(np.float32)
    b = np.random.randn(3).astype(np.float32)
    nodes = [helper.make_node('Gemm', ['x', 'W', 'b'], ['y'], transB=1)]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 5])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 3])
    inits = [_const('W', W), _const('b', b)]
    m = _make_onnx(nodes, [inp], [out], inits)

    layers, weights, in_name, in_shape, out_name = O.walk(m)
    types = [L.type for L in layers]
    assert 'FullyConnectedLayer' in types, f'expected FC, got {types}'


def test_walk_emits_fc_for_matmul_with_const_b():
    W = np.random.randn(5, 3).astype(np.float32)   # [K, M] - to become FC W [M, K]
    nodes = [helper.make_node('MatMul', ['x', 'W'], ['y'])]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 5])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 3])
    m = _make_onnx(nodes, [inp], [out], [_const('W', W)])
    layers, _, _, _, _ = O.walk(m)
    types = [L.type for L in layers]
    assert types.count('FullyConnectedLayer') == 1


def test_walk_handles_dynamic_matmul_as_placeholder():
    nodes = [helper.make_node('MatMul', ['x', 'y_in'], ['z'])]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 5])
    inp2 = helper.make_tensor_value_info('y_in', TensorProto.FLOAT, [5, 3])
    out = helper.make_tensor_value_info('z', TensorProto.FLOAT, [1, 3])
    m = _make_onnx(nodes, [inp, inp2], [out], [])
    # find_data_input picks the first non-initializer input
    layers, _, _, _, _ = O.walk(m)
    types = [L.type for L in layers]
    assert 'PlaceholderLayer' in types


def test_walk_handles_neg():
    nodes = [helper.make_node('Neg', ['x'], ['y'])]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 4])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 4])
    m = _make_onnx(nodes, [inp], [out], [])
    layers, _, _, _, _ = O.walk(m)
    types = [L.type for L in layers]
    assert 'ElementwiseAffineLayer' in types


def test_walk_handles_mul_with_init_as_affine():
    s = np.array([2.0], dtype=np.float32)
    nodes = [helper.make_node('Mul', ['x', 's'], ['y'])]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 4])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 4])
    m = _make_onnx(nodes, [inp], [out], [_const('s', s)])
    layers, _, _, _, _ = O.walk(m)
    types = [L.type for L in layers]
    assert 'ElementwiseAffineLayer' in types


def test_walk_handles_sub_with_init():
    bias = np.array([0.1, 0.2, 0.3, 0.4], dtype=np.float32)
    nodes = [helper.make_node('Sub', ['x', 'b'], ['y'])]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 4])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 4])
    m = _make_onnx(nodes, [inp], [out], [_const('b', bias)])
    layers, _, _, _, _ = O.walk(m)
    types = [L.type for L in layers]
    assert 'ElementwiseAffineLayer' in types


def test_walk_handles_dynamic_sub_with_addition_layer():
    nodes = [helper.make_node('Sub', ['a', 'b'], ['y'])]
    inp = helper.make_tensor_value_info('a', TensorProto.FLOAT, [1, 4])
    inp2 = helper.make_tensor_value_info('b', TensorProto.FLOAT, [1, 4])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 4])
    m = _make_onnx(nodes, [inp, inp2], [out], [])
    layers, _, _, _, _ = O.walk(m)
    types = [L.type for L in layers]
    # Dynamic Sub becomes neg-affine + Addition
    assert 'AdditionLayer' in types


def test_reshape_rank2_becomes_flatten():
    target = np.array([1, 12], dtype=np.int64)
    nodes = [helper.make_node('Reshape', ['x', 'shape'], ['y'])]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 3, 2, 2])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 12])
    sh = helper.make_tensor('shape', TensorProto.INT64, [2], target.tolist())
    m = _make_onnx(nodes, [inp], [out], [sh])
    layers, _, _, _, _ = O.walk(m)
    types = [L.type for L in layers]
    assert 'FlattenLayer' in types, f'expected Flatten, got {types}'


def test_reshape_rank4_becomes_HWC():
    target = np.array([1, 4, 5, 5], dtype=np.int64)
    nodes = [helper.make_node('Reshape', ['x', 'shape'], ['y'])]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 100])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 4, 5, 5])
    sh = helper.make_tensor('shape', TensorProto.INT64, [4], target.tolist())
    m = _make_onnx(nodes, [inp], [out], [sh])
    layers, _, _, _, _ = O.walk(m)
    reshape_specs = [L for L in layers if L.type == 'ReshapeLayer']
    assert len(reshape_specs) == 1
    # NNV target shape should be [H, W, C]
    assert reshape_specs[0].attrs['TargetShape'] == [5, 5, 4]


def test_constant_node_hoisted_into_initializers():
    val = np.array([0.5, 1.5], dtype=np.float32)
    cnode = helper.make_node('Constant', [], ['c'], value=helper.make_tensor(
        'c', TensorProto.FLOAT, [2], val.tolist()))
    mul = helper.make_node('Mul', ['x', 'c'], ['y'])
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 2])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 2])
    m = _make_onnx([cnode, mul], [inp], [out], [])
    layers, _, _, _, _ = O.walk(m)
    types = [L.type for L in layers]
    # Constant + Mul-with-constant = ElementwiseAffine (Constant hoisted, then Mul folds)
    assert 'ElementwiseAffineLayer' in types


def test_find_data_input_skips_initializers():
    # Some PyTorch exports list every initializer as a graph input too
    W = np.random.randn(3, 5).astype(np.float32)
    inputs = [
        helper.make_tensor_value_info('W', TensorProto.FLOAT, [3, 5]),
        helper.make_tensor_value_info('actual_input', TensorProto.FLOAT, [1, 5]),
    ]
    nodes = [helper.make_node('Gemm', ['actual_input', 'W'], ['y'], transB=1)]
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 3])
    m = _make_onnx(nodes, inputs, [out], [_const('W', W)])
    di = O.find_data_input(m)
    assert di.name == 'actual_input'


def test_matmul_rank1_rhs_becomes_fc_with_row_W():
    # ml4acopf MatMul with rank-1 B: A @ B where B is [K] (column vector).
    # Importer must emit FC W of shape [1, K], not [K, 1] from B.T.
    B = np.random.randn(8).astype(np.float32)
    nodes = [helper.make_node('MatMul', ['x', 'B'], ['y'])]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 8])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1])
    m = _make_onnx(nodes, [inp], [out], [_const('B', B)])
    layers, weights, _, _, _ = O.walk(m)
    fc_specs = [L for L in layers if L.type == 'FullyConnectedLayer']
    assert len(fc_specs) == 1, 'expected exactly one FC layer'
    spec = fc_specs[0]
    Wkey, bkey = spec.weight_keys
    Wm = weights[Wkey]
    bm = weights[bkey]
    assert Wm.shape == (1, 8), f'expected W shape (1,8), got {Wm.shape}'
    assert bm.shape == (1,) or bm.shape == (1, 1), f'expected b length 1, got {bm.shape}'


def test_slice_rank3_emits_selector_fc():
    # Rank-3 Slice on a [1,1,8] input selecting axis=2 [0:8] = identity-like
    # but should still emit an FC selector (general case).
    nodes = [
        helper.make_node('Slice', ['x', 'starts', 'ends', 'axes', 'steps'], ['y']),
    ]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 1, 8])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 1, 4])
    inits = [
        helper.make_tensor('starts', TensorProto.INT64, [1], [0]),
        helper.make_tensor('ends',   TensorProto.INT64, [1], [4]),
        helper.make_tensor('axes',   TensorProto.INT64, [1], [2]),
        helper.make_tensor('steps',  TensorProto.INT64, [1], [1]),
    ]
    m = _make_onnx(nodes, [inp], [out], inits)
    layers, weights, _, _, _ = O.walk(m)
    fc_specs = [L for L in layers if L.type == 'FullyConnectedLayer']
    assert len(fc_specs) >= 1, 'expected Slice to emit an FC selector'
    Wm = weights[fc_specs[0].weight_keys[0]]
    assert Wm.shape == (4, 8), f'expected selector W (4,8), got {Wm.shape}'
    # First 4 cols, identity selector (one-hot)
    expected = np.eye(8, dtype=np.float32)[:4]
    assert np.allclose(Wm, expected), 'selector should pick first 4 of 8'


def test_gather_constant_index_emits_selector_fc():
    # Gather on [1,1,8] with scalar index 3, axis=2 → output [1,1] with the
    # 3rd element. Should emit a 1×8 selector.
    nodes = [
        helper.make_node('Gather', ['x', 'idx'], ['y'], axis=2),
    ]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 1, 8])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 1])
    inits = [helper.make_tensor('idx', TensorProto.INT64, [], [3])]
    m = _make_onnx(nodes, [inp], [out], inits)
    layers, weights, _, _, _ = O.walk(m)
    fc_specs = [L for L in layers if L.type == 'FullyConnectedLayer']
    assert len(fc_specs) == 1, 'expected Gather to emit one FC selector'
    Wm = weights[fc_specs[0].weight_keys[0]]
    assert Wm.shape == (1, 8), f'expected selector W (1,8), got {Wm.shape}'
    expected = np.zeros((1, 8), dtype=np.float32); expected[0, 3] = 1.0
    assert np.allclose(Wm, expected), 'selector should pick index 3'


def test_no_spatial_op_input_is_feature_layer():
    # Pure-feature graph (no Conv/Pool) should produce FeatureInputLayer
    # even when the declared input rank is >2.
    W = np.random.randn(2, 8).astype(np.float32)
    nodes = [helper.make_node('Gemm', ['x', 'W'], ['y'], transB=1)]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 1, 8])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 2])
    m = _make_onnx(nodes, [inp], [out], [_const('W', W)])
    layers, _, _, _, _ = O.walk(m)
    assert layers[0].type == 'FeatureInputLayer', \
        f'expected FeatureInputLayer for no-spatial graph, got {layers[0].type}'


def test_concat_includes_inrank_attr():
    # Concat handler should stash InRank so the loader can pick the right
    # MATLAB cat dim from the ONNX axis.
    nodes = [helper.make_node('Concat', ['a', 'b'], ['y'], axis=-1)]
    a = helper.make_tensor_value_info('a', TensorProto.FLOAT, [1, 3])
    b = helper.make_tensor_value_info('b', TensorProto.FLOAT, [1, 5])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 8])
    # Use Add to force walk to start at concat's dependencies present as inputs
    m = _make_onnx(nodes, [a, b], [out], [])
    layers, _, _, _, _ = O.walk(m)
    cc = [L for L in layers if L.type == 'ConcatenationLayer']
    assert len(cc) == 1
    assert cc[0].attrs.get('InRank') == 2, f'expected InRank=2, got {cc[0].attrs}'


def test_reshape_rank2_feature_flow_emits_reshape_not_flatten():
    # When the upstream tensor is [1, F] (feature flow), Reshape to [1, F']
    # should produce ReshapeLayer with target [F', 1] (column vector), not
    # FlattenLayer (which adds singleton dims for image flattening). Add a
    # value_info entry so the walker sees the intermediate shape.
    W = np.random.randn(8, 5).astype(np.float32)
    nodes = [
        helper.make_node('Gemm', ['x', 'W'], ['z'], transB=1),
        helper.make_node('Reshape', ['z', 'shape'], ['y']),
    ]
    inp = helper.make_tensor_value_info('x', TensorProto.FLOAT, [1, 5])
    out = helper.make_tensor_value_info('y', TensorProto.FLOAT, [1, 8])
    z_info = helper.make_tensor_value_info('z', TensorProto.FLOAT, [1, 8])
    inits = [_const('W', W),
             helper.make_tensor('shape', TensorProto.INT64, [2], [1, -1])]
    g = helper.make_graph(nodes, 'm', [inp], [out],
                          initializer=inits, value_info=[z_info])
    m = helper.make_model(g, opset_imports=[helper.make_opsetid("", 14)])
    m.ir_version = 7
    layers, _, _, _, _ = O.walk(m)
    types = [L.type for L in layers]
    assert 'ReshapeLayer' in types, f'expected feature-flow Reshape, got {types}'
    assert 'FlattenLayer' not in types, f'should not emit FlattenLayer for feature flow, got {types}'


# Test runner: run all functions whose name starts with test_
def _run_all():
    fns = sorted(k for k in globals() if k.startswith('test_'))
    n_pass = 0
    n_fail = 0
    for fn_name in fns:
        try:
            globals()[fn_name]()
            print(f'  PASS  {fn_name}')
            n_pass += 1
        except Exception as e:
            print(f'  FAIL  {fn_name}: {e}')
            n_fail += 1
    print(f'\n{n_pass}/{n_pass + n_fail} tests passed')
    return 0 if n_fail == 0 else 1


if __name__ == '__main__':
    sys.exit(_run_all())
