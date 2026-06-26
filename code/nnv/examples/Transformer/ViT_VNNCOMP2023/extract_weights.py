"""Extract VNN-COMP 2023 ViT weights from ONNX into a .mat the NNV driver loads.

Mirrors n2v benchmarks/vit_vnncomp2023/model.py::load_from_onnx exactly (same
initializer names, same MatMul-by-order convention), but stores ONNX-native
matrices (in,out) and pre-folds BatchNorm to per-channel (scale, shift) so the
MATLAB side needs no torch/onnx. Also vendors the 100 inputs+labels per model
and the CIFAR-10 mean/std so the driver can build the normalized eps-box.

Usage:  python extract_weights.py [--benchmark-root PATH]
        The benchmark checkout (the n2v benchmarks/vit_vnncomp2023 directory,
        containing onnx/ and instances/) is located, in order, from:
          1) --benchmark-root PATH
          2) the $N2V_VIT_BENCHMARK environment variable
          3) ./vit_vnncomp2023 next to this script (vendored/symlinked checkout)
Outputs (next to this script):  pgd_2_3_16.mat, ibp_3_3_8.mat
"""
import argparse
import os
from pathlib import Path
import numpy as np
import onnx
from onnx import numpy_helper
from scipy.io import savemat

HERE = Path(__file__).resolve().parent


def resolve_benchmark_root(cli_root=None):
    """Locate the n2v vit_vnncomp2023 benchmark checkout (see module docstring)."""
    if cli_root is not None:
        root = Path(cli_root)
    elif os.environ.get("N2V_VIT_BENCHMARK"):
        root = Path(os.environ["N2V_VIT_BENCHMARK"])
    else:
        root = HERE / "vit_vnncomp2023"
    root = root.expanduser().resolve()
    if not (root / "onnx").is_dir() or not (root / "instances").is_dir():
        raise SystemExit(
            f"benchmark root {root} is missing onnx/ or instances/. Pass "
            f"--benchmark-root PATH or set $N2V_VIT_BENCHMARK to the n2v "
            f"benchmarks/vit_vnncomp2023 checkout.")
    return root

MEAN = np.array([0.4914, 0.4822, 0.4465], dtype=np.float64)
STD = np.array([0.2023, 0.1994, 0.2010], dtype=np.float64)

CONFIGS = {"pgd_2_3_16": dict(patch=16, depth=2),
           "ibp_3_3_8": dict(patch=8, depth=3)}
EPS = 1e-5


def arr(init):
    return numpy_helper.to_array(init).astype(np.float64)


def bn_affine(inits, prefix):
    w = arr(inits[prefix + ".weight"])
    b = arr(inits[prefix + ".bias"])
    rm = arr(inits[prefix + ".running_mean"])
    rv = arr(inits[prefix + ".running_var"])
    scale = w / np.sqrt(rv + EPS)
    shift = b - scale * rm
    return scale.astype(np.float64), shift.astype(np.float64)


def extract(name, patch, depth, onnx_dir, npz_dir):
    proto = onnx.load(str(onnx_dir / f"{name}.onnx"))
    inits = {i.name: i for i in proto.graph.initializer}

    out = {}
    out["name"] = name
    out["patch"] = patch
    out["depth"] = depth
    out["dim"] = 48
    out["heads"] = 3
    out["dim_head"] = 16
    out["mean"] = MEAN
    out["std"] = STD

    # patch embed conv: weight (48,3,p,p), bias (48,)
    out["proj_w"] = arr(inits["0.projection.weight"])  # (48,3,p,p)
    out["proj_b"] = arr(inits["0.projection.bias"])    # (48,)
    out["cls"] = arr(inits["0.cls_token"]).reshape(-1)  # (48,)
    out["pos"] = arr(inits["0.positions"])              # (N, 48)

    # weight-bearing MatMuls in order: per block Q,K,V,out,ff.l0,ff.l3
    matmul_nodes = [n for n in proto.graph.node
                    if n.op_type == "MatMul" and any(inp in inits for inp in n.input)]
    assert len(matmul_nodes) == 6 * depth, (len(matmul_nodes), 6 * depth)

    def mm_w(idx):
        node = matmul_nodes[idx]
        wname = next(n for n in node.input if n in inits)
        return arr(inits[wname])  # ONNX MatMul weight: (in, out)

    blocks = []
    for i in range(depth):
        blk = {}
        s, sh = bn_affine(inits, f"1.{i}.0.fn.0.norm")
        blk["bn1_scale"], blk["bn1_shift"] = s, sh
        blk["Mq"], blk["bq"] = mm_w(i*6+0), arr(inits[f"1.{i}.0.fn.1.query.bias"])
        blk["Mk"], blk["bk"] = mm_w(i*6+1), arr(inits[f"1.{i}.0.fn.1.key.bias"])
        blk["Mv"], blk["bv"] = mm_w(i*6+2), arr(inits[f"1.{i}.0.fn.1.value.bias"])
        blk["Mo"], blk["bo"] = mm_w(i*6+3), arr(inits[f"1.{i}.0.fn.1.out.bias"])
        s2, sh2 = bn_affine(inits, f"1.{i}.1.fn.0.norm")
        blk["bn2_scale"], blk["bn2_shift"] = s2, sh2
        blk["ff1_W"], blk["ff1_b"] = mm_w(i*6+4), arr(inits[f"1.{i}.1.fn.1.0.bias"])
        blk["ff2_W"], blk["ff2_b"] = mm_w(i*6+5), arr(inits[f"1.{i}.1.fn.1.3.bias"])
        blocks.append(blk)
    out["blocks"] = blocks

    hs, hsh = bn_affine(inits, "2.1")
    out["head_bn_scale"], out["head_bn_shift"] = hs, hsh
    out["head_W"] = arr(inits["2.2.weight"])  # (in=48, ... ) check orientation below
    out["head_b"] = arr(inits["2.2.bias"])

    # vendored inputs
    npz = np.load(npz_dir / f"{name}.npz")
    images = npz["images"].astype(np.float64)  # (100,3,32,32) in [0,1]
    out["images"] = images
    out["labels"] = npz["labels"].astype(np.float64)  # (100,)

    # onnxruntime reference logits on the first 5 NORMALISED images (parity oracle)
    import onnxruntime as ort
    sess = ort.InferenceSession(str(onnx_dir / f"{name}.onnx"),
                                providers=["CPUExecutionProvider"])
    iname = sess.get_inputs()[0].name
    norm = (images - MEAN.reshape(1, 3, 1, 1)) / STD.reshape(1, 3, 1, 1)
    ref = []
    for k in range(5):
        xin = norm[k:k+1].astype(np.float32)
        y = sess.run(None, {iname: xin})[0].reshape(-1)
        ref.append(y.astype(np.float64))
    out["ref_logits"] = np.array(ref)        # (5,10)
    out["ref_norm_images"] = norm[:5]        # (5,3,32,32) normalised

    savemat(str(HERE / f"{name}.mat"), out, do_compression=True)
    # report shapes for a sanity log
    print(f"[{name}] proj_w{out['proj_w'].shape} pos{out['pos'].shape} "
          f"Mq{blocks[0]['Mq'].shape} ff1_W{blocks[0]['ff1_W'].shape} "
          f"head_W{out['head_W'].shape} images{out['images'].shape}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description="Extract VNN-COMP 2023 ViT weights from ONNX into .mat files.")
    ap.add_argument("--benchmark-root", default=None,
                    help="n2v benchmarks/vit_vnncomp2023 checkout (with onnx/ and "
                         "instances/). Defaults to $N2V_VIT_BENCHMARK or "
                         "./vit_vnncomp2023 next to this script.")
    args = ap.parse_args()
    root = resolve_benchmark_root(args.benchmark_root)
    onnx_dir = root / "onnx"
    npz_dir = root / "instances"
    for nm, cfg in CONFIGS.items():
        extract(nm, onnx_dir=onnx_dir, npz_dir=npz_dir, **cfg)
    print("OK")
