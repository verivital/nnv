"""Dispatch-time single-ONNX -> weights-only .mat extractor for the VNN-COMP 2023 ViT.

Same parsing contract as extract_weights.py (initializer-by-name + weight-bearing
MatMul-by-graph-order + BatchNorm-fold to per-channel scale/shift), reduced to ONE
onnx in / ONE .mat out with NO images, NO labels, NO onnxruntime oracle: the vnnlib
supplies the verification box at dispatch, so the bundle only needs the weights +
architecture scalars that ViTReach.load / ViTCrown.toOps consume.

MATLAB's importNetworkFromONNX is unusable here (opset-9 Slice-v1 + the per-token
BatchNorm transpose make the importer fuse/refuse the graph -- see ViTReach.m header
and model.py), so a protobuf-level extractor is required; this is it.

Usage (dispatch):
    python extract_weights_one.py MODEL.onnx OUT.mat
       [--mean R G B] [--std R G B]   # CIFAR-10 defaults baked in

depth and patch are read straight from the initializers (no config dict needed),
so the SAME script handles pgd_2_3_16, ibp_3_3_8, or any ViT_BN-shaped export.
"""
import argparse
import numpy as np
import onnx
from onnx import numpy_helper
from scipy.io import savemat

# CIFAR-10 normalisation (the vnnlib box is in normalised space already in this
# benchmark, but ViTReach.epsBox/normImg keep mean/std so the .mat stays self-
# describing and a raw-[0,1] box can also be normalised MATLAB-side if needed).
MEAN = np.array([0.4914, 0.4822, 0.4465], dtype=np.float64)
STD = np.array([0.2023, 0.1994, 0.2010], dtype=np.float64)
BN_EPS = 1e-5


def arr(init):
    return numpy_helper.to_array(init).astype(np.float64)


def bn_affine(inits, prefix):
    """Fold a BatchNorm (weight/bias/running_mean/running_var) to per-channel
    (scale, shift) so y = scale*x + shift -- an exact affine map in eval mode."""
    w = arr(inits[prefix + ".weight"])
    b = arr(inits[prefix + ".bias"])
    rm = arr(inits[prefix + ".running_mean"])
    rv = arr(inits[prefix + ".running_var"])
    scale = w / np.sqrt(rv + BN_EPS)
    shift = b - scale * rm
    return scale, shift


def extract_one(onnx_path, mean=MEAN, std=STD):
    """Parse ONE benchmark ONNX into the weights-only field dict ViTReach.load wants."""
    proto = onnx.load(str(onnx_path))
    inits = {i.name: i for i in proto.graph.initializer}

    # architecture scalars, inferred from the initializers (no config dict)
    pos = arr(inits["0.positions"])                       # (N, E)
    proj_w = arr(inits["0.projection.weight"])            # (E, 3, p, p)
    N, E = pos.shape
    patch = proj_w.shape[-1]
    depth = len({k.split(".")[1] for k in inits
                 if k.startswith("1.") and ".0.fn.0.norm.weight" in k})

    out = {
        "patch": patch,
        "depth": depth,
        "dim": E,             # 48
        "heads": 3,           # fixed by the benchmark family (ViT_BN)
        "dim_head": 16,
        "mean": mean,
        "std": std,
        # patch-embed + token bookkeeping
        "proj_w": proj_w,                                 # (E,3,p,p)
        "proj_b": arr(inits["0.projection.bias"]),        # (E,)
        "cls": arr(inits["0.cls_token"]).reshape(-1),     # (E,)
        "pos": pos,                                        # (N,E)
    }

    # weight-bearing MatMuls in graph order: per block Q,K,V,out,ff.l0,ff.l3.
    # (QK^T and A*V have activation inputs on both sides -> no initializer ->
    # excluded.) ONNX MatMul weight is (in,out); stored ONNX-native (ViTReach
    # does Xt*Mw, so no transpose).
    matmul_nodes = [n for n in proto.graph.node
                    if n.op_type == "MatMul" and any(inp in inits for inp in n.input)]
    if len(matmul_nodes) != 6 * depth:
        raise ValueError(f"expected {6*depth} weight-matmuls, got {len(matmul_nodes)}")

    def mm_w(idx):
        node = matmul_nodes[idx]
        wname = next(n for n in node.input if n in inits)
        return arr(inits[wname])                          # (in, out)

    blocks = []
    for i in range(depth):
        blk = {}
        blk["bn1_scale"], blk["bn1_shift"] = bn_affine(inits, f"1.{i}.0.fn.0.norm")
        blk["Mq"], blk["bq"] = mm_w(i*6+0), arr(inits[f"1.{i}.0.fn.1.query.bias"])
        blk["Mk"], blk["bk"] = mm_w(i*6+1), arr(inits[f"1.{i}.0.fn.1.key.bias"])
        blk["Mv"], blk["bv"] = mm_w(i*6+2), arr(inits[f"1.{i}.0.fn.1.value.bias"])
        blk["Mo"], blk["bo"] = mm_w(i*6+3), arr(inits[f"1.{i}.0.fn.1.out.bias"])
        blk["bn2_scale"], blk["bn2_shift"] = bn_affine(inits, f"1.{i}.1.fn.0.norm")
        blk["ff1_W"], blk["ff1_b"] = mm_w(i*6+4), arr(inits[f"1.{i}.1.fn.1.0.bias"])
        blk["ff2_W"], blk["ff2_b"] = mm_w(i*6+5), arr(inits[f"1.{i}.1.fn.1.3.bias"])
        blocks.append(blk)
    out["blocks"] = blocks

    # head: BN (folded) + linear. head_W is (10,E) = torch (out,in), used directly.
    out["head_bn_scale"], out["head_bn_shift"] = bn_affine(inits, "2.1")
    out["head_W"] = arr(inits["2.2.weight"])              # (10,E)
    out["head_b"] = arr(inits["2.2.bias"])                # (10,)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("onnx", help="path to ONE benchmark .onnx")
    ap.add_argument("mat", help="output weights-only .mat path")
    ap.add_argument("--mean", type=float, nargs=3, default=None)
    ap.add_argument("--std", type=float, nargs=3, default=None)
    a = ap.parse_args()
    mean = np.array(a.mean, dtype=np.float64) if a.mean else MEAN
    std = np.array(a.std, dtype=np.float64) if a.std else STD
    out = extract_one(a.onnx, mean=mean, std=std)
    savemat(a.mat, out, do_compression=True)
    b0 = out["blocks"][0]
    print(f"[{a.onnx}] depth={out['depth']} patch={out['patch']} dim={out['dim']} "
          f"N={out['pos'].shape[0]} proj_w{out['proj_w'].shape} "
          f"Mq{b0['Mq'].shape} ff1_W{b0['ff1_W'].shape} head_W{out['head_W'].shape} "
          f"-> {a.mat}")


if __name__ == "__main__":
    main()
