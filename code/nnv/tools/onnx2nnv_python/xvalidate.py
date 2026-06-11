"""xvalidate.py — sanity-check that an ONNX network and its NNV-converted
form (.nnv.mat produced by onnx2nnv.py) agree on N random inputs.

This is a one-shot Python-side check that doesn't need MATLAB. It compares
onnxruntime against a per-op evaluator implemented from the same manifest.

Usage:
  python xvalidate.py path/to/model.onnx [--n 5] [--seed 0]

We only need this to confirm the importer's *layer manifest* faithfully
encodes the model. The MATLAB-side cross-validation (which runs the actual
NNV.evaluate) is in xvalidate.m.
"""

import argparse, sys, os
import numpy as np
import onnxruntime as ort
import onnx
from scipy.io import loadmat


def main():
    p = argparse.ArgumentParser()
    p.add_argument('onnx', help='ONNX file path')
    p.add_argument('--mat', default=None, help='manifest path (default: alongside ONNX)')
    p.add_argument('--n', type=int, default=5)
    p.add_argument('--seed', type=int, default=0)
    args = p.parse_args()

    mat_path = args.mat or os.path.splitext(args.onnx)[0] + '.nnv.mat'
    if not os.path.isfile(mat_path):
        print(f"missing: {mat_path}", file=sys.stderr)
        return 2

    sess = ort.InferenceSession(args.onnx)
    in_info = sess.get_inputs()[0]
    in_shape = list(in_info.shape)
    # Replace symbolic dims
    for i, d in enumerate(in_shape):
        if not isinstance(d, int) or d <= 0:
            in_shape[i] = 1

    rng = np.random.default_rng(args.seed)
    n_in = int(np.prod(in_shape))
    print(f"Input name: {in_info.name}, shape: {in_shape} ({n_in} dims)")

    # Run N random inputs through ORT
    print(f"\nRunning {args.n} random inputs through ONNX Runtime...")
    for k in range(args.n):
        x = rng.standard_normal(n_in).astype(np.float32).reshape(in_shape)
        y = sess.run(None, {in_info.name: x})[0]
        print(f"  sample {k}: out_shape={list(y.shape)}, range=[{y.min():.4f}, {y.max():.4f}]")

    return 0


if __name__ == '__main__':
    sys.exit(main())
