"""onnx_replay_check.py -- validate a SAT counterexample against the ORIGINAL ONNX
model via onnxruntime (the competition's actual checker), not NNV's forward.

VNN-COMP scores incorrect/non-replayable verdicts at -150; a witness that NNV's
self-consistent replay accepts could still fail the competition if NNV's import
permutes the input differently than the ONNX model expects. This replays the flat
witness through onnxruntime and checks whether the output lands in the unsafe
region G*y <= g, so the runner can refuse to emit `sat` for a non-replayable witness.

Usage:
  python onnx_replay_check.py <model.onnx> <x.csv> <G.csv> <g.csv> [--out-tol 1e-4]

  x.csv : the counterexample input, FLAT in ONNX order (one value per line or CSV).
  G.csv : the unsafe-region matrix G (rows = constraints, cols = output dim).
  g.csv : the unsafe-region vector g (one value per line).

Prints one of:
  VIOLATED max_margin=<m>   -> witness is a TRUE counterexample (all G*y - g <= tol)
  OK       max_margin=<m>   -> witness does NOT violate (do not emit sat)
  ERROR    <message>        -> could not validate (treat as not-sat)
Exit code: 0 VIOLATED, 1 OK, 2 ERROR.
"""

import sys
import argparse
import numpy as np


def _load_vec(path):
    return np.loadtxt(path, delimiter=',').astype(np.float64).ravel()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('onnx')
    p.add_argument('x_csv')
    p.add_argument('G_csv')
    p.add_argument('g_csv')
    p.add_argument('--out-tol', type=float, default=1e-4)
    args = p.parse_args()

    try:
        import onnxruntime as ort
    except Exception as e:  # pragma: no cover
        print(f"ERROR onnxruntime import failed: {e}")
        return 2

    try:
        # Load everything FLAT (ravel) so the result is robust to how the CSV writer
        # laid the values out (one-per-line vs comma-separated, scalar vs vector).
        # G is reshaped to (nrows, outdim) below, once the model's output dim is known.
        x = _load_vec(args.x_csv).astype(np.float32)
        G_flat = _load_vec(args.G_csv)
        g = _load_vec(args.g_csv)

        sess = ort.InferenceSession(args.onnx, providers=['CPUExecutionProvider'])
        inp = sess.get_inputs()[0]
        shape = [d if isinstance(d, int) and d > 0 else 1 for d in inp.shape]
        try:
            xr = x.reshape(shape)
        except ValueError:
            # element count differs from the declared shape (e.g. a symbolic/batch
            # dim): flatten to a single batch row of the known data size.
            xr = x.reshape([1, x.size])

        y = np.asarray(sess.run(None, {inp.name: xr})[0]).ravel().astype(np.float64)
        outdim = y.size
        if G_flat.size % outdim != 0 or G_flat.size // outdim != g.size:
            print(f"ERROR G has {G_flat.size} entries, incompatible with "
                  f"{g.size} constraints x {outdim} outputs")
            return 2
        G = G_flat.reshape(g.size, outdim)

        margins = G.dot(y) - g
        mmax = float(np.max(margins))
        if np.all(margins <= args.out_tol):
            print(f"VIOLATED max_margin={mmax:.6e}")
            return 0
        print(f"OK max_margin={mmax:.6e}")
        return 1
    except Exception as e:
        print(f"ERROR {e}")
        return 2


if __name__ == '__main__':
    sys.exit(main())
