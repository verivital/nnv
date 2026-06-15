# VNN-COMP 2026 — AWS box tooling (reproducible sweeps)

Scripts for setting up an AWS box and running the full NNV benchmark sweep. **Everything is
pulled from public repositories** — no box-to-box file transfer. The only out-of-band item is
the MATLAB license (it is a **secret**; never commit it).

## Source-of-truth repositories (all public)

| What | Repo |
|------|------|
| NNV (this code) | `ttj/nnv` (fork) → PRs to `verivital/nnv`; migration branch `feat/vnncomp2026-migration` |
| 2025 benchmarks | `https://github.com/VNN-COMP/vnncomp2025_benchmarks` |
| 2026 benchmarks | `https://github.com/VNN-COMP/vnncomp2026_benchmarks` |
| 2025 official results | `https://github.com/VNN-COMP/vnncomp2025_results` |
| Strategy / status | `https://github.com/ttj/nnv-vnncomp2026-status` |

All official VNN-COMP repos are under the **`VNN-COMP` org** (centralized in 2025 from legacy
`ChristopherBrix/*` locations).

## One-command fresh-box bootstrap

```bash
bash setup_aws_box.sh           # clone NNV + benchmarks + results + status, decompress,
                                # install ONNX/PyTorch converters (mpm) + NNV + python + matlab.engine
```
Idempotent (re-run safe). Disk budget ~60 GB. Tunables via env: `NNV_BRANCH` (default the 2026
migration branch — switch to `master` after merge), `RUN_BENCH_SETUP=0` to skip decompress,
`*_REPO` overrides. Provision the MATLAB license separately (file `.lic` into
`~/.matlab/R2026a_licenses/`, or `MLM_LICENSE_FILE=port@host`, or the AMI's online sign-in).

## Pre-sweep gate (run this FIRST, every time)

```bash
bash smoke_test_benchmarks.sh           # 1 instance/benchmark, short timeout; flags ERRORING benchmarks
```
Runs one representative instance per benchmark through the **same `run_vnncomp_instance` path the
sweep uses**, so it catches load/ONNX-import/dispatch breakage (and, once GPU-BaB is routed, any
GPU-BaB breakage) in minutes — before a multi-hour sweep is wasted. Exits non-zero if any benchmark
errors. Run it before every sweep and after any benchmark or NNV code change. Known errored
benchmarks to clear first are tracked in the status repo (`sweeps/.../ERRORS_TO_FIX.md`).

## Sweep workflow

| Script | Role |
|--------|------|
| `smoke_test_benchmarks.sh` | **Pre-sweep gate** (above): 1 instance/benchmark, reports erroring benchmarks. |
| `overnight_sweep.sh` | Full 2026 sweep: 11 concurrent `matlab -batch` workers over benchmark folder-groups, 300 s/instance, timestamped CSVs. Checks out `NNV_BRANCH` (default migration branch) and parse-guards `run_all_benchmarks.m`. |
| `disk_guard.sh` | Reclaims decompressed `.onnx` files (each has a `.gz` sibling) during the sweep so the disk doesn't fill. |
| `harvest_and_selfstop.sh` | Waits for the sweep to finish, tars results onto the persistent EBS, then **stops** (not terminates) the box to save cost. |
| `fix_validation.sh` | Re-runs the previously-erroring/timing-out benchmarks to confirm merged fixes unlocked them. |
| `preimport.m`, `dryrun_install.m`, `orchestrate_dryrun.sh` | Dry-run helpers (manifest pre-import, install timing, install→suite orchestration). |

## License & cost notes

- The MATLAB license is a **secret** — never in any repo. The MathWorks AWS AMI forces online
  licensing; the token is flaky across stop/start but **recovers on a clean restart**.
- An `EC2AutoShutdown` Lambda stops the box on an uptime timer — fine for cost, but it will stop
  an active sweep; `harvest_and_selfstop.sh` self-stops when done so the box isn't left running.
