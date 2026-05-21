# NNV 3.0 — Online (individual) MATLAB licence path

The standard `Dockerfile` at the repository root builds the artifact
image against a **network MATLAB licence** (port@host via
`--build-arg LICENSE_SERVER=`). That path is the canonical one for
ATVA 2026 artifact-evaluation reviewers whose institutions provide a
licence server.

This document covers the alternative `Dockerfile.online`, which
targets users with an **individual (signed-in) MathWorks account
licence** -- the kind issued for a single MathWorks profile via
`matlab.mathworks.com` and used by `mathworks/matlab:r2025b`'s browser
sign-in. Network-licence reviewers do **not** need this path; use the
standard `Dockerfile`.

## When to use this

Pick `Dockerfile.online` if **any** of:

- You have a MathWorks account but no `port@host` licence server.
- Your institution's licence server isn't reachable from your build
  host (firewall / VPN issues).
- You want to verify the artifact builds without sharing or hardcoding
  a network licence value.

Pick the standard `Dockerfile` if you have a network licence -- it's
faster overall (one build, no interactive sign-in, fully unattended
`run_all.sh`).

## Trade-offs vs the network-licence path

| | Standard `Dockerfile` | `Dockerfile.online` |
|---|---|---|
| MATLAB licence | network (`port@host`) | individual (MathWorks account) |
| Build-time MATLAB invocations | `install.m`, AIVL extract | none |
| `docker build` time | ~15-25 min | ~20 min |
| `docker build` needs licence? | yes (consumed at first `matlab` step) | no |
| First-run setup | none | sign in via browser, run `setup-online-license.sh` (~5 min) |
| Headless runs after setup | `docker run … bash run_all.sh` | same, with two volume mounts |
| Sign-in persistence | n/a | cached in named Docker volumes; may need periodic refresh |
| Multi-user fan-out (e.g. CI) | scales with licence-server capacity | one container per signed-in MathWorks user at a time |

## Step-by-step

### Step 1 — Build the image (~20 min, no licence needed)

```bash
docker build -t nnv3.0-online -f Dockerfile.online .
```

The build installs the same toolbox set the network-licence path uses,
**plus the Deep Learning Toolbox Verification Library** (the AIVL
Support Package), via `mpm`. `mpm` does not validate a MATLAB licence,
so this completes without any sign-in.

> **Licence requirement for AIVL.** This flow assumes the reviewer's
> MathWorks account entitlement includes the **Deep Learning Toolbox
> Verification Library**. It is normally free for any Deep Learning
> Toolbox licence holder, but some restricted seats may exclude it.
> If excluded, `mpm install` still succeeds but
> `verifyNetworkRobustness` / `estimateNetworkOutputBounds` fail at
> runtime with a licence-checkout error; the AIVL availability check
> in `setup_online.m` will flag this and ToolComparison will run
> NNV-only.

> **Windows / PowerShell users**: the same substitutions noted at the
> top of the root README apply to every `docker run` below (`${PWD}`,
> full-quoted `-v` values, no `\` line continuations).

### Step 2 — Sign in + one-time setup (~5 min, interactive, once)

Start the container in browser mode with `--gpus all` and two named
volumes that will hold the cached MATLAB activation and the NNV / AIVL
install state:

```bash
docker run -it --rm --gpus all \
    -p 8888:8888 --shm-size=512M \
    -v nnv3-matlab-prefs:/home/matlab/.matlab \
    -v nnv3-matlab-mw:/home/matlab/.MathWorks \
    --name nnv3-setup \
    nnv3.0-online -browser
```

Drop `--gpus all` on CPU-only hosts; ProbVer will auto-skip and
GNNV / VideoStar fall back to CPU (with proportional slowdowns).

> **Windows + Docker Desktop GPU passthrough**: `--gpus all` requires
> the **WSL2 backend** plus an **NVIDIA driver on the host** (not inside
> WSL2). Enable WSL2 in Docker Desktop *Settings → General* and confirm
> *Settings → Resources → WSL Integration* lists your distro. Sanity
> check with: `docker run --rm --gpus all nnv3.0-online nvidia-smi`
> — it should print your device name and driver version. If it errors
> with `could not select device driver`, install/enable the NVIDIA
> Container Toolkit (or, on Docker Desktop, update to a recent version
> that ships it bundled). From inside the MATLAB Command Window, verify
> with `gpuDeviceCount` (returns `1` or more if passthrough is working).

Open <http://localhost:8888> in a browser and sign in with your
MathWorks account. Once the MATLAB Command Window is ready, paste
this single line to install NNV paths and verify AIVL:

```matlab
run('/home/matlab/nnv/code/nnv/examples/NNV3.0/setup_online.m')
```

`setup_online.m` runs `install.m` (equivalent to the standard
Dockerfile's build-time NNV path setup, deferred here because online
sign-in is not available at build time) and then checks that AIVL is
on the path. AIVL itself is already installed by the Dockerfile build
via `mpm` -- no manual tarball staging needed. The `AIVL availability
check` at the end of the script confirms whether ToolComparison will
include the MathWorks-side rows. When the script finishes, do **not**
exit the container: the smoke / full runs below reuse this same
MATLAB session (and its cached licence).

### Step 3 — Run experiments in the same browser session

Online sign-in is only valid for the current MATLAB process; a fresh
`matlab -batch` launched from `docker exec` cannot use it. So unlike
the network-licence path, experiments run **inside the same browser
MATLAB** that owns the sign-in. Two single-line entry points:

Smoke (~20-25 min on CPU, faster with `--gpus all`):

```matlab
run('/home/matlab/nnv/code/nnv/examples/NNV3.0/run_smoke.m')
```

Full reproduction (~3-5 h, renders ATVA 2026 paper Tables 5, 6, 7):

```matlab
run('/home/matlab/nnv/code/nnv/examples/NNV3.0/run_full.m')
```

Both runners loop over the six experiments in this session, autodetect
GPU availability (ProbVer auto-skips on CPU-only hosts), consolidate
per-experiment outputs into `repeatability_logs/results/`, and write a
single index file `repeatability_logs/PAPER_COMPARISON.md` mapping each
experiment to its primary output and to the paper artefact it backs.

### Step 4 — Extract results to the host

From a separate PowerShell on the host:

```powershell
docker cp nnv3-setup:/home/matlab/nnv/code/nnv/examples/NNV3.0/repeatability_logs .\repeatability_logs
```

Then open `repeatability_logs/PAPER_COMPARISON.md` -- it is the single
file a repeatability committee can open to find every output, with
explicit pointers to `results/ToolComparison/tables/table_main.{tex,txt}`
(paper Tables 5, 6, 7) and to each experiment's CSV / `.mat` outputs.

## Troubleshooting

**`license checkout failed` mid-run.** MathWorks may invalidate cached
online activations periodically. Restart from Step 2: launch the
container in `-browser` mode, sign in again, and re-run `setup_online`
+ `run_smoke` / `run_full`. The named volumes preserve the saved
MATLAB path and AIVL install, so only the activation needs refreshing.

**MATLAB browser at localhost:8888 won't load.** Check Docker Desktop
has port 8888 free and `--shm-size=512M` was passed (insufficient
shared memory will hang the VNC server).

**`setup_online.m` reports `AIVL NOT found on the MATLAB path`.** The
`mpm` install in the Dockerfile did not pull
`Deep_Learning_Toolbox_Verification_Library`, or your MathWorks
account's Deep Learning Toolbox entitlement does not include it.
Check inside the container that the support package was extracted:

```bash
docker exec -it nnv3-setup ls /opt/matlab/R2025b/toolbox/nnet/supportpackages/
```

If `aivnv/` is missing, rebuild from a clean image
(`docker build --no-cache -t nnv3.0-online -f Dockerfile.online .`)
and watch the mpm install output for a "skipped: not licensed" line
for `Deep_Learning_Toolbox_Verification_Library`. If the entitlement
is missing, AIVL cannot be installed; ToolComparison will still run
NNV-only and the rest of the artefact is unaffected.

or from a fresh container with the same volume mounted, then re-run
`setup-online-license.sh`.

**Fresh start**: `docker volume rm nnv3-matlab-prefs nnv3-matlab-mw`
wipes the cached activation and NNV install. The next `docker run -it
… -browser` will require a fresh sign-in.
