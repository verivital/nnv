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
via `mpm`. `mpm` does not validate a MATLAB licence, so this completes
without any sign-in.

> **Windows / PowerShell users**: the same substitutions noted at the
> top of the root README apply to every `docker run` below (`${PWD}`,
> full-quoted `-v` values, no `\` line continuations).

### Step 2 — Sign in + one-time setup (~5 min, interactive, once)

Start the container in browser mode with two named volumes that will
hold the cached MATLAB activation and the NNV / AIVL install state:

```bash
docker run -it --rm \
    -p 8888:8888 --shm-size=512M \
    -v nnv3-matlab-prefs:/home/matlab/.matlab \
    -v nnv3-matlab-mw:/home/matlab/.MathWorks \
    --name nnv3-setup \
    nnv3.0-online -browser
```

Open <http://localhost:8888> in a browser, sign in with your MathWorks
account. Once MATLAB is running in the browser, drop to a second
terminal and run the setup script inside the same container:

```bash
docker exec -it nnv3-setup \
    bash /home/matlab/nnv/code/nnv/examples/NNV3.0/utils/setup-online-license.sh
```

That runs NNV's `install.m` and the AIVL `toolbox_install.m`, then
writes a marker so the next setup invocation is a no-op. When it
finishes, exit the first terminal's interactive container (Ctrl-C).

### Step 3 — Run experiments (headless, ~30 min smoke)

```bash
docker run --rm --gpus all \
    -v nnv3-matlab-prefs:/home/matlab/.matlab \
    -v nnv3-matlab-mw:/home/matlab/.MathWorks \
    -v "$PWD/results":/out \
    --entrypoint /bin/bash nnv3.0-online \
    -lc "bash run_all.sh && \
        cp -r /home/matlab/nnv/code/nnv/examples/NNV3.0/repeatability_logs /out/"
```

Drop `--gpus all` on CPU-only hosts (ProbVer auto-skips). The two
named volumes carry the cached activation and the NNV install from
Step 2, so `run_all.sh` runs headlessly with no sign-in prompt.

For the full ~5-7 h reproduction, add `-e TOOLCOMPARISON_MODE=full`
inside the run command, as in the standard path.

## Troubleshooting

**`license checkout failed` mid-run.** MathWorks may invalidate cached
online activations periodically. Re-run Step 2 (the named volumes will
be refreshed, the `setup-online-license.sh` marker keeps the NNV
install from re-running unnecessarily).

**MATLAB browser at localhost:8888 won't load.** Check Docker Desktop
has port 8888 free and `--shm-size=512M` was passed (insufficient
shared memory will hang the VNC server).

**`setup-online-license.sh` errors with "no licence"**. You're not yet
signed in. Switch to the first terminal, complete the browser sign-in,
wait for MATLAB to fully start, then re-run the `docker exec` step.

**Re-run the setup**: delete the marker inside the container:

```bash
docker exec -it nnv3-setup rm /home/matlab/.matlab/.nnv-online-setup-done
```

or from a fresh container with the same volume mounted, then re-run
`setup-online-license.sh`.

**Fresh start**: `docker volume rm nnv3-matlab-prefs nnv3-matlab-mw`
wipes the cached activation and NNV install. The next `docker run -it
… -browser` will require a fresh sign-in.
