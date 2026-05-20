# NNV 3.0 — Reviewer helper scripts

Small cross-platform utilities to make `docker build` setup easier
for ATVA 2026 artifact reviewers.

## `find-matlab-license.ps1` / `find-matlab-license.sh`

Locate the `port@host` value required by

```bash
docker build -t nnv3.0 --build-arg LICENSE_SERVER=<port>@<host> .
```

The scripts search every standard MATLAB license source on the host:

- `MLM_LICENSE_FILE` and `LM_LICENSE_FILE` environment variables
  (process, user, and machine scopes on Windows).
- `matlab` on `PATH`, resolved to its `MATLABROOT`.
- Windows registry (`HKLM`, `HKLM\WOW6432Node`, `HKCU` under
  `SOFTWARE\MathWorks\MATLAB`).
- Per-release install and preferences directories:
  - Windows: `C:\Program Files\MATLAB\R*`, `C:\Program Files (x86)\MATLAB\R*`,
    `C:\MATLAB\R*`, `%PROGRAMDATA%\MathWorks\MATLAB\R*`,
    `%APPDATA%\MathWorks\MATLAB\R*`, `%LOCALAPPDATA%\MathWorks\MATLAB\R*`.
  - macOS: `/Applications/MATLAB_R*.app`,
    `~/Library/Application Support/MathWorks/MATLAB/R*`.
  - Linux: `/usr/local/MATLAB/R*`, `/opt/matlab/R*`, `/opt/MATLAB/R*`, `~/.matlab/R*`.

For each MATLAB root the scripts parse:

- FlexLM-style `*.lic` files in `licenses/` for `SERVER <host> <hostid> <port>`
  lines, emitting `port@host` for each.
- Modern `license_info.xml` activation manifests, extracting `<licmode>` and any
  host/port fields (`lichost`/`licport`, `serverhost`/`serverport`, or a
  combined `<server>port@host</server>`).

The scripts print every search location they touched so users can verify
the search was thorough, then either:

- emit a table of `port@host` values discovered (success), or
- explain the declared license mode and what to do next — for example,
  an `onlinelicensing` (signed-in individual) license has no `port@host`
  and is incompatible with the build's `--build-arg LICENSE_SERVER=`
  flag.

### Usage

Windows (PowerShell):
```powershell
powershell -ExecutionPolicy Bypass -File code\nnv\examples\NNV3.0\utils\find-matlab-license.ps1
```

Linux / macOS / WSL / Git Bash:
```bash
bash code/nnv/examples/NNV3.0/utils/find-matlab-license.sh
```

Exit code is `0` if at least one `port@host` was found, `1` otherwise.
