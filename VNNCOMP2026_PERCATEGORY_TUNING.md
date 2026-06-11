# VNN-COMP 2026 -- Per-Category Falsification + Reachability Tuning Spec

*Design artifact produced 2026-06-11 from a 25-benchmark strategy sweep (5 parallel design agents +
synthesis), then IMPLEMENTED (PR #306 has merged). Two of its findings were independently verified against
the live code and fixed: the `reachOptionsList{1}`-overwrite bug (dist_shift run_vnncomp_instance.m:555/559,
linearize :568/572, cora :538/541, malbeware :597/601 -- the intended cheap-first method was silently
overwritten so only the wrong/exponential method ran) and the PGD dlnetwork-gate that starves the NN-manifest
categories (traffic_signs, lsnc_relu) of gradient falsification. STATUS: Sec 3.5 (overwrite fixes) and
Sec 3.1-3.3 (per-category `falsifyOpts` + threading) are implemented; the `nRand` column was tuned DOWN from
the table below after measuring ~15 ms/sample random sampling (PGD is the primary falsifier, so the random
fallback is kept light). Sec 3.4 (PGD on the NN-manifest path) and the onnxruntime SAT guard wiring remain as
follow-ups.*

---

# VNN-COMP 2026 Per-Category Tuning Spec for NNV

> **Scope.** This document is the per-category falsification + reachability tuning plan for NNV's VNN-COMP 2026 entry, derived from the 25-row strategy sweep. It is implementation-ready against `code/nnv/examples/Submission/VNN_COMP2025/run_vnncomp_instance.m` (the runner), `pgd_falsify.m`, `falsify_single`, and `load_vnncomp_network`. **Dependency:** this tuning builds on PR #306 ("needReshape-aware gradient falsification"), which added the `falsify_single -> pgd_falsify(...)` path the per-category `falsifyOpts` thread plugs into. If #306 is not present in your tree, apply/rebase onto it first. (See the intro above for current implementation status.)

---

## 1. Tuning philosophy

The scoring asymmetry drives everything: a sound `sat`/`unsat` is **+10**, an unknown is **0**, and an *incorrect* verdict is **-150**. So the policy is three-tiered by instance profile. **Falsify-first wherever SAT is plausible and the input box is cheap to attack** (acasxu, sat_relu, the control/Lyapunov nets, collins_rul, cora, safenlp, the dist_shift/linearizenn linearizers) -- gradient PGD on a low-dim box is nearly free and every validated counterexample is a free +10 that NNV missed in 2025 (354 SAT found vs ~1000 for the field). **Speed-first wherever reach is timeout-prone** (cifar100, tinyimagenet, vggnet16, and the big detection/transformer CNNs) -- there we shrink the PGD budget, cut `nRand`, and lead the reach ladder with the cheapest *sound* over-approximations (`approx-zono` -> `abs-dom` -> very-loose `relax-star-area`), never `approx-star`/`exact-star`, so the pass actually finishes. **Soundness-or-unknown everywhere**: `cp-star` is probabilistic and is demoted off every verdict-emitting path (kept only as a non-binding last-resort fallback); a SAT is emitted only after `validate_witness` (re-checked against onnxruntime on the hazardous benches), and on the four ERROR/import-blocked categories we deliberately return `unknown` rather than risk a -150. Net effect of the per-category ladders: every method that can *emit a verdict* is a sound over-approximation or a validated concrete witness, while the loosest-but-fastest sound method always leads so we trade precision for in-budget completion only where the net is too big to certify tightly.

---

## 2. Per-category configuration table

| category | sweep | falsify_first | pgd (restarts/steps/max_s) | nRand | reach_ladder | relaxFactor | rationale (1-line) |
|---|---|:--:|:--:|:--:|---|:--:|---|
| **acasxu** | unknown (45.2s) | yes | 40 / 60 / 5 | 1000 | approx-star -> exact-star | 0 | Tiny 5-in FC; falsify HARD (many SAT props score 0 today), then cheap approx->exact UNSAT. |
| **sat_relu** | SAT (4.5s) | yes | 50 / 80 / 5 | 1000 | approx-star -> exact-star | 0 | By-design SAT; max-out PGD for free CEs, sound ladder backs UNSAT. |
| **cersyve** | unknown (16.3s) | yes | 30 / 50 / 4 | 500 | approx-star -> relax-star-area | 0.8 | Small safety FC; drop probabilistic cp-star as verdict source, use sound ladder. |
| **lsnc_relu** | unknown (2.2s) | yes | 50 / 80 / 4 | 1000 | approx-star -> exact-star | 0 | Tiny [6]->[8] (xval'd 9.2e-7); PGD+exact both free. **NN-manifest path: see PGD-gate fix (Sec 3.4).** |
| **relusplitter** | unknown (12.2s) | yes | 30 / 50 / 3 | 500 | approx-star -> relax-star-area | 0.85 | Adversarial-split ReLU; tighten loose 1.0->0.85 so UNSAT actually closes; cap PGD (oval-conv heavier). |
| **dist_shift** | unknown (11.0s) | yes | 30 / 60 / 4 | 300 | approx-star -> relax-star-area | 0.7 | **Fix {1}-overwrite bug** (only exact-star runs today); drop exact-star, sound cheap->loose ladder. |
| **linearizenn** | unknown (12.1s) | yes | 30 / 60 / 4 | 300 | approx-star -> relax-star-area | 0.7 | Same {1}-overwrite bug; remove exact-star; cp-star only as non-path catch. |
| **tllverifybench** | unknown (15.4s) | yes | 40 / 80 / 5 | 400 | approx-star -> relax-star-area | 0.8 | TLL nets loose under over-approx -> reach rarely closes; invest max in PGD, lead with tighter approx-star. |
| **collins_rul** | SAT (3.3s) | yes | 30 / 60 / 4 | 200 | approx-star -> relax-star-area | 0.7 | Already SAT via sampling; PGD finds it faster/more; keep reach short. needReshape=2 preserved. |
| **nn4sys** | unknown (18.2s) | yes | 25 / 50 / 3 | 200 | approx-star -> relax-star-area | 0.8 | INVALID/all-SAT history -> soundness first: only validated PGD claims SAT; demote cp-star; mscn stays import-blocked. |
| **safenlp** | SAT (3.9s) | yes | 30 / 60 / 4 | 500 | approx-zono -> abs-dom -> approx-star | 0 | Small MLP, SAT win; push PGD; slow_cats already drops exact-star. |
| **malbeware** | UNSAT (9.4s) | yes | 20 / 40 / 3 | 200 | approx-star -> exact-star | 0 | Sound UNSAT already; **fix {1}-overwrite**, lead approx-star, exact only as small-net fallback; PGD moderate. |
| **ml4acopf** | unknown (15.3s) | yes | 25 / 50 / 3 | 300 | cp-star *(non-binding)* | 0 | No sound NNV net (nnvnet=""); only SAT is scorable -> push PGD on dlnetwork; cp-star never asserts a verdict. |
| **cora** | SAT (3.7s) | yes | 30 / 60 / 4 | 500 | approx-zono -> abs-dom -> relax-star-area -> approx-star | 0.9 | SAT win; **fix 0.9->0.7 {1}-overwrite**; clean cheap->precise ladder, all sound. |
| **metaroom** | UNSAT (15.6s) | yes | 15 / 30 / 2 | 100 | relax-star-area -> approx-star | 0.85 | Sound UNSAT on CNN near budget edge; keep PGD SMALL, lead loose-but-sound relax. needReshape=2 preserved. |
| **cifar100** | TIMEOUT (120s) | yes | 5 / 15 / 2 | 50 | approx-zono -> abs-dom -> relax-star-area | 0.95 | Big CNN timeout; cheapest sound ladder only, very loose relax; tiny PGD, cut nRand. |
| **tinyimagenet** | TIMEOUT (120s) | yes | 5 / 15 / 2 | 100 | approx-zono -> abs-dom -> relax-star-area | 0.95 | 200-class deep CNN; same cheap ladder; drop nRand 500->100 (sampling low-yield). |
| **vggnet16** | TIMEOUT (120s) | yes | 4 / 12 / 1.5 | 50 | approx-zono -> abs-dom | 0.9 | Heaviest net; restrict to two cheapest sound methods; minimal PGD/nRand. |
| **traffic_signs** | unknown (9.2s) | **no->yes (needs fix)** | 20 / 40 / 5 | 300 | approx-star -> relax-star-area | 0.5 | **NN-manifest path: PGD never fires today (Sec 3.4 fix).** Until then push nRand hard on small 30x30x3 net. |
| **vit** | unknown (26.2s) | yes | 10 / 25 / 3 | 100 | relax-star-area -> approx-star | 0.8 | Attention/softmax is the slow part; lead looser relax (0.5->0.8) for a verdict before tight approx. |
| **yolo** | unknown (18.8s) | yes | 20 / 40 / 3 | 200 | approx-zono -> abs-dom -> relax-star-area | 0.9 | Conv detector; cp-star needs GPU/unsound -> drop for verdicts; cheap sound ladder, capped PGD. |
| **cctsdb_yolo** | **ERROR** (import WIP) | yes | 30 / 50 / 5 | 300 | *(none -- blocked)* | 0 | Branch raises error() before any net built. **IMPORT FIX NEEDED**; until then PGD+sampling only, else unknown. |
| **cgan** | **ERROR** (ReshapeLayer) | yes | 20 / 40 / 3 | 200 | relax-star-area -> approx-star | 0.8 | ReshapeLayer import fails. **IMPORT FIX NEEDED**; ladder is ready once import works. |
| **collins_aerospace** | **ERROR** (Unsupported Layer) | yes | 20 / 40 / 3 | 200 | approx-star -> relax-star-area | 0.5 | Unsupported layer + documented INVALID-SAT history. **IMPORT FIX NEEDED**; only validated witnesses, needReshape=1 locked. |
| **soundnessbench** | **ERROR** (ReshapeLayer) | yes | 30 / 50 / 4 | 200 | approx-star -> relax-star-area | 0.5 | Purpose-built to catch unsound verifiers (max -150 hazard). **IMPORT FIX NEEDED**; SAT only via validate_witness+onnxruntime; never UNSAT from cp-star. |

---

## 3. Implementation plan

The design goal: thread a **per-category falsification options struct** from `load_vnncomp_network` through to `pgd_falsify`, **without touching the `reachOptionsList` logic at all** (the reach ladders are already produced by the existing per-category branches + Phase-1.5 + slow_cats post-processing). We add one new return value, one new argument, and a small table. The `relaxFactor`/`reach_ladder` columns in the table above are documentation of the *intended* reach behavior; most are already emitted by the existing branches and the handful of `{1}`-overwrite bugs are fixed in-place in those branches (Sec 3.5). The genuinely new machinery is the falsification thread.

### 3.1 Add a per-category falsify-opts table in `load_vnncomp_network`

Add an 8th return value `falsifyOpts` (a struct). Keep all existing returns and positions unchanged so no other caller breaks.

```matlab
function [net,nnvnet,needReshape,reachOptionsList,inputSize,inputFormat,nRand,falsifyOpts] = ...
        load_vnncomp_network(category, onnx, vnnlib)
```

At the **top** of the function (right after `nRand = 100;`, line ~419), seed a sound default, then override per category from a single lookup table. This is additive -- it does not interact with any `reachOptions` assignment.

```matlab
    % ---- Per-category PGD/falsification budget (VNNCOMP2026 tuning) ----
    % Defaults mirror the pre-tuning hardcoded struct('seed',0,'max_time',5).
    falsifyOpts = struct('n_restarts',20, 'n_steps',40, 'lr',0.1, ...
                         'fgsm',true, 'max_time',5, 'seed',0, 'nRand',nRand);

    % rows: {match, n_restarts, n_steps, max_time_s, nRand}
    ftab = {
        "acasxu",        40,60,5,  1000
        "sat_relu",      50,80,5,  1000
        "cersyve",       30,50,4,  500
        "lsnc_relu",     50,80,4,  1000
        "relusplitter",  30,50,3,  500
        "dist_shift",    30,60,4,  300
        "linearizenn",   30,60,4,  300
        "tllverifybench",40,80,5,  400
        "collins_rul",   30,60,4,  200
        "nn4sys",        25,50,3,  200
        "safenlp",       30,60,4,  500
        "malbeware",     20,40,3,  200
        "ml4acopf",      25,50,3,  300
        "cora",          30,60,4,  500
        "metaroom",      15,30,2,  100
        "cifar100",       5,15,2,  50
        "tinyimagenet",   5,15,2,  100
        "vggnet",         4,12,1.5,50      % matches vggnet16
        "traffic_signs", 20,40,5,  300
        "vit",           10,25,3,  100
        "yolo",          20,40,3,  200
        "cctsdb_yolo",   30,50,5,  300
        "cgan",          20,40,3,  200
        "collins_aerospace",20,40,3,200
        "soundnessbench",30,50,4,  200
    };
    for r = 1:size(ftab,1)
        if contains(category, ftab{r,1})
            falsifyOpts.n_restarts = ftab{r,2};
            falsifyOpts.n_steps    = ftab{r,3};
            falsifyOpts.max_time   = ftab{r,4};
            nRand                  = ftab{r,5};
            falsifyOpts.nRand      = nRand;
            break
        end
    end
```

> **Match-order caveat:** `contains` is a substring test. `"cctsdb_yolo"` also contains `"yolo"`, and `"vggnet"` is used (not `"vggnet16"`) to match the existing `slow_cats` convention. Put the **most specific keys first** (`cctsdb_yolo` before `yolo`, `collins_aerospace` and `collins_rul` are disjoint substrings so order-independent) and `break` on first hit. This mirrors how the existing branch `if/elseif` chain already disambiguates these categories.

### 3.2 Thread `falsifyOpts` through the runner header

At line 16, capture the new return:

```matlab
[net, nnvnet, needReshape, reachOptionsList, inputSize, inputFormat, nRand, falsifyOpts] = ...
    load_vnncomp_network(category, onnx, vnnlib);
```

`nRand` is still returned and used unchanged; `falsifyOpts.nRand` is set consistently with it for self-documentation.

### 3.3 Extend `falsify_single` to accept opts, and pass through to `pgd_falsify`

Add a trailing optional `opts` argument so the three existing call sites can be updated minimally and any other caller of `falsify_single` keeps working with the old default.

```matlab
function counterEx = falsify_single(net, lb, ub, inputSize, nRand, Hs, needReshape, inputFormat, opts)
    counterEx = nan;
    if nargin < 9 || isempty(opts)
        opts = struct('seed',0, 'max_time',5);   % preserve current default
    end
    if ~isfield(opts,'seed'), opts.seed = 0; end
    if isa(net, 'dlnetwork')
        try
            [cex, found] = pgd_falsify(net, lb, ub, Hs, inputSize, inputFormat, needReshape, opts);
            if found && validate_witness(net, cex{1}, lb, ub, Hs, inputSize, inputFormat, needReshape)
                counterEx = cex; return;
            end
        catch
            % fall through to random sampling
        end
    end
    ... % create_random_examples / sampling loop UNCHANGED
end
```

This is a **drop-in** change: `pgd_falsify` already reads `n_restarts`, `n_steps`, `lr`, `fgsm`, `max_time`, and `seed` from its `opts` struct via `getfielddef` (pgd_falsify.m lines 36-42), so no edit to `pgd_falsify` is required for the budget knobs.

Update the **three call sites** (lines 62, 65, 72) to forward `falsifyOpts`:

```matlab
% line 62
counterEx = falsify_single(net, lb, ub, inputSize, nRand, prop{1}.Hg, needReshape, inputFormat, falsifyOpts);
% line 65
counterEx = falsify_single(net, lb{spc}, ub{spc}, inputSize, nRand, prop{spc}.Hg, needReshape, inputFormat, falsifyOpts);
% line 72
counterEx = falsify_single(net, lb{arr}, ub{arr}, inputSize, nRand, prop{1}.Hg, needReshape, inputFormat, falsifyOpts);
```

Because `nRand` is already the per-category value returned by `load_vnncomp_network`, the random-sampling fallback automatically picks up the tuned counts; no extra wiring needed there.

### 3.4 Extend PGD to the NN-manifest path (the single highest-EV falsification fix)

Today `pgd_falsify` returns immediately for non-`dlnetwork` nets (`if ~isa(net,'dlnetwork'); return; end`, line 33), so **traffic_signs** and **lsnc_relu** (loaded as NNV `NN` via `load_manifest_net`) get *no gradient falsification* -- only random sampling. For traffic_signs (2025: largely SAT) this is the biggest single miss. Two options, lowest-risk first:

- **Option A (recommended, low-risk):** in `load_vnncomp_network`, for the manifest categories keep the MATLAB `dlnetwork` available (the importer has it pre-conversion) and pass *that* `net` to `falsify_single` while still using the `NN` `nnvnet` for reach. Then PGD fires on the dlnetwork and the validated witness is written back through the existing `needReshape`-aware mapping. No change to `pgd_falsify`.
- **Option B (fallback):** add a finite-difference / `NN.evaluate`-based FGSM branch inside `pgd_falsify` for `isa(net,'NN')` (NN.evaluate is differentiable-free, so estimate the margin gradient by central differences on the active half-space). Higher cost, only if the dlnetwork is genuinely unavailable.

Until A or B lands, traffic_signs stays `falsify_first=no` in practice and relies on `nRand=300` random sampling -- which is why its row carries the **"no->yes (needs fix)"** marker.

### 3.5 Fix the `reachOptionsList{1}`-overwrite bugs (in-place, no architecture change)

Four branches write two methods to the same index `{1}`, so only the second (often `exact-star`) ever runs. Fix each by appending to `{2}` (or `{end+1}`) and ordering cheap->precise per the table. These are localized edits inside the existing per-category branches and do **not** touch Phase-1.5 / slow_cats:

- **dist_shift** & **linearizenn**: stop writing `exact-star` to `{1}`; set `{1}=approx-star`, `{2}=relax-star-area (relaxFactor 0.7)`; **drop exact-star** entirely (blowup risk).
- **malbeware**: `{1}=approx-star`, `{2}=exact-star` (exact only as small-net precise fallback).
- **cora** (non-set branch): replace the `0.9->0.7` double-write with a single `relax-star-area, relaxFactor=0.9`; the slow_cats prepend already supplies `approx-zono`+`abs-dom` ahead of it.

These four fixes are independently valid and should be cherry-picked even if the falsification thread slips.

### 3.6 What explicitly does **not** change

- `reachOptionsList` construction, Phase-1.5 cp-star prepend (lines 870-879), and the `slow_cats` exact-star drop (lines 887-898) are **untouched** -- the reach ladders in the table are produced by code that already exists.
- The reach loop (`while ~isempty(reachOptionsList)`, lines 112+) and its sound "error->unknown" mapping are untouched.
- `validate_witness` gating is untouched; every SAT still passes through it.

---

## 4. Highest-EV changes first, and the ERROR-category triage

### 4.1 The 3-4 changes that move the score most (do these first)

1. **Land the per-category PGD budget + `falsifyOpts` thread (Sec 3.1-3.3).** This is the broad lever: it converts the already-SAT categories (**sat_relu, collins_rul, cora, safenlp**) from "sometimes found by random sampling" to "reliably found by gradient PGD," and--critically--turns **acasxu** from `unknown` into a stream of `sat` verdicts. ACAS-Xu has many SAT properties currently scored 0 purely because the 45s run only over-approximates; high-restart PGD on a 5-D box is essentially free. **This single change plausibly recovers the most points.**
2. **Extend PGD to the NN-manifest path for traffic_signs (Sec 3.4, Option A).** Highest *per-category* delta: traffic_signs was largely SAT in 2025 and currently gets **zero** gradient falsification because PGD is dlnetwork-gated. Re-enabling it (plus the dlnetwork being available pre-conversion) is a targeted, high-yield fix. lsnc_relu benefits from the same fix.
3. **Fix the four `{1}`-overwrite reach bugs (Sec 3.5).** Cheap, surgical, and directly converts `dist_shift`/`linearizenn` from "exact-star stalls to unknown" into "approx-star closes UNSAT in budget," and stops `cora`/`malbeware` from running an unintended method. Pure correctness wins with near-zero risk.
4. **Tighten the loose-relax verdict-killers: relusplitter 1.0->0.85, tllverifybench reorder approx-star-first, vit relax-first 0.5->0.8.** These three are `unknown` today specifically because the leading reach method is too loose to certify (relusplitter at relaxFactor 1.0) or in the wrong order (tll/vit lead with the loose method). Re-ordering/tightening lets genuine UNSAT proofs actually close without new machinery.

### 4.2 ERROR categories: importer-fix-needed vs. safe-as-unknown

| category | error | verdict | action |
|---|---|---|---|
| **cctsdb_yolo** | branch hard-`error()`s (import WIP) before any net is built | **IMPORTER FIX NEEDED** | ONNX forward-prop fails on this detection CNN; no reach possible until import lands. Until then: PGD+sampling on the dlnetwork only, else `unknown`. Reach ladder is empty by design. |
| **cgan** | ReshapeLayer / ONNXParams can't lower | **IMPORTER FIX NEEDED** | Patch the ReshapeLayer import; the sound ladder (`relax-star-area@0.8 -> approx-star`) is already correct and ready to enable once import works. |
| **collins_aerospace** | Unsupported Layer at conversion | **IMPORTER FIX NEEDED (high-caution)** | Has documented historical **INVALID-SAT** instances (the -150 trap). Fix the layer import, keep `needReshape=1` locked (2 and 0 returned invalid SAT), and accept witnesses *only* through validate_witness. Prefer `unknown` over any unvalidated verdict. |
| **soundnessbench** | ReshapeLayer import failure | **SAFE TO LEAVE AS UNKNOWN until import is sound** | Purpose-built to catch unsound verifiers -- the maximum -150 hazard. Do not rush the import. When fixed, emit SAT only via validate_witness + onnxruntime double-check, and **never** report UNSAT from cp-star or any method a sound star can't certify. When in doubt, `unknown`. |

**Triage summary:** all four need an importer fix to enable *reach*, but their EV differs. **cctsdb_yolo and cgan** are pure "blocked, no soundness risk" -- fixing the importer is upside-only and their ladders are pre-staged. **collins_aerospace and soundnessbench** are the danger pair: importer work there must be paired with strict validate_witness/onnxruntime gating, and **soundnessbench is the one category where staying `unknown` is the correct, safe outcome** until the import is provably sound -- a wrong verdict there costs -150 and there is no falsification shortcut worth that risk. In all four, gradient PGD on the imported dlnetwork plus random sampling is the only currently-safe points source.

---

## 5. Relevant files (absolute paths)

- Runner / `load_vnncomp_network` / `falsify_single` / call sites: `C:\Users\taylo\Dropbox\Research\talks\vu-isis-2025-11-15\vsc_matlab_mcp\nnv\code\nnv\examples\Submission\VNN_COMP2025\run_vnncomp_instance.m`
- PGD engine (already opts-driven; no change for budget knobs): `C:\Users\taylo\Dropbox\Research\talks\vu-isis-2025-11-15\vsc_matlab_mcp\nnv\code\nnv\examples\Submission\VNN_COMP2025\pgd_falsify.m`
- Witness validator (SAT gate): `C:\Users\taylo\Dropbox\Research\talks\vu-isis-2025-11-15\vsc_matlab_mcp\nnv\code\nnv\examples\Submission\VNN_COMP2025\validate_witness.m`
- PGD test (extend with per-category budget cases): `C:\Users\taylo\Dropbox\Research\talks\vu-isis-2025-11-15\vsc_matlab_mcp\nnv\code\nnv\tests\vnncomp25\test_pgd_falsify.m`

> **Gate status: PR #306 is pending merge (CI + conditional auto-merge armed). Once it lands, rebase onto `main` (which will then contain the `falsify_single -> pgd_falsify(opts)` path) and apply Sec 3.1-3.5 on top.**