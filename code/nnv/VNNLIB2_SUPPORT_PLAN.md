# VNN-LIB 2.0 Support Plan for NNV (VNN-COMP 2026 Extended Track)

Status: PLANNING / research deliverable. No production code changed by this document.
Date: 2026-06-12.
Author: research+planning pass (Claude).

> ## IMPLEMENTATION STATUS (2026-06-12) — Phase 0 + 1 + 3a SHIPPED
> The recommended scope (§5.1) is implemented and tested:
> - **`engine/utils/load_vnnlib2.m`** — single-network 2.0 parser. STREAMING design
>   (one S-expression statement at a time) so the 96 MB / 124 MB image benchmarks
>   parse without a whole-file token tree; multi-network/multimodal files gate from
>   the header alone (smart_turn 124 MB gated in 4 ms). Sound-or-unknown throughout.
> - **`engine/utils/load_vnnlib.m`** — `(vnnlib-version <2.0>)` sniff dispatches 2.0
>   files to `load_vnnlib2`; the 1.0 path is byte-for-byte unchanged.
> - **`examples/Submission/VNN_COMP2025/run_vnncomp_instance.m`** — `property.unsupported`
>   → emit `unknown` (Phase 3a) rather than risk a −150 wrong verdict.
> - **`tests/utils/test_load_vnnlib2.m`** — 13 tests (forms, multi-dim row-major flatten,
>   ==, var-var, OR/conjunct distribution, dispatch, and the four gates).
> - **Validation:** acasxu prop_1 parsed via 2.0 == via 1.0 EXACTLY (lb/ub/G/g diff 0);
>   the single-net 2.0 benchmarks parse to sound boxes; multi-net/nonlinear/multimodal
>   /mixed gate to `unknown`.
> - **Adversarial review (19-agent workflow):** found + fixed one real soundness bug —
>   a standalone output conjunct adjacent to an `(or ...)` must be distributed into
>   every disjunct (one prop cell), else the runner's `length(prop)==1` path silently
>   ignores it. Strict `<`/`>` are treated as `<=`/`>=` to match the battle-tested 1.0
>   parser (sound for the UNSAT proof; the boundary SAT case is within VNN-COMP's 1e-4
>   constraint tolerance and guarded by the onnxruntime witness replay).
> - **Deferred (still per §5.2):** multi-network `equal-to`/`isomorphic-to` (Phase 2),
>   `smart_turn` multimodal input set (Phase 3b), pointing the sweep driver at the
>   `2.0/` instances.csv, and PGD-on-the-box for nonlinear-gated instances.
Discipline: any 2.0 path proposed here must be **sound-or-unknown** (the −150 rule:
never emit `sat`/`unsat` we cannot defend; on any doubt return `unknown`).

Sources are cited inline as `[Sn]` and collected at the bottom. Local file evidence
is cited as `path:line`. Benchmark evidence is cited by `zcat` of the actual
`*.vnnlib.gz` under `../vnncomp2026_benchmarks/benchmarks/<bench>/2.0/`.

---

## 0. TL;DR / headline

- **34 benchmark folders ship a `2.0/` directory.** Of those, **4 are 1.0-syntax
  files mislabeled into a `2.0/` folder** (cersyve, cora_2024, malbeware,
  safenlp_2024 — they still use `(declare-const X_0 Real)` and have *no*
  `(vnnlib-version <2.0>)` header), so the existing 1.0 parser already reads them.
  **28 are genuine single-network 2.0** (one `declare-network`, box input + linear
  output — i.e. 1.0 semantics re-expressed in 2.0 syntax). **2 are genuine
  multi-network 2.0** (`monotonic_acasxu_2026`, `isomorphic_acasxu_2026` — two
  networks `f`,`g` with `equal-to`/`isomorphic-to` + cross-network constraints).
- **Tractable by June 30 with a parser extension + existing NNV reach: the 28
  single-network 2.0 benchmarks** (plus the 4 mislabeled-1.0). Two of these are
  extended-only and *only* reachable with a 2.0 parser:
  `adaptive_cruise_control_non_linear_2026` (single-net but **nonlinear** property →
  parse yes, sound-verify mostly no) and `smart_turn_multimodal_2026` (single-net
  but **2 input tensors** / multimodal → needs multi-input set construction).
- **Research-grade / likely out of scope for June 30: the 2 multi-network
  benchmarks.** They need a product/difference construction over two networks plus
  a runner that loads a *list* of `(name, onnx)` pairs. A sound construction exists
  for `equal-to` (Section 3) but the runner plumbing is a substantial change.
- **Recommended path: extend the MATLAB parser** (`load_vnnlib.m`) to a new
  `load_vnnlib2.m` dispatched on the `(vnnlib-version <2.0>)` header, scoped first to
  single-network linear properties; **keep the `vnnlib` PyPI package (dlshriver) as a
  cross-check oracle, not the primary path.** Multi-network is a separate, optional
  phase.

---

## 1. VNN-LIB 2.0 grammar summary (what matters for NNV)

VNN-LIB 2.0 ("Rigorous Foundations for Neural Network Verification", Roy, Antony,
Gimelli, Daggitt, CAV 2026) was released **2025-12-15** as `standard.pdf` v2.0.0
[S1][S2][S6]. The grammar is given in LBNF in `syntax/grammar.cf` [S3]. A 2.0 query
is: a `(vnnlib-version <2.0>)` header, then a list of typed `declare-network` blocks,
then a list of `(assert <bool-expr>)` statements [S3][S6].

### 1.1 Declarations

- **Version:** `(vnnlib-version <2.0>)` — first non-comment line; this is the
  dispatch key for "is this a 2.0 file?" [S3]. Confirmed present in every genuine
  2.0 file, e.g. `monotonic_acasxu_2026/2.0/vnnlib/instance_0.vnnlib.gz`.
- **Network:** `(declare-network <name> [<equivalence>] (declare-input ...)+
  (declare-output ...)+ (declare-hidden ...)* )` [S3][S5]. A network is a *named*
  scope; inputs/outputs are *namespaced by network name* in genuine multi-network
  files (`X_f`, `Y_g`).
- **Input/Output:** `(declare-input <name> <elementType> <shape>)`,
  `(declare-output <name> <elementType> <shape>)` [S5]. `<elementType>` is a type
  from the network theory (`float32`, `float64`, `int32`, ...) **or** the
  distinguished type `real` for "real-valued queries" [S6]. `<shape>` is a bracketed
  dim list: `[5]`, `[1,5]`, `[1,3,32,32]`, `[]` for scalar.
- **Hidden (new in 2.0):** `(declare-hidden <name> <elementType> <shape>
  <onnxNodeName>)` — lets a property constrain an *intermediate* ONNX node by its
  graph name [S4][S5]. Use case: encoder–decoder bottlenecks, attention. **No 2026
  benchmark uses this** (grep over all `2.0/` files: zero `declare-hidden`), so it is
  out of scope for NNV 2026.
- **`real` vs `float32` semantics:** declaring inputs as `real` is "explicit
  permission for the solver to interpret ONNX networks as real-valued functions and
  acknowledges that the result may be unsound due to numerical imprecision" [S6]. The
  2026 single-net benchmarks use `float32`; the multi-net acasxu ones use `real`.

### 1.2 Network equivalence relations (the genuinely new semantics)

- **`(equal-to <other>)`:** "each `declare-network` declaration is intended to map to
  the same network implementation" [S5]; formally the two networks "must actually be
  the *same* model file" and "both the element types and shapes of the input and
  output declarations for the two networks have to match" [S6]. Evidence:
  `monotonic_acasxu_2026/2.0/instances.csv` maps both `f` and `g` to the *same* onnx
  (`('f','onnx/ACASXU_..._2_2_...onnx'), ('g','onnx/ACASXU_..._2_2_...onnx')`).
- **`(isomorphic-to <other>)`:** "the same graph structure but with different weights"
  [S5]; "the shapes, but not the element types ... have to match" [S6]. Evidence:
  `isomorphic_acasxu_2026/2.0/instances.csv` maps `f`→`onnx/original/...onnx`,
  `g`→`onnx/perturbed/..._perturbed_0.onnx`.
- **No transitive chaining:** "an `equal-to` or `isomorphic-to` declaration cannot
  reference another network declaration that also contains an equivalence
  statement" [S6] → at most 2 networks effectively coupled in 2026.

### 1.3 Variable references and indexing

- **Multi-dimensional, 0-based:** `X[0,2]` = "row 0, column 2" [S6]. **Partial
  indexing is not allowed** — the number of indices must equal the tensor rank [S6].
  Confirmed: `adaptive_cruise.../instance_0.vnnlib.gz` uses `X[0,0]`, `X[0,1]`;
  `acasxu_2023/2.0` uses `X[0,0,0,0]`; `smart_turn` uses `X1[0,0,0]`.
- **Network-namespaced names:** in multi-net files variables are `X_f[0]`, `Y_g[3]`;
  in single-net files they are just `X`, `Y` (or `X1`,`X2` for multimodal). The
  trailing-`_f` is part of the *variable name*, not an index.

### 1.4 Operators

- **Comparisons (binary):** `<`, `<=`, `>`, `>=`, `==`, `!=` [S3][S5][S6]. (1.0 NNV
  supports only `<`,`<=`,`>`,`>=` — `load_vnnlib.m:16` "only >= or <= are
  supported".) `==` appears in the acasxu multi-net input equalities
  (`(== X_f[3] 0.227...)`) and `!=` in the nonlinear ACC property.
- **Boolean (n-ary):** `and`, `or` [S3][S6]. (1.0 NNV supports limited `and`/`or`
  shapes; see `load_vnnlib.m:13-15`.)
- **Arithmetic (n-ary, n≥2, plus unary negation):** `(+ a b ...)`, `(- a b ...)`,
  `(- a)`, `(* a b ...)` [S5][S6]. **Nonlinear products of variables are
  grammatically allowed** — `(* X[0,1] X[0,1])` appears in
  `adaptive_cruise.../instance_0.vnnlib.gz`. The spec delegates the actual arithmetic
  to the network theory and does **not** forbid variable·variable products [S6].

### 1.5 Satisfiability / sat–unsat / witness convention (vs 1.0)

- **Same polarity as 1.0.** A 2.0 query is **satisfiable iff there exists an input
  assignment under which all assertions evaluate true** [S6]: *"A query is
  satisfiable if such an assignment exists, and unsatisfiable otherwise."* The
  benchmark assertions encode the **negation of the safety property** (the unsafe /
  counterexample region), exactly as in 1.0. So for the verifier:
  **`sat` = counterexample found = property disproved; `unsat` = property proved;
  `unknown`/timeout otherwise** [S7]. This matches `run_vnncomp_instance.m`'s
  existing convention (`status 0=sat, 1=unsat, 2=unknown`,
  `run_vnncomp_instance.m:361-377`) — **no change to the result encoding is needed.**
- **Witness format:** a satisfying assignment is a mapping from each declared input
  variable to a concrete tensor [S6]. For single-network queries this is identical to
  1.0's `(X_i v) ... (Y_j v)` text dump (`run_vnncomp_instance.m:1059-1092`
  `write_counterexample`). For **multi-network** queries the spec does *not* pin a
  standardized witness text format [S6] — the witness must give *both* `X_f` and
  `X_g`; VNN-COMP's accompanying tooling decides the exact serialization. **This is an
  open risk for the multi-net path** (Section 3/5).

---

## 2. The 2026 2.0 benchmark catalog

Method: for each `<bench>/2.0/vnnlib/<first>.vnnlib.gz`, `zcat | grep` for
`vnnlib-version`, `declare-network`, `declare-input`, `declare-output`, `equal-to`,
`isomorphic-to`, `(or `, `(* `, `==`/`!=`, and counted `declare-input`. The `1.0/`
column is whether a sibling `1.0/` dir exists (NNV's existing path already covers it).

### 2.1 Classification table

| Benchmark | Has 1.0? | 2.0 syntax? | #nets | equiv? | cross-net? | #inputs | property shape | NNV tractability (2.0 path) |
|---|---|---|---|---|---|---|---|---|
| acasxu_2023 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| adaptive_cruise_control_non_linear_2026 | **no (2.0-only)** | **2.0** | 1 | – | – | 1 | box+**nonlinear** in, deeply nested or/and/`!=`/`*` out | **Parse yes; sound-verify mostly NO (nonlinear)** |
| cctsdb_yolo_2023 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | Single-net parse OK (net unsupported already) |
| cersyve | yes | **1.0 (mislabeled)** | (1.0 consts) | – | – | 1 | 1.0 file | **Existing 1.0 parser** |
| cgan2026 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| cgan_2023 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| challenging_certified_training_2026 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| cifar100_2024 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| collins_aerospace_benchmark | yes | **2.0** | 1 | – | – | 1 | box in, **`or`** out | **Single-net: tractable** (disjunctive out) |
| collins_rul_cnn_2022 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| cora_2024 | yes | **1.0 (mislabeled)** | (1.0 consts) | – | – | 1 | 1.0 file | **Existing 1.0 parser** |
| dist_shift_2023 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| isomorphic_acasxu_2026 | **no (2.0-only)** | **2.0** | **2** | **isomorphic-to** | **yes (in==, out ±ε)** | 1 ea | ε-equivalence of two nets | **Multi-net: research-grade** |
| linearizenn_2024 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| lsnc_relu | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** (manifest net) |
| malbeware | yes | **1.0 (mislabeled)** | (1.0 consts) | – | – | 1 | 1.0 file | **Existing 1.0 parser** |
| metaroom_2023 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| ml4acopf_2024 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** (cp-star already) |
| monotonic_acasxu_2026 | **no (2.0-only)** | **2.0** | **2** | **equal-to** | **yes (in≥, out<)** | 1 ea | monotonicity across 2 evals of *same* net | **Multi-net: research-grade (but `equal-to` ⇒ feasible, Section 3)** |
| nn4sys | yes | **2.0** | 1 | – | – | 1 | box in, **`or`** out | **Single-net: tractable** |
| relusplitter | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| relusplitter_2026 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| safenlp_2024 | yes | **1.0 (mislabeled)** | (1.0 consts) | – | – | 1 | 1.0 file | **Existing 1.0 parser** |
| sat_relu | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| smart_turn_multimodal_2026 | **no (2.0-only)** | **2.0** | 1 | – | – | **2 (X1,X2)** | box in (per modality), `(> Y[0,0] 0.5)` out | **Parse yes; needs multi-input-tensor set** |
| soundnessbench | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| soundnessbench_2026 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| test | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** (toy) |
| tinyimagenet_2024 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| tllverifybench_2023 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| traffic_signs_recognition_2023 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** (manifest net) |
| vggnet16_2022 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| vit_2023 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |
| yolo_2023 | yes | **2.0** | 1 | – | – | 1 | box in, linear out | **Single-net: tractable** |

### 2.2 Headline counts

- **34** folders with a `2.0/` directory.
- **4** are **1.0 syntax mislabeled** into `2.0/` (cersyve, cora_2024, malbeware,
  safenlp_2024) → the existing `load_vnnlib.m` already parses these; "supporting 2.0"
  is a *no-op* for them (only a folder-path/dispatch concern).
- **28** are **genuine single-network 2.0** (one `declare-network`, box input,
  linear or disjunctive output) → **tractable** with a 2.0 *parser extension* + the
  existing NNV single-network reach pipeline. (Network-import support is orthogonal
  and unchanged — whatever `load_vnncomp_network.m` does today.)
- **2** are **genuine multi-network 2.0** (monotonic, isomorphic acasxu) →
  **research-grade**; need product/difference reasoning + new runner plumbing.
- **The 4 extended-only (2.0-only) benchmarks** break down as:
  - `monotonic_acasxu_2026` — **multi-net, `equal-to`**, cross-net in/out.
  - `isomorphic_acasxu_2026` — **multi-net, `isomorphic-to`**, cross-net in/out (±ε).
  - `adaptive_cruise_control_non_linear_2026` — **single-net but nonlinear** property
    (variable·variable products, `!=`, deep nesting). Parser can read it; NNV's linear
    Star/ImageStar reach cannot soundly verify the nonlinear assertion in general.
  - `smart_turn_multimodal_2026` — **single-net but 2 input tensors** (multimodal:
    `X1 [1,80,800]`, `X2 [1,3,32,112,112]`); output `(> Y[0,0] 0.5)` is trivial. The
    blocker is constructing one input set over two tensors, not the property.

### 2.3 Concrete 2.0 examples (verbatim, abbreviated)

`monotonic_acasxu_2026/2.0/vnnlib/instance_0.vnnlib.gz` (multi-net, `equal-to`):
```
(vnnlib-version <2.0>)
(declare-network f (declare-input X_f real [5]) (declare-output Y_f real [5]))
(declare-network g (equal-to f) (declare-input X_g real [5]) (declare-output Y_g real [5]))
(assert (and (<= X_f[0] 0.667...) (>= X_f[0] -0.162...)))   ; box on f
(assert (== X_f[3] 0.227...))                               ; equality (point) on f
(assert (and (>= X_f[0] X_g[0]) (>= X_g[0] -0.162...)))     ; cross-net INPUT relation
(assert (== X_f[1] X_g[1])) ...                             ; X_g pinned to X_f elsewhere
(assert (< Y_f[3] Y_g[3]))                                  ; cross-net OUTPUT property
```

`isomorphic_acasxu_2026/2.0/vnnlib/instance_0.vnnlib.gz` (multi-net, `isomorphic-to`,
output ±ε band — and note the band is `> Y_f+0.05` AND `< Y_f-0.05`, an *empty*
region, i.e. the unsafe set is unsatisfiable by construction unless interpreted as OR):
```
(declare-network g (isomorphic-to f) ...)
(assert (== X_f[0] X_g[0]) ...)                             ; inputs identified
(assert (and (> Y_g[0] (+ Y_f[0] 0.05)) (< Y_g[0] (- Y_f[0] 0.05))))   ; ε-equivalence
```

`adaptive_cruise_control_non_linear_2026` (single-net, **nonlinear**):
```
(assert (and (> X[0,0] 0.0) (>= (* X[0,0] 200.0) (* X[0,1] X[0,1]))))   ; nonlinear input region
(assert (or (or (< Y[0,0] -100.001) ...) (and ... (* (+ X[0,0] (* X[0,1] 0.1)) 200.0) ...)))
```

`smart_turn_multimodal_2026` (single-net, **2 inputs**):
```
(declare-input X1 real [1, 80, 800])
(declare-input X2 real [1, 3, 32, 112, 112])
(declare-output Y real [1, 1])
(assert (and (>= X1[0,0,0] 0.894...) (<= X1[0,0,0] 0.994...)))   ; per-element box on X1
... (asserts on X2) ...
(assert (> Y[0,0] 0.5))                                          ; trivial linear output
```

---

## 3. Verification feasibility for multi-network properties

The interesting question: can NNV's **single-network reachability** decide a
cross-network property like `(X_f ≥ X_g) ⇒ (Y_f[3] < Y_g[3])` on two copies of a net?

### 3.1 `equal-to` (monotonic_acasxu) — *feasible in principle* via a difference/product net

Because `g` `equal-to` `f` means **literally the same ONNX** (confirmed in
`instances.csv`), the query is: over inputs `(X_f, X_g)` constrained by a box on `X_f`
and `X_g` plus cross-input linear relations (`X_f[0] ≥ X_g[0]`, `X_f[i] == X_g[i]`
for i≠0), is the unsafe output relation `Y_f[3] ≥ Y_g[3]` (negation of the property
`Y_f[3] < Y_g[3]`) reachable? Two sound NNV-compatible constructions:

1. **Product network + joint input Star.** Build a single network `h(X_f, X_g) =
   [f(X_f); f(X_g)]` (two parallel copies of `f` sharing weights). The cross-input
   linear constraints (`X_f[0] ≥ X_g[0]`, `X_f[i]=X_g[i]`) are exactly linear
   predicate constraints on a joint **Star** over the stacked input `[X_f; X_g]` —
   NNV `Star(V, C, d, pred_lb, pred_ub)` already supports arbitrary linear predicate
   constraints `C·a ≤ d`, so the coupled input region is representable *without
   approximation*. Reach `h` to a joint output Star `[Y_f; Y_g]`, then check the
   linear unsafe halfspace `Y_f[3] − Y_g[3] ≥ 0` with the existing
   `verify_specification`/HalfSpace machinery. **Soundness: yes** (over-approx reach +
   linear output halfspace is the standard NNV flow). The `X_f[i]==X_g[i]` equalities
   collapse shared predicates so the *effective* number of input predicates is ~5,
   not 10 — i.e. it is essentially a single acasxu reach with a wider/линейно-coupled
   box, well within NNV's acasxu exact-star capability.
   - **Catch (monotonic over-approx):** with `approx-star`, `Y_f` and `Y_g` are
     over-approximated *independently*, losing the correlation that they come from the
     *same* network on *correlated* inputs — so `Y_f[3] − Y_g[3]` bounds blow up and
     you likely get `unknown`, not `unsat`. The product-net construction above keeps
     the correlation **only if reach is run jointly on `h`** (shared predicates), which
     `exact-star` honors. For acasxu (small) **exact-star is already the configured
     method** (`run_vnncomp_instance.m:481`), so monotonic acasxu is plausibly
     *soundly verifiable* — promising but unproven; must be validated empirically.
2. **Difference net `d(X_f,X_g)=f(X_f)−f(X_g)` then check `d[3] ≥ 0`.** Equivalent;
   appends a final linear `(−I, I)` layer onto the product net. Same soundness.

**Bottom line for `equal-to`:** **sound construction exists and reuses NNV reach**;
the work is (a) parsing the cross-network constraints into a joint Star, (b) building
the weight-shared product net (a small `NN` assembly), (c) runner plumbing to load
*one* onnx but instantiate it twice. Tractable as a *research-grade* item, **not a
safe June-30 commitment** because of the empirical-tightness risk and runner rework.

### 3.2 `isomorphic-to` (isomorphic_acasxu) — *harder*

`f` and `g` are **different weights** (`onnx/original` vs `onnx/perturbed`), so the
product net stacks **two distinct** networks `h(X)=[f(X); g(X)]` (inputs identified
by `X_f==X_g`, so a *single* shared input `X`). The construction is the same (joint
reach, linear output check), and **still sound**, but:
- Two networks must be imported and composed into one `NN` — more plumbing than the
  weight-shared case.
- The output property is an **ε-equivalence band** written as `(and (> Y_g (+ Y_f ε))
  (< Y_g (- Y_f ε)))` per output — note this literal AND is the *empty* set (can't be
  both `> Y_f+ε` and `< Y_f−ε`), so either the benchmark intends per-clause OR
  semantics or the unsafe region is "outside the ±ε band" (`|Y_f−Y_g| > ε`), which is
  a **disjunction of two halfspaces per output** → NNV can represent it as an OR of
  HalfSpaces (the existing `prop` cell-of-HalfSpaces structure already encodes OR).
  This ambiguity must be resolved against the spec/VNN-COMP tooling before trusting
  any verdict (−150 risk).

### 3.3 Honest verdict

- **Single-network linear 2.0 (28 benches):** fully within NNV today — *parser only*.
- **`equal-to` multi-net (monotonic):** sound construction exists, likely verifiable
  with exact-star on acasxu; **research-grade, not guaranteed by June 30.**
- **`isomorphic-to` multi-net (isomorphic):** sound construction exists but more
  plumbing + an output-semantics ambiguity to resolve; **research-grade, lower
  priority.**
- **Nonlinear single-net (adaptive_cruise):** NNV's linear reach **cannot soundly
  verify** the variable·variable product constraints in general → **`unknown` is the
  honest answer** (parser can still *read* it; falsification/PGD can still find `sat`).
- **Multimodal single-net (smart_turn):** property is trivial; the only work is a
  two-tensor input set. Feasible but needs `create_input_set` to handle a *list* of
  input boxes feeding a multi-input ONNX — moderate effort.

---

## 4. Phased implementation plan (least-effort-highest-value first)

### Phase 0 — Dispatch + recognize 2.0 (tiny, do first)

- Add a one-line sniff: read the first non-comment line; if it matches
  `(vnnlib-version <2.0>)`, route to `load_vnnlib2`; else fall through to the existing
  `load_vnnlib` (which keeps handling the 4 *mislabeled-1.0* files in `2.0/` dirs —
  they have no version header, so they correctly take the 1.0 path).
- File: a small `load_vnnlib_dispatch.m` (or a guard at the top of a renamed entry
  point) that `run_vnncomp_instance.m:41` calls instead of `load_vnnlib` directly.
- **Risk: none** (pure addition; 1.0 behavior byte-for-byte unchanged).

### Phase 1 — Single-network 2.0 parser (`load_vnnlib2.m`) — THE high-value phase

Scope: one `declare-network`, `float32`/`real` inputs/outputs, box input + linear
(possibly disjunctive `or`) output. **Recommendation: extend in MATLAB, do NOT make
the Python `vnnlib` package the primary path** (rationale in §4.5).

**Output data structure (must be drop-in compatible with `run_vnncomp_instance.m`):**
```
property.lb   : single column vector (single)   OR cell{M} of vectors (disjunctive input)
property.ub   : same shape as lb
property.prop : cell{N} of structs, each with field .Hg = array of HalfSpace
                (HalfSpace.G·y <= HalfSpace.g encodes the UNSAFE output region;
                 array length > 1 = OR of disjuncts, exactly as 1.0 produces)
```
This is the **identical contract** the runner already consumes
(`run_vnncomp_instance.m:42-44, 79-97, 159-184`), so **no runner change is needed for
single-network 2.0** beyond Phase 0 dispatch. New parser responsibilities:
1. Parse `declare-network`/`declare-input`/`declare-output`; record per-network input
   shape (product of dims = flat dimension) and output dimension. Single-net ⇒ ignore
   the `_f` namespacing (variable is `X`/`Y`, or `X1`/`X2`).
2. **Multi-dim index → flat index.** `X[i,j,...]` with declared shape `[d0,d1,...]` →
   row-major flat `idx = ((i*d1)+j)*d2+...`. (Matches ONNX/vnnlib row-major; the
   runner already permutes per category via `needReshape`, so the parser should emit
   the **flat ONNX order** the 1.0 parser emits, keeping `needReshape` semantics
   intact.) **Soundness-critical:** an off-by-one or wrong-major-order here silently
   mislabels dimensions → wrong box → −150. Add a unit test against a hand-checked
   small instance (`test/2.0`).
3. Box input from `(<= X[..] u)` / `(>= X[..] l)` and `(and ...)` of these → fill
   `lb/ub`. Support the `(assert (and (>= ..) (<= ..)))` paired form used by 2026
   files (1.0 used two separate asserts; 2.0 pairs them in one `and`).
4. Output assertions → `HalfSpace` arrays exactly as `process_assertion`/`process_or`
   do in 1.0 (`load_vnnlib.m:190-362`); reuse that code. Handle `(or ...)` (nn4sys,
   collins_aerospace) as OR of HalfSpaces (already supported shape).
5. **`==`/`!=` handling (single-net):** `(== X[i] c)` ⇒ `lb[i]=ub[i]=c` (degenerate
   box; runner already special-cases `lb==ub` as trivially verified,
   `run_vnncomp_instance.m:149,205,281`). `(== Y[i] c)` in an unsafe region ⇒ two
   halfspaces `Y[i]≤c ∧ Y[i]≥c`. `!=` in an *unsafe* region is a disjunction
   (`Y<c ∨ Y>c`) → OR of halfspaces. **If a single-net file mixes `==`/`!=` into a
   form the HalfSpace encoding can't represent soundly, return a sentinel that makes
   the runner emit `unknown`** (never guess).
6. **Multimodal (smart_turn): two input tensors.** Emit `property.lb`/`ub` as the
   concatenation of the two flattened boxes **plus** a record of the per-tensor split
   (`property.inputSplit = [numel(X1), numel(X2)]`) so the runner can build a
   multi-input set. This is the one place the runner *does* need a small change
   (Phase 3b). If multimodal is descoped, the parser can still flag it `unknown`.
7. **Nonlinear (adaptive_cruise): detect and refuse soundly.** If any assertion
   contains a `(* <var> <var>)` (variable·variable) or nested nonlinear arithmetic the
   linear HalfSpace encoding can't represent, set `property.unsupported = true`. The
   runner then **skips reach (→ unknown)** but **may still run falsification/PGD** to
   find a `sat` witness (sound: a concrete violating point is always valid). Wire this
   as: parse the *input* box if it is linear (most ACC input bounds are simple boxes)
   so PGD has a region to sample, and mark the *output* property nonlinear→unknown.

**Deliverables of Phase 1:** `load_vnnlib2.m`, a `flatten_index.m` helper, and a
script-test under `engine/utils/` or the submission test dir that parses one file per
single-net benchmark and asserts `lb≤ub`, correct dims, and HalfSpace polarity against
a hand-checked oracle for `test/2.0` and `acasxu_2023/2.0`.

### Phase 2 — Multi-network handling (OPTIONAL, research-grade)

Only if Phase 1 lands with time to spare. Scope `equal-to` first (monotonic), then
`isomorphic-to` (isomorphic).

Parser additions (a *separate* output shape, since the runner contract differs):
```
property.networks   : struct array, one per declare-network:
                        .name, .inputName(s), .outputName, .inShape, .outShape,
                        .equivOf (name or ''), .equivKind ('equal'|'isomorphic'|'')
property.jointLb/jointUb : box over the STACKED input [X_f; X_g]
property.jointC, .jointd : linear predicate constraints coupling X_f,X_g
                           (from cross-network input asserts: X_f[0]>=X_g[0], ==, ...)
property.crossProp  : HalfSpace(s) over STACKED output [Y_f; Y_g] (unsafe region),
                       e.g. [0..,+1 at Y_f[3], -1 at Y_g[3]] <= 0  for  Y_f[3] >= Y_g[3]
```
Verification (new helper, e.g. `verify_multinet.m`):
1. Build product net `h` = stack of the per-network NNs (weight-shared for `equal-to`;
   two imports for `isomorphic-to`). A thin `NN` wrapper that runs each sub-net on its
   slice of the stacked input and concatenates outputs.
2. Build the joint input **Star** from `jointLb/jointUb` + `jointC/jointd`.
3. `h.reach(IS, exact-star)` → joint output Star; check `crossProp` with
   `verify_specification`. Sound by construction (§3.1).
4. Falsification: sample the joint box (respecting equality couplings), evaluate both
   sub-nets, check the cross-output halfspace → sound `sat` witnesses.

**Gate hard:** if the cross-network constraints are not expressible as linear joint
predicates, or the product-net assembly fails, **return `unknown`.**

### Phase 3 — Wiring into the sweep

- **3a (single-net, required):** in `run_vnncomp_instance.m:41`, call the Phase-0
  dispatch. **Nothing else changes** for the 28 single-net benches — same `lb/ub/prop`
  contract. Add the `category` strings for the 2.0-only single-net benches to
  `load_vnncomp_network.m` (adaptive_cruise → linear-input net but mark nonlinear-out
  unknown; smart_turn → multimodal). The sweep driver must point at the `2.0/`
  instances.csv (same `onnx,vnnlib,timeout` format — confirmed for all single-net
  2.0, including the 2 extended-only single-net ones).
- **3b (multimodal, small):** teach `create_input_set` (and `create_random_examples`)
  to consume `property.inputSplit` and build a set over two input tensors for
  smart_turn. Only needed if smart_turn is in scope.
- **3c (multi-net, optional):** a *separate* runner branch that, on
  `property.networks` present, loads the `(name,onnx)` **tuple list** from the
  multi-net `instances.csv` (Python-literal column 0 — needs a small CSV/tuple
  parser), builds the product net, and calls `verify_multinet`. This is the only place
  the runner's `(category, onnx, vnnlib)` signature must be generalized to a *list* of
  onnx paths.

### 4.5 Parser approach: extend MATLAB vs Python bridge — RECOMMENDATION

**Recommend: extend in MATLAB (`load_vnnlib2.m`), use the Python `vnnlib` package
(dlshriver) as an offline cross-check oracle only.**

- The Python `vnnlib` package *does* parse 2.0 — it released **v1.0.0 on 2025-12-15**,
  the same day as the 2.0 standard, and provides "a full VNN-LIB parser which
  generates an AST ... and a transformer class," with `declare-network` support
  [S8][S9]. So a bridge is *possible*.
- **But** the NNV runner is MATLAB and already has a battle-tested 1.0 parser whose
  helper functions (`process_assertion`, `process_or`, `process_constraint`,
  HalfSpace assembly) are **directly reusable** for single-net 2.0 — the only genuinely
  new parsing is the `declare-network` header and multi-dim index flattening. A MATLAB
  extension avoids a Python⇄MATLAB serialization layer on the hot path and keeps the
  −150-critical box/halfspace construction in one auditable place.
- The vendored `tools/onnx2nnv_python/` shows the team already runs Python for ONNX
  import; adding `pip install git+https://github.com/dlshriver/vnnlib.git` and a tiny
  script that parses a 2.0 file to a normalized JSON (`{lb, ub, prop_halfspaces}`)
  makes an excellent **differential oracle** for the MATLAB parser in tests (parse the
  same file both ways, assert identical lb/ub/HalfSpaces). Use it for *trust*, not for
  the live verdict.
- For the **multi-network** case specifically, the Python AST is more attractive
  (cross-network constraint extraction is fiddly), so Phase 2 may justify a Python
  preprocessing step that emits the `property.networks`/`jointC`/`crossProp` JSON,
  consumed by MATLAB. Decide at Phase 2 start.

---

## 5. Scope recommendation for June 30, risks, fallback

### 5.1 Recommended scope (commit)

1. **Phase 0 + Phase 1 + Phase 3a** — single-network 2.0 parser + dispatch + sweep
   wiring. This unlocks **all 28 genuine single-net 2.0 benchmarks** plus correctly
   routes the 4 mislabeled-1.0 ones, and reads `adaptive_cruise`/`smart_turn` even if
   it can only *falsify* them. **This is the high-value, low-risk core and should be
   the June-30 commitment.**
2. **Stretch:** Phase 3b (smart_turn multimodal input set) — small, self-contained.

### 5.2 Stretch / research (do not commit)

3. **Phase 2 + Phase 3c** — multi-network `equal-to` (monotonic acasxu) via the
   product-net construction. Sound in principle (§3.1), plausibly verifiable with
   exact-star, but needs runner generalization and empirical validation. Treat as
   upside, not a deliverable.
4. `isomorphic-to` (isomorphic acasxu) — lowest priority; output-band semantics
   ambiguity (§3.2) must be resolved first.

### 5.3 Risks

- **R1 (−150, highest): index flattening / major-order.** A wrong row-major vs
  column-major flatten silently mislabels input dims → wrong box → invalid sat/unsat.
  Mitigation: differential test vs the Python `vnnlib` oracle and vs the existing 1.0
  files (acasxu has both 1.0 and 2.0 — parse both, assert identical effective box).
- **R2 (−150): `==`/`!=` polarity.** Mis-encoding `!=`/`==` unsafe regions flips
  sat/unsat. Mitigation: encode `==` only as degenerate box (input) or paired
  halfspaces (output); encode `!=`/disjunctions as OR-of-HalfSpaces; **on any form not
  cleanly representable, return `unknown`.**
- **R3: multimodal split correctness.** Wrong per-tensor split or ordering for
  smart_turn → wrong input. Mitigation: descope smart_turn to `unknown` if not
  validated; it has trivial output so little is lost.
- **R4: nonlinear over-claim.** Never attempt to "linearize" `(* X X)` and claim
  soundness. Mitigation: detect nonlinear → `unknown` for reach; PGD-only for sat.
- **R5 (multi-net): witness format undefined by spec** [S6]. A multi-net `sat` must
  serialize *both* `X_f` and `X_g`; the exact text format is set by VNN-COMP tooling,
  not the standard. Mitigation: defer multi-net `sat` reporting until the format is
  confirmed against the 2026 rules repo; until then, multi-net emits only `unsat`/
  `unknown` (never an unverifiable `sat`).
- **R6: over-approx looseness on monotonic.** approx-star decorrelates `Y_f`,`Y_g` →
  `unknown`. Mitigation: use exact-star (already default for acasxu) and the joint
  product-net reach; accept `unknown` rather than a wrong verdict.

### 5.4 Fallback

**The 1.0 track is already fully covered** by `load_vnnlib.m` + `run_vnncomp_instance.m`
— every benchmark with a `1.0/` dir (30 of 34) is reachable through the existing,
hardened path. **2.0 support is pure upside.** If Phase 1 slips, NNV loses only the 4
extended-only 2.0-only benchmarks (and the 2.0-syntax variants of benches it already
scores on via 1.0). Therefore the safe posture is: ship Phase 0/1/3a, keep everything
**sound-or-unknown**, and treat multi-network as a post-deadline research item.

---

## Sources

- [S1] VNN-LIB official site — version, 2.0 release, links to spec PDF and GitHub:
  https://www.vnnlib.org/
- [S2] VNNLIB-Standard repo (standard.tex, syntax/grammar.cf, v2.0.0 released
  2025-12-15): https://github.com/VNNLIB/VNNLIB-Standard
- [S3] `syntax/grammar.cf` (LBNF grammar: vnnlib-version, declare-network,
  equal-to/isomorphic-to, declare-input/output/hidden, comparison/boolean/arithmetic
  operators, assertions): https://github.com/VNNLIB/VNNLIB-Standard (syntax/grammar.cf)
- [S4] `document/standard.tex` structure (chapters: models, query_language, logics,
  verifier_interface): https://github.com/VNNLIB/VNNLIB-Standard
- [S5] `document/chapters/query_language.tex` (verbatim definitions of declare-network,
  declare-input/output/hidden, equal-to, isomorphic-to, operators, satisfiability):
  https://raw.githubusercontent.com/VNNLIB/VNNLIB-Standard/main/document/chapters/query_language.tex
- [S6] VNN-LIB 2.0 paper, "Rigorous Foundations for Neural Network Verification"
  (Roy, Antony, Gimelli, Daggitt, CAV 2026) — satisfiability convention, equal-to vs
  isomorphic-to semantics, `real` type, nonlinear products allowed, multi-dim indexing,
  declare-hidden, witness mapping, no transitive chaining:
  https://arxiv.org/html/2605.07451
- [S7] VNN-COMP convention (unsat=proved, sat+counterexample=disproved,
  unknown/timeout otherwise): https://vnn-comp.github.io/
- [S8] dlshriver/vnnlib Python package (full VNN-LIB parser, AST + transformer,
  declare-network; install via pip from GitHub): https://github.com/dlshriver/vnnlib
- [S9] vnnlib on PyPI (VNNLIB-Python 1.0.0 released 2025-12-15, part of the official
  ecosystem): https://pypi.org/project/vnnlib/

### Local files cited
- `code/nnv/engine/utils/load_vnnlib.m` — 1.0 parser; output contract
  (`property.lb/ub/prop`), assumptions (`:10-18`), `>=`/`<=`-only (`:16`),
  HalfSpace assembly (`:190-362`).
- `code/nnv/examples/Submission/VNN_COMP2025/run_vnncomp_instance.m` — consumes
  `property.lb/ub/prop` (`:41-44`), sat/unsat/unknown status encoding (`:361-377`),
  per-category net loading (`load_vnncomp_network`, `:454-1011`), counterexample
  writer (`:1059-1092`).
- `code/nnv/tools/onnx2nnv_python/README.md` — vendored Python ONNX importer (basis
  for a Python `vnnlib` differential-oracle bridge).

### Benchmark files cited (under `../vnncomp2026_benchmarks/benchmarks/`)
- `monotonic_acasxu_2026/2.0/vnnlib/instance_0.vnnlib.gz`, `.../2.0/instances.csv`
  (multi-net `equal-to`, same onnx for f and g).
- `isomorphic_acasxu_2026/2.0/vnnlib/instance_0.vnnlib.gz`, `.../2.0/instances.csv`
  (multi-net `isomorphic-to`, original vs perturbed onnx).
- `adaptive_cruise_control_non_linear_2026/2.0/vnnlib/instance_0.vnnlib.gz`
  (single-net nonlinear).
- `smart_turn_multimodal_2026/2.0/vnnlib/instance_0.vnnlib.gz` (single-net, 2 inputs).
- `acasxu_2023/2.0`, `test/2.0`, `cgan2026/2.0`, `nn4sys/2.0`, `collins_aerospace_benchmark/2.0`,
  etc. (single-net 2.0, `float32`, box+linear/or).
- `cersyve/2.0`, `cora_2024/2.0`, `malbeware/2.0`, `safenlp_2024/2.0` (1.0 syntax
  mislabeled into `2.0/` — `(declare-const X_0 Real)`, no version header).
