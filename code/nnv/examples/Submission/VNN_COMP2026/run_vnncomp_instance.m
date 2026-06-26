function [status, tTime] = run_vnncomp_instance(category, onnx, vnnlib, outputfile)

% single script to run all instances (supported) from the vnncomp2025

t = tic;
status = 2; % unknown (to start with)

% CATEGORY-CONDITIONAL competition env (ANTI-DRIFT). These per-category tuning flags used to live ONLY in
% run_instance.sh, so dev harnesses that re-declare config (robust_runner.m, sweep_lambda.sh) silently
% drifted -- the 2026-06-24 cifar/safenlp class. Setting them HERE, keyed off the category, makes EVERY
% harness faithful automatically. Guarded with isempty so an explicitly-set env still wins, and the official
% run_instance.sh (which sets them first) is a no-op here. NOT the conv SOUNDNESS-policy opt-ins
% (NNV_CONV_TRUST_FP32 etc.) -- those stay deliberately harness-set. See BENCHMARK_RUN_MATRIX.md.
if contains(category, 'safenlp') && isempty(getenv('NNV_FALSIFY_MAXTIME'))
    setenv('NNV_FALSIFY_MAXTIME', '8');     % safenlp 20s timeout: default 30s PGD over-runs -> unsats lost
end
if contains(category, 'cifar100')
    if isempty(getenv('NNV_BAB_BETA_ITERS')),     setenv('NNV_BAB_BETA_ITERS', '3');     end
    if isempty(getenv('NNV_CONV_BETA_FRONTIER')), setenv('NNV_CONV_BETA_FRONTIER', '64'); end
    if isempty(getenv('NNV_AMORT_ALPHA')),        setenv('NNV_AMORT_ALPHA', '20');        end
end

% SOUNDNESS HARDENING (persistent/shared-session safety): the per-instance reach budget lives in GLOBALS
% (NNV_REACH_T0/_BUD) that ANOTHER verifier entry point (NN.verify_vnnlib / NN.verify_robustness) reads via
% verify_specification, where under exactReach a result==2 is promoted to 0 (sat). A budget left ARMED after
% this instance could make a LATER exactReach call in the SAME MATLAB session read a "gave-up -> 2" as a
% FALSE sat. So clear any stale budget on ENTRY and guarantee teardown on EXIT (normal return OR error) via
% onCleanup. No effect on the scored one-process-per-instance harness; matters for the persistent sweep.
clear global NNV_REACH_T0 NNV_REACH_BUD %#ok<GVMIS>
reachBudgetCleanup = onCleanup(@() i_arm_reach_budget(NaN)); %#ok<NASGU> % disarm budget globals on any exit

% CODEGEN PACKAGES (custom-layer +package, e.g. +mnist_concat / +<relusplitter-onnx>). importNetworkFromONNX
% writes a GENERATED custom-layer +package (named after the onnx) into the CURRENT FOLDER, and the per-onnx
% netcache (.netcache.mat, built at prepare) stores a dlnetwork that REFERENCES those custom-layer classes.
% So a cache load only reconstructs the real net when that +package is present on the MATLAB path -- if it is
% missing, MATLAB SILENTLY substitutes default layer objects ("Unable to load instances of class ... default
% objects will be substituted") and the net degrades to a useless 'unknown'. The packages are therefore
% REQUIRED ARTIFACTS: prepare_instance.sh pre-generates them into code/nnv (on startup_nnv's genpath) and they
% must persist (never delete; .gitignore'd). At a timed run the netcache HITS, so importNetworkFromONNX is NOT
% re-invoked -> no regeneration, no concurrent-write 'getExternalLayers' collision, and no per-instance cwd
% juggling is needed. (An earlier "give each instance a fresh scratch cwd" isolation was REMOVED: it sent any
% regeneration to an ephemeral tempdir that was then deleted, so the NEXT process's cache load could not find
% the package -> silent degrade -> 'unknown'. It was incompatible with the netcache and is the regression that
% broke dist_shift + relusplitter.) Competition runs one process per instance (no concurrent import), so a
% rare cache miss regenerates into code/nnv safely; the sweep relies on a serially pre-warmed cache.

% disp("We are running...")

%% 0) cctsdb_yolo: complete-enumeration path (bypasses net import + reach)
%
% Every cctsdb_yolo_2023 spec fixes ALL inputs (lb == ub) except the two patch
% position coordinates X_12288/X_12289 in [0,62], and the ONNX graph consumes
% those two ONLY through Cast(int64) truncation (Gather(12288/12289) -> Cast),
% so the network output is PIECEWISE CONSTANT on unit cells and enumerating the
% <= 3969 integer points is SOUND AND COMPLETE. cctsdb_enumerate.py re-verifies
% all of that structurally per instance (any deviation -> unknown), runs the
% enumeration through onnxruntime, and returns SAT-with-witness / UNSAT /
% UNKNOWN. MATLAB's importNetworkFromONNX cannot handle this model (the old
% stub errored "Working on supporting this one"), so this branch must route
% AROUND load_vnncomp_network / falsify / reach entirely and write the result
% file in the same format as section 4 below.
if contains(category, "cctsdb_yolo")
    [status, counterEx] = verify_cctsdb_enumeration(onnx, vnnlib);
    tTime = toc(t);
    disp("Verification result: " + string(status));
    disp("Total Time: " + string(tTime));
    fid = fopen(outputfile, 'w');
    if status == 0
        fprintf(fid, 'sat \n');
        fclose(fid);
        write_counterexample(outputfile, counterEx);
    elseif status == 1
        fprintf(fid, 'unsat \n');
        fclose(fid);
    else
        fprintf(fid, 'unknown \n');
        fclose(fid);
    end
    return;
end

%% 0a) nn4sys lindex / lindex_deep: huge OR over ULTRA-THIN X-boxes -> sound batched-IBP decider.
% A 1-in/1-out learned-index FC+ReLU net whose spec is a disjunction over 200-60000 ultra-thin X-boxes
% (~8e-8 wide), each requiring Y in a safe band. NNV's general star/LP reach walks the clauses one at a
% time and TIMES OUT (lindex_200 ran 444s). Route to nn4sys_lindex_decide.py: sound IBP + ADAPTIVE box
% subdivision (IBP is exact at a point, so subdividing the thin box is sound AND tightens it ~linearly)
% -> decides every lindex instance in <6s. SOUND-OR-UNKNOWN: unsat only when every box's sound interval
% is provably inside its band; sat only with an onnxruntime witness; else unknown. (Validated standalone:
% 24/24 lindex instances -> unsat, IBP-contains-forward self-check green.)
if contains(category, "nn4sys") && contains(onnx, "lindex")
    [status, counterEx] = verify_nn4sys_lindex(onnx, vnnlib);
    tTime = toc(t);
    disp("Verification result: " + string(status));
    fid = fopen(outputfile, 'w');
    if status == 0
        fprintf(fid, 'sat \n'); fclose(fid); write_counterexample(outputfile, counterEx);
    elseif status == 1
        fprintf(fid, 'unsat \n'); fclose(fid);
    else
        fprintf(fid, 'unknown \n'); fclose(fid);
    end
    return;
end

%% 0b) adaptive_cruise: NONLINEAR-property falsification path (bypasses NNV's gated load_vnnlib)
%
% adaptive_cruise_control_non_linear_2026 has NONLINEAR vnnlib (X*Y, X*X) which NNV deliberately
% GATES as unsupported (never linearize -> sound), so the normal path always returns unknown.
% adaptive_cruise_falsify.py finds a SAT witness by sampling the input box + forwarding the REAL
% onnx via onnxruntime, and accepts it ONLY if EVERY assertion holds under the AUTHORITATIVE vnnlib
% parser (the `vnnlib` pypi package: parse + to_dnf, exact arith eval on its AST) with a ROBUST
% margin (1e-2) -- which excludes the thin saturation artifact (net saturates ~100.0011, only 1.2e-4
% above the 100.001 threshold). SAT-or-FALL-THROUGH: a confirmed SAT is emitted + returns; anything
% else falls through to the normal (unknown) path -> ZERO regression. Never emits unsat (sampling
% can't prove it). The witness is authoritatively-parsed + ORT-forwarded + robust-margin -> sound.
if contains(category, "adaptive_cruise")
    [acStatus, acCounterEx] = verify_adaptive_cruise_falsify(onnx, vnnlib);
    if acStatus == 0   % SAT, authoritatively confirmed
        status = 0;    % FIX: set the RETURN status (was left at the default 2). The result FILE was
                       % written "sat" below, so the OFFICIAL path (execute.py reads the file) was fine,
                       % but the DEV SWEEP (robust_runner.m / run_all_benchmarks read the RETURN value)
                       % recorded `unknown` -> adaptive_cruise sats were silently undercounted.
        tTime = toc(t);
        disp("Verification result: 0");
        disp("Total Time: " + string(tTime));
        fid = fopen(outputfile, 'w');
        fprintf(fid, 'sat \n'); fclose(fid);
        write_counterexample(outputfile, acCounterEx);
        return;
    elseif acStatus == 1   % UNSAT: input region PROVABLY EMPTY (nonlinear input constraint infeasible
        status = 1;        % over the box) -> no input -> no counterexample -> vacuously UNSAT. Proven by
        % interval arithmetic over the AUTHORITATIVE parsed AST (+ dense-grid refutation backstop) in
        % adaptive_cruise_falsify.py, which emits 11 ONLY when proven; any uncertainty there is 12 ->
        % unknown -> falls through below. MUST return here (like SAT): otherwise the normal load path
        % overwrites this with the gated nonlinear-vnnlib 'unknown' (the 2026-06-25 first-cut bug).
        tTime = toc(t);
        disp("Verification result: 1");
        disp("Total Time: " + string(tTime));
        fid = fopen(outputfile, 'w');
        fprintf(fid, 'unsat \n'); fclose(fid);
        return;
    end
    % unknown -> fall through to the normal load/reach path (still unknown; preserves prior behaviour)
end

% collins_aerospace_benchmark: sound FALSIFICATION-ONLY path, fully outside MATLAB's
% importer. importNetworkFromONNX dies on the YOLOv5 Detect-head custom layers, the
% old import+matlab2nnv route produced invalid SAT instances (the vnnlib X order is
% HWC flat, not CHW -- see collins_falsify.py), and set-based reach at 640x640x3 is
% infeasible (multi-GB), so UNSAT is out of scope for this category. The Python
% falsifier validates the HWC mapping with a center-consistency gate, falsifies with
% finite-difference gradients via onnxruntime, and re-validates any candidate with
% margin before claiming SAT (exit 10). Anything else -> unknown. NEVER emits unsat.
if contains(category, 'collins_aerospace_benchmark')
    [status, tTime] = run_collins_falsifier(onnx, vnnlib, outputfile, t);
    return;
end

% smart_turn_multimodal_2026: a 692-node INT8-quantized MULTIMODAL transformer
% (199 DequantizeLinear + 119 QuantizeLinear + 22 Conv + 9 Erf + 5 Softmax, with a
% 5D video input). MATLAB's importer collapses it to an opaque fused blob NNV
% cannot soundly evaluate. Abstain (unknown) rather than risk feeding the fused
% net to cp-star/manifest and emitting an unsound verdict -- unknown is always a
% safe VNN-COMP verdict. A real fix (multimodal + quantization) is out of scope.
if contains(category, "smart_turn")
    status = 2;            % unknown
    tTime = toc(t);        % assign both outputs before the early return
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unknown \n');
    fclose(fid);
    return;
end

% nn4sys mscn: abstain ALL mscn models to unknown. Two distinct causes, one sound outcome:
%   (1) mscn_2048d_dual.onnx is corrupt upstream -- fails BOTH MATLAB's ONNX parser AND onnx 1.20's
%       ParseFromString ("Wire format was corrupt"); no tool can load it. [confirmed 2026-06-15]
%   (2) mscn_128d / mscn_2048d produce FALSE-SAT on the cardinality_* specs. [confirmed 2026-06-16:
%       10/10 NNV `sat` contradicted the abCROWN/NeuralSAT/PyRAT `unsat` consensus -- 10 x -150; NNV
%       had 0 correct mscn solves.] Root cause: these specs are MULTI-OUTPUT, so the falsifier takes
%       the ungated path and the onnxruntime SAT-witness replay (run_vnncomp_instance Pillar-2 gate)
%       never runs -- it only wraps the single-Hg path; and the box's system python3 lacks
%       onnxruntime regardless. So a spurious witness was emitted as a wrong verdict.
% Abstaining to unknown (0 points, never -150) is the sound verdict and costs no real solves.
% TODO durable fix: replay EVERY sat witness through the full vnnlib via the venv python before
% emitting (gate the multi-output path too), then this mscn abstain can be lifted. See diagnosis doc.
if contains(category, "nn4sys") && contains(onnx, "mscn")
    status = 2;            % unknown
    tTime = toc(t);
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unknown \n');
    fclose(fid);
    return;
end

%% 1) Load components

% Load networks

% LOAD GUARD: a failure at the load stage (importNetworkFromONNX / matlab2nnv / load_manifest_net /
% initialize -- e.g. a missing python-importer manifest, or an unsupported layer) would otherwise
% throw out of this function and abort the instance with NO output file written -- a crash the sweep
% cannot score. Convert it to a sound 'unknown' (0 points, never a wrong verdict), matching the
% sound-or-unknown discipline the reach stage already uses (the reach try/catch blocks below). During
% the untimed prepare phase (NNV_PREP_CACHE=1) just return without writing a result, since prepare
% ignores the verdict and the timed run re-attempts (and re-guards).
try
    [net, nnvnet, needReshape, reachOptionsList, inputSize, inputFormat, nRand, falsifyOpts] = i_load_vnncomp_network_cached(category, onnx, vnnlib);

    if isempty(inputSize)
        inputSize = net.Layers(1, 1).InputSize;
    end

    % Competition lever (NNV_FALSIFY_MAXTIME, seconds): cap the falsify-first PGD budget so the SOUND reach
    % ladder keeps its window within the OFFICIAL per-instance timeout. Several ftab rows set a PGD max_time
    % >= the timeout (safenlp 30s/to=20s, cora 30s/to=30s): on a ROBUST instance PGD burns the whole wall
    % finding nothing and the external kill fires BEFORE the (often <1s) reach proof -> a decidable unsat is
    % lost to timeout. Capping PGD leaves reach its window. STRICTLY SOUND: falsification is SAT-or-unknown
    % only, so less PGD only ever yields fewer SAT (-> unknown), never a wrong verdict; reach is untouched.
    % Applied HERE (after the cached load) so it covers BOTH the cache-hit and cache-miss paths -- the cached
    % falsifyOpts otherwise carries the un-capped ftab max_time. No-op when unset (tests/other callers
    % unaffected); the competition harness sets it, mirroring how execute.py drives NNV_REACH_BUDGET.
    fmt = getenv('NNV_FALSIFY_MAXTIME');
    if ~isempty(fmt) && isstruct(falsifyOpts) && isfield(falsifyOpts, 'max_time')
        fv = str2double(fmt);
        if isfinite(fv) && fv > 0
            falsifyOpts.max_time = min(falsifyOpts.max_time, fv);
        end
    end
catch loadME
    fprintf('LOAD FAILED for %s -> unknown: %s\n', category, loadME.message);
    status = 2; tTime = toc(t);
    if ~strcmp(getenv('NNV_PREP_CACHE'), '1')
        % Defensive: if the output file cannot be opened (fid == -1), fprintf/fclose would THROW
        % and re-break the "never crash" guarantee, so degrade to log-only on an fopen failure.
        fid = fopen(outputfile, 'w');
        if fid > 0
            fprintf(fid, 'unknown \n'); fclose(fid);
        else
            fprintf('WARN: could not open output file %s (fopen=-1) -> unknown logged only\n', outputfile);
        end
    end
    return;
end

% Prepare-phase net-cache pre-warm: prepare_instance.sh sets NNV_PREP_CACHE=1 and invokes this
% ONLY to build the net cache (done inside i_load_vnncomp_network_cached above), WITHOUT running
% the timed verification (VNN-COMP allows this ONNX->MATLAB conversion in the untimed prepare
% phase). Unset in normal runs -> no effect.
if strcmp(getenv('NNV_PREP_CACHE'), '1')
    status = 2; tTime = 0; return;   % cache built; result ignored by prepare
end

% Decide the reach input-set type from the NET's input-layer TYPE (not the input
% SHAPE) whenever the net exposes Layers -- dlnetwork AND SeriesNetwork/DAGNetwork:
% ImageInputLayer/Image3DInputLayer -> ImageStar; FeatureInputLayer -> Star. Leave
% [] (unknown) otherwise (e.g. NNV NN manifest nets, sequence inputs) so
% create_input_set falls back to its shape heuristic. Fixes acasxu (a [1 5 1]
% ImageInputLayer): the shape heuristic built a Star and every acasxu reach errored
% ("Input is not an ImageStar"). Using [] for the unknown case (not false) avoids
% forcing a Star on a non-dlnetwork image net (which would reintroduce that bug).
useImageStar = [];
% Guard the net.Layers access: a malformed imported dlnetwork (e.g. some quantized
% traffic_signs nets under R2026a) throws "Undefined function 'getExternalLayers' for
% input arguments of type 'double'" from dlnetwork.get.Layers. Catching it here keeps
% useImageStar = [] (the sound shape-heuristic fallback) and lets the normal reach path
% run -> a sound 'unknown' instead of an uncaught crash that writes NO result file.
try
    if isprop(net, 'Layers') && ~isempty(net.Layers)
        L1 = net.Layers(1);
        if isa(L1, 'nnet.cnn.layer.ImageInputLayer') || isa(L1, 'nnet.cnn.layer.Image3DInputLayer')
            useImageStar = true;
        elseif isa(L1, 'nnet.cnn.layer.FeatureInputLayer') || isa(L1, 'nnet.onnx.layer.FeatureInputLayer')
            useImageStar = false;
        end
    end
catch ME
    fprintf('net.Layers inspection failed (%s) -> useImageStar=[] (shape heuristic)\n', ME.message);
end

% Load property to verify
property = load_vnnlib(vnnlib);

% VNN-LIB 2.0 gate: load_vnnlib2 (dispatched for 2.0 files) flags any construct NNV
% cannot soundly verify -- 3+ networks / isomorphic-to, nonlinear/arithmetic
% output, multimodal (>1 input tensor), declare-hidden, mixed input/output disjunction.
% Emit `unknown` (0 points) rather than parse it unsoundly and risk a -150 wrong
% verdict. (1.0 properties never set this field, so this is a no-op for them.)
if isfield(property, 'unsupported') && property.unsupported
    if isfield(property, 'reason') && ~isempty(property.reason)
        fprintf('vnnlib 2.0 unsupported -> unknown: %s\n', property.reason);
    end
    status = 2;
    tTime = toc(t);
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unknown \n');
    fclose(fid);
    return;
end

% Phase 3c: a clean multi-network `equal-to` pair parses into property.multinet
% (unsupported = false) with the single-network fields lb/ub/prop intentionally
% EMPTY. `equal-to` means both declared networks are the SAME model (e.g.
% monotonic_acasxu: g equal-to f), evaluated on two inputs coupled by the parsed
% jointC/jointd. verify_multinet builds the product net [f(x_f); f(x_g)], falsifies
% the cross-network unsafe region FIRST (sound, witness validated by concrete
% evaluation) for sat, then reaches the joint input Star for unsat; ANY doubt ->
% unknown (sound-or-unknown; the -150 rule holds). isomorphic-to (different g) is
% flagged property.unsupported above and never reaches here.
if isfield(property, 'multinet')
    [status, counterEx] = verify_multinet(nnvnet, property);
    % AUTHORITATIVE two-network witness gate (isomorphic/monotonic_acasxu, vnnlib-2.0 equal-to).
    % verify_multinet is a COMPLETE, sound relational verifier (falsify-first on the REAL net +
    % product-net joint-Star reach), and its forward matches onnxruntime to ~1e-9 (verified
    % 2026-06-24: the earlier "issue #7 shape mismatch -> spurious instant sat" worry was FALSE --
    % monotonic_acasxu sats are REAL). So instead of a blanket sat->unknown downgrade (which
    % discarded ~50 real sats), re-validate the sat witness through the TWO-network onnxruntime gate
    % (validate_witness_multinet.py: forwards BOTH x_f and x_g, checks the LITERAL two-network spec):
    %   valid    -> trust the sat (ORT-confirmed counterexample)
    %   spurious -> downgrade to unknown (a would-be -150 -> 0 points)
    %   cant     -> downgrade to unknown (sound default: never ship an un-validated multinet sat)
    % The unsat (joint-Star reach) verdict has no replayable witness, so it stays a conservative
    % `unknown` until separately audited. NNV_TRUST_MULTINET=1 bypasses the gate (trusts both
    % directions -- for testing / once the unsat direction is independently validated).
    trust = strcmp(strtrim(getenv('NNV_TRUST_MULTINET')), '1');
    if status == 0
        % stacked witness [x_f; x_g] -> [f(x_f); f(x_g)] (write_counterexample's X_i/Y_i format)
        xf = counterEx{1}; xg = counterEx{2};
        yf = nnvnet.evaluate(xf); yg = nnvnet.evaluate(xg);
        fid = fopen(outputfile, 'w'); fprintf(fid, 'sat \n'); fclose(fid);
        write_counterexample(outputfile, {[xf(:); xg(:)], [yf(:); yg(:)]});
        if ~trust
            gate = authoritative_witness_gate(onnx, vnnlib, outputfile, 'validate_witness_multinet.py');
            if ~strcmp(gate, 'valid')
                fprintf('vnnlib 2.0 multi-network: sat witness gate=%s -> DOWNGRADE to unknown\n', gate);
                status = 2;
                fid = fopen(outputfile, 'w'); fprintf(fid, 'unknown \n'); fclose(fid);
            else
                fprintf('vnnlib 2.0 multi-network: sat witness gate=VALID (onnxruntime-confirmed) -> trust sat\n');
            end
        end
    elseif status == 1
        if trust
            fid = fopen(outputfile, 'w'); fprintf(fid, 'unsat \n'); fclose(fid);
        else
            fprintf(['vnnlib 2.0 multi-network: unsat (joint-Star reach) has no replayable witness ' ...
                     '-> sound unknown (NNV_TRUST_MULTINET=1 to trust)\n']);
            status = 2;
            fid = fopen(outputfile, 'w'); fprintf(fid, 'unknown \n'); fclose(fid);
        end
    else
        fid = fopen(outputfile, 'w'); fprintf(fid, 'unknown \n'); fclose(fid);
    end
    fprintf('vnnlib 2.0 multi-network (equal-to) -> status=%d (0 sat / 1 unsat / 2 unknown)\n', status);
    tTime = toc(t);
    return;
end

lb = property.lb; % input lower bounds
ub = property.ub; % input upper bounds
prop = property.prop; % output spec to verify

% Optional smoke-test epsilon shrink (env NNV_EPS_SHRINK = fraction in (0,1)):
% tighten the input perturbation box toward its center so a compute-bound net can
% reach a verdict at all on an easier property -- a smoke test of tractability, not
% a competition run. Default (unset/empty) is a NO-OP; the property is untouched.
eps_shrink = str2double(getenv('NNV_EPS_SHRINK'));
if ~isnan(eps_shrink) && eps_shrink > 0 && eps_shrink < 1
    if iscell(lb)
        for ii = 1:numel(lb)
            c = (lb{ii} + ub{ii})/2;
            lb{ii} = c - eps_shrink*(c - lb{ii});
            ub{ii} = c + eps_shrink*(ub{ii} - c);
        end
    else
        c = (lb + ub)/2;
        lb = c - eps_shrink*(c - lb);
        ub = c + eps_shrink*(ub - c);
    end
end

% fid = fopen(outputfile, 'w');
% fprintf(fid, 'unknown \n');
% fclose(fid);


%% 2) SAT?

% nRand = 100; % number of random inputs (this can be changed)
% We got some penalties last year, why?
% Wrong vnnlib parsing? Wrong counterrexample writing? Do we need to reshape it?
% Let's test last years properties and make sure those errors/bugs are
% fixed before this year's submission

% Choose how to falsify based on vnnlib file
hitSpc = 0;   % >=1 -> index of the multi-output spec region the witness hit (used by the gate below)
if ~isa(lb, "cell") && length(prop) == 1 % one input, one output
    counterEx = falsify_single(net, lb, ub, inputSize, nRand, prop{1}.Hg, needReshape, inputFormat, falsifyOpts);
elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
    for spc = 1:length(lb) % try parfeval, parfor does not work for early return
        counterEx = falsify_single(net, lb{spc}, ub{spc}, inputSize, nRand, prop{spc}.Hg, needReshape, inputFormat, falsifyOpts);
        if iscell(counterEx)
            hitSpc = spc;
            break
        end
    end
elseif isa(lb, "cell") && length(prop) == 1 % can violate the output property from any of the input sets
    for arr = 1:length(lb) % try parfeval, parfor does not work for early return
        counterEx = falsify_single(net, lb{arr}, ub{arr}, inputSize, nRand, prop{1}.Hg, needReshape, inputFormat, falsifyOpts);
        if iscell(counterEx)
            break
        end
    end
else
    warning("Working on adding support to other vnnlib properties");
end

cEX_time = toc(t);
if ~isempty(getenv('NNV_PHASE_LOG')), fpl=fopen(getenv('NNV_PHASE_LOG'),'a'); if fpl>0, fprintf(fpl,'PHASE import+falsify cEX_time=%.2f\n', cEX_time); fclose(fpl); end, end


%% 3) UNSAT?

% Check if property was violated earlier
if iscell(counterEx)
    status = 0;

    % Pillar-2 final gate (config-robust, ALL falsify paths): replay the SAT witness through
    % onnxruntime on the ORIGINAL ONNX -- the competition's own checker -- before committing to
    % `sat`. validate_witness_onnx finds an ort-capable python whether it lives in a venv or a
    % direct/system interpreter (see ort_python). This catches an import/encoding divergence between
    % NNV's net and the real ONNX, which is the nn4sys mscn false-SAT class: NNV's manifest-imported
    % net found a "counterexample" the real ONNX does NOT have (3 gold tools prove unsat). Unlike the
    % old gate, the MULTI-output path is gated too: use the Hg region the witness actually hit
    % (prop{hitSpc}); the single-output paths fall back to prop{1}.
    if hitSpc >= 1, gateHs = prop{hitSpc}.Hg; else, gateHs = prop{1}.Hg; end
    % "riskyNet" = NNV's net is NOT a faithful independent stand-in for the real ONNX, so when
    % onnxruntime cannot POSITIVELY confirm the witness we must NOT keep `sat`. True when:
    %   - the net is a python/manifest import (isa(net,'NN')) -- the in-memory class is authoritative
    %     and robust to the .netcache.mat cache-hit path, where no .nnv.mat is on disk (a file-only
    %     check would miss it and leave a spurious manifest SAT standing); OR
    %   - the witness depended on a non-identity input reshape/permute (needReshape ~= 0) -- then
    %     validate_witness applies the SAME permute as the falsifier, so it CANNOT see an import/
    %     encoding divergence between NNV's net and the real ONNX; only the onnxruntime replay can
    %     (this is the cifar100/challenging/malbeware/metaroom -150 class). The on-disk .nnv.mat is
    %     OR'd in as a belt-and-suspenders signal.
    % A needReshape==0 standard MATLAB-imported flat-feature net (acasxu-style) is NORMALLY not risky:
    % validate_witness is faithful there, so a `sat` it accepted may stand even without onnxruntime.
    % EXCEPTION (Phase 1): the batched falsifier manufactures THIN, near-boundary acasxu witnesses
    % (2e6 samples + coordinate-descent polish targeting margin -> 0). Such a witness is the most
    % divergence-prone class: an NNV-fp32-vs-real-ONNX sign flip (NNV says -eps=violated, ONNX says
    % +eps=safe) would be a false `sat`. Given the -150 asymmetry, treat ALL acasxu as risky -> REQUIRE
    % onnxruntime confirmation. No effect when ORT is present (every acasxu SAT here is ORT-confirmed);
    % only the no-ORT fallback changes: acasxu SAT -> unknown (sound, loses +10) rather than trusting
    % validate_witness on a boundary witness.
    riskyNet = isa(net, 'NN') ...
            || (~isempty(needReshape) && any(needReshape(:) ~= 0)) ...
            || isfile(regexprep(char(onnx), '\.onnx$', '.nnv.mat')) ...
            || contains(category, "acasxu") ...
            || contains(category, "lsnc_relu");   % SAME thin near-boundary batched+polish witnesses as
                                                  % acasxu -> the most divergence-prone class -> REQUIRE
                                                  % onnxruntime confirmation (don't lean on the .nnv.mat
                                                  % file coupling; principled match to the acasxu reason).
    try
        [orVio, orAvail] = validate_witness_onnx(onnx, counterEx{1}, gateHs);
        if orAvail && ~orVio
            % onnxruntime ran and the witness violates NOTHING on the real ONNX -> spurious -> drop.
            fprintf('onnxruntime replay rejected the SAT witness -> unknown\n');
            status = 2; counterEx = nan;
        elseif ~orAvail && riskyNet
            % No independent onnxruntime confirmation AND NNV's net can diverge from the real ONNX
            % (manifest import or non-identity reshape) -> cannot stand behind this `sat`. Sound-or-
            % unknown rather than risk a -150. (onnxruntime is a stated requirement -- requirements.txt;
            % this is the safety net if it is somehow missing on the run host.)
            fprintf('SAT witness unconfirmed (no onnxruntime) on a divergence-risk net -> unknown\n');
            status = 2; counterEx = nan;
        end
    catch
        % Gate itself errored. Divergence-risk nets -> prefer unknown; a needReshape==0 MATLAB-imported
        % net falls back to validate_witness (which already accepted the witness on NNV's own forward).
        if riskyNet, status = 2; counterEx = nan; end
    end
end

vT = tic;

quickRun = false;
% 
% if quickRun 
%     tTime = toc(t);
%     disp("Quiting early...")
%     return
% end

if status == 2 && ~quickRun % no counterexample found and supported for reachability (otherwise, skip step 3 and write results)

% Choose how to verify based on vnnlib file
    if ~isa(lb, "cell") && length(prop) == 1 % one input, one output 

        if ~nnz(lb-ub) % lb == ub, not a set

            status = 1; % verified, since  we already tested this

        else

            % GPU-BaB sound additive unsat pre-check (FC argmax-robustness timeout
            % categories): try the batched GPU-BaB before Star; a 'robust' verdict is a
            % sound unsat -> emit + skip Star. Anything else falls through unchanged.
            [status, reachOptionsList] = i_gpu_bab_precheck(category, nnvnet, lb, ub, prop, status, reachOptionsList, needReshape, inputSize);

            % Per-instance reach time budget (opt-in NNV_REACH_BUDGET seconds): arms the SOUND cap in
            % PosLin.reach_star_exact so a slow exact-star aborts -> unknown instead of hanging. Reset
            % here per instance so a persistent-session runner gets a fresh budget each instance.
            clear global NNV_REACH_T0 NNV_REACH_BUD; global NNV_REACH_T0 NNV_REACH_BUD %#ok<GVMIS>
            i_reachbud = str2double(getenv('NNV_REACH_BUDGET'));   % reset EVERY instance (no budget leak in a persistent session); arm only when valid
            if isfinite(i_reachbud) && i_reachbud > 0
                NNV_REACH_T0 = tic; NNV_REACH_BUD = i_reachbud;
            end

            while ~isempty(reachOptionsList)

                reachOptions = reachOptionsList{1};
                rt0 = tic;

                IS = create_input_set(lb, ub, inputSize, needReshape, useImageStar);

                % Compute reachability. Reach may FAIL LOUD by design (layers
                % refuse to return unsound sets); an error is mapped to
                % "unknown" for this method and we try the next one -- that is
                % sound (claims nothing), unlike swallowing a wrong set.
                try
                    if ~strcmp(reachOptions.reachMethod, "cp-star")
                        if ~is_nnvnet_valid(nnvnet)
                            % matlab2nnv conversion failed earlier; report unknown.
                            % PRINT it: this exact silent skip masked the dist_shift
                            % regression (50 provable unsats looked like loose bounds).
                            fprintf('reach skipped: matlab2nnv conversion failed earlier -> unknown\n');
                            status = 2; break;
                        end
                        if ~isempty(getenv('NNV_PHASE_LOG')), nc=1; if isfield(reachOptions,'numCores'), nc=reachOptions.numCores; end; fpl=fopen(getenv('NNV_PHASE_LOG'),'a'); if fpl>0, fprintf(fpl,'PHASE reachcall %s START t=%.2f numCores=%g\n', reachOptions.reachMethod, toc(rt0), nc); fclose(fpl); end, end
                        ySet = nnvnet.reach(IS, reachOptions);
                    else
                        ySet = Prob_reach(net, IS, reachOptions);
                    end
                catch ME
                    fprintf('reach (%s) errored: %s -> unknown\n', ...
                        reachOptions.reachMethod, ME.message);
                    status = 2;
                    reachOptionsList = reachOptionsList(2:end);
                    continue;
                end

                % Verify property
                status = verify_specification(ySet, prop);
                if ~isempty(getenv('NNV_PHASE_LOG')), fpl=fopen(getenv('NNV_PHASE_LOG'),'a'); if fpl>0, fprintf(fpl,'PHASE reach %s %.2fs status=%d\n', reachOptions.reachMethod, toc(rt0), status); fclose(fpl); end, end
                if status == 1 && i_is_probabilistic(reachOptions.reachMethod)
                    % cp-star 'unsat' is PROBABILISTIC (conformal), not a sound proof -> tracked.
                    % (When NNV_QUARANTINE_CPSTAR is set, cp-star is stripped upstream so this
                    % branch is never reached; default OFF preserves the original behavior.)
                    i_log_cpstar(onnx);
                end

                if status == 1 % verified, then stop
                    break
                else
                    reachOptionsList = reachOptionsList(2:end);
                end

            end

        end

    elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
        
        local_status = 2*ones(length(lb),1); % track status for each specification in the vnnlib
        i_reachbud = str2double(getenv('NNV_REACH_BUDGET'));  % broadcast the per-instance reach budget into the parfor
                                                             % (globals do not cross to workers; arms the SOUND PosLin cap per worker)

        parfor spc = 1:length(lb) % We can compute these in parallel for faster computation

            i_arm_reach_budget(i_reachbud);  % arm NNV_REACH_T0/_BUD on THIS worker (global decls are illegal in a parfor body -> helper)
            lb_spc = lb{spc};
            ub_spc = ub{spc};

            if ~nnz(lb_spc-ub_spc) % lb == ub, not a set

                local_status(spc) = 1; % verified, since we already tested this
                
            else

                reachOptPar = reachOptionsList;
                
                tempStatus = 2;
                while ~isempty(reachOptPar)

                    reachOptions = reachOptPar{1};

                    IS = create_input_set(lb_spc, ub_spc, inputSize, needReshape, useImageStar);

                    % Compute reachability
                    if ~strcmp(reachOptions.reachMethod, "cp-star")
                        if ~is_nnvnet_valid(nnvnet)
                            % same observability as the single-spec path: never skip silently
                            fprintf('reach skipped: matlab2nnv conversion failed earlier -> unknown\n');
                            tempStatus = 2; break;
                        end
                    end
                    % [29] Layers are now fail-loud by design (a sound refusal,
                    % not a crash). An uncaught reach error here would abort the
                    % whole instance/parfor; instead record UNKNOWN for this reach
                    % option and fall through to the next one.
                    try
                        if ~strcmp(reachOptions.reachMethod, "cp-star")
                            ySet = nnvnet.reach(IS, reachOptions);
                        else
                            ySet = Prob_reach(net, IS, reachOptions);
                        end

                        % Verify property
                        if isempty(ySet.C)
                            dd = ySet.V; DD = ySet.V;
                            ySet = Star(dd,DD);
                        end

                        % Add verification status
                        tempStatus = verify_specification(ySet, prop(spc));
                        if tempStatus == 1 && i_is_probabilistic(reachOptions.reachMethod)
                            i_log_cpstar(onnx);   % parity: parfor cp-star verdicts were previously unlogged
                        end
                    catch
                        tempStatus = 2;   % unknown; try the next reach option
                    end

                    if tempStatus ~= 2 % verified, then stop (or falsified)
                        break
                    else
                        reachOptPar = reachOptPar(2:end);
                    end

                end

                local_status(spc) = tempStatus;

            end
            
        end

        % Check for the global verification result
        if all(local_status == 1)
            status = 1;
        else
            status = 2;
        end

    elseif isa(lb, "cell") && length(prop) == 1 % one specification, multiple input definitions 

        local_status = 2*ones(length(lb),1); % track status for each specification in the vnnlib, initialize as unknown
        i_reachbud = str2double(getenv('NNV_REACH_BUDGET'));  % broadcast the per-instance reach budget into the parfor:
                                                             % globals set on the client do NOT cross to parfor workers,
                                                             % so without this the SOUND cap in PosLin.reach_star_exact is a
                                                             % no-op on every multi-input-box prop (disjunctive-input overruns).

        parfor spc = 1:length(lb) % We can compute these in parallel for faster computation

            i_arm_reach_budget(i_reachbud);  % arm NNV_REACH_T0/_BUD on THIS worker (global decls are illegal in a parfor body -> helper)
            reachOptPar = reachOptionsList;
            
            lb_spc = lb{spc};
            ub_spc = ub{spc};

            if ~nnz(lb_spc-ub_spc) % lb == ub, not a set

                local_status(spc) = 1; % verified, since we already tested this earlier
                
            else

                tempStatus = 2;
                while ~isempty(reachOptPar)

                    reachOptions = reachOptPar{1};

                    IS = create_input_set(lb_spc, ub_spc, inputSize, needReshape, useImageStar);

                    % Compute reachability
                    if ~strcmp(reachOptions.reachMethod, "cp-star")
                        if ~is_nnvnet_valid(nnvnet)
                            % same observability as the single-spec path: never skip silently
                            fprintf('reach skipped: matlab2nnv conversion failed earlier -> unknown\n');
                            tempStatus = 2; break;
                        end
                    end
                    % [29] fail-loud reach must not abort the parfor (see above).
                    try
                        if ~strcmp(reachOptions.reachMethod, "cp-star")
                            ySet = nnvnet.reach(IS, reachOptions);
                        else
                            ySet = Prob_reach(net, IS, reachOptions);
                        end

                        % Add verification status.
                        % This branch is "one output spec, MULTIPLE input boxes"
                        % (length(prop)==1, lb is a cell): each input box is verified
                        % against the SINGLE shared spec prop{1}. The old prop(spc) made
                        % parfor try to SLICE prop by spc=1..length(lb), but prop has only
                        % one element -> "Index exceeds the number of array elements" thrown
                        % at the parfor supply stage (OUTSIDE this try/catch), surfacing as
                        % an uncaught ERR for disjunctive-input props (e.g. acasxu prop_6).
                        % prop(1) is a constant index -> parfor broadcasts prop, no slicing.
                        tempStatus = verify_specification(ySet, prop(1));
                        if tempStatus == 1 && i_is_probabilistic(reachOptions.reachMethod)
                            i_log_cpstar(onnx);   % parity: parfor cp-star verdicts were previously unlogged
                        end
                    catch
                        tempStatus = 2;   % unknown; try the next reach option
                    end

                    if tempStatus ~= 2 % verified, then stop (or falsified)
                        break
                    else
                        reachOptPar = reachOptPar(2:end);
                    end

                    local_status(spc) = tempStatus;

                end
                local_status(spc) = tempStatus;

            end

        end

        % Check for the global verification result
        if all(local_status == 1)
            status = 1;
        else
            status = 2;
        end

    else
        warning("Working on adding support to other vnnlib properties")
    end

end

vT = toc(vT);



%% 4) Process results

tTime = toc(t); % save total computation time

% (removed) a dead, commented-out "status==2 && exact-star -> status=0" promotion used to live here.
% Deleted as a soundness trap: re-enabling it would turn any exact-star UNKNOWN (incl. a sound budget/cap
% abort) into a FALSE sat. Unknown must never be promoted to a verdict.

disp("Verification result: " + string(status));
disp("Counterexample search time: " + string(cEX_time));
disp("Reachability time: " + string(vT));
disp("Total Time: "+ string(tTime));
disp( " ");

% Write results to output file
if status == 0
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'sat \n');
    % fprintf(fid, '%f \n', tTime); % remove this line when running on submission site
    fclose(fid);
    write_counterexample(outputfile, counterEx)
elseif status == 1
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unsat \n');
    % fprintf(fid, '%f \n', tTime); % remove this line when running on submission site
    fclose(fid);
elseif status == 2
    fid = fopen(outputfile, 'w');
    fprintf(fid, 'unknown \n');
    % fprintf(fid, '%f \n', tTime); % remove this line when running on submission site
    fclose(fid);
end

end


%% Helper functions

% collins_aerospace_benchmark: delegate the whole instance to collins_falsify.py
% (sat-or-unknown ONLY; reach/UNSAT is out of scope for this category, see the
% dispatch comment at the top of run_vnncomp_instance). The script:
%   * parses the vnnlib (1.0 subset used by these specs) itself -- MATLAB never
%     touches the 1.2M-variable spec;
%   * validates the HWC-flat -> CHW input mapping with a center-consistency gate
%     (a 'violated' box center means the mapping is wrong -> it refuses, exit 2);
%   * falsifies with finite-difference gradients through onnxruntime and
%     RE-VALIDATES the exact candidate with >1e-6 margin before claiming SAT.
% Exit 10 -> sat + witness CSV (flat vnnlib order = HWC, which is exactly the
% order write_counterexample emits and the official checker replays); anything
% else -> unknown. NEVER unsat.
function [status, tTime] = run_collins_falsifier(onnx, vnnlib, outputfile, t)
    status = 2;
    script = fullfile(fileparts(mfilename('fullpath')), 'collins_falsify.py');
    witness = [tempname '.csv'];
    % Per-instance timeout is 3600s; leave slack for parsing + this wrapper.
    % NNV_COLLINS_BUDGET (seconds) overrides for tests/smoke runs.
    budget = str2double(getenv('NNV_COLLINS_BUDGET'));
    if isnan(budget) || budget <= 0
        budget = 3300;
    end
    cmd = sprintf('%s "%s" "%s" "%s" "%s" %g', python_exe(), script, ...
        char(onnx), char(vnnlib), witness, budget);
    [st, out] = system(cmd);
    disp(out);
    % Trust the SAT verdict ONLY on the exact contract: exit code 10 AND a witness
    % file AND the explicit SAT marker (guards against a shell mangling the code).
    if st == 10 && isfile(witness) && contains(out, 'SAT')
        v = readmatrix(witness);
        nx = round(v(1)); ny = round(v(2));
        if numel(v) == 2 + nx + ny
            x = v(3:2+nx);          % flat vnnlib order (HWC) -- written as X_i as-is
            y = v(3+nx:2+nx+ny);    % flat row-major [25200x11] = vnnlib Y order
            status = 0;
            fid = fopen(outputfile, 'w');
            fprintf(fid, 'sat \n');
            fclose(fid);
            write_counterexample(outputfile, {x, y});
        end
    end
    if status ~= 0   % anything else (incl. malformed witness) -> unknown, never unsat
        fid = fopen(outputfile, 'w');
        fprintf(fid, 'unknown \n');
        fclose(fid);
    end
    if isfile(witness), delete(witness); end
    tTime = toc(t);
end

% Interpreter resolution BASED ON validate_witness_onnx.m's idiom (pyenv first),
% but deliberately DIFFERENT in its PATH fallback: 'python3' on Linux (bare
% 'python' often does not exist on the competition VM; run_instance.sh itself
% uses python3), 'python' on Windows dev boxes.
function py = python_exe()
%PYTHON_EXE  Return a quoted python interpreter that can import the verification stack
%   (onnx + onnxruntime), discovered ROBUSTLY across machine configs -- a venv OR a
%   direct/system python (lambda dev box has the stack in ~/taylor_venv; install_tool.sh
%   provisions it into system python3 on the eval box). Used by the python verifiers that
%   need onnxruntime: the cctsdb_yolo complete-enumeration shell-out (cctsdb_enumerate.py)
%   and the witness-replay helpers. The old version returned bare python3/pyenv, which on a
%   box whose system python3 lacks onnx/ort made cctsdb_enumerate.py fail -> all-unknown
%   (the cctsdb_yolo 39->0 regression). Probe order: $NNV_ORT_PYTHON -> MATLAB pyenv ->
%   common venv -> python3 -> python; pick the first that imports onnx+onnxruntime; if none
%   has the stack, fall back to a plain interpreter so non-ort callers still run. Cached.
    persistent cached
    if ~isempty(cached), py = cached{1}; return; end
    cands = {};
    e = strtrim(getenv('NNV_ORT_PYTHON'));
    if ~isempty(e), cands{end+1} = e; end
    % Prefer the project venv (~/taylor_venv) BEFORE MATLAB's pyenv: the witness-replay gate
    % and the falsifiers (adaptive_cruise) must use the SAME onnxruntime BUILD the witnesses
    % were validated against. MATLAB pyenv can resolve a python with onnx+ort of a DIFFERENT
    % ort version that then REJECTS real witnesses (the 2026-06-24 adaptive_cruise all-unknown:
    % 50/50, incl. the validated instance_2). NNV_ORT_PYTHON still wins for an explicit
    % eval-box override; on the eval box (no taylor_venv) this falls through to pyenv/system
    % unchanged, so cctsdb_enumerate.py keeps its onnx+ort interpreter (no 39->0 regression).
    home = getenv('HOME');
    if ~isempty(home), cands{end+1} = fullfile(home, 'taylor_venv', 'bin', 'python'); end
    try
        pe = pyenv;
        if ~isempty(pe.Executable) && isfile(pe.Executable), cands{end+1} = char(pe.Executable); end
    catch
    end
    if ispc
        cands = [cands, {'python', 'python3'}];
    else
        cands = [cands, {'python3', 'python'}];
    end
    found = '';
    for i = 1:numel(cands)
        c = cands{i};
        if isempty(c), continue; end
        [st, ~] = system(sprintf('"%s" -c "import onnx, onnxruntime"', c));
        if st == 0, found = ['"' c '"']; break; end
    end
    if isempty(found)
        if ispc, found = 'python'; else, found = 'python3'; end   % no stack found -> plain fallback
    end
    cached = {found};
    py = found;
end

function IS = create_input_set(lb, ub, inputSize, needReshape, useImageStar)

    % Choose the set type from the NET's input-layer TYPE when the caller knows it
    % (useImageStar), not from the input SHAPE. A net imported as an ImageInputLayer
    % (e.g. acasxu, inputSize [1 5 1]) needs an ImageStar even though its spatial dims
    % are singleton; the old shape-only heuristic built a Star and EVERY acasxu reach
    % errored ("Input is not an ImageStar"), losing all robust/UNSAT acasxu verdicts.
    % Star is still used for flat feature / NN-manifest inputs, because some downstream
    % NNV layer reach() implementations read Set.dim, which Star has but ImageStar lacks.
    if nargin >= 5 && ~isempty(useImageStar)
        is_feature_input = ~useImageStar;
    else
        is_feature_input = isscalar(inputSize) || ...
            (numel(inputSize) <= 3 && nnz(inputSize > 1) <= 1);
    end
    if is_feature_input
        IS = Star(double(lb(:)), double(ub(:)));
        return;
    end

    % Image input: original behavior
    if ~isscalar(inputSize)
        lb = reshape(lb, inputSize);
        ub = reshape(ub, inputSize);
    end

    % Format bounds into correct dimensions
    if needReshape == 1
        lb = permute(lb, [2 1 3]);
        ub = permute(ub, [2 1 3]);
    elseif needReshape == 2
        newSize = [inputSize(2) inputSize(1) inputSize(3:end)];
        lb = reshape(lb, newSize);
        lb = permute(lb, [2 1 3 4]);
        ub = reshape(ub, newSize);
        ub = permute(ub, [2 1 3 4]);
    elseif needReshape == 3
        % ONNX row-major NHWC flat (C fastest, then W, then H): the reshape
        % above already produced [C W H] (inputSize must be given as
        % [C W H]); permute to NNV's [H W C] array convention.
        lb = permute(lb, [3 2 1]);
        ub = permute(ub, [3 2 1]);
    end

    % Create input set
    IS = ImageStar(lb, ub);

    % Delete constraints for non-interval dimensions
    try
        xxx = find((lb-ub)); % do this for now as it is easier, but it can get created using the (ImageStar(V,C,d,lb,ub) way)
        xxx = reshape(xxx, [], 1);
        if numel(lb) ~= length(xxx)
            IS.C = IS.C(:,xxx);
            IS.pred_lb = IS.pred_lb(xxx);
            IS.pred_ub = IS.pred_ub(xxx);
            xxx = xxx + 1;
            IS.V = IS.V(:,:,:,[1;xxx]);
            IS.numPred = length(xxx);
        end
    end

end

function ok = is_nnvnet_valid(nnvnet)
% Sanity-check that nnvnet was built (some dispatchers set it to "" when
% matlab2nnv conversion fails inside a try/catch). We need an NN object,
% not a string sentinel.
    ok = ~isempty(nnvnet) && ~ischar(nnvnet) && ~isstring(nnvnet);
end

function i_arm_reach_budget(bud)
% Arm (or clear) the SOUND per-instance reach time budget on the CURRENT MATLAB
% workspace by setting the globals PosLin.reach_star_exact / verify_specification read.
% Routed through a helper because `global` declarations are ILLEGAL inside a parfor
% body -- calling this from a parfor iteration sets the globals on that WORKER (where
% the reach actually runs), which a client-side assignment cannot do (globals do not
% cross the parfor boundary). Sound: an exceeded budget only ever aborts to unknown.
    global NNV_REACH_T0 NNV_REACH_BUD %#ok<GVMIS>
    if isfinite(bud) && bud > 0
        NNV_REACH_T0 = tic; NNV_REACH_BUD = bud;
    else
        NNV_REACH_T0 = []; NNV_REACH_BUD = [];   % no budget configured -> cap stays a no-op
    end
end

function i_log_cpstar(onnx)
% A cp-star (probabilistic, conformal) reach produced the accepted verdict -- NOT a sound proof.
% Logged in EVERY reach branch (single-spec AND both parfor branches) so the per-benchmark -150
% exposure from a probabilistic verdict is tracked uniformly (previously only the single-spec
% branch logged it, hiding parfor cp-star verdicts). i_is_probabilistic + i_finalize_reach_options
% live in their OWN files (same folder) so the routing logic is unit-testable -- see
% tests/nn/vnncomp/test_run_vnncomp_routing.m.
    fprintf('VERDICT VIA CP-STAR (probabilistic, not a sound proof): %s\n', onnx);
end

function [status, reachOptionsList] = i_gpu_bab_precheck(category, nnvnet, lb, ub, prop, status, reachOptionsList, needReshape, inputSize)
% Sound, ADDITIVE GPU-BaB unsat pre-check for argmax-robustness timeout categories: the FC nets
% (safenlp / sat_relu / relusplitter) AND the conv resnets (cifar100 / tinyimagenet). Runs the
% batched GPU-BaB (sound CROWN + ReLU-split BaB, DOUBLE precision -- FP-sound) BEFORE Star. A
% 'robust' verdict is a sound UNSAT proof -> emit it (status=1) and skip Star. Any other outcome
% (unknown / skip / unsafe / error) leaves status + reachOptionsList UNCHANGED, so Star runs
% exactly as before. gpu_bab_try_verify verifies the spec is argmax-robustness for the net's own
% center prediction (target derived + checked inside) AND requires the op-list evaluation to
% match net.evaluate (orientation guard), so it can only ADD a fast sound unsat -- never a wrong
% verdict (non-matching spec / unsupported net / failed guard -> 'skip' -> Star). For the CONV
% nets the flat vnnlib box is first remapped to the NET's input order via the SAME curated
% needReshape Star uses (so gpu-bab bounds the SAME function reach verifies); a wrong order still
% fails the guard -> skip, so soundness holds either way (validated: 10/10 cifar100, 0 contra).
    if status ~= 2, return; end                          % already decided -> nothing to do
    if nargin < 8 || isempty(needReshape), needReshape = 0; end
    if nargin < 9, inputSize = []; end
    isFC   = contains(category,"safenlp") || contains(category,"sat_relu") || contains(category,"relusplitter");
    % NOTE: challenging_certified is deliberately NOT added here. Routing it to the conv GPU-BaB
    % precheck CRASHED MATLAB on the box (hard std::exception "MATLAB process cannot be terminated" --
    % the orientation guard did NOT catch it, so it is NOT a safe sound-or-skip). Needs the conv
    % GPU-BaB path hardened against challenging's net before it can be enabled. Left on Star for now.
    % malbeware: the malimg CNNs (4-25 / 16-25) are plain ImageInput->Conv2D->ReLU->Flatten->FC
    % nets that nn_to_ops fully supports (probe: op-list builds for every flatten order, orientation
    % guard matches net.evaluate, gold-unsat instances certify 'robust' in 1 node on GPU-single AND
    % CPU-double). Routing malbeware through the SAME sound conv GPU-BaB precheck adds fast sound
    % unsat certs for the ~62 CNN instances that otherwise time out in exact-star. STRICTLY ADDITIVE:
    % a non-certifying instance leaves status/reachOptionsList unchanged -> the existing approx/exact-
    % star ladder still runs (the linear-25 FC net is decided there as before). The NNV_CONV_NO_STAR
    % "fast unknown" strip below is deliberately NOT applied to malbeware (unlike the resnets, whose
    % Star never certifies) so malbeware keeps its sound Star fallback -- see the guard at line ~1071.
    isConv = contains(category,"cifar100") || contains(category,"tinyimagenet") || contains(category,"vggnet") || contains(category,"malbeware");
    % General-halfspace control benchmarks (NOT argmax): prove the output avoids every unsafe
    % disjunct in prop.Hg. gpu_bab_halfspace_verify is sound-or-skip (FP64 CROWN + orientation
    % guard), so it can only ADD a sound unsat. Routed separately below (own predicate).
    % tllverifybench (TwoLevelLattice control nets, 2-D input, single-output >= c specs) is added
    % here too: its spec is a general halfspace, so the same sound-or-skip halfspace path applies.
    isHalfspace = contains(category,"cersyve") || contains(category,"lsnc_relu") || contains(category,"linearizenn") || contains(category,"tllverifybench");
    % INPUT-BISECTION BaB BACKSTOP categories: low-input-dim (<=4) control nets whose general-
    % halfspace spec the ROOT CROWN bound (gpu_bab_halfspace_verify) cannot decide, but the FP64
    % input-bisection CROWN BaB (gpu_bab_halfspace_input_bab, the acas-proven routine) CAN -- halving
    % a tiny input box quickly shrinks the relaxation gap (empirically: tll 0.6s/69 nodes, cersyve
    % 0.05-5.9s). Sound-or-unknown by construction (children partition the parent box): validated
    % 0/5 gold-SAT falsely certified (they hit the minWidth barrier -> 'unknown' -> upstream PGD
    % finds the sat), gold-UNSAT certified via partition-covering FP64. Backstop is gated to these.
    isInputBab = contains(category,"tllverifybench") || contains(category,"cersyve") || contains(category,"linearize");
    % ACAS-Xu (Phase 2): low-dim (5-input) pure-ReLU, general-halfspace specs. Certified via FP64
    % input-bisection CROWN BaB (gpu_bab_halfspace_input_bab). That routine targets product/bilinear
    % nets, but for the 5-D acas boxes input bisection DOES converge (validated: prop_1 8/8 certified,
    % 55-4473 nodes, 0.2-13.9s) and is sound-or-unknown by construction. Routed separately below.
    isAcas = contains(category,"acasxu");
    if ~(isFC || isConv || isHalfspace || isAcas), return; end     % only the gated timeout categories
    if ~is_nnvnet_valid(nnvnet), return; end             % need a valid NNV net for nn_to_ops + evaluate
    if isAcas
        % FP64-ONLY additive unsat precheck. SOUNDNESS: gpu_bab_halfspace_input_bab returns 'robust'
        % only when input bisection covers the box with FP64-certified-avoided sub-boxes (children
        % partition the parent box); anything else (unknown / barrier / timeCap / maxNodes) leaves
        % status + reachOptionsList UNCHANGED so Star runs as before. The op-list is built with the
        % SAME orientation guard (degenerate-IBP == nnvnet.evaluate) used by the other gpu-bab paths,
        % so a wrong flatten order -> guard fails -> skip (never bounds the wrong function). No FP32
        % is trusted (acas nets are tiny; FP64 runs in <15s). Falsification/SAT is handled upstream.
        try
            ops = i_acas_build_ops(nnvnet, lb, ub);       % [] if no flatten order matches the net
            if isempty(ops)
                fprintf('acas BaB pre-check: orientation guard skip -> Star reach\n');
                return;
            end
            % Default BaB timeCap raised 30 -> 60s (per-category acasxu, NOT per-instance): the 3 hardest
            % prop_1 nets (3_9/4_7/4_9) certify robust in 19-43s but were cut off at 30s -> raising to 60s
            % yields 45/45 prop_1 (solo/competition). Sound: more FP64 BaB time only adds certs, never a
            % wrong verdict. 60s leaves room under the ~116s per-instance budget (a cert skips Star).
            tcap = 60;     ev = str2double(getenv('NNV_ACAS_BAB_TIMECAP'));  if isfinite(ev) && ev > 0,  tcap = ev; end
            mnod = 300000; ev = str2double(getenv('NNV_ACAS_BAB_MAXNODES')); if isfinite(ev) && ev >= 1, mnod = ev; end
            babOpts = struct('precision','double','maxNodes',mnod,'timeCap',tcap);
            if iscell(lb)                                 % defensive: each input set must be robust (single box is the usual acas call path)
                allRobust = ~isempty(lb);
                for s = 1:numel(lb)
                    [Gd, gd] = i_acas_halfspaces(prop{min(s, numel(prop))});
                    if isempty(Gd), allRobust = false; break; end   % no halfspaces -> never vacuously certify; defer to Star
                    bt = tic; v = gpu_bab_halfspace_input_bab(ops, lb{s}, ub{s}, Gd, gd, babOpts);
                    if ~strcmp(v, 'robust'), allRobust = false; break; end
                end
                if allRobust
                    status = 1; reachOptionsList = {};
                    fprintf('acas BaB pre-check: robust/unsat (all %d input sets) -> skip Star\n', numel(lb));
                end
            else
                [Gd, gd] = i_acas_halfspaces(prop{1});
                if isempty(Gd), return; end   % no halfspaces -> never vacuously certify; defer to Star (sound)
                bt = tic; [v, info] = gpu_bab_halfspace_input_bab(ops, lb, ub, Gd, gd, babOpts);
                if strcmp(v, 'robust')
                    status = 1; reachOptionsList = {};
                    fprintf('acas BaB pre-check: robust/unsat (%d nodes, %.1fs) -> skip Star\n', info.nodes, toc(bt));
                else
                    fprintf('acas BaB pre-check: %s (%s, %.1fs) -> Star reach\n', v, info.reason, toc(bt));
                end
            end
        catch ME
            fprintf('acas BaB pre-check errored (%s) -> Star reach\n', ME.message);
        end
        return;                                           % acas path is terminal (no argmax fall-through)
    end
    if isHalfspace
        try
            if iscell(lb)                                % multiple input sets -> every set must be safe
                allRobust = ~isempty(lb);
                for s = 1:numel(lb)
                    pv = prop{min(s, numel(prop))};
                    hv = gpu_bab_halfspace_verify(nnvnet, lb{s}, ub{s}, pv);
                    if ~strcmp(hv, 'robust'), allRobust = false; break; end
                end
                if allRobust
                    status = 1; reachOptionsList = {};
                    fprintf('halfspace pre-check: robust/unsat (all %d input sets) -> skip Star\n', numel(lb));
                end
            else
                [hv, hinfo] = gpu_bab_halfspace_verify(nnvnet, lb, ub, prop{1});
                if strcmp(hv, 'robust')
                    status = 1; reachOptionsList = {};
                    fprintf('halfspace pre-check: robust/unsat (%d disjuncts) -> skip Star\n', hinfo.nDisjuncts);
                else
                    fprintf('halfspace pre-check: %s (%s) -> Star reach\n', hv, hinfo.reason);
                end
            end
        catch ME
            fprintf('halfspace pre-check errored (%s) -> Star reach\n', ME.message);
        end
        % INPUT-BISECTION BaB BACKSTOP (tll + cersyve): if the root-only CROWN bound above did not
        % certify (status still 2), try the acas-proven FP64 input-bisection CROWN BaB, which decides
        % the low-dim tll/cersyve boxes the single root bound cannot. STRICTLY ADDITIVE + sound: it
        % emits status=1 ONLY on a partition-covering 'robust' proof (every sub-box FP64-certified
        % avoided); anything else (unknown / barrier / [] / error) leaves status + reachOptionsList
        % UNCHANGED so the existing Star ladder runs exactly as before.
        if status == 2 && isInputBab
            [status, reachOptionsList] = i_halfspace_input_bab(nnvnet, lb, ub, prop, status, reachOptionsList);
        end
        return;                                          % halfspace path is terminal (no argmax fall-through)
    end
    if isConv && needReshape ~= 0 && ~isempty(inputSize) && ~isscalar(inputSize)
        try
            lb = i_remap_box_to_net(lb, inputSize, needReshape);
            ub = i_remap_box_to_net(ub, inputSize, needReshape);
        catch
            return;                                      % remap failed -> leave to Star (sound)
        end
    end
    % Conv resnets are ~195s/instance in CPU-double (the bottleneck) but ~10s on GPU-single. So
    % for conv use a 2-STAGE check: a fast GPU-SINGLE screen, then a CPU-DOUBLE re-confirm ONLY on
    % a 'robust' screen (rare -- conv certifies ~0). SOUNDNESS: the EMITTED unsat comes from the
    % DOUBLE-precision confirm; GPU-single is just a fast filter that can never emit a verdict.
    % FC nets are cheap in double -> single-stage CPU-double (maxNodes 5000). No GPU -> conv also
    % falls back to CPU-double (slow but sound).
    cMaxNodes = 2000; cFrontier = 32;                        % conv BaB budget (env-tunable: NNV_CONV_MAXNODES)
    % NOTE (2026-06-18): default raised 64 -> 2000. The S6 sweep diagnostic showed certified cifar
    % instances need 9-1039 BaB nodes (the FP64 confirm), so the old default of 64 silently truncated
    % most certs to 'unknown'. 2000 covers the observed cert depth with margin. SOUND: a node cap only
    % turns a verdict into 'unknown' (early stop), never a wrong verdict; the cost is wall-clock on
    % non-certs (bounded by the per-instance cap). Lower it via the env for fast dev/smoke runs.
    % env overrides only when finite + sane (str2double of unset/non-numeric is NaN -> keep the
    % default; a NaN budget would make info.nodes>maxNodes false forever and silently un-bound the BaB)
    ev = str2double(getenv('NNV_CONV_MAXNODES')); if isfinite(ev) && ev >= 1, cMaxNodes = ev; end
    ev = str2double(getenv('NNV_CONV_FRONTIER')); if isfinite(ev) && ev >= 1, cFrontier = ev; end
    try
        if isConv && gpuDeviceCount >= 1
            % CONV: an optional fast GPU-single SCREEN, then the sound CPU-DOUBLE pass that decides
            % the emit (single is never trusted). The screen filters cheaply, but DOUBLE bounds are
            % tighter (no FP32 rounding loosening) so the double pass crosses the convex barrier in
            % far fewer nodes -- on robust conv the screen grinds the launch-bound tail much longer
            % than the double pass takes. NNV_CONV_GPU_SCREEN=0 skips the screen -> straight to double.
            % NNV_CONV_TRUST_FP32 emits FROM the GPU-single screen, so it FORCES the screen on (else
            % NNV_CONV_GPU_SCREEN=0 from FIX A1 would skip the screen -> straight to the FP64-CPU
            % confirm and trust-FP32 could never fire / the GPU would stay idle).
            useScreen = ~isequal(getenv('NNV_CONV_GPU_SCREEN'), '0') || ~isempty(getenv('NNV_CONV_TRUST_FP32'));
            soundEmit = ~isempty(getenv('NNV_SOUND_FP32_TIGHT'));     % sound-FP32 fast-emit attempt enabled
            screenPass = true;
            if useScreen
                % FAST UNSOUND-FP32 SCREEN: the cheap, TIGHT candidate filter. soundFP32=false FORCES it
                % unsound even when NNV_SOUND_FP32_TIGHT is set, so a HARD cert (which the widened sound
                % bound cannot certify in budget) still passes the filter -> reaches the FP64 confirm.
                [gv, ginfo] = gpu_bab_try_verify(nnvnet, lb, ub, prop, ...
                    struct('engine','batched','maxNodes',cMaxNodes,'maxFrontier',cFrontier,'device','gpu','allowUnsoundSingle',true,'soundFP32',false));
                screenPass = strcmp(gv, 'robust');
                if ~screenPass
                    fprintf('GPU-BaB pre-check: %s (gpu-screen, %s) -> Star reach\n', gv, ginfo.reason);
                end
            end
            if screenPass
                % TRUST-FP32 EMIT (abCROWN/NeuralSAT model -- opt-in NNV_CONV_TRUST_FP32, default OFF):
                % the GPU-SINGLE batched BaB above (the screen) certified 'robust' from the raw CROWN
                % bound run ENTIRELY on the GPU. Falsification (PGD) already ran FALSIFY-FIRST -- this
                % precheck is only reached when NO counterexample was found -- so a screen-robust that is
                % ALSO PGD-clean is the SAME two-mechanism soundness the GPU-winning tools rely on (an
                % FP32 GPU lower bound + an independent adversarial attack; abCROWN runs stock FP32 with
                % PGD, NeuralSAT likewise). EMIT directly here, skipping the ~200s FP64-CPU reconfirm that
                % leaves the GPU idle and caps the conv decide rate. SOUNDNESS POLICY: this TRUSTS FP32
                % rounding (~1e-6, negligible vs real robustness margins) instead of the rigorous FP64
                % confirm; gated OFF by default and to be validated 0 false-robust vs the alpha-beta-CROWN
                % gold set before competition use. PGD is the backstop for the residual FP32 gap.
                if useScreen && ~isempty(getenv('NNV_CONV_TRUST_FP32'))
                    status = 1; reachOptionsList = {};
                    fprintf('GPU-BaB pre-check: robust/unsat (TRUSTED gpu-single screen, %d nodes; PGD falsify-first clean) -> skip Star\n', ginfo.nodes);
                    return;
                end
                % SOUND-FP32 EMIT (M3b): on a screen-robust candidate, try the SOUND-FP32 BaB FIRST -- a
                % 'robust' carrying soundFP32=true is a provably-sound unsat (every CROWN bound outward-
                % widened; frontier-G1 0/1584 + root-G2 0/200 + 2-round adversarial review) -> emit
                % directly, skipping the ~284s FP64 confirm. If it returns unknown (the widened bound is
                % too loose for this hard instance), FALL BACK to FP64 below -- so the emit only ADDS fast
                % certs, never loses the ones FP64 gets. Gated on NNV_SOUND_FP32_TIGHT (default off).
                if soundEmit
                    [gvE, giE] = gpu_bab_try_verify(nnvnet, lb, ub, prop, ...
                        struct('engine','batched','maxNodes',cMaxNodes,'maxFrontier',cFrontier,'device','gpu','allowUnsoundSingle',true,'soundFP32',true));
                    if strcmp(gvE,'robust') && isfield(giE,'soundFP32') && giE.soundFP32
                        status = 1; reachOptionsList = {};
                        fprintf('GPU-BaB pre-check: robust/unsat (sound-FP32 emit, %d nodes) -> skip Star\n', giE.nodes);
                        return;
                    end
                end
                [gv2, gi2] = gpu_bab_try_verify(nnvnet, lb, ub, prop, struct('engine','batched','maxNodes',cMaxNodes,'maxFrontier',cFrontier));  % CPU double = sound emit
                if strcmp(gv2, 'robust')
                    status = 1; reachOptionsList = {};
                    if useScreen, src = 'gpu-screen + '; else, src = ''; end
                    fprintf('GPU-BaB pre-check: robust/unsat (%s%d-node double) -> skip Star\n', src, gi2.nodes);
                else
                    fprintf('GPU-BaB pre-check: double=%s -> Star reach\n', gv2);
                end
            end
            % MEASUREMENT: cifar/tinyimagenet Star reach never certifies (GPU-BaB exists precisely
            % because Star cannot do these resnets) and burns the whole per-instance cap. With
            % NNV_CONV_NO_STAR set, a non-certifying conv pre-check emits a FAST sound 'unknown'
            % instead of grinding Star to the cap. SOUND: 'unknown' is always safe, and the
            % pre-check's unsat certs are unchanged. Default-OFF -> no behaviour change.
            % malbeware is EXCLUDED from the strip: unlike the resnets (cifar100/tinyimagenet/vggnet,
            % whose Star provably never certifies), malbeware's approx/exact-star ladder DOES decide
            % instances (esp. the linear-25 FC net). Stripping it would regress the 88 baseline. So a
            % non-certified malbeware instance keeps its sound Star fallback (strictly additive).
            if status == 2 && ~isempty(getenv('NNV_CONV_NO_STAR')) && ~contains(category,"malbeware")
                reachOptionsList = {};
                fprintf('GPU-BaB pre-check: conv not certified + NNV_CONV_NO_STAR -> skip Star (fast unknown)\n');
            end
        else
            mn = 5000; if isConv, mn = 64; end                 % conv-without-GPU: bounded (slow but sound)
            [gv, ginfo] = gpu_bab_try_verify(nnvnet, lb, ub, prop, struct('engine','batched','maxNodes',mn));
            if strcmp(gv, 'robust')
                status = 1; reachOptionsList = {};
                fprintf('GPU-BaB pre-check: robust/unsat (%d nodes, %s) -> skip Star reach\n', ginfo.nodes, ginfo.reason);
            else
                fprintf('GPU-BaB pre-check: %s (%s) -> Star reach\n', gv, ginfo.reason);
            end
        end
    catch ME
        fprintf('GPU-BaB pre-check errored (%s) -> Star reach\n', ME.message);
    end
end

function ops = i_acas_build_ops(nnvnet, lb, ub)
% Build the gpu_bab op-list for a low-dim acas net, choosing the flatten order whose degenerate-IBP
% matches nnvnet.evaluate at several sample points (the orientation guard shared with the other gpu-bab
% paths). Returns [] if no order matches -> caller skips to Star (sound). A wrong order produces grossly
% different outputs and fails the guard, so the op-list that passes bounds the SAME function reach verifies.
    if iscell(lb), lb1 = lb{1}; ub1 = ub{1}; else, lb1 = lb; ub1 = ub; end
    lb1 = double(lb1(:)); ub1 = double(ub1(:)); n = numel(lb1);
    idx = (1:n)';                                          % center + 3 low-discrepancy interior points
    f1 = mod(idx*0.6180339887498949, 1);
    f2 = mod(idx*1.3247179572447460 + 0.37, 1);
    f3 = mod(idx*0.7548776662466927 + 0.11, 1);
    pb = [(lb1+ub1)/2, lb1+f1.*(ub1-lb1), lb1+f2.*(ub1-lb1), lb1+f3.*(ub1-lb1)];
    orders = {'colmajor', 'chw_rowmajor', 'hwc_rowmajor'}; ops = [];
    for oi = 1:numel(orders)
        try, cand = nn_to_ops(nnvnet, orders{oi}, n); catch, continue; end
        ok = true;
        for pp = 1:size(pb, 2)
            cp = pb(:, pp);
            yo = gpu_bab_ibp(cand, cp, cp, 'double'); yo = yo(:);
            yn = nnvnet.evaluate(reshape(cp, [n 1]));  yn = yn(:);
            if numel(yo) ~= numel(yn) || max(abs(yo - yn)) > 1e-4*max(1, max(abs(yn))), ok = false; break; end
        end
        if ok, ops = cand; return; end
    end
end

function [Gd, gd] = i_acas_halfspaces(pv)
% Extract the unsafe-region halfspace disjuncts (G, g) from a loaded vnnlib prop for gpu_bab_*.
    Hg = pv.Hg; Gd = cell(1, numel(Hg)); gd = cell(1, numel(Hg));
    for i = 1:numel(Hg), Gd{i} = double(Hg(i).G); gd{i} = double(Hg(i).g(:)); end
end

function [status, reachOptionsList] = i_halfspace_input_bab(nnvnet, lb, ub, prop, status, reachOptionsList)
% FP64 input-bisection CROWN BaB backstop for low-input-dim general-halfspace control nets
% (tllverifybench + cersyve). Mirrors the acas BaB path: an orientation-guarded op-list
% (i_acas_build_ops) + gpu_bab_halfspace_input_bab over the unsafe-region halfspaces
% (i_acas_halfspaces). SOUND/ADDITIVE: returns status=1 ONLY on a 'robust' verdict (the input
% box is fully covered by FP64-certified-avoided sub-boxes -- the two children of every bisection
% PARTITION the parent box, so 'robust' is a sound UNSAT proof); any other outcome (unknown /
% barrier / no-halfspaces / [] op-list / error) leaves status + reachOptionsList UNCHANGED so the
% existing Star ladder runs exactly as before. A wrong flatten order fails the shared orientation
% guard -> [] -> skip, so it never bounds the wrong function. Validated (FP64): gold-UNSAT certified
% (tll 69 nodes/0.6s, cersyve 159-1805 nodes); gold-SAT correctly NOT certified (minWidth barrier
% -> 'unknown', upstream PGD then emits the sat). No FP32 is trusted (these nets are tiny).
    if status ~= 2, return; end
    try
        ops = i_acas_build_ops(nnvnet, lb, ub);          % [] if no flatten order matches the net
        if isempty(ops)
            fprintf('halfspace input-bisection BaB: orientation guard skip -> Star reach\n');
            return;
        end
        % Same FP64 BaB budget as the acas path (env-tunable via the shared NNV_ACAS_BAB_* knobs,
        % which execute.py scales to the official per-instance timeout). 60s leaves room for Star.
        tcap = 60;     ev = str2double(getenv('NNV_ACAS_BAB_TIMECAP'));  if isfinite(ev) && ev > 0,  tcap = ev; end
        mnod = 300000; ev = str2double(getenv('NNV_ACAS_BAB_MAXNODES')); if isfinite(ev) && ev >= 1, mnod = ev; end
        babOpts = struct('precision','double','maxNodes',mnod,'timeCap',tcap);
        if iscell(lb)                                    % defensive: every input set must be robust
            allRobust = ~isempty(lb);
            for s = 1:numel(lb)
                [Gd, gd] = i_acas_halfspaces(prop{min(s, numel(prop))});
                if isempty(Gd), allRobust = false; break; end   % no halfspaces -> never vacuously certify
                v = gpu_bab_halfspace_input_bab(ops, lb{s}, ub{s}, Gd, gd, babOpts);
                if ~strcmp(v, 'robust'), allRobust = false; break; end
            end
            if allRobust
                status = 1; reachOptionsList = {};
                fprintf('halfspace input-bisection BaB: robust/unsat (all %d input sets) -> skip Star\n', numel(lb));
            end
        else
            [Gd, gd] = i_acas_halfspaces(prop{1});
            if isempty(Gd), return; end                  % no halfspaces -> defer to Star (sound)
            bt = tic; [v, info] = gpu_bab_halfspace_input_bab(ops, lb, ub, Gd, gd, babOpts);
            if strcmp(v, 'robust')
                status = 1; reachOptionsList = {};
                fprintf('halfspace input-bisection BaB: robust/unsat (%d nodes, %.1fs) -> skip Star\n', info.nodes, toc(bt));
            else
                fprintf('halfspace input-bisection BaB: %s (%s, %.1fs) -> Star reach\n', v, info.reason, toc(bt));
            end
        end
    catch ME
        fprintf('halfspace input-bisection BaB errored (%s) -> Star reach\n', ME.message);
    end
end

function v = i_remap_box_to_net(x, inputSize, needReshape)
% Map a flat vnnlib box vector to the NET's input order (returned flat, column-major), using the
% SAME reshape+permute create_input_set applies for Star. gpu_bab_try_verify reshapes the flat
% box column-major into the net's [H W C], so feeding it this net-ordered box makes it bound the
% correctly-oriented image. Mirrors create_input_set's needReshape cases EXACTLY (1/2/3).
    x = double(x(:));
    if needReshape == 2
        img = reshape(x, [inputSize(2) inputSize(1) inputSize(3:end)]);
        img = permute(img, [2 1 3 4]);
    else
        img = reshape(x, inputSize);
        if needReshape == 1
            img = permute(img, [2 1 3]);
        elseif needReshape == 3
            img = permute(img, [3 2 1]);
        end
    end
    v = img(:);
end

function L = relax_ladder(ks)
% Build an ordered cell array of sound relax-star-area reachOptions.
%   ks : vector of relaxFactor values, ordered cheapest(high k) -> tightest(low k).
% Every k in [0,1] is a SOUND over-approximation (PosLin.reach_relaxed_star_area only
% REPLACES an LP-tight bound with an estimateRanges superset, never tightens below the
% true reachable set), so the whole ladder is sound-or-unknown regardless of order; the
% verify loop takes the FIRST rung that returns a decisive status. Used to (a) lead with
% a fast high-k pass so reach COMPLETES within budget (the verify loop has no per-method
% timeout, so a timeout kills the instance before later rungs run) and (b) follow with
% tighter rungs that can certify loose-unknowns.
    L = cell(1, numel(ks));
    for i = 1:numel(ks)
        o = struct();
        o.reachMethod = 'relax-star-area';
        o.relaxFactor = ks(i);
        L{i} = o;
    end
end

function [net,nnvnet,needReshape,reachOptionsList,inputSize,inputFormat,nRand,falsifyOpts] = load_vnncomp_network(category, onnx, vnnlib)
% load participating vnncomp 2025 benchmark NNs
% Not yet supported:
% - cctsdb (some errrors when forward propagating)
% - lsnc_relu
% - traffic_signs_recognition (last year all instances were sat, maybe we are not wrong?)
% Handled OUTSIDE this dispatcher (early routes at the top of run_vnncomp_instance):
% - cctsdb_yolo (complete-enumeration path)
% - collins aerospace (the invalid SAT instances were a CHW/HWC input-order mix-up;
%   now routed to collins_falsify.py and never reaches this function)


    needReshape = 0; % default is to use MATLAB reshape, otherwise use the python reshape
    % reachOptions = struct;
    % reachOptions.reachMethod = 'approx-star'; % default parameters
    numCores = feature('numcores'); % in case we select exact method
    inputSize = [];
    inputFormat = "default"; % no need to change for most of them, but needed for some ("UU")
    nRand = 100; % default from previous competitions
    % Per-category falsification budget (VNNCOMP2026 tuning). Default raised to 15s (was 5s)
    % so an unlisted SAT-capable category still gets a meaningful PGD window before the sound
    % reach ladder runs; the ftab at the END of this function overrides per category. Sound:
    % falsification is SAT-or-unknown only (witness gated by validate_witness + onnxruntime
    % replay), so a larger PGD budget can never manufacture a false verdict -- it only reclaims
    % wall-clock that reach was otherwise spending timing out.
    falsifyOpts = struct('n_restarts',20, 'n_steps',40, 'lr',0.1, 'fgsm',true, ...
                         'max_time',15, 'seed',0);

    if contains(category, "acasxu")
        % acasxu: MUST import as BCSS -- BC (FeatureInputLayer) makes matlab2nnv fail with
        % "Unsupported Class of Layer" (tested on the box). BCSS yields an ImageInputLayer that
        % matlab2nnv converts, but create_input_set then builds an ImageStar whose per-layer spatial
        % overhead makes exact-star slow on these props. NOTE: acasxu is NOT a cheap win -- BCSS exact
        % times out on the hard props even multi-core, and approx-star (DeepZ-level) is too loose to
        % decide them. The approx-first + multi-core ladder below is SOUND and helps any approx-
        % decidable prop, but the CROWN-needing props need a tighter method (separate, deeper work).
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS");
        nnvnet = matlab2nnv(net);
        % SOUND-FIRST ladder for ALL props: fast approx-star (tightest non-exact; ~1-2s on these
        % tiny 5->5 nets, settles the robust/UNSAT props) then MULTI-CORE exact-star as the complete
        % closer, run only if approx returns unknown. Fixes two regressions: (1) prop_1/2/5-10
        % previously ran exact-star ONLY (no approx rung) -> timed out the easy props; (2) numCores
        % was pinned to 1 -> single-core exact-star blew past the 116s budget on the hard props.
        % Restored to feature('numcores') (the 2024 strong config). SAFE: NN.start_pool returns early
        % inside a parfeval worker (NN.m:1373 getCurrentTask) and the official harness runs ONE
        % instance per process, so no nested parpool. approx-star is sound (no-intersection = proof);
        % no cp-star; PGD (nRand) still finds SAT. (relax-star rungs omitted: relax is LOOSER than
        % approx, so on these tiny nets it can never decide what approx couldn't.)
        reachOptions = struct; reachOptions.reachMethod = 'approx-star';
        reachOptionsList{1} = reachOptions;
        reachOptions = struct; reachOptions.reachMethod = 'exact-star';
        % acasxu nets are tiny (5->5): parpool gives NO speedup (1 input set) and is BROKEN on some boxes
        % (worker status-1 crash -> exact-star errored->unknown). Run serial; speed comes from the prefilter.
        reachOptions.numCores = 1;
        reachOptionsList{2} = reachOptions;
        nRand = 500;

    elseif contains(category, "adaptive_cruise")
        % adaptive_cruise_control_non_linear_2026: a pure FFNN (8x Gemm + 5x Relu,
        % [1,2]->[1,1]) that imports cleanly; it only erred because no dispatcher
        % branch matched and control fell to the catch-all "ONNX model not
        % supported". Route it through the standard FFNN reach path.
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BC");
        nnvnet = matlab2nnv(net);
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star';
        reachOptionsList{1} = reachOptions;
        reachOptions = struct;
        reachOptions.reachMethod = 'exact-star'; % tiny net -> sound+complete fallback
        reachOptions.numCores = 1;
        reachOptionsList{2} = reachOptions;

    elseif contains(category, "cctsdb_yolo")
        % Handled by the complete-enumeration path at the TOP of
        % run_vnncomp_instance (verify_cctsdb_enumeration -> cctsdb_enumerate.py);
        % control never reaches this dispatcher for cctsdb_yolo. Fail loud if it
        % somehow does: MATLAB's importer mis-handles this model (Cast/Gather
        % truncation on a flat 12296 input).
        error("cctsdb_yolo is handled by the enumeration path in run_vnncomp_instance; load_vnncomp_network should not be reached for it");

    elseif contains(category, "cersyve")
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BC");
        nnvnet = matlab2nnv(net);
        % SOUND-FIRST ladder (was cp-star only). cp-star is probabilistic, so its
        % 'holds'/unsat is not a proof -> run the SOUND methods first and keep cp-star
        % as a LOGGED last-resort. approx-star/relax-star soundly decide the 'con'
        % instances; cp-star fires only if both return unknown.
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star';
        reachOptionsList{1} = reachOptions;
        reachOptions = struct;
        reachOptions.reachMethod = 'relax-star-area';
        reachOptions.relaxFactor = 0.5;
        reachOptionsList{2} = reachOptions;
        reachOptions = struct;
        reachOptions.reachMethod = 'cp-star';   % probabilistic last-resort (logged at verify)
        reachOptionsList{3} = reachOptions;

    elseif contains(category, "cgan")
        % cgan: route through the Python-importer manifest. MATLAB's matlab2nnv fails on
        % the custom ReshapeLayer (shape in .Vars, no ONNXParams). The manifest forward
        % pass is cross-validated vs onnxruntime (2026-06-11: 6/7 variants xval < 5e-7).
        % The _upsample variant uses an unsupported ONNX Resize op, so the importer REFUSES
        % identity evaluation -> that instance errors -> unknown (never unsound). FeatureInput [5].
        if ~contains(onnx, 'transformer')
            nnvnet = load_manifest_net(onnx);
            net = nnvnet;   % falsify_single dispatches NN.evaluate for NNV nets
            inputSize = nnvnet.Layers{1}.InputSize;
            % cgan models are GENERATORS with SAT-expected specs; approx-star reach
            % NEVER decides them (0 unsat across all 57 sweep instances, and the one
            % instance whose reach completed at 69.5 s returned unknown) -- it only
            % burns the whole per-instance budget timing out. Drop reach so the budget
            % goes to PGD falsification (where cgan's SAT verdicts actually come from,
            % boosted in the ftab below): a no-counterexample instance now returns
            % unknown after the falsification budget (ftab max_time, ~25 s) instead of
            % a much longer reach timeout. Sound: unknown is always safe.
            % (The 'transformer' cgan variant in the else-branch is a different model
            % and method path -- cp-star, not approx-star -- so it is left unchanged
            % here; the dead-reach finding was specifically the approx-star path.)
            reachOptionsList = {};
        else
            net = importNetworkFromONNX(onnx,"InputDataFormats","BC");
            nnvnet = "";
            reachOptions = struct;
            reachOptions.train_epochs = 100;
            reachOptions.train_lr = 0.001;
            reachOptions.dims = [-1 -1];
            reachOptions.coverage = 0.999;
            reachOptions.confidence = 0.999;
            reachOptions.train_mode = 'Linear';
            reachOptions.surrogate_dim = [10, 10];
            reachOptions.threshold_normal = 1e-5;
            reachOptions.dlarrayType = 'CB';
            reachOptions.reachMethod = "cp-star";
            reachOptionsList{1} = reachOptions;
        end


    elseif contains(category, "soundnessbench")
        % soundnessbench: MATLAB's importer produces a CustomInputLayerMultiOutput that NNV
        % cannot handle; route through the Python-importer manifest instead. Forward pass
        % cross-validated vs onnxruntime (2026-06-11: xval 7.7e-6). FeatureInput [128].
        % This benchmark is purpose-built to catch UNSOUND verifiers: NNV's sound
        % over-approximate reach yields unknown-or-correct, never an unsound verdict, and
        % any SAT witness is replayed through onnxruntime before being emitted.
        nnvnet = load_manifest_net(onnx);
        net = nnvnet;   % falsify_single dispatches NN.evaluate for NNV nets
        inputSize = nnvnet.Layers{1}.InputSize;
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star';
        reachOptionsList{1} = reachOptions;


    elseif contains(category, "cifar100")
        % cifar100: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        needReshape = 1;
        reachOptions = struct;
        reachOptions.reachMethod = 'cp-star';
        reachOptionsList{1} = reachOptions;

    elseif contains(category, "challenging")
        % challenging_certified_training_2026: cifar10 CNNs (Conv/BN/ReLU/FC),
        % same BCSS image import + CHW-flat vnnlib order as cifar100 -> same
        % needReshape. Center-image output cross-checked vs onnxruntime
        % (2026-06-12, max|diff| < 1e-5) so the input-box orientation is right.
        % Deterministic sound ladder (NOT probabilistic cp-star): approx-star
        % first, looser relax-star fallback; anything unresolved -> unknown.
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        needReshape = 1;
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star';
        reachOptionsList{1} = reachOptions;
        reachOptions = struct;
        reachOptions.reachMethod = 'relax-star-area';
        reachOptions.relaxFactor = 0.7;
        reachOptionsList{2} = reachOptions;

    elseif contains(category, 'collins_aerospace_benchmark')
        % Routed EARLY in run_vnncomp_instance to the sound falsification-only
        % path (collins_falsify.py): MATLAB's importer cannot handle the YOLOv5
        % Detect-head custom layers, and the old import+matlab2nnv route emitted
        % invalid SAT instances (vnnlib X order is HWC flat, not CHW). This branch
        % is unreachable; fail loud if anything ever re-enters it.
        error('run_vnncomp_instance:collinsRouted', ...
            'collins_aerospace_benchmark is handled by the collins_falsify.py path');

    elseif contains(category, 'collins_rul')
        net = importNetworkFromONNX(onnx);
        nnvnet = matlab2nnv(net);
        needReshape = 2;
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star'; % default parameters
        reachOptionsList{1} = reachOptions;

    elseif contains(category, "cora")
        net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        if contains(onnx, '-set')
            reachOptionsList = relax_ladder([0.5, 0.2]);
            oo = struct(); oo.reachMethod = 'approx-star';
            reachOptionsList{end+1} = oo;
        else
            % Adaptive relax-star ladder for the 133 cora timeouts: the verify loop has
            % NO per-method timeout, so a slow method kills the instance before any later
            % rung runs -> LEAD with the fastest (high-k) pass so reach COMPLETES in
            % budget, then tighter rungs to certify loose-unknowns. Every k is a sound
            % over-approx (PosLin estimateRanges); slow_cats still prepends abs-dom.
            reachOptionsList = relax_ladder([0.95, 0.85, 0.5]);
        end
        nRand = 500;

    elseif contains(category, "dist_shift")
        % dist_shift: onnx to matlab, , matlab to nnv?
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape
        try 
            nnvnet = matlab2nnv(net);
        catch
            nnvnet = "";
        end
        % Cheap-to-precise SOUND ladder: approx-star, looser relax-star, then a GATED
        % exact-star closer (added below). (Was approx-star at {1} then exact-star
        % OVERWRITING {1}, so only the exponential exact-star ran and stalled to
        % 'unknown'. exact-star is re-added as {3} -- now tractable WITH predicate
        % contraction -- running only after approx/relax both return unknown.)
        % approx-star -> relax-star tightening ladder (cheap, for the 36 unknowns) ->
        % GATED exact-star closer LAST (only if everything above returns unknown). Each
        % rung is a sound over-approx; the loop takes the first decisive status.
        reachOptionsList = {};
        o1 = struct(); o1.reachMethod = 'approx-star'; reachOptionsList{1} = o1;
        reachOptionsList = [reachOptionsList, relax_ladder([0.7, 0.4, 0.15])];
        % exact-star + predicate contraction (acasxu Step 1+2) as the precise closer --
        % tractable WITH contraction on these small nets (validated: dist_shift probe ->
        % sound unsat @120s); runs only after approx-star and the relax ladder return unknown.
        oe = struct(); oe.reachMethod = 'exact-star'; oe.numCores = 1;
        reachOptionsList{end+1} = oe;

    elseif contains(category, "linearize")
        % matlab2nnv cannot parse linearizenn's custom static SliceLayer (AllInOne_*.SliceLayer ->
        % "Unsupported Class of Layer"), so the old try/catch fell through to cp-star (probabilistic)
        % on a GENERAL-HALFSPACE control benchmark -- a -150 exposure (a cp-star 'unsat' is not a proof).
        % Route through the Python importer instead: onnx2nnv.py lowers the STATIC Slice to an EXACT
        % sparse-selector FullyConnectedLayer (y = W_sel*x, b=0), so nnvnet is a valid NN and the SOUND
        % approx-star/relax-star ladder (and the gpu_bab halfspace pre-check) run -- cp-star is gone.
        % The manifest forward pass is cross-validated vs onnxruntime at generation; a missing/failed
        % manifest degrades to 'unknown' via the load guard (sound). The .nnv.mat is built untimed at
        % prepare (prepare_instance.sh gen_manifest *linearize*).
        try
            nnvnet = load_manifest_net(onnx);
            net = nnvnet;   % falsify_single dispatches NN.evaluate for NNV nets
            inputSize = nnvnet.Layers{1}.InputSize;
            % sound tightness ladder: relax 0.5 -> 0.25 -> approx-star (tightest). First decisive wins.
            reachOptionsList = relax_ladder([0.5, 0.25]);
            oo = struct(); oo.reachMethod = 'approx-star';
            reachOptionsList{end+1} = oo;
        catch
            % Manifest unavailable (not generated at prepare) -> FALSIFICATION-ONLY: keep the raw net
            % for PGD (sat-or-unknown), never cp-star. Sound, and no regression vs the PGD sat path.
            net = importNetworkFromONNX(onnx, "InputDataFormats", "BC", 'OutputDataFormats',"BC");
            nnvnet = "";
            reachOptionsList = {};
        end

    elseif contains(category, "lsnc_relu")
        % lsnc_relu: flat [6] feature input, [8] output. SPLIT net usage:
        %   REACH (UNSAT side) -> the Python-importer manifest (nnvnet), cross-validated
        %     vs onnxruntime (max diff 9.2e-07, 2026-06-09); reach is byte-identical to
        %     before -> no unsat regression.
        %   FALSIFY (SAT side) -> a MATLAB-imported dlnetwork (net) so the acasxu-proven
        %     batched-random falsifier + autodiff PGD can run. The manifest's one-at-a-time
        %     numerical-gradient PGD MISSED the gold-sat that a huge batched dlnetwork
        %     forward finds (probe: state_1 violated at margin -3.2e-3). Falsify uses
        %     `net`, reach uses `nnvnet` (confirmed: reach=nnvnet, falsify=net).
        % SOUNDNESS: every lsnc SAT witness is onnxruntime-replayed (riskyNet, below) on
        % the ORIGINAL onnx before `sat` stands -> SAT-or-unknown, independent of any
        % dlnetwork-vs-ONNX import divergence. Import is best-effort: on failure fall back
        % to the manifest-as-net (prior numerical-gradient behaviour, still sound).
        % The older "importNetworkFromONNX cannot parse this model" note is STALE: the
        % import succeeds + is netcached; the manifest is kept only as the reach net.
        nnvnet = load_manifest_net(onnx);
        inputSize = 6;
        inputFormat = "CB";   % FeatureInputLayer wants channel-first: data [6 N] = 'CB' (NOT 'BC')
        try
            net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");
        catch
            net = nnvnet;        % fall back: falsify_single dispatches NN.evaluate (numerical gradient)
            inputFormat = "default";
        end
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star';
        reachOptionsList{1} = reachOptions;

    elseif contains(category, "malbeware")
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BCSS");
        nnvnet = matlab2nnv(net);
        reachOptions = struct;                    % initialize before field assignment
        reachOptions.reachMethod = 'approx-star'; % cheap sound method first
        reachOptionsList{1} = reachOptions;
        reachOptions = struct;
        reachOptions.reachMethod = 'exact-star';   % exact UNSAT fallback (was
        reachOptions.numCores = 1;                 % overwriting {1} so approx never ran)
        reachOptionsList{2} = reachOptions;
        needReshape = 2;

    elseif contains(category, "metaroom")
        % metaroom: onnx to matlab
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        needReshape = 2;
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star'; % default parameters
        reachOptionsList{1} = reachOptions;
        % exact-star + predicate contraction (acasxu Step 1+2) as the precise closer.
        % metaroom 4cnn/6cnn are tiny; WITH contraction exact-star decides the
        % approx-star 'unknown's within budget (validated: ~12/14 probe -> sound unsat
        % @7-26s; master timed out). Runs only after approx-star returns unknown.
        reachOptions = struct;
        reachOptions.reachMethod = 'exact-star';
        reachOptions.numCores = 1;
        reachOptionsList{2} = reachOptions;

    elseif contains(category, "ml4acopf")
        % ml4acopf: matlab2nnv cannot parse the ACOPF net (custom fused Slice_To_MatMul/Mul/Transpose +
        % Gather_To_Sub). Route through the Python importer manifest (onnx2nnv.py lowers those static
        % ops); the net's nonlinear ElementwiseProduct layers ARE soundly OVER-APPROXIMATED by
        % approx-star (validated: approx-star reaches in ~8s with no refusal), so run the SOUND
        % approx-star reach -- an upgrade from the prior falsification-only/cp-star (cp-star was a -150
        % exposure). PGD still provides SAT. Robust: if the manifest is unavailable, fall back to
        % FALSIFICATION-ONLY (sat-or-unknown), never cp-star. Manifest built untimed at prepare
        % (prepare_instance.sh gen_manifest *ml4acopf*).
        try
            nnvnet = load_manifest_net(onnx);
            net = nnvnet;   % falsify_single dispatches NN.evaluate for NNV nets
            inputSize = nnvnet.Layers{1}.InputSize;
            reachOptions = struct; reachOptions.reachMethod = 'approx-star';
            reachOptionsList{1} = reachOptions;
        catch
            net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");
            nnvnet = "";
            reachOptionsList = {};
        end

    elseif contains(category, "nn4sys")
        if contains(onnx, "mscn")
            % nn4sys mscn: MATLAB's importNetworkFromONNX cannot parse these models
            % (Gather/embedding + elementwise product/division). Route via the
            % Python-importer manifest (forward pass cross-validated vs onnxruntime).
            % The net has NONLINEAR ElementwiseProduct/Division layers with no sound
            % approx/cp-star reach, so run FALSIFICATION-ONLY (reachOptionsList = {}): a
            % validated SAT witness (replayed through onnxruntime in falsify_single) or
            % unknown -- never an unsound unsat. mscn_2048d_dual.onnx is corrupt upstream
            % (fails MATLAB AND onnx 1.20 parsers) and abstains to unknown earlier in
            % run_vnncomp_instance.
            nnvnet = load_manifest_net(onnx);
            net = nnvnet;   % falsify_single dispatches NN.evaluate for NNV nets
            inputSize = nnvnet.Layers{1}.InputSize;
            reachOptionsList = {};
        elseif contains(onnx, "lindex")
            net = importNetworkFromONNX(onnx, "InputDataFormats", "BC");
            nnvnet = matlab2nnv(net);
            reachOptions = struct;
            reachOptions.reachMethod = 'approx-star'; % default parameters
            reachOptionsList{1} = reachOptions;
        else
            % pensieve: matlab2nnv cannot parse it (custom Reshape_To_GemmLayer + a generic
            % nnet InputLayer). Route through the Python importer manifest -> the net is FC/Relu/
            % Reshape/Concat only -> SOUND approx-star (validated: approx-star -> unsat in ~1.6s),
            % an upgrade from the prior probabilistic cp-star. PGD provides SAT. Robust: if the
            % manifest is unavailable, fall back to FALSIFICATION-ONLY (sat-or-unknown), never cp-star.
            try
                nnvnet = load_manifest_net(onnx);
                net = nnvnet;
                inputSize = nnvnet.Layers{1}.InputSize;
                reachOptions = struct; reachOptions.reachMethod = 'approx-star';
                reachOptionsList{1} = reachOptions;
            catch
                net = importNetworkFromONNX(onnx);
                if contains(onnx, "pensieve_big_parallel")
                    inputSize = [12,8];
                    inputFormat = "UU";
                    X = dlarray(rand(12,8), inputFormat);
                elseif contains(onnx, "pensieve_small_parallel")
                    inputSize = [12,8];
                    inputFormat = "UU";
                    X = dlarray(rand(12,8), inputFormat);
                    needReshape = 1;
                else
                    inputSize = [1,6,8];
                    inputFormat = "UUU";
                    X = dlarray(rand(1,6,8), inputFormat);
                end
                net = initialize(net, X);
                nnvnet = "";
                reachOptionsList = {};   % falsify-only (was cp-star -- probabilistic)
            end
        end
        % Somehow, some of these networks have discrepancies  (all sat (invalid))
        

    elseif contains(category, "relusplitter")
        if contains(onnx, "mnist")
            net = importNetworkFromONNX(onnx, "InputDataFormats", "BCT");
            inputFormat = "BCT";
            inputSize = [1 784];
        elseif contains(onnx, "oval")
            net = importNetworkFromONNX(onnx, "InputDataFormats", "BCSS");
            needReshape = 1;
        else
            net = importNetworkFromONNX(onnx, "InputDataFormats", "BC");
            if contains(onnx, "base")
                needReshape = 1;
            end
        end
        nnvnet = matlab2nnv(net);
        reachOptions = struct;
        reachOptions.reachMethod = 'relax-star-area';
        reachOptions.relaxFactor = 1;
        reachOptions.inputFormat = inputFormat;
        reachOptionsList{1} = reachOptions;
        % reachOptions.reachMethod = 'relax-star-area';
        % reachOptions.relaxFactor = 0.5;
        % reachOptionsList{2} = reachOptions;

    elseif contains(category, "safenlp")
        % safeNLP: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        % needReshape = ?
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star'; % default parameters
        reachOptionsList{1} = reachOptions;
        reachOptions = struct;
        reachOptions.reachMethod = 'exact-star'; % default parameters
        reachOptions.numCores = 1;  % Phase 1.3: avoid nested-parpool errors when called inside parfeval workers
        reachOptionsList{2} = reachOptions;
        nRand = 500;

    elseif contains(category, "sat_relu")
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BC");
        nnvnet = matlab2nnv(net);
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star'; % default parameters
        reachOptionsList{1} = reachOptions;
        reachOptions = struct;
        reachOptions.reachMethod = 'exact-star'; % default parameters
        reachOptions.numCores = 1;  % Phase 1.3: avoid nested-parpool errors when called inside parfeval workers
        reachOptionsList{2} = reachOptions;

    elseif contains(category, "soundness")
        
        net = importNetworkFromONNX(onnx, "InputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        reachOptions = struct;
        reachOptions.train_epochs = 500;
        reachOptions.train_lr = 0.001;
        reachOptions.coverage = 0.999;
        reachOptions.confidence = 0.999;
        reachOptions.train_mode = 'Linear';
        reachOptions.surrogate_dim = [-1,-1];
        reachOptions.threshold_normal = 1e-5;
        reachOptions.reachMethod = 'cp-star';
        reachOptionsList{1} = reachOptions;

    elseif contains(category, "tinyimagenet")
        % tinyimagenet: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS", "OutputDataFormats","BC");
        nnvnet = matlab2nnv(net);
        reachOptions = struct;
        reachOptions.train_epochs = 150;
        reachOptions.train_lr = 0.0001;
        reachOptions.coverage = 0.999;
        reachOptions.confidence = 0.999;
        reachOptions.train_mode = 'Linear';
        reachOptions.surrogate_dim = [-1,-1];
        reachOptions.threshold_normal = 1e-5;
        reachOptions.reachMethod = "cp-star";
        reachOptionsList{1} = reachOptions;
        needReshape = 1;
        nRand = 500;

    elseif contains(category, "tllverify")
        % tllverify: onnx to nnv
        net = importNetworkFromONNX(onnx,"InputDataFormats", "BC", 'OutputDataFormats',"BC");
        nnvnet = matlab2nnv(net);
        % 0.95 fast pass for the 18 big-N timeouts; tighter rungs for the 12 unknowns;
        % approx-star (full LP) as the tightest fallback. All sound (any k is over-approx).
        reachOptionsList = relax_ladder([0.95, 0.8, 0.5]);
        oo = struct(); oo.reachMethod = 'approx-star';
        reachOptionsList{end+1} = oo;

    elseif contains(category, "traffic")
        % traffic_signs_recognition: MATLAB's importNetworkFromONNX cannot
        % parse this binarized (Sign) model, so load via the Python-importer
        % manifest. Cross-validated against onnxruntime: max diff 4.3e-19,
        % argmax 8/8 over random inputs (2026-06-09; required the regenerated
        % manifest -- SignLayer + proper transpose handling -- plus the
        % BCHW<->BHWC identity-by-convention placeholders).
        % Input is ONNX [1,30,30,3] (NHWC); vnnlib X order is ONNX row-major
        % (C fastest), so unflatten via needReshape=3:
        %   img = permute(reshape(x, [3 30 30]), [3 2 1])  -> [30 30 3] HWC
        nnvnet = load_manifest_net(onnx);
        net = nnvnet;   % falsify_single dispatches NN.evaluate for NNV nets
        % reshape size BEFORE the [3 2 1] permute. The old hardcoded [3 30 30]
        % only fit the 30x30 models and threw "Number of elements must not
        % change" on the 48x48 / 64x64 traffic-sign models (6912 / 12288 != 2700).
        % Derive it from the net's own input layer: InputSize [H W C] -> [C W H].
        inputSize = nnvnet.Layers{1}.InputSize([3 2 1]);
        needReshape = 3;
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star';
        reachOptionsList{1} = reachOptions;

    elseif contains(category, "vggnet")
        % error("TODO: add support")
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BCSS", 'OutputDataFormats',"BC");
        try
            nnvnet = matlab2nnv(net);
        catch
            nnvnet = "";
        end
        needReshape = 1; %?
        reachOptions = struct;
        reachOptions.train_epochs = 150;
        reachOptions.train_lr = 0.0001;
        reachOptions.coverage = 0.999;
        reachOptions.confidence = 0.999;
        reachOptions.train_mode = 'Linear';
        reachOptions.surrogate_dim = [-1,-1];
        reachOptions.threshold_normal = 1e-5;
        reachOptions.reachMethod = 'cp-star';
        reachOptionsList{1} = reachOptions;
    
    elseif contains(category, "vit")
        % vit: MATLAB's importer collapses the 133-node ViT into ONE opaque
        % fused custom layer (Shape_To_ReduceMeanLayer) -> matlab2nnv refuses ->
        % the old path fell to cp-star (probabilistic + env-dependent) -> 200/200
        % unknown. Route through the PYTHON manifest instead: the layout-model
        % fix in onnx2nnv.py makes the manifest evaluate match onnxruntime to
        % ~5e-7 (MATLAB-validated 2026-06-13), so the NN-net numerical-gradient
        % PGD falsifier produces REAL sat witnesses (replayed through the
        % onnxruntime gate before emission). Reach through DynamicMatmulLayer
        % refuses by design -> unsat is out of scope -> sat-or-unknown only.
        nnvnet = load_manifest_net(onnx);
        net = nnvnet;   % falsify_single dispatches NN.evaluate for NNV nets
        inputSize = nnvnet.Layers{1}.InputSize;
        needReshape = 1;   % CHW-flat vnnlib -> [H,W,C] image, same as cifar-class
        reachOptionsList = {};   % no sound reach available -> skip to unknown
        % (legacy MATLAB-import + cp-star path removed; see git history pre-#345)

    elseif contains(category, "yolo")
        % yolo: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS"); % padlayer
        try
            nnvnet = matlab2nnv(net);
        catch
            nnvnet = "";
        end
        needReshape = 2; % ?
        reachOptions = struct;
        reachOptions.train_epochs = 200;
        reachOptions.train_lr = 0.0001;
        reachOptions.coverage = 0.999;
        reachOptions.confidence = 0.999;
        reachOptions.train_mode = 'Linear';
        reachOptions.surrogate_dim = [-1,-1];
        reachOptions.threshold_normal = 1e-5;
        reachOptions.reachMethod = 'cp-star'; % default parameters
        reachOptionsList{1} = reachOptions;

    elseif contains(category, "test")
        % VNN-COMP 'test' sanity benchmark: tiny FC/ReLU nets that MATLAB's
        % importNetworkFromONNX rejects ("ONNX model not supported"). Route via the
        % Python-importer manifest (forward pass cross-validated vs onnxruntime); the
        % nets are small enough that approx-star decides them soundly.
        nnvnet = load_manifest_net(onnx);
        net = nnvnet;   % falsify_single dispatches NN.evaluate for NNV nets
        inputSize = nnvnet.Layers{1}.InputSize;
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star';
        reachOptionsList{1} = reachOptions;

    else % all other benchmarks
        error("ONNX model not supported")
    end

    % Method-selection POST-PROCESSING. Extracted to a PURE, unit-testable helper so the routing
    % (which method actually runs per category, and whether a probabilistic method survives) is
    % guarded by tests instead of being implicit in a 1500-line function. Does, in order:
    %   (1) Phase-1.5 sound prepend: for a cp-star-led list, prepend {approx-star, relax-star 0.5}
    %       so the dispatcher tries the SOUND methods first (cp-star stays as a fallback);
    %   (2) slow_cats (cifar100/cora/safenlp/sat_relu/tinyimagenet/vggnet): drop exact-star and lead
    %       with fast still-SOUND over-approximations {approx-zono, abs-dom} to avoid timing out;
    %   (3) cp-star QUARANTINE: gated by env NNV_QUARANTINE_CPSTAR (default OFF -> behavior UNCHANGED).
    %       When set, strip cp-star entirely so a probabilistic 'unsat' can never be emitted as a
    %       sound verdict in ANY reach branch (single-spec OR the two parfor branches). Strictly
    %       sound: removing cp-star can only turn an unsat/unknown into unknown, never a wrong verdict.
    % See tests/nn/vnncomp/test_run_vnncomp_routing.m for the per-category routing assertions.
    reachOptionsList = i_finalize_reach_options(reachOptionsList, category, is_nnvnet_valid(nnvnet));

    % ---- Per-category PGD/falsification budget (VNNCOMP2026 tuning) ----------------
    % Applied LAST so it is the single source of truth for falsification effort,
    % overriding any per-category nRand set above. rows: {key, n_restarts, n_steps,
    % max_time_s, nRand}. `contains` is a substring test, so the MOST SPECIFIC key
    % must come first when one key is a substring of another: cctsdb_yolo BEFORE yolo.
    % (collins_aerospace vs collins_rul are disjoint; "vggnet" matches vggnet16.)
    %
    % nRand sizing: gradient PGD is the PRIMARY falsifier, so the random sampling is
    % only a light fallback for the dlnetwork categories -- measured ~15 ms/sample
    % (one-at-a-time predict), so nRand=1000 cost ~15 s/instance for ~zero marginal
    % SAT over PGD. Keep it small there. The two MANIFEST categories that get NO PGD
    % (pgd_falsify is dlnetwork-gated): lsnc_relu and traffic_signs -- rely on random
    % sampling, so keep their nRand higher until PGD reaches the NN path (a follow-up).
    ftab = {
        "acasxu",            20, 40, 8,   100   % trimmed PGD (was 60,80,30,200): the 30s budget was
                                                % wasted in full on the many SAFE acasxu props (exact/approx
                                                % is complete, so PGD is only a SAT accelerator). 8s still
                                                % finds easy counterexamples. Sound (falsify = SAT-or-unknown).
        "sat_relu",          60, 100,30,  200
        "cersyve",           30, 50, 4,   150
        "lsnc_relu",         50, 80, 4,   500   % manifest: no PGD -> random is primary
        "relusplitter",      30, 50, 3,   150
        "dist_shift",        30, 60, 4,   150
        "linearize",         30, 60, 4,   150
        "tllverifybench",    40, 80, 5,   200
        "collins_rul",       30, 60, 4,   150
        "collins_aerospace", 20, 40, 3,   150
        "nn4sys",            25, 50, 3,   150
        "safenlp",           40, 80, 30,  300
        "malbeware",         40, 80, 60,  150
        "ml4acopf",          25, 50, 3,   200
        "cora",              40, 80, 30,  200
        "metaroom",          15, 30, 2,   100
        "cifar100",          15, 40, 20,  60
        "tinyimagenet",      15, 40, 20,  60
        "vggnet",            4,  12, 1.5, 20
        "traffic_signs",     20, 40, 5,   1000  % manifest + BINARIZED (Sign): PGD is gradient-blind, so
                                                 % the per-coord box VERTICES in create_random_examples
                                                 % (~(nRand-2)/2 = 499, > the ~400 the seeded probe needs)
                                                 % are the PRIMARY finder. Reliably decides the small-net
                                                 % (30x30) sats; the 64x64 net's per-sample forward is too
                                                 % slow to exhaust enough vertices in-budget -> needs a
                                                 % BATCHED-vertex forward (+ bit-flip local search, §4.5).
        "cctsdb_yolo",       30, 50, 5,   150
        "yolo",              20, 40, 3,   100
        "vit",               20, 50, 25,  50
        "cgan",              80, 120,40,  2000   % reach dropped for cgan -> spend the whole budget on PGD (5-dim input, cheap)
        "soundnessbench",    30, 50, 4,   100
    };
    for r = 1:size(ftab,1)
        if contains(category, ftab{r,1})
            falsifyOpts.n_restarts = ftab{r,2};
            falsifyOpts.n_steps    = ftab{r,3};
            falsifyOpts.max_time   = ftab{r,4};
            nRand                  = ftab{r,5};
            break
        end
    end

    % Phase 1 -- batched-random + coordinate-descent polish falsifier (low-dim flat-feature nets, e.g.
    % acasxu). A single high-N vectorized predict finds most thin counterexamples the trimmed PGD/random
    % miss (the 15-instance gold-sat "falsify gap"); a wall-capped local polish recovers the rest. Pure
    % SAT-or-unknown: every returned witness is validate_witness-checked AND onnxruntime-replayed by the
    % caller before 'sat' stands, so this can only ADD a sound SAT, never a wrong verdict. acasxu opts-in.
    if contains(category, "acasxu")
        falsifyOpts.batched_N    = 2e6;   % one batched predict over uniform samples (validated 15/15 recovery)
        falsifyOpts.near_miss    = 0.5;   % polish only when best sample margin is within this of 0 (else clearly unsat -> skip)
        falsifyOpts.polish       = true;  % coordinate-descent polish from the top-K lowest-margin samples
        falsifyOpts.polish_max_t = 3;     % wall cap (s) on the whole polish stage -> bounds cost on unsat instances
    end

    % lsnc_relu opts-in to the SAME batched-random falsifier (scalar [6] inputSize, 'CB' dlnetwork). The
    % batched dlnetwork forward finds the gold-sat the manifest's numerical-gradient PGD missed (probe:
    % state_1 at margin -3.2e-3). lsnc is currently 0-solved (all unknown) -> pure-gain category. Every
    % witness is validate_witness-checked AND onnxruntime-replayed (riskyNet) -> SAT-or-unknown.
    if contains(category, "lsnc_relu")
        falsifyOpts.batched_N    = 2e6;   % d=6 -> cheap; probe found the violation by 3e5 samples
        falsifyOpts.near_miss    = 0.5;
        falsifyOpts.polish       = true;
        falsifyOpts.polish_max_t = 3;
    end

end

% Create an array of random examples from input set and reshape if necessary
% We use dlnetwork for simulation (MATLAB data structure)
function [net,nnvnet,needReshape,reachOptionsList,inputSize,inputFormat,nRand,falsifyOpts] = i_load_vnncomp_network_cached(category, onnx, vnnlib)
% Cache wrapper for load_vnncomp_network. The ONNX import + matlab2nnv is paid PER INSTANCE but
% the net + category options are ONNX-determined, so import ONCE per onnx and reuse across that
% onnx's many vnnlib instances (VNN-COMP allows this format conversion in the prepare phase;
% prepare_instance.sh can pre-build the cache, else the first run builds it). Keyed PER-ONNX
% (path) on the onnx file's size+mtime + the NNV version -- a category with MULTIPLE networks
% caches each separately, and a changed onnx OR a new NNV re-imports. SOUNDNESS: the net used to
% verify MUST be identical to the onnx, so a stale cache would be a -150 -> the size+mtime+version
% key guarantees a hit only for the SAME onnx+NNV (matlab2nnv is deterministic). Best-effort: any
% cache error falls back to a fresh import (never fails the instance).
    p = [char(onnx) '.netcache.mat'];
    d = dir(char(onnx));
    if isfile(p) && ~isempty(d)
        try
            % SELF-HEAL against a missing custom-layer +package. The cache stores a dlnetwork that
            % REFERENCES importer-generated custom-layer classes (+<onnx>, e.g. +mnist_concat). If that
            % package is absent from the path, load() does NOT error -- it SILENTLY substitutes default
            % layer objects (warning "...Default objects will be substituted") and the net degrades to a
            % useless 'unknown'. Force warnings ON so lastwarn captures the substitution, then treat a
            % degraded load as a cache MISS: the fresh import below regenerates the +package AND re-warms
            % the cache. Sound either way (worst case = an extra re-import); makes a deleted/absent package
            % auto-recover instead of silently producing unknowns in a sweep or the competition.
            owarn = warning('on','all');
            restoreWarn = onCleanup(@() warning(owarn));  % restore warning state even if load() THROWS
            lastwarn('');
            S = load(p, 'C'); C = S.C;
            wmsg = lastwarn;
            clear restoreWarn;                            % success path: restore immediately
            if contains(wmsg,'Default objects will be substituted') || contains(wmsg,'Unable to load instances of class')
                error('netcache:degraded','custom-layer +package missing -> re-import (%s)', wmsg);
            end
            % TYPE-VALIDITY guard (soundness): a silent dlnetwork.loadobj failure can substitute a plain
            % `double` for C.net while the metadata (bytes/datenum/ver) still matches -> the HIT is accepted
            % and a later dot-index throws "Dot indexing is not supported for variables of type double"
            % (cifar100, 5 instances) -> unknown; worse, a partially-substituted net could reach and emit a
            % CONFIDENT WRONG verdict. Reject any cache whose net is not a real network object (or whose NNV
            % net is invalid) and fall through to a fresh re-import. Can only help: a valid cache has a
            % non-numeric C.net + a valid C.nnvnet, so genuine HITs are unaffected.
            cacheNetOK = isfield(C,'net') && ~isnumeric(C.net) ...
                      && isfield(C,'nnvnet') && ~isnumeric(C.nnvnet) && is_nnvnet_valid(C.nnvnet);
            if isfield(C,'onnx_bytes') && isequal(C.onnx_bytes, d.bytes) ...
                    && abs(C.onnx_datenum - d.datenum) < 1e-9 ...
                    && isfield(C,'nnv_ver') && strcmp(C.nnv_ver, i_nnv_ver()) ...
                    && isfield(C,'cfg_ver') && strcmp(C.cfg_ver, i_cfg_ver()) ...   % invalidate on config change
                    && cacheNetOK                                                    % reject substituted-double / invalid net
                net=C.net; nnvnet=C.nnvnet; needReshape=C.needReshape; reachOptionsList=C.reachOptionsList;
                inputSize=C.inputSize; inputFormat=C.inputFormat; nRand=C.nRand; falsifyOpts=C.falsifyOpts;
                fprintf('net cache HIT (%s)\n', p);
                return;
            end
        catch
            % corrupt/foreign cache -> fall through to a fresh import
        end
    end
    [net,nnvnet,needReshape,reachOptionsList,inputSize,inputFormat,nRand,falsifyOpts] = load_vnncomp_network(category, onnx, vnnlib);
    if ~isempty(d) && is_nnvnet_valid(nnvnet)
        try
            C = struct();
            C.net=net; C.nnvnet=nnvnet; C.needReshape=needReshape; C.reachOptionsList=reachOptionsList;
            C.inputSize=inputSize; C.inputFormat=inputFormat; C.nRand=nRand; C.falsifyOpts=falsifyOpts;
            C.onnx_bytes=d.bytes; C.onnx_datenum=d.datenum; C.nnv_ver=i_nnv_ver(); C.cfg_ver=i_cfg_ver();
            tmp = [p '.tmp'];   % same dir as p; no two workers write the same onnx's cache concurrently
            save(tmp, 'C', '-v7.3');
            movefile(tmp, p, 'f');
            fprintf('net cache BUILT (%s)\n', p);
        catch
            % best-effort: a failed save just means the next instance re-imports
        end
    end
end

function v = i_nnv_ver()
    try, v = char(NNVVERSION()); catch, v = 'na'; end
end

function v = i_cfg_ver()
% Config version for the net cache: the .netcache.mat stores reachOptionsList/nRand/falsifyOpts (the
% verification CONFIG), but the cache validity guard otherwise only checks the onnx+NNV version -- so a
% CODE change to load_vnncomp_network's config (e.g. acasxu reach ladder / numCores / PGD budget) would
% be masked by a stale cache (served the OLD config -> wrong/slow verdicts). BUMP this string whenever
% the config logic in load_vnncomp_network changes, to invalidate all stale caches.
    v = '2026-06-25.pgd-full-disjunction-min-margin';
end

function xRand = create_random_examples(net, lb, ub, nR, inputSize, needReshape,inputFormat)
    % Draw half per-coord BERNOULLI VERTICES (each coord independently at lb OR ub) and half uniform
    % interior. Vertices are the ONLY way to falsify GRADIENT-BLIND nets -- binarized / Sign / quantized
    % models (e.g. traffic_signs): a Sign activation has zero gradient a.e., so PGD cannot descend and
    % the adversarials sit at L-inf box CORNERS. Vertices are valid points in [lb,ub], so this can only
    % ADD candidate sats (each is property-checked here + the witness is validate_witness/ORT-replayed
    % before emit) -- harmless to smooth nets, where PGD is the primary finder and this is only a fallback.
    % DETERMINISTIC via a LOCAL RNG stream so the falsifier is REPRODUCIBLE (a random falsifier is a
    % defect: the official run could draw a different state than dev) WITHOUT resetting the GLOBAL rng --
    % which would change the caller's stream for the rest of the session (esp. the persistent-session
    % runner). 'twister' Seed 0 matches the default generator, so the draws are identical to a `rng(0)`
    % version but with no global side effect. (Plan §4.5.)
    rs = RandStream('twister', 'Seed', 0);
    nrand = max(0, nR - 2);
    nv = floor(nrand/2);                  % => (nR-2)/2 vertices (the rest uniform); lb,ub added below
    ncont = nrand - nv;
    if nv > 0, Xvert = lb + (ub - lb) .* (rand(rs, numel(lb), nv) < 0.5); else, Xvert = []; end
    if ncont > 0, Xcont = lb + (ub - lb) .* rand(rs, numel(lb), ncont); else, Xcont = []; end
    xRand = [lb, ub, Xvert, Xcont];
    if needReshape
        if needReshape ==2 % for collins only (full_window_40) and metaroom
            newSize = [inputSize(2) inputSize(1) inputSize(3:end)];
            xRand = reshape(xRand, [newSize nR]);
            xRand = permute(xRand, [2 1 3 4]);
        elseif needReshape == 3
            % ONNX row-major NHWC flat with inputSize = [C W H]:
            % unflatten then permute to NNV's [H W C] arrays.
            xRand = reshape(xRand, [inputSize nR]);
            xRand = permute(xRand, [3 2 1 4]);
        else
            xRand = reshape(xRand, [inputSize nR]);
            xRand = permute(xRand, [2 1 3 4]);
        end
%         xRand = python_reshape(xRand, [inputSize nR]);
    else
        xRand = reshape(xRand,[inputSize nR]); % reshape vectors into net input size
    end
    if isa(net, 'dlnetwork') % need to convert to dlarray
        if strcmp(inputFormat, "default")
            if isa(net.Layers(1, 1), 'nnet.cnn.layer.ImageInputLayer')
                xRand = dlarray(xRand, "SSCB");
            elseif isa(net.Layers(1, 1), 'nnet.cnn.layer.FeatureInputLayer') || isa(net.Layers(1, 1), 'nnet.onnx.layer.FeatureInputLayer')
                xRand = dlarray(xRand, "CB");
            else
                disp(net.Layers(1,1));
                error("Unknown input format");
            end
        else
            if contains(inputFormat, "U")
                xRand = dlarray(xRand, inputFormat+"U");
            else
                xRand = dlarray(xRand, inputFormat);
            end
        end
        
    end
end

% Write counterexample to output file
function write_counterexample(outputfile, counterEx)
    % First line - > sat
    % after that, write the variables for each input dimension  of the counterexample
    %
    % Example:
    %  ( (X_0 0.12132)
    %    (X_1 3.45454)
    %    ( .... )
    %    (Y_0 2.32342)
    %    (Y_1 3.24355)
    %    ( ... )
    %    (Y_N 0.02456))
    %

    precision = '%.16g'; % set the precision for all variables written to txt file
    % open file and start writing counterexamples
    fid = fopen(outputfile, 'a+');
    x = counterEx{1};
    x = reshape(x, [], 1);
    % begin specifying value for input example
    fprintf(fid,'(');
    for i = 1:length(x)
        fprintf(fid, "(X_" + string(i-1) + " " + num2str(x(i), precision)+ ")\n");
    end
    y = counterEx{2};
    y = reshape(y, [], 1);
    % specify values for output example
    for j =1:length(y)
        fprintf(fid, "(Y_" + string(j-1) + " " + num2str(y(j), precision)+ ")\n");
    end
    fprintf(fid, ')');
    % close and save file
    fclose(fid);

end

% cctsdb_yolo complete-enumeration verifier: shell out to cctsdb_enumerate.py
% (same folder as this runner). The script is sound-or-unknown by construction:
% it structurally re-verifies, per instance, that (a) the only free inputs are
% X_12288/X_12289, (b) the ONNX consumes them only through Cast(int64)
% truncation (=> output piecewise-constant on unit cells), (c) the output spec
% is a single half-space on Y_0 -- and enumerates all <= 3969 cells through
% onnxruntime. Exit protocol: 10 = SAT (witness CSV written + "SAT p q y" on
% stdout), 11 = UNSAT, anything else = unknown. ANY parse/replay irregularity
% on the MATLAB side also degrades to unknown -- never an unsound verdict.
function [status, counterEx] = verify_nn4sys_lindex(onnx, vnnlib)
    % Sound batched-IBP + adaptive box-subdivision decider for nn4sys lindex (nn4sys_lindex_decide.py).
    % Crash-safe exit-code map: 10=sat, 11=unsat, ANYTHING ELSE -> unknown (fail-open). A Python crash
    % (uncaught exception -> exit 1) therefore degrades to a sound `unknown`, never a wrong unsat (-150).
    status = 2; counterEx = nan;
    here = fileparts(mfilename('fullpath'));
    script = fullfile(here, 'nn4sys_lindex_decide.py');
    if ~isfile(script), fprintf('nn4sys_lindex_decide.py not found -> unknown\n'); return; end
    py = python_exe();   % needs onnx+numpy+onnxruntime (it parses the vnnlib itself; no vnnlib pkg)
    % python_exe() falls back to a PLAIN interpreter when the stack is absent, so isempty alone never
    % trips. Probe the imports explicitly -> fail fast to a sound `unknown` instead of spawning a doomed
    % subprocess. (A bad stack would still fail-open via the crash-safe exit map, but this is cleaner.)
    [probe_rc, ~] = system(sprintf('%s -c "import onnx, onnxruntime, numpy"', py));
    if isempty(py) || probe_rc ~= 0
        fprintf('no onnx+onnxruntime+numpy python for the lindex decider -> unknown\n'); return;
    end
    if ispc
        cmd = sprintf('%s "%s" "%s" "%s"', py, script, char(onnx), char(vnnlib));
    else
        cmd = sprintf('timeout 60 %s "%s" "%s" "%s"', py, script, char(onnx), char(vnnlib));
    end
    [rc, out] = system(cmd);
    disp(strtrim(out));
    switch rc
        case 11, status = 1;   % unsat: every box's sound interval is provably inside its band
        case 10                % sat: a concrete onnxruntime witness
            xv = regexp(out, '\(X_0\s+([-+0-9.eE]+)\)', 'tokens', 'once');
            yv = regexp(out, '\(Y_0\s+([-+0-9.eE]+)\)', 'tokens', 'once');
            if ~isempty(xv) && ~isempty(yv)
                status = 0; counterEx = {str2double(xv{1}); str2double(yv{1})};
            end                % else: claimed sat but unparseable witness -> stay unknown
        otherwise, status = 2; % 12 (cant/unknown) / 124 (timeout) / 1 (crash) / anything -> fail-open
    end
end

function [status, counterEx] = verify_cctsdb_enumeration(onnx, vnnlib)
    status = 2; counterEx = nan;
    here = fileparts(mfilename('fullpath'));
    script = fullfile(here, 'cctsdb_enumerate.py');
    if ~isfile(script)
        fprintf('cctsdb_enumerate.py not found next to the runner -> unknown\n');
        return;
    end
    witness_csv = [tempname '.csv'];
    cleanup = onCleanup(@() delete_if_exists(witness_csv)); %#ok<NASGU>
    py = python_exe();
    % The script self-limits via --timeout: 300 s leaves margin inside the 350 s
    % per-instance benchmark budget (the enumeration itself measures 2-25 s).
    cmd = sprintf('%s "%s" "%s" "%s" "%s" --timeout 300', ...
        py, script, char(onnx), char(vnnlib), witness_csv);
    [st, out] = system(cmd);
    disp(strtrim(out));
    if st == 10                  % SAT with witness
        % stdout carries "SAT p q y"; the CSV carries the FULL input vector in
        % FLAT vnnlib order X_0..X_{N-1} (one value per line).
        tok = regexp(out, '\<SAT\s+(-?\d+)\s+(-?\d+)\s+([-+0-9.eE]+)', 'tokens', 'once');
        try
            x = readmatrix(witness_csv);
        catch
            x = [];
        end
        if isempty(tok) || isempty(x) || ~all(isfinite(x(:)))
            fprintf('malformed SAT report from cctsdb_enumerate.py -> unknown\n');
            return;
        end
        y = str2double(tok{3});
        if ~isfinite(y)
            fprintf('non-finite witness output from cctsdb_enumerate.py -> unknown\n');
            return;
        end
        % {flat input; output} -- the exact cell write_counterexample expects
        % (same shape falsify_single produces).
        counterEx = {reshape(x, [], 1); y};
        status = 0;
    elseif st == 11              % UNSAT: complete enumeration, all cells safe
        status = 1;
    end                          % anything else: guard violation/timeout -> unknown
end

function delete_if_exists(f)
    if exist(f, 'file'), delete(f); end
end

% adaptive_cruise NONLINEAR falsifier: shell out to adaptive_cruise_falsify.py, which samples the input
% box + forwards the REAL onnx via onnxruntime and accepts a witness ONLY if EVERY assertion holds under
% the authoritative `vnnlib` parser with a robust margin. Exit 10 = SAT (+ flat witness csv); anything
% else (12 unknown / missing deps / error) -> unknown. Never returns unsat. Mirrors verify_cctsdb_enumeration.
function [status, counterEx] = verify_adaptive_cruise_falsify(onnx, vnnlib)
    status = 2; counterEx = nan;
    here = fileparts(mfilename('fullpath'));
    script = fullfile(here, 'adaptive_cruise_falsify.py');
    if ~isfile(script)
        fprintf('adaptive_cruise_falsify.py not found next to the runner -> unknown\n');
        return;
    end
    witness_csv = [tempname '.csv'];
    cleanup = onCleanup(@() delete_if_exists(witness_csv)); %#ok<NASGU>
    % python_exe() resolves an onnx+onnxruntime interpreter; this helper ALSO needs `vnnlib`. If the
    % resolved one lacks it, prefer a candidate that imports all three (else the path silently degrades
    % to unknown even when vnnlib is installed elsewhere). Sound either way: no vnnlib -> exit 12 -> unknown.
    py = python_exe();
    [stv, ~] = system(sprintf('%s -c "import vnnlib"', py));
    if stv ~= 0
        extra = {strtrim(getenv('NNV_ORT_PYTHON')), ...
                 fullfile(getenv('HOME'), 'taylor_venv', 'bin', 'python'), 'python3', 'python'};
        for k = 1:numel(extra)
            c = extra{k};
            if isempty(c), continue; end
            [s2, ~] = system(sprintf('"%s" -c "import onnx, onnxruntime, vnnlib"', c));
            if s2 == 0, py = ['"' c '"']; break; end
        end
    end
    % WALL BUDGET: bound the sampler to the OFFICIAL per-instance timeout so an UNFINDABLE instance
    % returns a graceful `unknown` (exit 12) IN-TIME instead of grinding --n 1.5M samples past the
    % official timeout and being hard-killed -> EMPTY result. execute.py exports NNV_REACH_BUDGET =
    % 0.95*official_timeout; we leave a further ~15s of headroom for python+onnxruntime startup, the
    % onnx load, and the harness's own slack. Falls back to a generous fixed cap for dev/test callers
    % that run the runner without execute.py (NNV_REACH_BUDGET unset). Sound either way: the cap only
    % ever turns a would-be timeout into `unknown`, never affects a found witness (which exits fast).
    reachBudget = str2double(getenv('NNV_REACH_BUDGET'));
    if isfinite(reachBudget) && reachBudget > 0
        wallSec = max(5, reachBudget - 15);
    else
        wallSec = 90;   % dev/test default (< the uniform 100s official adaptive_cruise timeout)
    end
    cmd = sprintf('%s "%s" "%s" "%s" "%s" --n 1500000 --tol 0.01 --max-seconds %.1f', ...
        py, script, char(onnx), char(vnnlib), witness_csv, wallSec);
    [st, out] = system(cmd);
    disp(strtrim(out));
    if st == 10                          % SAT: witness csv carries the flat input; stdout "SAT y0 y1 ..."
        tok = regexp(out, '\<SAT\s+([-+0-9.eE\s]+)', 'tokens', 'once');
        try
            x = readmatrix(witness_csv);
        catch
            x = [];
        end
        if isempty(tok) || isempty(x) || ~all(isfinite(x(:)))
            fprintf('malformed SAT report from adaptive_cruise_falsify.py -> unknown\n');
            return;
        end
        y = sscanf(strtrim(tok{1}), '%f');
        if isempty(y) || ~all(isfinite(y))
            fprintf('non-finite witness output from adaptive_cruise_falsify.py -> unknown\n');
            return;
        end
        counterEx = {reshape(x, [], 1); reshape(y, [], 1)};
        status = 0;
    elseif st == 11                      % UNSAT: input region PROVABLY EMPTY (the nonlinear input
        % constraint is infeasible over the box) -> no input -> no counterexample -> vacuously UNSAT.
        % adaptive_cruise_falsify.py proves emptiness by interval arithmetic over the AUTHORITATIVE
        % parsed AST (+ a dense-grid refutation backstop) and emits 11 ONLY when proven; any
        % uncertainty there exits 12 (unknown). Sound: an empty region cannot contain a witness.
        status = 1;
    end                                  % anything else (incl. 12) -> unknown (sound; never unsat)
end

% (python_exe is defined once above -- the Linux-aware version that
% defaults to python3 on the competition VM.)

% Load an NNV net from the Python-importer manifest written alongside the ONNX
% (tools/onnx2nnv_python/onnx2nnv.py). Used for models MATLAB's
% importNetworkFromONNX cannot parse (lsnc_relu, traffic_signs_recognition).
function nnvnet = load_manifest_net(onnx)
    manifest = regexprep(char(onnx), '\.onnx$', '.nnv.mat');
    if ~isfile(manifest)
        error('run_vnncomp_instance:noManifest', ...
            ['NNV manifest not found: %s\nGenerate it with:\n' ...
             '  python tools/onnx2nnv_python/onnx2nnv.py "%s" --vnnlib <spec.vnnlib>'], ...
            manifest, char(onnx));
    end
    nnvnet = load_nnv_from_mat(manifest);
end

% Falsification function (random simulation looking for counterexamples)
function counterEx = falsify_single(net, lb, ub, inputSize, nRand, Hs, needReshape, inputFormat, opts)
    counterEx = nan;
    % Per-category PGD budget (load_vnncomp_network's falsifyOpts). Default preserves
    % the previous hardcoded seed/max_time so any other caller keeps working unchanged.
    if nargin < 9 || isempty(opts), opts = struct('seed', 0, 'max_time', 5); end
    if ~isfield(opts, 'seed'),     opts.seed = 0;     end
    if ~isfield(opts, 'max_time'), opts.max_time = 5; end
    rng(opts.seed, 'twister');   % B: APPLY the seed -> reproducible falsification (PGD restarts + Box.sample); sound (witness still ort-validated)

    % Phase 1 -- BATCHED-RANDOM + coordinate-descent polish, tried FIRST for low-dim flat-feature
    % dlnetwork nets (acasxu footprint: d<=20, 3-elem image inputSize OR scalar flat-feature inputSize,
    % identity reshape). One vectorized predict over a huge uniform sample finds thin counterexamples that
    % the trimmed PGD/random miss (the gold-sat "falsify gap"); a wall-capped local polish recovers the
    % rest. SOUND: returns a witness ONLY if it passes validate_witness here (and the caller ORT-replays
    % it), so it can only ADD a sound SAT or fall through to PGD/random. Gated by falsifyOpts.batched_N
    % (set for acasxu (3-elem inputSize, SSCB) and lsnc_relu (scalar inputSize, CB) in the ftab block).
    if isfield(opts,'batched_N') && opts.batched_N > 0 && isa(net,'dlnetwork') ...
            && numel(lb) <= 20 && (numel(inputSize) == 3 || isscalar(inputSize)) ...
            && (isempty(needReshape) || all(needReshape(:) == 0))
        try
            nearMiss = 0.5; if isfield(opts,'near_miss'),    nearMiss = opts.near_miss;   end
            doPol    = isfield(opts,'polish') && opts.polish;
            polMaxT  = 3;   if isfield(opts,'polish_max_t'), polMaxT  = opts.polish_max_t; end
            cex = falsify_batched_random(net, lb, ub, Hs, inputSize, opts.batched_N, doPol, polMaxT, nearMiss);
            if iscell(cex) && validate_witness(net, cex{1}, lb, ub, Hs, inputSize, inputFormat, needReshape)
                if ~isempty(getenv('NNV_DEBUG_BATCHED')), fprintf('[batched-random] validated SAT witness (d=%d, inputSize=%s)\n', numel(lb), mat2str(inputSize)); end
                counterEx = cex; return;
            end
        catch
            % fall through to PGD / random sampling
        end
    end

    % Gradient-directed falsification FIRST (FGSM warm-start + PGD). pgd_falsify maps
    % the flat input to the network-input layout with the SAME reshape+permute the
    % runner uses (needReshape), so it works for image/permuted inputs too. NNV found
    % 354 SAT vs ~1000 for the field in 2025; this targets that gap. dlnetwork nets use
    % autodiff; NNV NN nets (the Python-importer manifest path, e.g. traffic_signs --
    % all-SAT in 2025 -- and lsnc_relu) use a numerical gradient inside pgd_falsify, so
    % they too get gradient falsification instead of random sampling alone. Wrapped in
    % try/catch and VALIDATED before acceptance, so it can only ADD a sound SAT or fall
    % through to the existing random sampling. [VNNCOMP2026_STRATEGY Pillar 1/2]
    if isa(net, 'dlnetwork') || isa(net, 'NN')
        try
            [cex, found] = pgd_falsify(net, lb, ub, Hs, inputSize, inputFormat, needReshape, opts);
            if found && validate_witness(net, cex{1}, lb, ub, Hs, inputSize, inputFormat, needReshape)
                counterEx = cex; return;
            end
        catch
            % fall through to random sampling
        end
    end
    try
        xRand = create_random_examples(net, lb, ub, nRand, inputSize, needReshape, inputFormat);
    catch
        % Random-sampling falsification unavailable for this instance -- e.g. inputSize/vnnlib
        % element-count mismatch makes the reshape in create_random_examples throw ("Number of
        % elements must not change"), as seen on nn4sys pensieve_*_parallel under the manifest
        % route. Skip it: sound (no counterexample found -> the reach verdict / 'unknown' stands;
        % falsification can only ever ADD a sound SAT, never decide UNSAT), and it stops one bad
        % input shape from crashing the whole instance to NO output.
        return;
    end
    s = size(xRand);
    n = length(s);
    %  look for counterexamples
    for i=1:s(n)
        x = get_example(xRand, i);
        try
            if isa(net, 'NN')   % NNV net (Python-importer manifest path)
                yPred = net.evaluate(x);
            else
                yPred = predict(net, x);
            end
            if isa(x, 'dlarray') % if net is a dlnetwork
                x = extractdata(x);
                yPred = extractdata(yPred);
            end
            % check if property violated
            yPred = reshape(yPred, [], 1); % convert to column vector (if needed)
            % disp([x;yPred']);
            for h=1:length(Hs)
                if Hs(h).contains(double(yPred)) % property violated
                    % check if the counter example needs to be reshaped
                    n = numel(x);
                    if needReshape == 2
                        % x = reshape(x, [n 1]);
                        x = permute(x, [2 1 3]);
                    elseif needReshape == 1
                        if ndims(x) == 3 % RGB  image
                            x = permute(x, [2 1 3]);
                        end
                    elseif needReshape == 3
                        % [28] create_random_examples built x via
                        % reshape([C W H]) then permute([3 2 1 4]) -> [H W C].
                        % Invert that permute so write_counterexample flattens the
                        % witness back in the ORIGINAL ONNX NHWC flat order
                        % (C fastest); otherwise the SAT counterexample is written
                        % H-fastest -> invalid / competition penalty. permute([3 2 1])
                        % is its own inverse. (Was missing: only 1/2 were handled.)
                        if ndims(x) == 3
                            x = permute(x, [3 2 1]);
                        end
                    end
                    counterEx = {x; yPred}; % save input/output of countex-example
                    break;
                end
            end
        end
    end
end

% Phase 1 -- batched-random + coordinate-descent polish falsifier for low-dim flat-feature dlnetwork nets
% (acasxu). Draws N uniform samples over [lb,ub], evaluates them in ONE batched predict, and computes each
% sample's worst margin = min over the OR-halfspaces Hs of (max row of G*y - g) (<=0 means that sample
% violates the unsafe region = a counterexample). If a raw sample violates, returns it. Otherwise, when the
% best margin is a near-miss (< nearMiss) and doPolish is set, runs a wall-capped coordinate descent from
% the top-K lowest-margin samples to push a witness past 0. Returns {x; y} (x in flat vnnlib order) or [].
% SOUND: only ever returns a concrete point that violates on the net's own forward pass; the CALLER
% re-validates it (validate_witness + onnxruntime replay) before any 'sat' stands. Validated 15/15 (all
% ORT-confirmed) on the acasxu falsify gap -- see status-repo/research/acas_next_run_evidence/test_sound15.m.
function cex = falsify_batched_random(net, lb, ub, Hs, inputSize, N, doPolish, polMaxT, nearMiss)
    cex = [];
    lb = double(lb(:)); ub = double(ub(:)); d = numel(lb); sp = ub - lb;
    X = lb + sp .* rand(d, N);                                   % uniform samples in flat vnnlib order
    if isscalar(inputSize)                                       % flat feature net (lsnc): [d N] = 'CB'
        Y = reshape(extractdata(predict(net, dlarray(single(reshape(X, [inputSize N])), 'CB'))), [], N);
    else                                                         % image net (acasxu): [inputSize N] = 'SSCB'
        Y = reshape(extractdata(predict(net, dlarray(single(reshape(X, [inputSize N])), 'SSCB'))), [], N);
    end
    nH = numel(Hs); M = inf(1, N); bH = ones(1, N);
    for h = 1:nH                                                 % worst margin per sample over the OR-halfspaces
        mh = max(Hs(h).G * Y - Hs(h).g(:), [], 1);
        b = mh < M; M(b) = mh(b); bH(b) = h;
    end
    [Msort, ord] = sort(M, 'ascend');
    if Msort(1) <= 0                                             % a raw sample already violates -> witness
        k = ord(1); cex = {X(:,k); reshape(Y(:,k), [], 1)}; return;
    end
    if ~doPolish || Msort(1) >= nearMiss, return; end           % clearly unsat region (no near point) -> skip polish
    topK = min(30, N); target = -1e-3; tp = tic;
    for kk = 1:topK
        if toc(tp) > polMaxT, break; end                        % wall cap bounds cost on (true-unsat) instances
        x = X(:, ord(kk)); h = bH(ord(kk)); G = Hs(h).G; g = Hs(h).g(:);
        m0 = max(G * i_eval_flat(net, x, inputSize) - g); step = sp * 0.05;
        for it = 1:400
            improved = false;
            for dd = 1:d
                for s = [-1 1]
                    xt = x; xt(dd) = min(max(x(dd) + s*step(dd), lb(dd)), ub(dd));
                    mt = max(G * i_eval_flat(net, xt, inputSize) - g);
                    if mt < m0 - 1e-15, x = xt; m0 = mt; improved = true; end
                end
            end
            if m0 <= target, break; end
            if ~improved, step = step * 0.5; if max(step ./ max(sp, eps)) < 1e-10, break; end; end
            if toc(tp) > polMaxT, break; end
        end
        if m0 <= 0                                              % polished a sample past the unsafe boundary -> witness
            cex = {x; reshape(i_eval_flat(net, x, inputSize), [], 1)}; return;
        end
    end
end

% Single forward pass of a flat input through a low-dim dlnetwork; returns a column.
% Scalar inputSize -> flat feature net ('CB', e.g. lsnc [6]); 3-elem -> image net ('SSCB', acasxu).
function y = i_eval_flat(net, x, inputSize)
    if isscalar(inputSize)
        y = double(reshape(extractdata(predict(net, dlarray(single(reshape(x, [inputSize 1])), 'CB'))), [], 1));
    else
        y = double(reshape(extractdata(predict(net, dlarray(single(reshape(x, inputSize)), 'SSCB'))), [], 1));
    end
end

% Get random example from input set
function x = get_example(xRand,i)
    s = size(xRand);
    n = length(s);
    if n == 4
        x = xRand(:,:,:,i);
    elseif n == 3
        x = xRand(:,:,i);
    elseif n == 2
        x = xRand(:,i);
        xsize = size(x);
        if xsize(1) ~= 1 && ~isa(x,"dlarray")
            x = x';
        end
    else
        error("InputSize = "+string(s));
    end
end