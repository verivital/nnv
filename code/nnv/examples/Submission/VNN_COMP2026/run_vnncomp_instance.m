function [status, tTime] = run_vnncomp_instance(category, onnx, vnnlib, outputfile)

% single script to run all instances (supported) from the vnncomp2025

t = tic;
status = 2; % unknown (to start with)

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

[net, nnvnet, needReshape, reachOptionsList, inputSize, inputFormat, nRand, falsifyOpts] = i_load_vnncomp_network_cached(category, onnx, vnnlib);

if isempty(inputSize)
    inputSize = net.Layers(1, 1).InputSize;
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
if isprop(net, 'Layers') && ~isempty(net.Layers)
    L1 = net.Layers(1);
    if isa(L1, 'nnet.cnn.layer.ImageInputLayer') || isa(L1, 'nnet.cnn.layer.Image3DInputLayer')
        useImageStar = true;
    elseif isa(L1, 'nnet.cnn.layer.FeatureInputLayer') || isa(L1, 'nnet.onnx.layer.FeatureInputLayer')
        useImageStar = false;
    end
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
    fprintf('vnnlib 2.0 multi-network (equal-to) -> status=%d (0 sat / 1 unsat / 2 unknown)\n', status);
    tTime = toc(t);
    fid = fopen(outputfile, 'w');
    if status == 0
        fprintf(fid, 'sat \n');
        fclose(fid);
        % witness = stacked product-net I/O [x_f; x_g] -> [f(x_f); f(x_g)] (the
        % order verify_multinet's product net uses); sound -- the verdict, not the
        % witness format, is what the sweep scores, and the CE is already validated.
        xf = counterEx{1}; xg = counterEx{2};
        yf = nnvnet.evaluate(xf); yg = nnvnet.evaluate(xg);
        write_counterexample(outputfile, {[xf(:); xg(:)], [yf(:); yg(:)]});
    elseif status == 1
        fprintf(fid, 'unsat \n');
        fclose(fid);
    else
        fprintf(fid, 'unknown \n');
        fclose(fid);
    end
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
    isManifest = isfile(regexprep(char(onnx), '\.onnx$', '.nnv.mat'));   % python-imported -> divergence risk
    try
        [orVio, orAvail] = validate_witness_onnx(onnx, counterEx{1}, gateHs);
        if orAvail && ~orVio
            % onnxruntime ran and the witness violates NOTHING on the real ONNX -> spurious -> drop.
            fprintf('onnxruntime replay rejected the SAT witness -> unknown\n');
            status = 2; counterEx = nan;
        elseif ~orAvail && isManifest
            % No independent onnxruntime confirmation AND the net came from a (python) manifest import
            % that can diverge from the real ONNX -> we cannot stand behind this `sat`. Sound-or-unknown
            % rather than risk a -150. (onnxruntime is a stated requirement; this is the safety net if it
            % is somehow missing.) Standard MATLAB-imported nets keep `sat` (validate_witness is reliable).
            fprintf('SAT witness unconfirmed (no onnxruntime) on a manifest-imported net -> unknown\n');
            status = 2; counterEx = nan;
        end
    catch
        % Gate itself errored. Manifest nets are the divergence risk -> prefer unknown; standard
        % MATLAB-imported nets fall back to validate_witness (which already accepted the witness).
        if isManifest, status = 2; counterEx = nan; end
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

            while ~isempty(reachOptionsList)

                reachOptions = reachOptionsList{1};

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
                if status == 1 && strcmp(reachOptions.reachMethod, "cp-star")
                    % cp-star is PROBABILISTIC (conformal prediction): this 'holds'/unsat
                    % is NOT a sound proof. It is the primary reach method for some
                    % categories (cifar100, ml4acopf, ...) and a sound-methods-first last
                    % resort for others (cersyve). Either way, LOG it so the per-benchmark
                    % -150 exposure from a probabilistic verdict is tracked (see
                    % PROGRESS_LOG cp-star policy).
                    fprintf('VERDICT VIA CP-STAR (probabilistic, not a sound proof): %s\n', onnx);
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
        
        parfor spc = 1:length(lb) % We can compute these in parallel for faster computation

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
        
        parfor spc = 1:length(lb) % We can compute these in parallel for faster computation

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

                        % Add verification status
                        tempStatus = verify_specification(ySet, prop(spc));
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

% if status == 2 && strcmp(reachOptions.reachMethod, 'exact-star')
%     status = 0;
% end

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
    if ispc
        py = 'python';
    else
        py = 'python3';
    end
    try
        pe = pyenv;
        if ~isempty(pe.Executable) && isfile(pe.Executable)
            py = ['"' char(pe.Executable) '"'];
        end
    catch
    end
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
    isConv = contains(category,"cifar100") || contains(category,"tinyimagenet");
    if ~(isFC || isConv), return; end                    % only the gated timeout categories
    if ~is_nnvnet_valid(nnvnet), return; end             % need a valid NNV net for nn_to_ops + evaluate
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
    try
        if isConv && gpuDeviceCount >= 1
            [gv, ginfo] = gpu_bab_try_verify(nnvnet, lb, ub, prop, ...
                struct('engine','batched','maxNodes',64,'device','gpu','allowUnsoundSingle',true));
            if strcmp(gv, 'robust')
                [gv2, gi2] = gpu_bab_try_verify(nnvnet, lb, ub, prop, struct('engine','batched','maxNodes',64));  % CPU double = sound emit
                if strcmp(gv2, 'robust')
                    status = 1; reachOptionsList = {};
                    fprintf('GPU-BaB pre-check: robust/unsat (gpu-screen + %d-node double-confirm) -> skip Star\n', gi2.nodes);
                else
                    fprintf('GPU-BaB pre-check: gpu-screen robust but double-confirm=%s -> Star reach\n', gv2);
                end
            else
                fprintf('GPU-BaB pre-check: %s (gpu-screen, %s) -> Star reach\n', gv, ginfo.reason);
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
        % acasxu: onnx to nnv
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS");
        nnvnet = matlab2nnv(net);
        if ~contains(vnnlib, "prop_3.") && ~contains(vnnlib, "prop_4.")
            reachOptions.reachMethod = 'exact-star';
            reachOptions.numCores = 1;  % Phase 1.3: avoid nested-parpool errors when called inside parfeval workers
            reachOptionsList{1} = reachOptions;
            nRand = 500;
        else
            reachOptions = struct;
            reachOptions.reachMethod = 'approx-star'; % default parameters
            reachOptionsList{1} = reachOptions;
            reachOptions.reachMethod = 'exact-star';
            reachOptions.numCores = 1;  % Phase 1.3: avoid nested-parpool errors when called inside parfeval workers
            reachOptionsList{2} = reachOptions;
        end

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
        % 
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape
        try 
            nnvnet = matlab2nnv(net);
            % 120/120 unknown today = bounds too loose (reach completes <2s, not a
            % timeout). Tightness ladder relax 0.5 -> 0.25 -> approx-star (full LP,
            % tightest). All sound; any rung that certifies is a pure +10.
            reachOptionsList = relax_ladder([0.5, 0.25]);
            oo = struct(); oo.reachMethod = 'approx-star';
            reachOptionsList{end+1} = oo;
        catch
            nnvnet = "";
            reachOptions = struct;
            reachOptions.reachMethod = 'cp-star';
            reachOptionsList{1} = reachOptions;
        end

    elseif contains(category, "lsnc_relu")
        % lsnc_relu: MATLAB's importNetworkFromONNX cannot parse this model
        % (IR/opset), so load via the Python-importer manifest
        % (tools/onnx2nnv_python/onnx2nnv.py writes <model>.nnv.mat alongside
        % the ONNX). Cross-validated against onnxruntime: max diff 9.2e-07
        % over random inputs (2026-06-09). Flat [6] feature input, [8] output.
        nnvnet = load_manifest_net(onnx);
        net = nnvnet;   % falsify_single dispatches NN.evaluate for NNV nets
        inputSize = 6;
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
        % ml4acopf: onnx to matlab
        net = importNetworkFromONNX(onnx, "InputDataFormats","BC", "OutputDataFormats","BC");
        nnvnet = "";
        reachOptions = struct;
        reachOptions.train_epochs = 500;
        reachOptions.train_lr = 0.0001;
        reachOptions.coverage = 0.999;
        reachOptions.confidence = 0.999;
        reachOptions.train_mode = 'Linear';
        reachOptions.surrogate_dim = [-1,-1];
        reachOptions.threshold_normal = 1e-5;
        reachOptions.reachMethod = "cp-star";
        reachOptionsList{1} = reachOptions;
        % inputFormat = "BC";
        % error("Not supported");

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
            reachOptions = struct;
            reachOptions.inputFormat = inputFormat;
            reachOptions.reachMethod = 'cp-star'; % default parameters
            reachOptionsList{1} = reachOptions;
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

    % Phase 1.5 (TODO_VNNCOMP25_V01): for any category whose first reach
    % method is "cp-star" (probabilistic, requires GPU on this machine),
    % prepend sound CPU-only methods so the dispatcher tries them first.
    % cp-star remains as fallback when the sound methods can't decide.
    if ~isempty(reachOptionsList) && isfield(reachOptionsList{1}, 'reachMethod') ...
            && strcmp(reachOptionsList{1}.reachMethod, 'cp-star') ...
            && is_nnvnet_valid(nnvnet)
        sound_opts = cell(1,2);
        o1 = struct(); o1.reachMethod = 'approx-star';
        sound_opts{1} = o1;
        o2 = struct(); o2.reachMethod = 'relax-star-area'; o2.relaxFactor = 0.5;
        sound_opts{2} = o2;
        reachOptionsList = [sound_opts, reachOptionsList];
    end

    % Compute-bound categories (large CNNs / NLP / a SAT-encoding) time out under
    % exact-star -- and even approx-star -- within a realistic per-instance budget.
    % Lead with fast, looser, but STILL-SOUND over-approximations (zonotope, then
    % abstract-domain) and DROP exact-star, so they at least produce a sound verdict
    % instead of timing out; the tighter per-category methods stay as fallbacks.
    % (Per the timeout-model guidance: never exact-star for these.)
    slow_cats = ["cifar100","cora","safenlp","sat_relu","tinyimagenet","vggnet"];
    if any(contains(category, slow_cats)) && is_nnvnet_valid(nnvnet)
        kept = {};
        for k = 1:numel(reachOptionsList)
            if ~strcmp(reachOptionsList{k}.reachMethod, 'exact-star')
                kept{end+1} = reachOptionsList{k}; %#ok<AGROW>
            end
        end
        zo = struct(); zo.reachMethod = 'approx-zono';
        ad = struct(); ad.reachMethod = 'abs-dom';
        reachOptionsList = [{zo, ad}, kept];
    end

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
        "acasxu",            60, 80, 30,  200
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
        "traffic_signs",     20, 40, 5,   400   % manifest: no PGD -> random is primary
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
            S = load(p, 'C'); C = S.C;
            if isfield(C,'onnx_bytes') && isequal(C.onnx_bytes, d.bytes) ...
                    && abs(C.onnx_datenum - d.datenum) < 1e-9 ...
                    && isfield(C,'nnv_ver') && strcmp(C.nnv_ver, i_nnv_ver())
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
            C.onnx_bytes=d.bytes; C.onnx_datenum=d.datenum; C.nnv_ver=i_nnv_ver();
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

function xRand = create_random_examples(net, lb, ub, nR, inputSize, needReshape,inputFormat)
    xB = Box(lb, ub); % lb, ub must be vectors
    xRand = xB.sample(nR-2);
    xRand = [lb, ub, xRand];
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
    xRand = create_random_examples(net, lb, ub, nRand, inputSize, needReshape, inputFormat);
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