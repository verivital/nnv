function [status, tTime] = run_vnncomp_instance(category, onnx, vnnlib, outputfile)

% single script to run all instances (supported) from the vnncomp2025

t = tic;
status = 2; % unknown (to start with)

% disp("We are running...")



%% 1) Load components

% Load networks

[net, nnvnet, needReshape, reachOptionsList, inputSize, inputFormat, nRand, falsifyOpts] = load_vnncomp_network(category, onnx, vnnlib);

if isempty(inputSize)
    inputSize = net.Layers(1, 1).InputSize;
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
% cannot soundly verify -- multi-network (equal-to/isomorphic-to), nonlinear/arithmetic
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
if ~isa(lb, "cell") && length(prop) == 1 % one input, one output 
    counterEx = falsify_single(net, lb, ub, inputSize, nRand, prop{1}.Hg, needReshape, inputFormat, falsifyOpts);
elseif isa(lb, "cell") && length(lb) == length(prop) % multiple inputs, multiple outputs
    for spc = 1:length(lb) % try parfeval, parfor does not work for early return
        counterEx = falsify_single(net, lb{spc}, ub{spc}, inputSize, nRand, prop{spc}.Hg, needReshape, inputFormat, falsifyOpts);
        if iscell(counterEx)
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

    % Pillar-2 final gate: replay the SAT witness through onnxruntime on the ORIGINAL
    % ONNX model (the competition's own checker) before committing to `sat`. A witness
    % that NNV's self-consistent validate_witness accepts could still fail the
    % competition if NNV's import permutes the input differently than the ONNX expects
    % (the suspected source of several of the 19 -150 penalties in 2025). If onnxruntime
    % is AVAILABLE and definitively says the witness violates NOTHING, downgrade to
    % `unknown` rather than risk a -150. If onnxruntime is UNAVAILABLE (not installed, or
    % a replay error), TRUST validate_witness and keep `sat` -- never suppress a SAT we
    % cannot disprove. Single-output-spec instances only (length(prop)==1) -- covers
    % BOTH the non-cell and the cell-lb (multiple input sets, single output spec) SAT
    % search paths above, which both falsify against prop{1}.Hg, the region the witness
    % hit. (Multi-output-spec instances aren't gated: which spec the witness hit isn't
    % tracked here; they still pass validate_witness.)
    if length(prop) == 1
        try
            [orVio, orAvail] = validate_witness_onnx(onnx, counterEx{1}, prop{1}.Hg);
            if orAvail && ~orVio
                fprintf('onnxruntime replay rejected the SAT witness -> unknown\n');
                status = 2; counterEx = nan;
            end
        catch
            % onnx guard could not run -> trust validate_witness, keep sat
        end
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

function [net,nnvnet,needReshape,reachOptionsList,inputSize,inputFormat,nRand,falsifyOpts] = load_vnncomp_network(category, onnx, vnnlib)
% load participating vnncomp 2025 benchmark NNs 
% Not yet supported:
% - cctsdb (some errrors when forward propagating)
% - lsnc_relu
% - traffic_signs_recognition (last year all instances were sat, maybe we are not wrong?)
% - collins aerospace (unsure of what is wrong, but we get invalid SAT instances)


    needReshape = 0; % default is to use MATLAB reshape, otherwise use the python reshape
    % reachOptions = struct;
    % reachOptions.reachMethod = 'approx-star'; % default parameters
    numCores = feature('numcores'); % in case we select exact method
    inputSize = [];
    inputFormat = "default"; % no need to change for most of them, but needed for some ("UU")
    nRand = 100; % default from previous competitions
    % Per-category falsification budget (VNNCOMP2026 tuning). Default preserves the old
    % hardcoded struct('seed',0,'max_time',5); the ftab at the END of this function
    % overrides n_restarts/n_steps/max_time and nRand per category.
    falsifyOpts = struct('n_restarts',20, 'n_steps',40, 'lr',0.1, 'fgsm',true, ...
                         'max_time',5, 'seed',0);

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
        

    elseif contains(category, "cctsdb_yolo")
        % net = importNetworkFromONNX(onnx);
        % nnvnet = "";
        % inputSize = [12296, 1];
        % inputFormat = "UU";
        % X = dlarray(rand(12296, 1), inputFormat);
        % net = initialize(net, X);
        % reachOptions = struct;
        % reachOptions.reachMethod = 'cp-star';
        % reachOptions.inputFormat = inputFormat;
        % reachOptionsList{1} = reachOptions;
        error("Working on supporting this one");

    elseif contains(category, "cersyve")
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BC");
        nnvnet = matlab2nnv(net);
        % reachOptions = struct;
        % reachOptions.reachMethod = 'approx-star'; % default parameters
        % reachOptionsList{1} = reachOptions;
        reachOptions.reachMethod = 'cp-star';
        % reachOptions.numCores = numCores;
        reachOptionsList{1} = reachOptions;

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
            reachOptions = struct;
            reachOptions.reachMethod = 'approx-star';
            reachOptionsList{1} = reachOptions;
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
        net = importNetworkFromONNX(onnx, "InputDataFormats","BCSS");
        nnvnet = matlab2nnv(net);
        reachOptions.train_epochs = 100;
        reachOptions.train_lr = 0.001;
        reachOptions.coverage = 0.999;
        reachOptions.confidence = 0.999;
        reachOptions.train_mode = 'Linear';
        reachOptions.surrogate_dim = [];
        reachOptions.threshold_normal = 1e-5;
        reachOptions.reachMethod = 'cp-star';
        reachOptionsList{1} = reachOptions;
        needReshape = 1; % 2 and 0 return invalid sat instances

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
            reachOptions = struct;
            reachOptions.reachMethod = 'relax-star-area';
            reachOptions.relaxFactor = 0.5;
            reachOptionsList{1} = reachOptions;
            reachOptions = struct;
            reachOptions.reachMethod = 'approx-star'; % default parameters
            reachOptionsList{2} = reachOptions;
        else
            % Was a 0.9-then-0.7 double-write to {1} (only 0.7 ever ran). Keep the
            % single sound relax-star-area@0.9; slow_cats prepends approx-zono+abs-dom.
            reachOptions = struct;
            reachOptions.reachMethod = 'relax-star-area';
            reachOptions.relaxFactor = 0.9;
            reachOptionsList{1} = reachOptions;
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
        % Cheap-to-precise SOUND ladder: approx-star first, looser relax-star fallback.
        % (Was approx-star at {1} then exact-star OVERWRITING {1}, so only the
        % exponential exact-star ran and stalled to 'unknown'. Drop exact-star: it
        % blows up on these nets and never finished in the per-instance budget.)
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star';
        reachOptionsList{1} = reachOptions;
        reachOptions = struct;
        reachOptions.reachMethod = 'relax-star-area';
        reachOptions.relaxFactor = 0.7;
        reachOptionsList{2} = reachOptions;

    elseif contains(category, "linearize")
        % 
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape
        try 
            nnvnet = matlab2nnv(net);
            % Sound approx-star -> looser relax-star fallback (was approx-star at {1}
            % then exact-star OVERWRITING {1}; only the exponential exact-star ran).
            reachOptions = struct;
            reachOptions.reachMethod = 'approx-star';
            reachOptionsList{1} = reachOptions;
            reachOptions = struct;
            reachOptions.reachMethod = 'relax-star-area';
            reachOptions.relaxFactor = 0.7;
            reachOptionsList{2} = reachOptions;
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
        if contains(onnx, "lindex")
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
            elseif contains(onnx, "mscn")
                % if contains(onnx, "dual")
                %     inputSize = [1,22,14];
                %     inputFormat = "UUU";
                %     X = dlarray(rand(1,22,14), inputFormat);
                % else
                %     needReshape = 2;
                %     inputSize = [1,11,14];
                %     inputFormat = "UUU";
                %     X = dlarray(rand(1,11,14), inputFormat);
                % end
                error("These are not supported yet.")
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
        reachOptions = struct;
        reachOptions.reachMethod = 'relax-star-area';
        reachOptions.relaxFactor = 0.9;
        reachOptionsList{1} = reachOptions;
        reachOptions = struct;
        reachOptions.reachMethod = 'approx-star'; % default parameters
        reachOptionsList{2} = reachOptions;

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
        inputSize = [3 30 30];   % reshape size BEFORE the [3 2 1] permute
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
        % vit: onnx to matlab. Try sound Star-set reach first (approx-star
        % then relax-star fallback), since cp-star is probabilistic and
        % doesn't exercise the actual reach() chain through attention layers.
        % If matlab2nnv conversion fails (e.g., unsupported layer in the ONNX
        % import), fall back to cp-star to get *some* result.
        net = importNetworkFromONNX(onnx, "InputDataFormats", "BCSS", 'OutputDataFormats',"BC");
        try
            nnvnet = matlab2nnv(net);
            needReshape = 1;
            reachOptionsList = {};
            % First: approx-star (tightest sound bound NNV offers without LP)
            reachOpts1 = struct;
            reachOpts1.reachMethod = 'approx-star';
            reachOptionsList{end+1} = reachOpts1;
            % Second: relax-star (faster fallback for deeper nets)
            reachOpts2 = struct;
            reachOpts2.reachMethod = 'relax-star-area';
            reachOpts2.relaxFactor = 0.5;
            reachOptionsList{end+1} = reachOpts2;
        catch
            nnvnet = "";
            % cp-star fallback ONLY when sound conversion failed
            reachOpts = struct;
            reachOpts.train_epochs = 100;
            reachOpts.train_lr = 0.001;
            reachOpts.coverage = 0.999;
            reachOpts.confidence = 0.999;
            reachOpts.train_mode = 'Linear';
            reachOpts.surrogate_dim = [-1,-1];
            reachOpts.threshold_normal = 1e-5;
            reachOpts.reachMethod = 'cp-star';
            reachOptionsList = {reachOpts};
            needReshape = 1;
        end

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
        "acasxu",            40, 60, 5,   200
        "sat_relu",          50, 80, 5,   200
        "cersyve",           30, 50, 4,   150
        "lsnc_relu",         50, 80, 4,   500   % manifest: no PGD -> random is primary
        "relusplitter",      30, 50, 3,   150
        "dist_shift",        30, 60, 4,   150
        "linearize",         30, 60, 4,   150
        "tllverifybench",    40, 80, 5,   200
        "collins_rul",       30, 60, 4,   150
        "collins_aerospace", 20, 40, 3,   150
        "nn4sys",            25, 50, 3,   150
        "safenlp",           30, 60, 4,   200
        "malbeware",         20, 40, 3,   150
        "ml4acopf",          25, 50, 3,   200
        "cora",              30, 60, 4,   200
        "metaroom",          15, 30, 2,   100
        "cifar100",          5,  15, 2,   30
        "tinyimagenet",      5,  15, 2,   30
        "vggnet",            4,  12, 1.5, 20
        "traffic_signs",     20, 40, 5,   400   % manifest: no PGD -> random is primary
        "cctsdb_yolo",       30, 50, 5,   150
        "yolo",              20, 40, 3,   100
        "vit",               10, 25, 3,   50
        "cgan",              20, 40, 3,   100
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