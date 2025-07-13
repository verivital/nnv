clear

clc


instances = {
"onnx/cGAN_imgSz32_nCh_1.onnx",	"vnnlib/cGAN_imgSz32_nCh_1_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_1.onnx",	"vnnlib/cGAN_imgSz32_nCh_1_prop_1_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_1.onnx",	"vnnlib/cGAN_imgSz32_nCh_1_prop_2_input_eps_0.020_output_eps_0.025.vnnlib";
"onnx/cGAN_imgSz32_nCh_1.onnx",	"vnnlib/cGAN_imgSz32_nCh_1_prop_3_input_eps_0.020_output_eps_0.025.vnnlib";
"onnx/cGAN_imgSz32_nCh_3.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";
"onnx/cGAN_imgSz32_nCh_3.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_prop_1_input_eps_0.015_output_eps_0.020.vnnlib";
"onnx/cGAN_imgSz32_nCh_3.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_prop_2_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_3.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_prop_3_input_eps_0.015_output_eps_0.020.vnnlib";
"onnx/cGAN_imgSz64_nCh_1.onnx",	"vnnlib/cGAN_imgSz64_nCh_1_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz64_nCh_1.onnx",	"vnnlib/cGAN_imgSz64_nCh_1_prop_1_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz64_nCh_1.onnx",	"vnnlib/cGAN_imgSz64_nCh_1_prop_2_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz64_nCh_1.onnx",	"vnnlib/cGAN_imgSz64_nCh_1_prop_3_input_eps_0.005_output_eps_0.010.vnnlib";
"onnx/cGAN_imgSz64_nCh_3.onnx",	"vnnlib/cGAN_imgSz64_nCh_3_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz64_nCh_3.onnx",	"vnnlib/cGAN_imgSz64_nCh_3_prop_1_input_eps_0.005_output_eps_0.010.vnnlib";
"onnx/cGAN_imgSz64_nCh_3.onnx",	"vnnlib/cGAN_imgSz64_nCh_3_prop_2_input_eps_0.005_output_eps_0.010.vnnlib";
"onnx/cGAN_imgSz64_nCh_3.onnx",	"vnnlib/cGAN_imgSz64_nCh_3_prop_3_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_3_nonlinear_activations.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_nonlinear_activations_prop_0_input_eps_0.015_output_eps_0.020.vnnlib";
"onnx/cGAN_imgSz32_nCh_1_transposedConvPadding_1.onnx",	"vnnlib/cGAN_imgSz32_nCh_1_transposedConvPadding_1_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_3_upsample.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_upsample_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_3_small_transformer.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_small_transformer_prop_0_input_eps_0.010_output_eps_0.015.vnnlib";
"onnx/cGAN_imgSz32_nCh_3_small_transformer.onnx",	"vnnlib/cGAN_imgSz32_nCh_3_small_transformer_prop_1_input_eps_0.010_output_eps_0.015.vnnlib"    
};


for i = 20:21
    
    if i<20
        method = 'd';
    else
        method = 'p';
    end
    
    status = -3;

    % Load specification
    property = load_vnnlib(instances{i, 2});
    lb = property.lb;     % Input lower bounds (5x1)
    ub = property.ub;     % Input upper bounds (5x1)
    prop = property.prop; % Output specification

    % Counterexample search via sampling
    t = tic;
    nR = 100;
    IB = Box(lb, ub);
    xRand = IB.sample(nR - 2);
    xRand = [lb, ub, xRand];          % Include bounds as edge cases
    xRand = dlarray(xRand, "CB");     % Features x Batch format

    if strcmp(method, 'd')
        net = importNetworkFromONNX( instances{i,1}, "InputDataFormats", "BC", "OutputDataFormats", "BC" );
        nnvnet = matlab2nnv(net);
    else
        netFcn = importONNXFunction(instances{i,1}, "MyNet");
        nnvnet = netFcn;
    end

    for k = 1:nR
        x = xRand(:, k);
        yPred = predict(net, x);
        yPred = extractdata(yPred);
        yPred = reshape(yPred, [], 1);  % Ensure column vector
        Hs = prop{1}.Hg;
        if Hs.contains(double(yPred))
            status = 0;
            counterEx = {extractdata(x); yPred};
            break;
        end
    end
    t = toc(t);

    % Construct input set for reachability
    IS = Star(lb, ub);

    % Reachability analysis
    if status ~= 0
        if strcmp( method , 'd')
            reachOptions = struct;
            reachOptions.reachMethod = 'approx-star';
            ySet = nnvnet.reach(IS, reachOptions);
            
        elseif strcmp( method, 'p')
            train_epochs = 200;
            train_lr = 0.0001;
            reachOptions.train_epochs = 200;
            reachOptions.train_lr = 0.0001;
            reachOptions.dims = [-1 -1];
            reachOptions.coverage = 0.9999;
            reachOptions.confidence = 0.9995;
            reachOptions.train_mode = 'Linear';
            reachOptions.surrogate_dim = [];
            reachOptions.threshold_normal = 1e-5;

            ySet = Prob_reach(nnvnet, IS, reachOptions);

        end

    end

    % Verify output against specification
    if status ~= 0
        status = verify_specification(ySet, prop);
    end

    % Fallback to exact method if inconclusive
    if status == 2
        reachOptions = struct;
        reachOptions.reachMethod = 'exact-star';
        reachOptions.numCores = 4;
        ySet = nnvnet.reach(IS, reachOptions);
        status = verify_specification(ySet, prop);
    end

    % Print results
    disp("Instance (" + string(i) + ") finished!! Result: " + string(status) + ...
        " , verification time: " + string(t) + " seconds");

end




% % Load instances
% instances = readmatrix('instances.csv', 'OutputType','string');
% 
% 
% for i=1:height(instances)
% 
%     status = -3;
%     net = importNetworkFromONNX(fullfile(onnxfiles(i).folder,onnxfiles(i).name),     "InputDataFormats","BC","OutputDataFormats","BC");
%     nnvnet = matlab2nnv(net);
% 
%     % Load specification
%     property = load_vnnlib(instances(i,2));
%     lb = property.lb; % input lower bounds
%     ub = property.ub; % input upper bounds
%     prop = property.prop; % output spec to verify
% 
%     % Counterexample search
%     t = tic;
%     nR = 100; % number of samples
%     IB = Box(lb,ub);
%     xRand = IB.sample(nR-2);
%     xRand = [lb, ub, xRand];
%     xRand = dlarray(xRand, "CB");
%     for k=1:nR
%         x = xRand(:,k);
%         yPred = predict(net, x);
%         yPred = extractdata(yPred);
%         yPred = reshape(yPred, [], 1); % convert to column vector (if needed)
%         Hs = prop{1}.Hg;
%         if Hs.contains(double(yPred)) % property violated
%             status = 0;
%             counterEx = {extractdata(x); yPred};
%             break;
%         end
%     end
% 
%     t = toc(t);
%     % disp("Counterexample search finished in "+string(t)+" seconds");
% 
%     % if status == 0
%         % disp("-------     Counterexample found");
%         % disp(counterEx{1});
%         % disp(counterEx{2});
%     % end
% 
%     % Create Input Set
%     IS = Star(lb,ub);
% 
%     % Define reachability options
%     reachOptions = struct;
%     % reachOptions.reachMethod = 'relax-star-area';
%     % reachOptions.relaxFactor = 0.0;
%     % reachOptions.reachMethod = 'exact-star';
%     reachOptions.reachMethod = 'approx-star';
% 
% 
%     % Compute reachability
%     t = tic;
%     if status~=0
%         ySet = nnvnet.reach(IS, reachOptions);
%     end
%     t = toc(t);
% 
%     % Verify property
%     if status~=0
%         status = verify_specification(ySet, prop);
%     end
% 
%     if status == 2
%         reachOptions = struct;
%         reachOptions.reachMethod = 'exact-star';
%         reachOptions.numCores = 4;
%         ySet = nnvnet.reach(IS, reachOptions);
%         status = verify_specification(ySet, prop);
%     end
% 
% 
%     % Print results
%     disp("Instance ("+string(i)+") finished!! Result: "+string(status)+" , verification time: "+string(t));
%     % disp("========================")
% 
% end

% Notes:
% some counterexamples
% most unknowns with relax-star (1)
% 