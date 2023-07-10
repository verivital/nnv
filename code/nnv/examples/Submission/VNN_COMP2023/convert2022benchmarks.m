%% Run some verification instances from the competition

% path to benchmarks
% vnncomp_path = "/home/manzand/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/";
vnncomp_path = "/home/dieman95/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/";

% Go through some of the instances for every benchmark
% It looks like unless the networks are fully fully supported, we get errors...
% Same problem we ra into when running with no GUI ...
% 
% Error using matlab.internal.indentcodeLegacy
%This feature is not supported because:
%Swing is not currently available.
%
%Error in indentcode (line 42)
%        indentedText = matlab.internal.indentcodeLegacy(text, language);
%
%Error in nnet.internal.cnn.onnx.fcn.ModelTranslation/genOutputFcn (line 100)
%            code = indentcode(code);                                                    % Indent the code
%
%Error in nnet.internal.cnn.onnx.fcn.ModelTranslation (line 39)
%            this.ModelFunctionCode        = genOutputFcn(this, outputFcnName, this.ToplevelGraphTranslation, isGeneratingCustomLayerOPTIONAL);
%
%Error in nnet.internal.cnn.onnx.CustomLayerManager/genCustomLayerFunction (line 499)
%            modelTranslation = nnet.internal.cnn.onnx.fcn.ModelTranslation(modelProto, modelFcnPath, modelFcnName, maxNameLength, true);
%
%Error in nnet.internal.cnn.onnx.CustomLayerManager (line 77)
%                this.ModelTranslations(i) = genCustomLayerFunction(this, this.ModelProtos(i), nameStem);
%
%Error in nnet.internal.cnn.onnx.ModelTranslationIntoLayers (line 228)
%                this.CustomLayerManager = %nnet.internal.cnn.onnx.CustomLayerManager(this.GraphProtoManager, ...
%
%Error in nnet.internal.cnn.onnx.importONNXNetwork (line 32)
%modelTranslationIntoLayers = nnet.internal.cnn.onnx.ModelTranslationIntoLayers(modelProto, ...
%
%Error in importONNXNetwork (line 113)
%Network = nnet.internal.cnn.onnx.importONNXNetwork(modelfile, varargin{:});
%
%Error in run_vnncomp2023_instance>load_vnncomp_network (line 127)
%        net = importONNXNetwork(onnx, "InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape


%% cifar_biasfield

disp("Running cifar_biasfield...")

biasfield_path = vnncomp_path + "cifar_biasfield/onnx/";

onnx = dir(biasfield_path);

for i=3:length(onnx)
    onnxFile = biasfield_path + onnx(i).name;
    net = load_vnncomp_network('cifar_biasfield', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name], "net", "-v7.3");
end

%% rl_benchmarks

disp("Running rl_benchmarks...")

rl_path = vnncomp_path + "rl_benchmarks/onnx/";

onnx = dir(rl_path);

for i=3:length(onnx)
    onnxFile = rl_path + onnx(i).name;
    net = load_vnncomp_network('rl_benchmarks', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name], "net", "-v7.3");
end

%% reach_prob_density

disp("Running reach_prob_density...")

reach_prob_path = vnncomp_path + "reach_prob_density/onnx/";

onnx = dir(reach_prob_path);

for i=3:length(onnx)
    onnxFile = reach_prob_path + onnx(i).name;
    net = load_vnncomp_network('reach_prob_density', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name], "net", "-v7.3");
end

%% mnist_fc

disp("Running mnist_fc...")

mnist_path = vnncomp_path + "mnist_fc/onnx/";

onnx = dir(mnist_path);

for i=3:length(onnx)
    onnxFile = mnist_path + onnx(i).name;
    net = load_vnncomp_network('mnist_fc', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name], "net", "-v7.3");
end

%% sri_resnet_a

disp("Running sri_resnets..")

sri_path = vnncomp_path + "sri_resnet_a/onnx/";

onnx = dir(sri_path);

for i=3:length(onnx)
    onnxFile = sri_path + onnx(i).name;
    net = load_vnncomp_network('sri_resnet_a', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name, '.mat'], "net", "-v7.3");
end

%% sri_resnet_b

sri_path = vnncomp_path + "sri_resnet_b/onnx/";

onnx = dir(sri_path);

for i=3:length(onnx)
    onnxFile = sri_path + onnx(i).name;
    net = load_vnncomp_network('sri_resnet_b', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name, '.mat'], "net", "-v7.3");
end

%% cifar2020

cifar2020_path = vnncomp_path + "cifar2020/onnx/";

onnx = dir(cifar2020_path);

for i=3:length(onnx)
    onnxFile = cifar2020_path + onnx(i).name;
    if ~contains(onnxFile, "PGD")
        net = load_vnncomp_network('cifar2020', onnxFile);
        name = onnx(i).name;
        name = name(1:end-5);
        save(['networks2023/', name], "net", "-v7.3");
    end
end

%% oval21

oval_path = vnncomp_path + "oval21/onnx/";

onnx = dir(oval_path);

for i=3:length(onnx)
    onnxFile = oval_path + onnx(i).name;
    net = load_vnncomp_network('oval21', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name], "net", "-v7.3");
end

%% cifar100_tinyimagenet_resnet 

cifar_path = vnncomp_path + "cifar100_tinyimagenet_resnet/onnx/";

onnx = dir(cifar_path);

for i=3:length(onnx)
    onnxFile = cifar_path + onnx(i).name;
    net = load_vnncomp_network('cifar100_tinyimagenet_resnet', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name], "net", "-v7.3");
end


%% nn4sys

disp("Running nn4sys...")

nn4sys_path = vnncomp_path + "nn4sys/onnx/";

onnx = dir(nn4sys_path);

for i=3:length(onnx)
    onnxFile = nn4sys_path + onnx(i).name;
    if contains(onnxFile, "lindex")
        net = load_vnncomp_network('nn4sys', onnxFile);
        name = onnx(i).name;
        name = name(1:end-5);
%         save(['networks2023/', name], "net", "-v7.3");
    end
end


%% tllverify

disp("Running tllverify..")

tll_path = vnncomp_path + "tllverifybench/onnx/";

onnx = dir(tll_path);

for i=3:length(onnx)
    onnxFile = tll_path + onnx(i).name;
    net = load_vnncomp_network('tllverifybench', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
%     save(['networks2023/', name], "net", "-v7.3");
end


%% Loading function
function net = load_vnncomp_network(category, onnx)
% load vnncomp 2022 benchmark NNs (subset support)

    if contains(category, "nn4sys")
        % nn4sys: onnx to matlab:
        net = importONNXNetwork(onnx, "OutputDataFormats", "BC"); % lindex (no reshape)
        
    elseif contains(category, "tllverify")
        % tllverify: onnx to matlab
        net = importONNXNetwork(onnx,"InputDataFormats", "BC", 'OutputDataFormats',"BC"); % no reshape
        
    elseif contains(category, 'biasfield')
        % onnx to nnv
        if contains(onnx, "base")
            net = importONNXNetwork(onnx); % needs reshape
        else
            net = importONNXNetwork(onnx,  "InputDataFormats", "BC"); % does not need reshape (feature input)
        end

    elseif contains(category, 'rl_bench')
        % rl: onnx to nnv
        net = importONNXNetwork(onnx, "InputDataFormats", "BC", "OutputDataFormats", "BC"); % no need to reshape (feature input)

    elseif contains(category, 'reach_prob')
        % reach_prob: onnx to nnv
        net = importONNXNetwork(onnx, "InputDataFormats", "BC"); % no reshape

    elseif contains(category, 'mnist_fc')
        % mnist: onnx to nnv
        net = importONNXNetwork(onnx, 'InputDataFormats',"SSC", 'OutputDataFormats',"BC"); % no reshape

    elseif contains(category, 'sri_resnet')
        net = importONNXNetwork(onnx, 'OutputDataFormats',"BC"); % reshape

    elseif contains(category, 'cifar100')
        net = importONNXNetwork(onnx,  'OutputDataFormats',"BC"); %reshape

    elseif contains(category, 'oval')
        net = importONNXNetwork(onnx, 'OutputDataFormats',"BC"); %reshape

    elseif contains(category, 'cifar2020')
        % onnx to nnv
        if ~contains(onnxFile, "PGD")
            net = importONNXNetwork(onnx, 'OutputDataFormats',"BC"); %reshape
        else
            net = importONNXNetwork(onnx, "InputDataFormats", "BCSS", "OutputDataFormats", "BC"); % reshape?
        end
        
    else % all other benchmarks
        % traffic: onnx to matlab: opset15 issues
        error("ONNX model not supported")
    end

end