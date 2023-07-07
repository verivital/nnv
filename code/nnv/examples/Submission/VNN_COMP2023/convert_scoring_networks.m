%% Run some verification instances from the competition

% path to benchmarks
% vnncomp_path = "/home/manzand/Documents/MATLAB/vnncomp2023_benchmarks/benchmarks/";
vnncomp_path = "/home/dieman95/Documents/MATLAB/vnncomp2023_benchmarks/benchmarks/";

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



%% acasxu

disp("Running acas xu...")

acas_path = vnncomp_path + "acasxu/onnx/";

onnx = dir(acas_path);

for i=3:length(onnx)
    onnxFile = acas_path + onnx(i).name;
    net = load_vnncomp_network('acasxu', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name], "net", "-v7.3");
end

%% ViT

disp("Running ViT...")

vit_path = vnncomp_path + "vit/onnx/";

onnx = dir(vit_path);

for i=3:length(onnx)
    onnxFile = vit_path + onnx(i).name;
    net = load_vnncomp_network('vit', onnxFile);
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
        save(['networks2023/', name], "net", "-v7.3");
    end
end


%% dist_shift

disp("Running dist_shift...")

dist_path = vnncomp_path + "dist_shift/onnx/";

onnxFile = dist_path + "mnist_concat.onnx";
net = load_vnncomp_network('dist_shift', onnxFile);
name = 'mnist_concat.onnx';
name = name(1:end-5);
save(['networks2023/', name], "net", "-v7.3");


%% traffic_sign
% we should be able to support it, but matlab does not, so we would have to create the models manually...
% can load most info with Keras importer

%% collins_rul

disp("Running collins_rul..")

rul_path = vnncomp_path + "collins_rul_cnn/onnx/";

onnx = dir(rul_path);

for i=3:length(onnx)
    onnxFile = rul_path + onnx(i).name;
    net = load_vnncomp_network('collins_rul_cnn', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name], "net", "-v7.3");
end

%% cgan

disp("Running cgan..")

cgan_path = vnncomp_path + "cgan/onnx/";

onnx = dir(cgan_path);

for i=3:length(onnx)
    onnxFile = cgan_path + onnx(i).name;
    if ~contains(onnxFile, 'transformer')
        net = load_vnncomp_network('cgan', onnxFile);
        name = onnx(i).name;
        name = name(1:end-5);
        save(['networks2023/', name], "net", "-v7.3");
    end
end

%% vggnet16

disp("Running vggnet16...")

vgg_path = vnncomp_path + "vggnet16/onnx/";

vgg_instances = ["onnx/vgg16-7.onnx","vnnlib/spec0_screw.vnnlib"];

onnxFile = vgg_path + "vgg16-7.onnx";
net = load_vnncomp_network('vggnet16', onnxFile);
name = 'vgg16-7.onnx';
name = name(1:end-5);
save(['networks2023/', name], "net", "-v7.3");


%% ml4acopf

disp("Running ml4acopf..")

ml4_path = vnncomp_path + "ml4acopf/onnx/";

onnx = dir(ml4_path);

for i=3:length(onnx)
    onnxFile = ml4_path + onnx(i).name;
    net = load_vnncomp_network('ml4acopf', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name], "net", "-v7.3");
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
    save(['networks2023/', name], "net", "-v7.3");
end

%% Other 2023 benchmarks

%% test 

disp("Running test examples...");

test_path = vnncomp_path + "test/";

onnxmodels = ["test_nano.onnx", "test_sat.onnx", "test_small.onnx", "test_tiny.onnx", "test_unsat.onnx"];

for i = 1:length(onnxmodels)
    if contains(onnxmodels(i), "sat")
        onnxFile = test_path + onnxmodels(i);
        net =  importONNXNetwork(onnxFile, "InputDataFormats", "BCSS");
        name = char(onnxmodels(i));
        name = name(1:end-5);
        save(['networks2023/', name], "net", "-v7.3");
    end
end


%% yolo

disp("Running yolo...");

yolo_path = vnncomp_path + "yolo/onnx/";

onnx = dir(yolo_path);

for i = 3:length(onnx)
    onnxFile = yolo_path + onnx(i).name;
    net = load_vnncomp_network('yolo', onnxFile);
    name = onnx(i).name;
    name = name(1:end-5);
    save(['networks2023/', name], "net", "-v7.3");
end


%% Loading function
function net = load_vnncomp_network(category, onnx)
% load vnncomp 2023 benchmark NNs (subset support)

    % collins_rul: onnx to nnvnet
    % collins_nets = load_collins_NNs;
    if contains(category, 'collins_rul')
        net = importONNXNetwork(onnx);

    elseif contains(category, "nn4sys")
        % nn4sys: onnx to matlab:
        net = importONNXNetowrk(onnx, "OutputDataFormats", "BC"); % lindex
        
    elseif contains(category, "dist_shift")
        % dist_shift: onnx to matlab:
        net = importONNXNetwork(onnx, "InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape
        
    elseif contains(category, "cgan")
        % cgan
        net = importONNXNetwork(onnx,"InputDataFormats", "BC", 'OutputDataFormats',"BC"); %reshape
        
    elseif contains(category, "vgg")
        % vgg16: onnx to matlab
        net = importONNXNetwork(onnx); % flattenlayer
        
    elseif contains(category, "tllverify")
        % tllverify: onnx to matlab
        net = importONNXNetwork(onnx,"InputDataFormats", "BC", 'OutputDataFormats',"BC");
        
    elseif contains(category, "vit")
        % vit: onnx to matlab
        net = importONNXNetwork(onnx, "TargetNetwork","dlnetwork" );
        
    elseif contains(category, "cctsdb_yolo")
        % cctsdb_yolo: onnx to matnet
        net = importONNXNetwork(onnx, "TargetNetwork","dlnetwork" );
        
    elseif contains(category, "collins_yolo")
        % collins_yolo: onnx to matlab:
        net = importONNXNetwork(onnx, "TargetNetwork","dlnetwork" );

    elseif contains(category, "yolo")
        % yolo: onnx to matlab
        net = importONNXNetwork(onnx); % padlayer

    elseif contains(category, "acasxu")
        % acasxu: onnx to nnv
        net = importONNXNetwork(onnx, "InputDataFormats","BCSS");

    elseif contains(category, "ml4acopf")
        % ml4acopf: ?
        net = importONNXNetwork(onnx, "InputDataFormats", "BC");
        
    else % all other benchmarks
        % traffic: onnx to matlab: opset15 issues
        error("ONNX model not supported")
    end

end