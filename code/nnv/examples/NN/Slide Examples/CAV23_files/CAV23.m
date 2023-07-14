%% Example 1 - CAV'23 - UNSAT
% RL Benchmarks (VNN-COMP'22)

disp("Running examples 1...")

t = tic;

% Load Options (all networks follow a similar architecture)
loadOptions.InputDataFormat = "BC";

% input files
onnxFile = "CAV23_files/cartpole.onnx";
vnnlibFile = "CAV23_files/cartpole_54.vnnlib";

% Load network
nn = onnx2nnv(onnxFile, loadOptions);

% Verify property
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
result1 = nn.verify_vnnlib(vnnlibFile, reachOptions);

toc(t);

disp("Verification result is UNSAT, property holds");
disp(" ");


%% Example 2 - CAV'23 - UNKNOWN
% RL Benchmarks (VNN-COMP'22)

disp("Running example 2...")

t = tic;

% input files
onnxFile = "CAV23_files/lunarlander.onnx";
vnnlibFile = "CAV23_files/lunarlander_40.vnnlib";

% Load network
nn = onnx2nnv(onnxFile, loadOptions);

% Verify property
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
result2 = nn.verify_vnnlib(vnnlibFile, reachOptions);

toc(t);

disp("Verification result is UNKNOWN");
disp(" ");


%% Example 3 - CAV'23 - SAT
% RL Benchmarks (VNN-COMP'22)

disp("Running example 3...")

t = tic;

% input files
onnxFile = "CAV23_files/lunarlander.onnx";
vnnlibFile = "CAV23_files/lunarlander_43.vnnlib";

% Load network
nn = onnx2nnv(onnxFile, loadOptions);

% Verify property
reachOptions = struct;
reachOptions.reachMethod = 'exact-star';
result3 = nn.verify_vnnlib(vnnlibFile, reachOptions);

toc(t);

disp("Verification result is SAT, property is violated");
