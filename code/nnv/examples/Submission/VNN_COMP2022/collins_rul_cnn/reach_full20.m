%% Reachability analysis of collins benchmark (NN_rul_full_window_20.onnx)

%% 1) Load network
vnnFolder = "/home/manzand/Documents/MATLAB/vnncomp2022_benchmarks/benchmarks/";
rulFull20File = "collins_rul_cnn/onnx/NN_rul_full_window_20.onnx";
t = tic;
rulFull20 = onnx2nnv(vnnFolder+rulFull20File);
t = toc(t);

%% 2) Set reachability parameters
% This is the idea
% lb = zeros(20,20,1);
% ub = zeros(20,20,1);
% X = ImageStar(lb,ub);
% Y_rul_a = rulFull20.reach(X, reachOptions); % Seems to be working
% rulOptions = reachOptions;
% rulOptions.reachMethod = 'exact-star';
% Y_rul_e = rulFull20.reach(X, reachOptions); % Seems to be working

% Property robustness_2perturbations_delta5_epsilon10_w20.vnnlib
[lb_x, ub_x] = load_vnnlib(vnnFolder);
X = ImageStar(lb_x',ub_x');
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
rT = tic;
Y = rulFull20.reach(X, reachOptions); % Seems to be working
rT = toc(rT);
% Evaluate property
[y_lb,y_ub] = Y.getRanges;
disp(' ');
disp('RESULTS')
disp(' ')
if (y_lb >= 469.74957275390625 && y_ub <= 574.1383666992188)
    disp(' ');
    disp('Property is satisfied');
else
    disp(' ');
    disp('Verification failed');
end
disp("Reachability computation time = "+string(rT) + " seconds")
disp("Time to load network into NNV = "+string(t) + " seconds")

%% Helper functions

function [lb_input, ub_input] = load_vnnlib(vnnFolder)
    vnnlib_file = vnnFolder + 'collins_rul_cnn/vnnlib/robustness_2perturbations_delta5_epsilon10_w20.vnnlib';
    fileID = fopen(vnnlib_file,'r');
    tline = fgetl(fileID);
    % output_start = false;
    input_start = false;
    lb_input = zeros(20,20);
    ub_input = zeros(20,20);
    while ischar(tline)
        if isempty(tline)
            % Go to next line
            tline = fgetl(fileID);
            continue;
        end
        if contains(tline, 'Input constraints')
            input_start = true;
            tline = fgetl(fileID);
            i = 1;
            continue;
        elseif contains(tline, 'Output constraints')
            input_start = false;
            break
        end
        if input_start
            s = split(tline);
            value = split(s{4},')');
            value = str2double(value{1});
            if contains(s{2}, '>=') % lower bound assignment
                lb_input(i) = value;
            else
                ub_input(i) = value; % upper bound assignment
                i = i+1;
            end
        end
        % Go to next line
        tline = fgetl(fileID);
        disp(tline);
    end
    fclose(fileID);
end