% Verification of mnist network

% Load network
load images.mat;
load('Large_ConvNet.mat');
nnvNet = matlab2nnv(net);

% Note: label = 1 --> digit 0
%       label = 2 --> digit 1
%       ...
%       label = 10 --> digit 9

% Adversarial attack variables
delta = [0.005 0.01 0.015];
d = [250 245 240]; % threshold for brightening attack

% Data to consider
M = length(delta);
P = length(d);
% N = 100; % number of test images used for testing robustness (original)
N = 10; % much faster 

%% construct input sets with different values of d and delta

inputSetStar = cell(P, M);
correct_labels = cell(P, M);
n_images = length(IM_labels);

for i=1:P
    for j=1:M
        inputStar = [];
        count = 0;
        labels = zeros(1, N);
        for k=1:n_images
            % Create adversarial attack
            IM = IM_data(:,:, k);
            lb = IM;
            ub = IM;
            for l=1:numel(IM)
                if IM(l) >= d(i)
                    lb(l) = 0;
                    ub(l) = delta(j)*IM(l);
                end
            end
            % Create input set
            lb = reshape(lb, [28 28 1]);
            ub = reshape(ub, [28 28 1]);
            S = ImageStar(lb, ub);
            if ~isempty(S.C) % ensure it is a set
                inputStar = [inputStar S];
                count = count + 1;
                labels(count) = IM_labels(k);
            end
            % When we get to N, stop 
            if count == N
                break;
            end
        end        
        inputSetStar{i, j} = inputStar;
        correct_labels{i, j} = labels;
    end
end



%% evaluate robustness

% Initialize result variables
VT_star = zeros(P, M); % verification time of the approx-star method
VT_absdom = zeros(P, M); % verification time of the DeepPoly abstract domain method
r_star = zeros(P, M); % robustness percentage on an array of N tested input sets obtained by the approx-star method
r_absdom = zeros(P, M); % robustness percentage on an array of N tested input sets obtained by the DeepPoly abstract domain method

% Reachability analysis using approx-star
reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
disp('Computing robustness using approx-star...');
for i=1:P
    for j=1:M
        t = tic;
        tot_combo = 0;
        for k=1:N
            res = nnvNet.verify_robustness(inputSetStar{i, j}(k), reachOptions, correct_labels{i, j}(k));
            tot_combo = tot_combo + res == 1;% add only those that are verified to be robust
        end
        r_star(i, j) = tot_combo/N; % save as percentage
        VT_star(i, j) = toc(t);
    end
end

% Reachability analysis using abstrac-domain method
reachOptions = struct;
reachOptions.reachMethod = 'abs-dom';
for i=1:P
    for j=1:M
        % absdom run into memory issue (set as 1 hour for delta = 0.01 and 0.015)
        if  ~((j== 3) || (i==3 && j==2))
            t = tic;
            tot_combo = 0;
            for k=1:N
                res = nnvNet.verify_robustness(inputSetStar{i, j}(k), reachOptions, correct_labels{i, j}(k));
                tot_combo = tot_combo + res == 1;% add only those that are verified to be robust
            end
            r_absdom(i, j) = tot_combo/N; % save as percentage
            VT_absdom(i, j) = toc(t);
        end
    end
end
       
% save Large_ConvNet_Results.mat r_star VT_star r_absdom VT_absdom;


%% Results - print to screen

fprintf('\n========================================================================================');
fprintf('\n          ROBUSTNESS VERIFICATION RESULTS (IN PERCENT) OF LARGE_CONVNET                 ');
fprintf('\n========================================================================================\n\n');

for j=1:M
    fprintf("             delta = %.5f", delta(j));
end
fprintf("\n");
for j=1:M
    fprintf("         Polytope   ImageStar");
end
fprintf("\n");
for i=1:P
    fprintf("d = %d", d(i));
    for j=1:M
        fprintf("    %.2f        %.2f      ", N*r_absdom(i, j), N*r_star(i, j));
    end
    fprintf("\n");
end

fprintf('\n========================================================================================');
fprintf('\n                VERIFICATION TIMES (IN SECONDS) OF LARGE_CONVNET                        ');
fprintf('\n========================================================================================\n\n');

for j=1:M
    fprintf("             delta = %.5f", delta(j));
end
fprintf("\n");
for j=1:M
    fprintf("         Polytope   ImageStar");
end
fprintf("\n");
for i=1:P
    fprintf("d = %d", d(i));
    for j=1:M
        fprintf("    %.2f       %.2f        ", VT_absdom(i, j), VT_star(i, j));
    end
    fprintf("\n");
end