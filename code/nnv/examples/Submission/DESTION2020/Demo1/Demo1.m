%% Using NNV to verify a FNN (AcasXu)
% Verify one of the AcasXu NNs with several input sets and properties

%% Step 1. Load the neural network
%  Load info
netInfo = load('acasxu_1_1.mat');

% Load network NNV format
net1 = Load_nn('acasxu_1_1.mat');

% Load network using nnvmw from Reluplex
net2 = Load_nn('AcasXu_1_1.nnet');

%% Step 2. Create input set
% Property 2 -> COC (index 1) will not be the maximal score
lb = [55947.69; -3.14; -3.14; 1145; 0];
ub = [60760; 3.14; 3.14; 1200; 60];
mean_scale = [19791.0910000000;0;0;650;600];
range_scale = [60261;6.28318530718000;6.28318530718000;1100;1200];
lb = (lb-mean_scale)./range_scale;
ub = (ub-mean_scale)./range_scale;
inp1 = Star(lb,ub);

%% Step 3. Reachability analysis (Over-approximate)
tic;
ApproxSet = net2.reach(inp1,'approx-star');
toc;

%% Step 4. Evaluate results
[m, M] = ApproxSet.getRanges;
disp(' ');

if all(m) < M(1)
    disp('Property 2 is satisfied. COC (index 1) is not the maximal score');
    p2sat = true;
else
    disp('Property 2 is unknown. COC (index 1) may be the maximal score');
    p2sat = false;
end

%% Repeat process (exact-star)

%% Step 2. Create input set
% Property 2 -> COC (index 1) will not be the maximal score
lb = [250; 0.2; -3.14; 100; 0];
ub = [400; 0.4; 3.14; 400; 400];
mean_scale = [19791.0910000000;0;0;650;600];
range_scale = [60261;6.28318530718000;6.28318530718000;1100;1200];
lb = (lb-mean_scale)./range_scale;
ub = (ub-mean_scale)./range_scale;
inp2 = Star(lb,ub);

%% Step 3. Reachability analysis
tic;
ExactSet = net2.reach(inp2,'exact-star');
toc;

%% Step 4. Evaluate results
[m, M] = ExactSet.getRanges;
disp('');

if all(m) < M(1)
    disp('Property 5 is satisfied. COC (index 1) is not the maximal score');
    p2sat = true;
else
    disp('Property 5 is unknown. COC (index 1) may be the maximal score');
    p2sat = false;
end