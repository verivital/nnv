load ACASXU_run2a_1_1_batch_2000.mat;
Layers = [];
n = length(b);
for i=1:n - 1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = Layer(Wi, bi, 'ReLU');
    Layers = [Layers Li];
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = Layer(Wn, bn, 'Linear');

Layers = [Layers Ln];
F = FFNN(Layers);

% Input Constraints
% 55947.69 <= i1(\rho) <= 60760,
% -3.14 <= i2 (\theta) <= 3.14,
%-3.14 <= i3 (\shi) <= -3.14
% 1145 <= i4 (\v_own) <= 1200, 
% 0 <= i5 (\v_in) <= 60

lb = [55947.69; -3.14; -3.14; 1145; 0];
ub = [60760; 3.14; 3.14; 1200; 60];

% normalize input
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end

V = Reduction.getVertices(lb, ub);
% Input Set
I = Polyhedron('V', V');

desired_csv = 100; % conservativeness of 15%
k_max = 1; % maximum division level, k=1 -> divide each x[i] into 2 segments, k = 2 -> futher divide 2 segments into 4 ..
numOfCores = 2; % number of cores used in computation
n_samples = 5000; % use 5000 samples to estimate the output range


[R1, t1] = F.reach_approx_with_CSV_guarantee(I, desired_csv, k_max, numOfCores, n_samples); % exact scheme
[csv_vec, r, computed_range, est_range] = F.estimate_CSV(I, 5000);


% normalize estimated output range
norm_est_range = zeros(5,2);
for i=1:5   
    norm_est_range(i, 1) = est_range(i, 1) * range_for_scaling(6) + means_for_scaling(6);
    norm_est_range(i,2) =  est_range(i, 2) * range_for_scaling(6) + means_for_scaling(6);        
end

