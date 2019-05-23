load ACASXU_run2a_4_9_batch_2000.mat;
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
% 1500 <= i1(\rho) <= 1800,
% -0.06 <= i2 (\theta) <= 0.06,
% 3.1 <= i3 (\shi) <= 3.14
% 1000 <= i4 (\v_own) <= 1200, 
% 700 <= i5 (\v_in) <= 800

lb = [1500; -0.06; 3.1; 1000; 700];
ub = [1800; 0.06; 3.14; 1200; 800];

% normalize input
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end

B = Box(lb, ub);

I = B.toStar;

numCores = [1, 10, 20, 30, 40]; % number of cores used in computation
N = length(numCores);
reachTime = zeros(N, 1);

for i=1:N
    [~, t] = F.reach(I, 'exact', numCores(i), []); % exact reach set
    reachTime(i) = t; 
end

save reachTime.mat reachTime;
save numCores.mat numCores;

