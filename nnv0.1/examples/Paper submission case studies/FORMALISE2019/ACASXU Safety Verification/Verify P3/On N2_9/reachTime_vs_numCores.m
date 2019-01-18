% This script runs reachability analysis for N2_9 with different number of
% cores, save the reachTime and then plot the figue of time reduction.

load ACASXU_run2a_2_9_batch_2000.mat;
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
% 980 <= i4 (\v_own) <= 1200, 
% 960 <= i5 (\v_in) <= 1200

lb = [1500; -0.06; 3.1; 980; 960];
ub = [1800; 0.06; 3.14; 1200; 1200];

% normalize input
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end

I = Polyhedron('lb', lb, 'ub', ub);

numCores = [4; 12; 16; 24; 28; 32];
reachTime = zeros(length(numCores), 1);
for i=1:length(numCores)
    [R, t] = F.reach(I, 'exact', numCores(i), []); % exact reach set
    reachTime(i) = t;
end

save numCores.mat numCores;
save reachTime.mat reachTime;

load reachTime.mat;
load numCores.mat;

fig = figure;
plot(numCores, reachTime, '--*' );
xlabel('N', 'Fontsize', 16);
ylabel('RT', 'Fontsize', 16);
set(gca, 'Fontsize', 16);
saveas(gcf,'time_reduction.pdf');

