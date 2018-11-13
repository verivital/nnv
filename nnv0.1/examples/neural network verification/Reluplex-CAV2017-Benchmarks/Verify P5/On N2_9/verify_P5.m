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
% 250 <= i1(\rho) <= 400,
% 0.2 <= i2 (\theta) <= 0.4,
%-3.141592 <= i3 (\shi) <= -3.141592 + 0.005
% 100 <= i4 (\v_own) <= 400, 
% 0 <= i5 (\v_in) <= 400

lb = [250; 0.2; -3.141592; 100; 0.0];
ub = [400; 0.4; -3.141592 + 0.005; 400; 400];

% normalize input
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end

V = Reduction.getVertices(lb, ub);

I = Polyhedron('V', V');

[R, t] = F.reach(I, 'exact', 4, []); % exact reach set
save result.mat; % save the verified network
F.print('F.info'); % print all information to a file
