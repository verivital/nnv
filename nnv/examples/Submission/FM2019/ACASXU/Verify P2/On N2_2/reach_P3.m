load ACASXU_run2a_2_2_batch_2000.mat;
Layers = [];
n = length(b);
for i=1:n - 1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers = [Layers Li];
end
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');

Layers = [Layers Ln];
F = FFNNS(Layers);

% Input Constraints
% 1500 <= i1(\rho) <= 1800,
% -0.06 <= i2 (\theta) <= 0.06,
% 3.1 <= i3 (\shi) <= 3.14
% 980 <= i4 (\v_own) <= 1200, 
% 960 <= i5 (\v_in) <= 1000

lb = [1500; -0.06; 3.1; 980; 960];
ub = [1800; 0.06; 3.14; 1200; 1000];

% normalize input
for i=1:5
    lb(i) = (lb(i) - means_for_scaling(i))/range_for_scaling(i);
    ub(i) = (ub(i) - means_for_scaling(i))/range_for_scaling(i);   
end

I = Star(lb, ub);

numCores = 6;

[R1, ~] = F.reach(I, 'exact-star', numCores); % exact reach set using polyhdedron
F.print('F_exact_star.info'); % print all information to a file

[R2, ~] = F.reach(I, 'approx-star'); % approximate reach set using star
F.print('F_approx_star.info'); % print all information to a file

[R3, ~] = F.reach(I.getZono, 'approx-zono'); % approximate reach set using zonotope
F.print('F_approx_zono.info'); % print all information to a file

[R4, ~] = F.reach(I, 'abs-dom'); % approximate reach set using abstract domain
F.print('F_abs_dom.info'); % print all information to a file

save outputSet.mat R1 R2 R3 R4;

