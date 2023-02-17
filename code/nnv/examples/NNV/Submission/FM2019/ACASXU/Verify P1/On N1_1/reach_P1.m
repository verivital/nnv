load ACASXU_run2a_1_1_batch_2000.mat;
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

I = Star(lb, ub);

numCores = 90;

[R1, ~] = F.reach(I, 'exact-star', numCores); % exact reach set using polyhdedron
F.print('F_exact_star.info'); % print all information to a file

[R2, ~] = F.reach(I, 'approx-star'); % approximate reach set using star
F.print('F_approx_star.info'); % print all information to a file

[R3, ~] = F.reach(I.getZono, 'approx-zono'); % approximate reach set using zonotope
F.print('F_approx_zono.info'); % print all information to a file

[R4, ~] = F.reach(I, 'abs-dom'); % approximate reach set using abstract domain
F.print('F_abs_dom.info'); % print all information to a file


% NOTE ********
% we just found out that I made a mistake of copying a wrong the outputSet.mat file for this property. 
% the Zonotope, abstract-domain and approx-star cannot prove the safety
% property for this network which is wrongly stated in the FM2019 paper. We are sorry about this mistake
% Dung Tran: 10/11/2019

save outputSet.mat R1 R2 R3 R4;

