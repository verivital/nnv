load NeuralNetwork7_3.mat;
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
C = [1 0 0; -1 0 0; 0 1 0; 0 -1 0; 0 0 1; 0 0 -1];
d = [1; 1; 1; 1; 1; 1];

V = [0 0 0; 1 0 0; 0 1 0; 0 0 1];
pred_lb = [-1;-1;-1];
pred_ub = [1;1;1];
I = Star(V', C, d, pred_lb, pred_ub); % input set as a Star set

% select option for reachability algorithm

%[R, t] = F.reach(I, 'exact-star', 1); % exact reach set using stars


F.flatten('exact-star');