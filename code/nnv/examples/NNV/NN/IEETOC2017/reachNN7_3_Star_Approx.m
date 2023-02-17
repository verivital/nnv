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

lb = [-1; -1; -1];
ub = -lb;

I = Star(lb, ub);
%[R1, t1] = F.reach(I, 'exact-star',4); % exact reachable set
[R2, t2] = F.reach(I, 'approx-star'); % over-approximate reach set using stars




