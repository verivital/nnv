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

% Input Constraints (see paper Reluplex: An Efficient SMT Solver for Verifying Deep Neural Networks, CAV 2017)
% 12000 <= i1 (\rho) <= 62000,  
% 0.7 <= i2 (\theta) <= 3.141592 or (-3.141592 <= i2 (\theta) <= -0.7)  
% -3.141592 <= i3 (\shi) <= -3.141592 + 0.005
% 100 <= i4 (\v_own) <= 1200,
% 0 <= i5 (\v_in) <= 1200

lb1 = [12000; 0.7; -3.141592; 100; 0.0];
ub1 = [62000; 3.141592; -3.141592 + 0.005; 1200; 1200];
I1 = Polyhedron('lb', lb1, 'ub', ub1);

lb2 = [12000; -3.141592; -3.141592; 100; 0.0];
ub2 = [62000; -0.7; -3.141592 + 0.005; 1200; 1200];
I2 = Polyhedron('lb', lb2, 'ub', ub2);

I = [I1 I2];

[R, t] = F.reach(I, 'exact', 4, []); % exact reach set
save F.mat F; % save the verified network
F.print('F.info'); % print all information to a file


