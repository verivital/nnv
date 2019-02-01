load ACASXU_run2a_1_1_batch_2000.mat;

% construct the network

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

% normalized inputs 
% for example:
% (12000 - 19791)/60261  <= i1 <= (62000 - 19791)/60261

for i=1:5
    lb1(i) = (lb1(i) - means_for_scaling(i))/range_for_scaling(i);
    ub1(i) = (ub1(i) - means_for_scaling(i))/range_for_scaling(i);   
end

load ACASXU_run2a_2_1_batch_2000.mat;
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

V1 = Reduction.getVertices(lb1, ub1);

I1 = Polyhedron('V', V1');

lb2 = [12000; -3.141592; -3.141592; 100; 0.0];
ub2 = [62000; -0.7; -3.141592 + 0.005; 1200; 1200];
% normalized inputs 
for i=1:5
    lb2(i) = (lb2(i) - means_for_scaling(i))/range_for_scaling(i);
    ub2(i) = (ub2(i) - means_for_scaling(i))/range_for_scaling(i);   
end

load ACASXU_run2a_2_1_batch_2000.mat;
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

V2 = Reduction.getVertices(lb2, ub2);

I2 = Polyhedron('V', V2');

I = [I1 I2];

[R, t] = F.reach(I, 'exact', 4, []); % exact reach set
save F.mat F; % save the verified network
F.print('F.info'); % print all information to a file


