I = ExamplePoly.randVrep;   % input set

W1 = [1.5 1; 0 0.5];  
b1 = [0.5; 0.5];
L1 = Layer(W1, b1, 'ReLU'); % first layer

W2 = [-1.5 1; 0.1 -1];  
b2 = [0.1; 0.2];
L2 = Layer(W2, b2, 'ReLU'); % second layer

Layers = [L1 L2]; % construct an array of layers

F = FFNN(Layers); % construct a feedforward neural network

% plot figure of input set and some sampled points of the input set
fig = figure;
I.plot;
Y = I.V';
hold on; 
plot(Y(1, :), Y(2, :), '*');
title('Input Set')

% plot figure of the output reach set, i.e., equal to R2, using FFNN.reach
% command, and some sampled points of the output from the input's sampled
% points to see if the reach set actually contained these points.

[R1, rn, t] = L1.reach_exact(I, 'single');
fig = figure;
R1.plot;
title('Layer 1 exact set')

%[R2, rn, t] = L2.reach_exact(R1, 'parallel');
%fig = figure;
%R2.plot;
%title('Layer 2 exact set')


[R11, t] = L1.reach_approx_polyhedron(I, 'single');
fig = figure;
R11.plot;
title('Layer 1 over-approximate reach set using polyhedron');

%[R22, t] = L2.reach_approx_polyhedron(R11, 'parallel');
%fig = figure;
%R22.plot;
%title('Layer 2 over-approximate reach set using polyhedron');

%[R1, rn, t] = L1.reach_exact(I, 'parallel');
%fig = figure;
%R1.plot;
%title('Layer 1 exact set')

%[R2, rn, t] = L2.reach_exact(R1, 'parallel');
%fig = figure;
%R2.plot;
%title('Layer 2 exact set')


%[R11, t] = L1.reach_approx_polyhedron(I, 'parallel');
%fig = figure;
%R11.plot;
%title('Layer 1 over-approximate reach set using polyhedron');

%[R22, t] = L2.reach_approx_polyhedron(R11, 'parallel');
%fig = figure;
%R22.plot;
%title('Layer 2 over-approximate reach set using polyhedron');


%[R, t] = L1.reach_approx_box(I, 'single');
%fig = figure;
%R.plot;
%title('Layer 1 box-set');

%[R, t] = L2.reach_approx_box(R, 'single');
%fig = figure;
%R.plot;
%title('Layer 2 box-set');

%R = F.reach(I, 'exact');
%Y3 = F.sample(I.V');
%Y32 = cell2mat(Y3(2));

%fig = figure; 
%plot(Y32(1, :), Y32(2, :), '*');
%hold on;
%R.plot;
%title('Exact output set');

%R1 = F.reach(I, 'approx-oriented-box');
%fig = figure; 
%R1.plot;
%hold on;
%plot(Y32(1, :), Y32(2, :), '*');
%title('Over-approximate output set using oriented box');

%R2 = F.reach(I, 'approx-box');
%fig = figure;
%R2.plot;
%hold on;
%plot(Y32(1, :), Y32(2, :), '*');
%title('Over-approximate output set using box');







