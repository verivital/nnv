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

% plot figure of the reach set of the first layer output and the sampled
% points of the output corresponding to the input's sampled points
R1 = L1.reach(I, 'exact');
Y1 = L1.sample(I.V');
fig = figure;
R1.plot;
hold on;
plot(Y1(1, :), Y1(2, :), '*');

% plot figure of the reach set of the second layer output, i.e., the output
% layer, and the sampled points of the output corresponding to the input's
% sampled points
R2 = L2.reach(R1, 'exact');
Y2 = L2.sample(Y1);
fig = figure; 
plot(Y2(1, :), Y2(2, :), '*');
hold on;
R2.plot;


% plot figure of the output reach set, i.e., equal to R2, using FFNN.reach
% command, and some sampled points of the output from the input's sampled
% points to see if the reach set actually contained these points.
R = F.reach(I, 'exact');
Y3 = F.sample(I.V');
Y32 = cell2mat(Y3(2));

fig = figure; 
plot(Y32(1, :), Y32(2, :), '*');
hold on;
R.plot;



