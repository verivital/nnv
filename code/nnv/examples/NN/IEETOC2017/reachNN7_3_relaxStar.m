%% Example of reachability and verification for a randomly generated neural network
% Smaller set, faster example

% Load the network (parameters)
load NeuralNetwork7_3.mat;

% Create NNV model
n = length(b);
Layers = cell(n,1);
% hidden layers are all ReLUs
for i=1:n - 1
    bi = cell2mat(b(i));
    Wi = cell2mat(W(i));
    Li = LayerS(Wi, bi, 'poslin');
    Layers{i} = Li;
end
% output layer is linear
bn = cell2mat(b(n));
Wn = cell2mat(W(n));
Ln = LayerS(Wn, bn, 'purelin');
Layers{end} = Ln;
% Create NN model
F = NN(Layers);

% Define input set for the network
lb = [0; 0 ; 0];
ub = [1; 1; 1];
I = Star(lb, ub);

% select option for reachability algorithm
reachOptions = struct; % initialize
reachOptions.reachMethod = 'relax-star-range'; % reachability method
reachOptions.relaxFactor = 0.0;

% Comute reachability (approx)
tr = tic;
Rr = F.reach(I, reachOptions); % exact reach set using stars
tr = toc(tr);

% generate some input to test the output (evaluate network)
e = 0.2;
x = [];
y = [];
for x1 = lb(1):e:ub(1)
    for x2 = lb(2):e:ub(2)
        for x3=lb(3):e:ub(3)
            xi = [x1; x2; x3];
            yi = F.evaluate(xi);
            x = [x, xi];
            y = [y, yi];
        end
    end
end


%% Visualize results

% Plot exact and approx sets 
fig = figure;
Star.plots(Rr,'r');
hold on;

% Plot some of the evaluated inputs
plot(y(1, :), y(2, :), '.', 'Color', 'k');

% Evaluate upper, lower bounds
y1 = F.evaluate(lb);
y2 = F.evaluate(ub);
y3 = F.evaluate((lb+ub)/2);

% Plot upper and lower bound results
plot(y1(1,:), y1(2,:), 'x', 'Color', 'r');
hold on;
plot(y2(1,:), y2(2,:), 'x', 'Color', 'r');
hold on;
plot(y3(1,:), y3(2,:), 'x', 'Color', 'r');
