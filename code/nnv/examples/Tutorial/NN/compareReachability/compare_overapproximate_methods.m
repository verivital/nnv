%% Example of reachability for a randomly generated neural network

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

%% Reachability 1

% select option for reachability algorithm
reachOptions = struct; % initialize
reachOptions.reachMethod = 'approx-star';

% Comute reachability 
t1 = tic;
R1 = F.reach(I, reachOptions); % exact reach set using stars
t1 = toc(t1);

%% Reachability 2

% select option for reachability algorithm
reachOptions = struct; % initialize
reachOptions.reachMethod = 'relax-star-range'; % other option (-area, -random, -bound)
reachOptions.relaxFactor = 0.1; % must declare a relaxFactor for any relax method
% relaxFactor must be between 0 and 1 (.1 -> do 90% of LP, skip the other 90%)

% Comute reachability 
t2 = tic;
R2 = F.reach(I, reachOptions); % exact reach set using stars
t2 = toc(t2);

%% Reachability 3

% select option for reachability algorithm
reachOptions = struct; % initialize
reachOptions.reachMethod = 'relax-star-range'; % other option (-area, -random, -bound)
reachOptions.relaxFactor = 0.25; % must declare a relaxFactor for any relax method
% relaxFactor must be between 0 and 1 (.25 -> do 75% of LP, skip the other 25%)

% Comute reachability 
t3 = tic;
R3 = F.reach(I, reachOptions); % exact reach set using stars
t3 = toc(t3);

%% Simulations

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

% Plot computed sets 
fig = figure;
Star.plot(R3,'b');
hold on;
Star.plot(R2,'r');
hold on;
Star.plot(R1,'c');
hold on;

% Plot some of the evaluated inputs
plot(y(1, :), y(2, :), '.', 'Color', 'k');

% Evaluate upper, lower bounds
y1 = F.evaluate(lb);
y2 = F.evaluate(ub);
y3 = F.evaluate((lb+ub)/2);

% Plot upper and lower bound results
plot(y1(1,:), y1(2,:), 'x', 'Color', 'k');
hold on;
plot(y2(1,:), y2(2,:), 'x', 'Color', 'k');
hold on;
plot(y3(1,:), y3(2,:), 'x', 'Color', 'k');

