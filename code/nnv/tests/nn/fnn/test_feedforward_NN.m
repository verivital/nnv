%% test 1:  falsify
W = [1 1; 0 1];
b = [0; 0.5];
L = LayerS(W, b, 'poslin');
Layers = {L};


F = NN(Layers);


lb = [-1; -1];
ub = [1; 1];

I = Star(lb, ub);

R = F.reach(I);

G = [-1 0];
g = [-1.5];

U = HalfSpace(G, g);

n_samples = 1000;

counter_inputs = F.falsify(I, U, n_samples);
counter_outputs = F.evaluate(counter_inputs);

figure;
subplot(1, 2, 1);
Star.plot(I);
hold on; 
plot(counter_inputs(1, :), counter_inputs(2, :), 'o');
title('Input Set and counter inputs');
subplot(1, 2, 2);
Star.plots(R);
hold on;
plot(counter_outputs(1, :), counter_outputs(2, :), 'o');
title('Output set and counter outputs');


%% test 2:  isRobust
W = [1 1; 0 1];
b = [0; 0.5];
L = LayerS(W, b, 'poslin');
Layers = {L};

F = NN(Layers);


lb = [0.5; 0.5];
ub = [1.5;1.5];
I = Star(lb,ub);

G = [-1 0];
g = [-1.5];

U = HalfSpace(G, g); % unsafe robust region

reachOptions = struct;
reachOptions.reachMethod = 'exact-star';

res = F.verify_robustness(I, reachOptions, U);
figure;
Star.plots(F.reachSet{end});
hold on;
U.plot;
title('Output set and unsafe region');



%% test 3: safety
W = [1 1; 0 1];
b = [0; 0.5];
L = LayerS(W, b, 'poslin');
Layers = {L};


F = NN(Layers);


lb = [-1; -1];
ub = [1; 1];

I = Star(lb, ub); % input set

reachOptions = struct;
reachOptions.reachMethod = 'exact-star';
R = F.reach(I, reachOptions);

G = [-1 0];
g = [-1.5];

U = HalfSpace(G, g); % unsafe region

n_samples = 100;

reachOptions = struct;
reachOptions.reachMethod = 'approx-star';
[safe, counter_inputs] = F.verify_safety(I, U, reachOptions, n_samples);
counter_outputs = F.evaluate(counter_inputs);

figure;
subplot(1, 2, 1);
Star.plot(I);
hold on;
plot(counter_inputs(1, :), counter_inputs(2, :), 'o');
title('Input Set and counter input set');
subplot(1, 2, 2);
Star.plots(R);
hold on;
plot(counter_outputs(1, :), counter_outputs(2, :), 'o');

title('Output set and unsafe region');
