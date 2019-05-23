load neural_network_information_10.mat
n = length(W);
L = [];
for i=1:n
    L1 = Layer(W{1, i}, b{1, i}, 'ReLU');
    L = [L L1];
end

F = FFNN(L); % neural network

lb = [-0.1; -0.1; -0.1];
ub = [0.1; 0.1; 0.1];

I = Box(lb, ub); % input set = Sherlock
t = tic;
[R1, t1] = F.reach(I.toStar, 'exact', 4, []);
exact_range = Star.get_hypercube_hull(R1);
reach_time = toc(t);
save exact_range exact_range;
save reach_time reach_time;
save outputSet.mat R1; % save the verified network
F.print('F.info'); % print all information to a file
