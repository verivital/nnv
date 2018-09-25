load SingleCarPlant.mat;

Layers = [];
for i=1:6
    W = double(SingleCarPlant.W{1, i}');
    b = double(SingleCarPlant.b{1, i}');
    Layers = [Layers Layer(W, b, 'ReLU')];
end

F = FFNN(Layers);

lb = SingleCarPlant.min';
ub = SingleCarPlant.max';

I = Polyhedron('lb', lb, 'ub', ub); % Input set

%[R, t] = F.reach(I, 'exact', 1, []); % exact reach set
[R, t] = F.reach(I, 'approx', 1, 500); % exact reach set
save F.mat F; % save the verified network
F.print('F.info'); % print all information to a file

fig = figure; 
R.plot;

