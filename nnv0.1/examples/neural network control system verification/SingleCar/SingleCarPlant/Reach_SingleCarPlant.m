load SingleCarPlant.mat;

Layers = [];
for i=1:6
    W = double(SingleCarPlant.W{1, i});
    b = double(SingleCarPlant.b{1, i}');
    Layers = [Layers Layer(W, b, 'ReLU')];
end

F = FFNN(Layers);

A1 = eye(7);
A2 = -eye(7);

b1 = SingleCarPlant.max;
b2 = -SingleCarPlant.min;

A1(4, :) = [];
b1(4) = [];
A2(4, :) = [];
b2(4) = [];


A = vertcat(A1, A2);
b = vertcat(b1, b2);
Ae = [0 0 0 1 0 0 0];

I1 = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', [0]); % corresponding to gear = 0

I2 = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', [1]); % corresponding to gear = 1

I3 = Polyhedron('A', A, 'b', b, 'Ae', Ae, 'be', [2]); % corresponding to gear = 2

I = [I1 I2 I3]; % input set


[R, t] = F.reach(I, 'exact', 4, []); % exact reach set
%[R, t] = F.reach(I, 'approx', 1, 500); % exact reach set
save F.mat F; % save the verified network
F.print('F.info'); % print all information to a file

fig = figure; 
R.plot;

