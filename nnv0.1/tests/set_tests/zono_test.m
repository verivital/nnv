c1 = [0; 0];
V1 = [1 -1; 1 1];
Z1 = Zono(c1, V1);

B1 = Z1.getBox();

c2 = [1; 1];
V2 = [2 1; -1 1];
Z2 = Zono(c2, V2);
fig = figure;
Z2.plot;
hold on;
Z1.plot;

W = [3 1; 1 0; 2 1];
b = [0.5; 1; 0];

Z3 = Z1.affineMap(W, b);
fig = figure;
Z3.plot;

B3 = Z3.getBox();
fig = figure;
B3.plot;

Z4 = Z1.MinkowskiSum(Z2);
fig = figure;
Z4.plot;