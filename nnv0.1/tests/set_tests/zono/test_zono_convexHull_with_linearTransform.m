c = [1; 1];
V = rand(2, 3);
Z1 = Zono(c, V);

W = [2 1; 0 -1];
b = [];

Z2 = Z1.affineMap(W, b);

Z12 = Z1.convexHull(Z2);


figure;
Z12.plot;
hold on;
Z1.plot;
hold on;
Z2.plot;

Z3 = Z1.convexHull_with_linearTransform(W);
figure;
Z3.plot;
hold on;
Z1.plot;
hold on;
Z2.plot;



