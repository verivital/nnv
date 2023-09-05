c = [1; 1];
V = rand(2, 3);
Z1 = Zono(c, V);

W = [2 1; 0 -1];
b = [];

Z2 = Z1.affineMap(W, b);

Z12 = Z1.convexHull(Z2);


figure;
Zono.plot(Z12);
hold on;
Zono.plot(Z1);
hold on;
Zono.plot(Z2);

Z3 = Z1.convexHull_with_linearTransform(W);
figure;
Zono.plot(Z3);
hold on;
Zono.plot(Z1);
hold on;
Zono.plot(Z2);



