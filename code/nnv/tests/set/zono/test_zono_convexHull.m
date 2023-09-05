c1 = [0; 0];
V1 = [1 0 -1; 1 1 1];
Z1 = Zono(c1, V1);


c2 = [1; 1];
V2 = [2 1 0; -1 1 0];
Z2 = Zono(c2, V2);

Z3 = Z2.convexHull(Z1);

figure;
Zono.plot(Z3);
hold on;
Zono.plot(Z2);
hold on; 
Zono.plot(Z1);