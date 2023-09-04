c1 = [0; 0];
V1 = [1 -1; 1 1; 0.5 0; -1 0.5];
Z1 = Zono(c1, V1');

B1 = Z1.getOrientedBox();

figure;
Box.plot(B1);
hold on;
Zono.plot(Z1);