c1 = [0; 0];
V1 = [1 -1; 1 1];
Z1 = Zono(c1, V1);
r1 = Z1.getSupInfinityNorm();

c2 = [1; 1];
V2 = [2 1; -1 1];
Z2 = Zono(c2, V2);
r2 = Z2.getSupInfinityNorm();

figure;
Z1.plot;
hold on;
Z2.plot;
