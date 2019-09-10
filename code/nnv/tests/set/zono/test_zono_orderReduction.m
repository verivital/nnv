c1 = [0; 0];
V1 = [1 -1; 1 1; 0.5 1; -1.2 1];
Z1 = Zono(c1, V1');

Z2 = Z1.orderReduction_box(3);

figure;
Z2.plot;
hold on;
Z1.plot;

