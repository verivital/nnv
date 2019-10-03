
c1 = [0; 0];
V1 = [1 -1; 1 1; 0.5 0; -1 0.5];
Z1 = Zono(c1, V1');
I1 = Z1.toStar();

I2 = I1.getOrientedBox();
I3 = I1.getBox();

figure;
I2.plot;
hold on;
I1.plot;

figure;
I3.plot;
hold on;
I1.plot;

