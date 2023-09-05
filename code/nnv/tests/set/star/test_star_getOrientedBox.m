
c1 = [0; 0];
V1 = [1 -1; 1 1; 0.5 0; -1 0.5];
Z1 = Zono(c1, V1');
I1 = Z1.toStar();

I2 = I1.getOrientedBox();
I3 = I1.getBox();

figure;
Box.plot(I2);
hold on;
Star.plot(I1);

figure;
Box.plot(I3);
hold on;
Star.plot(I1);

