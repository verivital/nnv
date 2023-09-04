I1 = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I1 = Star(V', I1.A, I1.b); % star set 1

W = [2 1; 1 -1];

I2 = I1.convexHull_with_linearTransform(W);

figure;
Star.plot(I2);
hold on;
Star.plot(I1);