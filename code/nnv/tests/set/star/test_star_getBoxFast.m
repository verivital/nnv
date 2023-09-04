I1 = ExamplePoly.randVrep;   
V = [0 0; 1 1; 0.5 1];
I1 = Star(V', I1.A, I1.b);

B1 = I1.getBox();
B2 = I1.getBoxFast();

figure;
Box.plot(B2);
hold on;
Box.plot(B1);

