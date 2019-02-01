center = [1; 1; 1];
V = [1 0 1; 0 1 1; 1 0 0];
P = ExamplePoly.randHrep('d',3);
S = Star([center V], P.A, P.b);

figure;
Star.plotBoxes_3D(S, 1, 2, 3, 'green');