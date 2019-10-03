
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star
S = SatLin.stepReach(I, 1);

I = I.toPolyhedron;
P = SatLin.stepReachPolyhedronExact(I, 1);
figure;
P.plot;

figure;
Star.plots(S);