
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

figure;
I.plot;
S = SatLin.stepReachStarApprox(I, 1);
figure;
S.plot;