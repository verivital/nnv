
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

figure;
I.plot;
S = SatLins.stepReachAbstractDomain(I, 1);
S1 = SatLins.stepReach(I,1);
figure;
S.plot; % over-approximate reach set
hold on;
Star.plots(S1); % exact reach set