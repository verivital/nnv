
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

figure;
I.plot;
S = PosLin.reach_star_approx(I); % over-approximate reach set
S1 = PosLin.reach(I); % exach reach set
figure;
S.plot;
hold on;
Star.plots(S1);