I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

S = LogSig.reach_star(I);

figure;
I.plot;
hold on;
S.plot;