% Reachability of HardSig
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

S = HardSig.reach_star_exact(I, 'single');