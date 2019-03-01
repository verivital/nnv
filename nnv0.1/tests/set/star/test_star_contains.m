I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
S = Star(V', I.A, I.b); % star set 1

s1 = [-0.5; 1];
s2 = [1.5; 0.5];

b1 = S.contains(s1);
b2 = S.contains(s2);

S.plot;