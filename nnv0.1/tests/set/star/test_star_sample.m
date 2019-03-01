I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
S = Star(V', I.A, I.b); % star set 1

V = S.sample(100);

figure;
S.plot;
hold on;
plot(V(1, :), V(2,:), '*');