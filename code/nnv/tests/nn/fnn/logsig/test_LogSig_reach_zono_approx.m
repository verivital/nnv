
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

B = I.getBox; 
I = B.toZono;


W = [0.5 1; -1 1];

I = I.affineMap(W, []);

I_star = I.toStar;
S1 = LogSig.reach_star_approx_split(I_star);
S2 = LogSig.reach_star_approx_no_split(I_star);

X = I_star.sample(10);
Y = LogSig.evaluate(X);

Z = LogSig.reach_zono_approx(I);

figure;
I.plot; % input set
hold on;
Z.plot; % zonotope reach set
hold on;
S2.plot; % star reach set = zonotope reach set
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

