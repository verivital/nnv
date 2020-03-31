
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

B = I.getBox; 
I = B.toZono;

W = [0.5 1; -1 1];

I = I.affineMap(W, []);

I_star = I.toStar;
S1 = TanSig.reach(I_star, 'approx-star');
S2 = TanSig.reach(I_star, 'abs-dom');

X = I_star.sample(10);
Y = TanSig.evaluate(X);

Z = TanSig.reach(I, 'approx-zono');

figure;
I.plot; % input set
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

figure;
Z.plot; % zonotope reach set

figure;
Star.plots(S2, 'blue'); % abs-dom
hold on;
Star.plots(S1, 'red'); % approx-star 
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs

