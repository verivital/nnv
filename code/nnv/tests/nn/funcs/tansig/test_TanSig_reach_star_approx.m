
I = ExamplePoly.randVrep;
I.outerApprox;
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub); % input star

S = TanSig.reach_star_approx(I);
relaxFactor = 0.5;
S2 = TanSig.reach_star_approx(I, 'approx-star', relaxFactor);
X = I.sample(10);
Y = TanSig.evaluate(X);



figure;
I.plot; % input set
hold on;
S.plot; % reach set
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

figure;
S2.plot;
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs