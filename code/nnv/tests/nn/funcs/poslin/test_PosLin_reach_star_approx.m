
I0 = ExamplePoly.randVrep; 
I0.outerApprox;
V = [0 0; 1 1;1 0];
I = Star(V', I0.A, I0.b, I0.Internal.lb, I0.Internal.ub); % input star
X = I.sample(100);

figure;
I.plot;
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

t = tic;
S = PosLin.reach_star_approx(I); % over-approximate reach set
t1 = toc(t);
S1 = PosLin.reach(I); % exach reach set
t = tic;
S2 = PosLin.reach_star_approx2(I); % new-over-approximate method
t2 = toc(t);

Y = PosLin.evaluate(X);

figure;
S.plot;
hold on;
Star.plots(S1);
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs

figure;
Star.plots(S2);
