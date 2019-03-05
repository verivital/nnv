
I = ExamplePoly.randVrep;   
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b); % input star

S = LogSig.reach_star_approx(I);
X = I.sample(10);
Y = LogSig.evaluate(X);



figure;
I.plot; % input set
hold on;
S.plot; % reach set
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

