
lb = [-1; 1];
ub = [1; 2];


I = Box(lb, ub); 
I = I.toZono;

W = [0.5 1; -1 1];

I = I.affineMap(W, []);

I_star = I.toStar;
S = TanSig.reach_star_approx(I_star);
X = I_star.sample(10);
Y = TanSig.evaluate(X);

Z = TanSig.reach_zono_approx(I);

figure;
I.plot; % input set
hold on;
Z.plot; % zonotope reach set
hold on;
Star.plots(S, 'yellow'); % star reach set = zonotope reach set
hold on;
plot(Y(1, :), Y(2, :), '*'); % sampled outputs
hold on;
plot(X(1, :), X(2, :), 'ob'); % sampled inputs

