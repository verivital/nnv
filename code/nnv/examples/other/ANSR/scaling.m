% Another one of the affine Mapping that we can do is scaling

% scaling + rotation + translation

rng(3);

% Create random set
I = ExamplePoly.randVrep;
I.outerApprox;
V = [0 0; 1 0; 0 1];
I = Star(V', I.A, I.b, I.Internal.lb, I.Internal.ub); % input star

%% Can we make if twice biiger/smaller?

% x2 bigger
R2 = I;
R2.V(:,2:end) = R2.V(:,2:end) * 2;

% x0.50 (half size)
R = I;
R.V(:,2:end) = R.V(:,2:end) * 0.5;

% Visualize results
f = figure;
Star.plot(R2,'b');
hold on;
Star.plot(I, 'r');
hold on;
Star.plot(R, 'c');
saveas(f, "scale_randSet.png");


%% Let's get a square to demonstrate this as well

% center at [0, 0]
% side lengths = 1
V = [0 0.5 0; 0 0 0.5];
C = [0,0];
d = 0;
pred_lb = [-1;-1];
pred_ub = [1;1];

X = Star(V,C,d,pred_lb, pred_ub);

% Do the scaling now

% x2 bigger
X1 = X;
X1.V(:,2:end) = X1.V(:,2:end) * 2;

% x0.50 (half size)
X2 = X;
X2.V(:,2:end) = X2.V(:,2:end) * 0.5;

% Visualize results
f = figure;
Star.plot(X1,'b');
hold on;
Star.plot(X, 'r');
hold on;
Star.plot(X2, 'c');
saveas(f, "scale_square.png");
