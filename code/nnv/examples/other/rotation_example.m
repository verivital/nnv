%% Show rotation examples

% Let's define a 2-dimensional shape (rectangle)

% center at [0.5, 0.5]
% side lengths = 1
V = [0.5 0.5 0; 0.5 0 0.5];
C = [0,0];
d = 0;
pred_lb = [-1;-1];
pred_ub = [1;1];

X = Star(V,C,d,pred_lb, pred_ub);

figure;
Star.plots(X,'b');

% Now, if we want to rotate this set 45 degrees, what would it look like?

% the center needs to be the same
V = [0, sqrt(2)/4, -sqrt(2)/4;
sqrt(2)/2, sqrt(2)/4, sqrt(2)/4];
C = [0 0];
d = [0];
pred_lb = [-1;-1];
pred_ub = [ 1; 1];

X2 = Star(V,C,d,pred_lb, pred_ub);

figure;
Star.plots(X2,'b');
% xlim([-3 3])
% ylim([-3 3])


%% Let's try the rotation here as well

% rotate the rectangle by an angle
angle = 45;
angle = deg2rad(angle);

rot_mat = [cos(angle), -sin(angle);
    sin(angle), cos(angle)];

X_r = X.affineMap(rot_mat, []);

figure;
Star.plots(X_r,'b');


%% Let's add constraints to the C and d variables rather than the predicate lower and uper bounds

% the center needs to be the same
V = [0, sqrt(2)/4, -sqrt(2)/4;
sqrt(2)/2, sqrt(2)/4, sqrt(2)/4];
C = [1 0;
    0 1;
    -1 0;
    0 -1];
d = [1; 1; 1; 1];

X_cd = Star(V,C,d);

figure;
Star.plot(X_cd, 'r');

% Let's rotate the same way we did above

X_cd_r = X_cd.affineMap(rot_mat, []);

figure;
Star.plot(X_cd_r, 'r');


%% Rotation examples with more basis vectors

% What if we have several basis vectors and constraints?

% Center should always be [0 sqrt(2)/2], only the vectors and constraints
% should change

V = [0, sqrt(2)/4, -sqrt(2)/4 1;
sqrt(2)/2, sqrt(2)/4, sqrt(2)/4 1];
C = [1 0 0;
    0 1 0;
    -1 0 0;
    0 -1 0;
    0 0 1;
    0 0 -1];
d = [1; 1; 1; 1; 1; 1];

R = Star(V,C,d);

figure;
Star.plot(R, 'c');
hold on;
scatter(R.V(1,1), R.V(2,1), 'xr');
xlim([-3 3])
ylim([-3 3])


% Rotate the set above 45 degrees

R_45 = R.affineMap(rot_mat, []);

figure;
Star.plot(R_45, 'c');
hold on;
scatter(R_45.V(1,1), R_45.V(2,1), 'xr');
xlim([-3 3])
ylim([-3 3])

% Rotate another 45 degrees

R_90 = R_45.affineMap(rot_mat, []);

% figure;
% Star.plot(R_90, 'c');
% hold on;
% scatter(R_90.V(1,1), R_90.V(2,1), 'xr');

% Rotate another 45 degrees

R_135 = R_90.affineMap(rot_mat, []);

% figure;
% Star.plot(R_135, 'c');
% hold on;
% scatter(R_135.V(1,1), R_135.V(2,1), 'xr');

% Rotate another 45 degrees

R_180 = R_135.affineMap(rot_mat, []);

% figure;
% Star.plot(R_180, 'c');
% hold on;
% scatter(R_180.V(1,1), R_180.V(2,1), 'xr');

% Rotate another 45 degrees

R_225 = R_180.affineMap(rot_mat, []);

% figure;
% Star.plot(R_225, 'c');
% hold on;
% scatter(R_225.V(1,1), R_225.V(2,1), 'xr');

% Rotate another 45 degrees

R_270 = R_225.affineMap(rot_mat, []);

% figure;
% Star.plot(R_270, 'c');
% hold on;
% scatter(R_270.V(1,1), R_270.V(2,1), 'xr');

% Rotate another 45 degrees

R_315 = R_270.affineMap(rot_mat, []);

% figure;
% Star.plot(R_315, 'c');
% hold on;
% scatter(R_315.V(1,1), R_315.V(2,1), 'xr');

% Rotate another 45 degrees

R_360 = R_315.affineMap(rot_mat, []);

% figure;
% Star.plot(R_360, 'c');
% hold on;
% scatter(R_360.V(1,1), R_360.V(2,1), 'xr');


% The sets are rotating, but are also getting translated, or at least that
% what it looks like...

% This could be fixed by adding a translating vector (offset) to the star,
% which could be calculated afterwards and applied in two steps

% Target center
center = [R.V(1,1), R.V(2,1)];

% Rotated set (45 degrees)
center_45 = [R_45.V(1,1), R_45.V(2,1)];
center_offset = center_45 - center;

R_45_fixed = R_45;
R_45_fixed.V(:,1) = R_45.V(:,1) - center_offset';
figure;
Star.plot(R_45_fixed, 'c');
hold on;
scatter(R_45_fixed.V(1,1), R_45_fixed.V(2,1), 'xr');
xlim([-3 3])
ylim([-3 3])


% We can also do it in a single step by eaving the center (obj.V(:,1)) unchanged

% make copy of set
R45 = R;
% matrix multiplication part of affineMapping
R45.V(:,2:end) = rot_mat * R.V(:, 2:end);

figure
Star.plot(R45, 'c');
hold on;
scatter(R45.V(1,1), R45.V(2,1), 'xr');
xlim([-3 3])
ylim([-3 3])

% The second approach is better and faster
% With some more testing, we can add some rotating function to the Star
% method

% One of the things we need to do is convert the predicate bounds to the
% C,d variables for the rotation to work

%% Example 3: rotation without C,d, only predicate bounds

V = [0, sqrt(2)/4, -sqrt(2)/4 1;
sqrt(2)/2, sqrt(2)/4, sqrt(2)/4 1];
% Instead of this, we define the constraints as predicate bounds
C = [0 0 0];
d = [0];
pred_lb = [-1; -1; -1];
pred_ub = [1; 1; 1];

S = Star(V,C,d, pred_lb, pred_ub);

figure;
Star.plot(S, 'g');
hold on;
scatter(S.V(1,1), S.V(2,1), 'xr');

S45 = S;
S45.V(:,2:end) = rot_mat * S.V(:, 2:end);

figure
Star.plot(S45, 'g');
hold on;
scatter(S45.V(1,1), S45.V(2,1), 'xr');
xlim([-3 3])
ylim([-3 3])

% In this example, this also works, but does it always?

%% Example 4: rotation with both C,d and predicate bounds

V = [0, sqrt(2)/4, -sqrt(2)/4 1;
sqrt(2)/2, sqrt(2)/4, sqrt(2)/4 1];
C = [1 -0.5 0;
    0.5 1 0;
    -1 0 0;
    0 -1 0;
    0 0 1;
    0 0 -1];
d = [1; 1; 1; 1; 1; 1];
pred_lb = [0.5; -1; -1];
pred_ub = [2; 1; 1];

I = Star(V,C,d, pred_lb, pred_ub);

figure;
Star.plot(I, 'm');
hold on;
scatter(I.V(1,1), I.V(2,1), 'xk');
xlim([-3 3])
ylim([-3 3])

% rotate 45 degrees

I45 = I;
I45.V(:,2:end) = rot_mat * I.V(:, 2:end);

figure
Star.plot(I45, 'm');
hold on;
scatter(I45.V(1,1), I45.V(2,1), 'xk');
xlim([-3 3])
ylim([-3 3])

%% Example 5: exact reachable sets

load("example_stars.mat");

E = Re(11);

figure;
Star.plot(E,'m');
hold on;
scatter(E.V(1,1), E.V(2,1), 'xk');
xlim([-13 -6])
ylim([6 12])

% rotate 45 degrees
angle = 45;
angle = deg2rad(angle);
% rotation matrix
rot_mat = [cos(angle), -sin(angle);
    sin(angle), cos(angle)];

E45 = E;
% rotate
E45.V(:,2:end) = rot_mat * E45.V(:,2:end);
figure
Star.plot(E45, 'm');
hold on;
scatter(E45.V(1,1), E45.V(2,1), 'xk');
xlim([-12 -5])
ylim([6 12])

% Looks good too, so let's rotate it back -45 degrees
% rotate 45 degrees
angle = -45;
angle = deg2rad(angle);
% rotation matrix
rot_mat = [cos(angle), -sin(angle);
    sin(angle), cos(angle)];

E_ = E45;
% rotate
E_.V(:,2:end) = rot_mat * E_.V(:,2:end);
figure
Star.plot(E_, 'm');
hold on;
scatter(E_.V(1,1), E_.V(2,1), 'xk');
xlim([-13 -6])
ylim([6 12])


% TODO: update affineMap to allow scalar multiplication (W) or empty so
% that we it can be used for translation only as well (b)

% QUESTION: how can we rotate a union of stars? (e.g. compute a set using
% exact reachability, which outputs N Stars, and the output set is the
% union of these sets, and then rotate that set by X degrees?
% - Exact method -> ?
% - Approx -> compute convex hull and rotate