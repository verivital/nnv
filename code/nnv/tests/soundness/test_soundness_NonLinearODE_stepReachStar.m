% test_soundness_NonLinearODE_stepReachStar
% Regression test: NonLinearODE.stepReachStar must over-approximate the
% control set. Star.getZono returns [] when an LP in Star.getBox fails, and
% stepReachStar used to replace the control set by the single point u = 0
% in that case, so the plant reach set ignored the actual control inputs.
% It now falls back to the LP-free Star.getBoxFast, and raises an error if
% that box is unbounded.
% To run: results = runtests('test_soundness_NonLinearODE_stepReachStar')

%% Test 1: input set whose getBox fails is over-approximated, not replaced by u = 0
rng(42);

Tr = 0.001; % reachability time step for the plant
Tc = 0.1; % control period of the plant
C = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % output matrix

lb = [49; 25; 9; 20; 1; 0];
ub = [51; 25.2; 11; 20.2; 1; 0];
B0 = Box(lb, ub);
init_set = B0.toStar;

% non-empty control set u = 0.95 + 0.05*a, a in [-1, 1], i.e. u in [0.9, 1.0],
% with no stored zonotope and a getBox that reports an LP failure
u_lb = 0.9;
u_ub = 1.0;
input_set = StarFailingGetBox([0.95 0.05], [1; -1], [1; 1], -1, 1);
assert(isempty(input_set.getZono), 'getZono must fail for this input set');

car = NonLinearODE(6, 1, @car_dynamics, Tr, Tc, C);
s = warning('off', 'NonLinearODE:inputSetBoxFallback');
R = car.stepReachStar(init_set, input_set);
warning(s);
assert(~isempty(R), 'reach set must not be empty');

R_lb = inf(6, 1);
R_ub = -inf(6, 1);
for k = 1:numel(R)
    Bk = R(k).getBoxFast;
    R_lb = min(R_lb, Bk.lb);
    R_ub = max(R_ub, Bk.ub);
end

% every simulated state at t = Tc, for any u in [0.9, 1.0], must be in the
% reach set (with u = 0 the ego car's internal state x6 stays near 0, while
% u in [0.9, 1.0] drives it to about 0.16..0.18)
sim_car = NonLinearODE(6, 1, @car_dynamics, Tr, Tc, C);
tol = 1e-6;
n_samples = 10;
X0 = [lb, ub, lb + rand(6, n_samples) .* (ub - lb)];
U0 = [u_lb, u_ub, u_lb + rand(1, n_samples) * (u_ub - u_lb)];
for k = 1:size(X0, 2)
    [~, y] = sim_car.evaluate(X0(:, k), U0(k));
    x_end = y(end, :)';
    assert(all(x_end >= R_lb - tol) && all(x_end <= R_ub + tol), ...
        'simulated state with u = %g is outside the reach set', U0(k));
end

% the fallback must be reported (promote its warning to an error to detect it)
car = NonLinearODE(6, 1, @car_dynamics, Tr, Tc, C);
s = warning('error', 'NonLinearODE:inputSetBoxFallback'); %#ok<CTPCT>
err_id = '';
try
    car.stepReachStar(init_set, input_set);
catch ME
    err_id = ME.identifier;
end
warning(s);
assert(strcmp(err_id, 'NonLinearODE:inputSetBoxFallback'), ...
    'stepReachStar must warn when it falls back to getBoxFast');

%% Test 2: input set with a stored zonotope uses getZono directly (no fallback)
rng(42);

Tr = 0.001;
Tc = 0.1;
C = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0];

lb = [49; 25; 9; 20; 1; 0];
ub = [51; 25.2; 11; 20.2; 1; 0];
B0 = Box(lb, ub);
init_set = B0.toStar;

B1 = Box(0.9, 1.0);
input_set = B1.toStar;

car = NonLinearODE(6, 1, @car_dynamics, Tr, Tc, C);
s = warning('error', 'NonLinearODE:inputSetBoxFallback'); %#ok<CTPCT>
R = car.stepReachStar(init_set, input_set);
warning(s);
assert(~isempty(R), 'reach set must not be empty');

%% Test 3: unbounded fallback box raises an error instead of returning a reach set
rng(42);

Tr = 0.001;
Tc = 0.1;
C = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0];

lb = [49; 25; 9; 20; 1; 0];
ub = [51; 25.2; 11; 20.2; 1; 0];
B0 = Box(lb, ub);
init_set = B0.toStar;

% predicate variable has no finite bounds, so getBoxFast is unbounded
input_set = StarFailingGetBox([0.95 0.05], [1; -1], [1; 1], -Inf, Inf);

car = NonLinearODE(6, 1, @car_dynamics, Tr, Tc, C);
err_id = '';
try
    car.stepReachStar(init_set, input_set);
catch ME
    err_id = ME.identifier;
end
assert(strcmp(err_id, 'NonLinearODE:unboundedInputSet'), ...
    'stepReachStar must reject an input set it cannot bound');
