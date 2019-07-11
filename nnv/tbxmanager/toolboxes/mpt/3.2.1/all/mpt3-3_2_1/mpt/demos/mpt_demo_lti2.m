function mpt_demo_lti2
%
% mpt_demo_lti2
%

% demonstrates most basic control synthesis
close all

% Define an LTI dynamics x^+ = A x + B u
% (note that we can ommit the output equation)
A = [1 1; 0 1];
B = [1; 0.5];
lti = LTISystem('A', A, 'B', B);

% Define an MPC controller using "lti" as the prediction model
ctrl = MPCController(lti);

% Set the prediction horizon
ctrl.N = 10;

% Add constraints on predicted states
ctrl.model.x.min = [-5; -5];
ctrl.model.x.max = [5; 5];

% Add constraints on predicted control inputs
ctrl.model.u.min = -1;
ctrl.model.u.max = 1;

% Use quadratic state penalty with identity weighting matrix
W = eye(2);
ctrl.model.x.penalty = QuadFunction(W);

% Set quadratic input penalty with identity weighting matrix
W = 1;
ctrl.model.u.penalty = QuadFunction(W);

% Obtain the optimal control input for a given initial condition
x0 = [-4; 0];
u = ctrl.evaluate(x0);

% We can also ask for the open-loop predictions:
[u, feasible, openloop] = ctrl.evaluate(x0);
openloop

% Plot the open-loop trajectories
ctrl.model.plot()

end
