function mpt_demo_lti3
%
% mpt_demo_lti3
%

% This demo shows how to simulate a closed-loop system consisting of a
% controlled system and an MPC controller

% First we create the prediction model:
A = [1 1; 0 1];
B = [1; 0.5];
C = [1 0];
D = 0;
model = LTISystem('A', A, 'B', B, 'C', C, 'D', D);

% Next we define an MPC controller:
horizon = 5;
ctrl = MPCController(model, horizon);

% Specify the MPC problem (constraints and penalties):
ctrl.model.x.min = [-5; -5];
ctrl.model.x.max = [5; 5];
ctrl.model.u.min = -1;
ctrl.model.u.max = 1;
ctrl.model.x.penalty = OneNormFunction(eye(model.nx)); % 1-norm type penalty
ctrl.model.u.penalty = InfNormFunction(eye(model.nu)); % Inf-norm type penalty

% Create the closed-loop system:
loop = ClosedLoop(ctrl, model);

% Simulate the closed loop from a given initial condition
x0 = [-4; 0];
N_sim = 20;
data = loop.simulate(x0, N_sim);

% Plot the simulated state trajectories
plot(0:N_sim, data.X);

end
