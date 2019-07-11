function mpt_demo_lti4
%
% mpt_demo_lti4
%

% This demo demonstrates explicit MPC

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
% Use quadratic cost function
ctrl.model.x.penalty = QuadFunction(eye(model.nx));
ctrl.model.u.penalty = QuadFunction(eye(model.nu));

% Finally, we can convert the controller to an explicit form:
disp('Generating the explicit solution:');
expctrl = ctrl.toExplicit()

% Compare optimal solutions
x0 = [-4; 0];
disp('Optimal control input obtained by evaluating the explicit solution:');
Uexp = expctrl.evaluate(x0)

disp('Optimal control input obtained by solving the optimization problem on-line:');
Uonl = ctrl.evaluate(x0)

close all
% plot the explicit optimizer
expctrl.feedback.fplot();
title('PWA representation of the optimal control input')

% plot the value function
figure
expctrl.cost.fplot();
title('Explicit cost function');

% plot the regions
figure
expctrl.partition.plot()
title('Regions of the polyhedral partition');


end
