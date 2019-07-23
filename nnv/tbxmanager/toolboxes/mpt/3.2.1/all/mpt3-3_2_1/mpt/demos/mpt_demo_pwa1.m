function mpt_demo_pwa1
% This demo shows how to model PWA systems in MPT3.
%
% As an example, we will assume the following dynamics:
%
%        { A1*x + B*u if x_1 <= 0
%  x^+ = {
%        { A2*x + B*u if x_1 >= 0
%
%    y = C*x + D*u
%

% PWA systems are created by defining each dynamics as an LTI system:
B = [0; 1]; C = [1 0]; D = 0;

% First dynamics:
alpha = -pi/3;
A1 = 0.8*[cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];

dyn1 = LTISystem('A', A1, 'B', B, 'C', C, 'D', D);
% We need to tell that dynamics #1 should be active if x_1<=0:
dyn1.setDomain('x', Polyhedron([1 0], 0));

% Second dynamics:
alpha = pi/3;
A2 = 0.8*[cos(alpha) -sin(alpha); sin(alpha) cos(alpha)];

dyn2 = LTISystem('A', A2, 'B', B, 'C', C, 'D', D);
% Region of validity of the dynamics (x_1>=0):
dyn2.setDomain('x', Polyhedron([-1 0], 0));

% Create the PWA description using an array of LTI systems:
pwa = PWASystem([dyn1 dyn2]);

% Optionally we can set constraints:
pwa.x.min = [-10; -10];
pwa.x.max = [10; 10];
pwa.y.min = -10;
pwa.y.max = 10;
pwa.u.min = -1;
pwa.u.max = 1;

% Define an on-line MPC controller for such a system
horizon = 2;
onl_ctrl = MPCController(pwa, horizon);
% Set panalties used in the cost function:
onl_ctrl.model.x.penalty = OneNormFunction(10*eye(2));
onl_ctrl.model.u.penalty = OneNormFunction(1);

% Construct the explicit solution
exp_ctrl = onl_ctrl.toExplicit();

% Obtain the closed-loop optimizers for a particular initial condition
x0 = [-4; 0];
Uonl = onl_ctrl.evaluate(x0)
Uexp = exp_ctrl.evaluate(x0)

end
