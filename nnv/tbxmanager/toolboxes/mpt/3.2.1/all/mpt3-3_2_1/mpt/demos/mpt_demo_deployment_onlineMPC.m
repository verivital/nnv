function mpt_demo_deployment_onlineMPC
% demostration how to deploy online MPC controller in real-time
% using RTW

% oscillator model

A=[   0.5403   -0.8415; 0.8415    0.5403];
B=[ -0.4597; 0.8415];
C=[1 0];
D=0;

% linear discrete-time model with sample time 1
sys = ss(A,B,C,D,1);

model = LTISystem(sys);

% set constraints on output
model.y.min = -10;
model.y.max = 10;

% set constraints on input
model.u.min = -1;
model.u.max = 1;

% weights on states/inputs
model.x.penalty = QuadFunction(eye(2));
model.u.penalty = QuadFunction(1);

% terminal set
Tset = model.LQRSet;

% terminal weight
PN = model.LQRPenalty;

% add terminal set and terminal penalty
model.x.with('terminalSet');
model.x.terminalSet = Tset;
model.x.with('terminalPenalty');
model.x.terminalPenalty = PN;


% MPC controller
ctrl = MPCController(model,5);

% export variables to YALMIP
v = ctrl.toYALMIP;

% initial condition is the parameter theta
theta = v.variables.x(:,1); 

% extract remaining variables
states = v.variables.x(:,2:end);
inputs = v.variables.u;
outputs = v.variables.y;

% transform the problem to parametric LCP
problem = Opt(v.constraints,v.objective,theta,[inputs(:); states(:); outputs(:)]);
problem.qp2lcp;

% open the simulink scheme
mpt_demo_rtw_onlinempc

% export the problem data to workspace for simulation to run
assignin('base','problem',problem);

end