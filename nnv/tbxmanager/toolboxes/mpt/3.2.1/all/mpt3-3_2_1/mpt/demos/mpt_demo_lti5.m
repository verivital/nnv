function mpt_demo_lti5
%
% mpt_demo_lti5
%


% demonstrates usage of additional properties
close all

% define an LTI prediction model
lti = LTISystem('A', [1 1; 0 1], 'B', [1; 0.5]);

% define the MPC controller
horizon = 5;
ctrl = MPCController(lti, horizon);

% define quadratic penalties
ctrl.model.x.penalty = QuadFunction(eye(2));
ctrl.model.u.penalty = QuadFunction(1);

% add a terminal set constraint (see help SystemSignal/filter_terminalSet)
ctrl.model.x.with('terminalSet');
ctrl.model.x.terminalSet = Polyhedron('lb', [-1; -1], 'ub', [1; 1]);

% add an LQR terminal penalty (see help SystemSignal/filter_terminalPenalty)
lqr_penalty = ctrl.model.LQRPenalty();
ctrl.model.x.with('terminalPenalty');
ctrl.model.x.terminalPenalty = lqr_penalty;

% add a move-blocking constraint (the last 3 moves are to be constant)
ctrl.model.u.with('block');
ctrl.model.u.block.from = ctrl.N-2;
ctrl.model.u.block.to = ctrl.N;

% obtain the optimal control input
x0 = [-4; 0];
[Uonl, feasible] = ctrl.evaluate(x0)

% we can also ask for full open-loop predictions:
[~, ~, openloop] = ctrl.evaluate(x0)

% plot the open-loop predictions
ctrl.model.plot()

end
