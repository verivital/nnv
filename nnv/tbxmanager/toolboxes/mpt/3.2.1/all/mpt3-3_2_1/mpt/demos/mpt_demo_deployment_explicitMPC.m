function mpt_demo_deployment_explicitMPC
% demostration how to deploy explicit MPC controller in real-time
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

% explicit controller
ectrl = ctrl.toExplicit;

% export explicit controller to C
dir_name = 'rtw_explicitMPC';
file_name = 'EMPCcontroller';
ectrl.exportToC(file_name,dir_name);

% compile the S-function
p = [pwd,filesep,dir_name,filesep];
mex(['-LargeArrayDims -I',dir_name],[p,file_name,'_sfunc.c'],[p,file_name,'.c']);

% open the simulink scheme
mpt_demo_rtw_explicitmpc

% set the source files in the Simulink scheme
str = sprintf('"%s%s_sfunc.c" "%s%s.c"',p,file_name,p,file_name);
set_param('mpt_demo_rtw_explicitmpc','CustomSource',str);
set_param('mpt_demo_rtw_explicitmpc/Controller','FunctionName',[file_name,'_sfunc']);

end