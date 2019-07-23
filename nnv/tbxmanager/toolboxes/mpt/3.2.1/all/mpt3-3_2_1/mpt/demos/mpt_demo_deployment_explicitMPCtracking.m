function mpt_demo_deployment_explicitMPCtracking
% demostration how to deploy explicit MPC controller in real-time
% using RTW for tracking of a time-varying signal

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
model.u.min = -5;
model.u.max = 5;

% The objective is to track the time varying output while satisfying the
% constraints. To incorporate a time varying signal one needs to active the
% appropriate reference filter and mark it as "free"  
model.y.with('reference');
model.y.reference = 'free';

% objective function is to penalize ||y-ref||_1 + ||\Delta u||_1
model.y.penalty = OneNormFunction( 3 );
model.u.without('penalty');
model.u.with('deltaPenalty');
model.u.deltaPenalty = OneNormFunction( 0.5 );

% online controller with the horizon 4
ctrl = MPCController(model, 4);

% explicit controller
ectrl = ctrl.toExplicit();

% export explicit controller to C
dir_name = 'rtw_explicitMPCtracking';
file_name = 'EMPCcontroller';
ectrl.exportToC(file_name,dir_name);

% compile the S-function
p = [pwd,filesep,dir_name,filesep];
mex(['-LargeArrayDims -I',dir_name],[p,file_name,'_sfunc.c'],[p,file_name,'.c']);

% open the simulink scheme
mpt_demo_rtw_explicitmpctracking

% set the source files in the Simulink scheme
str = sprintf('"%s%s_sfunc.c" "%s%s.c"',p,file_name,p,file_name);
set_param('mpt_demo_rtw_explicitmpctracking','CustomSource',str);
set_param('mpt_demo_rtw_explicitmpctracking/Controller','FunctionName',[file_name,'_sfunc']);


end