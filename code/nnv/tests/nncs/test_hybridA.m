% Test hybridA Construct and functions
% HybridA construct
test_HA = HybridA(2,0,bball,1);
dim = 2;
inputs = 0;
modes = 1;
dynamicsHA = bball();
test_HA1 = HybridA(dim,inputs,dynamicsHA,modes);

% Test simple functions to modify HybridA properties
test_HA1.set_errorOrder(1);
test_HA1.set_taylorTerms(8);
test_HA1.set_polytopeOrder(9);
test_HA1.set_tFinal(1.7);

% Reachability functions
% Zonotope
inp_set = Zono; % input set (0 inputs)
init_set = Zono([1;0],[0.05 0;0 0.05]); % initial set
[R,rt] = test_HA.reach_zono(init_set, inp_set, 0.1,1.7); % R = Reach set (2D), rt = reachability time
% Star set
inp_setS = Star; % input set (0 inputs)
init_setS = init_set.toStar; % initial set
[Rs] = test_HA.stepReachStar(init_setS, inp_setS);

% Simulation
x0 = [1; 0];
inp = 0; % No inputs
y = test_HA.evaluate(inp,x0);

