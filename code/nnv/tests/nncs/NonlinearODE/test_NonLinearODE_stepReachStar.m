
Tr = 0.001; % reachability time step for the plant
Tc = 0.1; % control period of the plant
% output matrix
C = [0 0 0 0 1 0;1 0 0 -1 0 0; 0 1 0 0 -1 0]; % output matrix
car = NonLinearODE(6, 1, @car_dynamics, Tr, Tc, C);

lb = [49; 25; 9; 20; 1; 0];
ub = [51; 25.2; 11; 20.2; 1; 0];

B1 = Box(lb, ub);

init_set = B1.toStar; % initial set of state

lb = 0.9;
ub = 1.0;
B1 = Box(lb, ub);

input_set = B1.toStar; % input set

R = car.stepReachStar(init_set, input_set); % return reach set at t = tFinal