%% Test 1: construct and reachability

car = NonLinearODE(6, 1, @car_dynamics, 0.01, 0.1, eye(6));

lb = [49; 25; 9; 20; 1; 0];
ub = [51; 25.2; 11; 20.2; 1; 0];

B1 = Box(lb, ub);

init_set = B1.toZono; % initial set of state

lb = 0.9;
ub = 1.0;
B1 = Box(lb, ub);

input_set = B1.toZono; % input set

timeStep = 0.01;
tFinal = 0.1;

[R, reachTime] = car.reach_zono(init_set, input_set, timeStep, tFinal);
