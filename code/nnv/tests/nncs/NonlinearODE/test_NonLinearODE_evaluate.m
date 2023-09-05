%% Test 1: constructor and evaluate

car = NonLinearODE(6, 1, @car_dynamics, 0.1, 1, eye(6));
% initial state
x0 = [49; 25; 9; 20; 1; 0];
% input
u = 1;
% evaluate / simulate
[t, y] = car.evaluate(x0, u);