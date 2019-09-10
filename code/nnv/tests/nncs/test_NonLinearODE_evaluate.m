car = NonLinearODE(4, 1, @car_dynamics);

tspan = [0 1];
x0 = [49; 25; 9; 20];
u = 1;
[t, y] = car.evaluate(tspan, x0, u);