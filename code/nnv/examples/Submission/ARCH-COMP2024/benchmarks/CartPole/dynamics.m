function dx = dynamics(x, f)

% Cartpole Swingup, 4 states x and the input f (action of the controller)
% The controller takes (x(1), x(2), x(3), x(4)) as input, its output can
% be used directly.

dx(1,1) = x(2);
dx(2,1) = 2*f;
dx(3,1) = x(4);
dx(4,1) = (0.08*0.41*(9.8*sin(x(3))-2*f*cos(x(3)))-0.0021*x(4))/0.0105;

end