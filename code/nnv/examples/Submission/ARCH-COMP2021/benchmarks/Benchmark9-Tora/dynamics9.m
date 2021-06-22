function [dx] = dynamics9(x,u)
% Ex_tora
% Initial state range:
% [0.6, 0.7] � [-0.7, -0.6] � [-0.4, -0.3] � [0.5, 0.6]
%
% The output of the neural network f(x) needs to be normalized 
% in order to obtain u.
% Namely, u = f(x) - 10 
%
dx(1,1) = x(2);
dx(2,1) = - x(1) + 0.1*sin(x(3));
dx(3,1) = x(4);
dx(4,1) = u;

end