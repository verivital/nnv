function [dx] = dynamicsMC(t,x,u)
% Mountain car dynamics
% x(1) = position
% x(2) = velocity

dx(1,1) = x(2);
dx(2,1) = 0.0015*u-0.0025*cos(3*x(1));
end

