function [dx]=dynamics(x,u)

dx(1,1) = x(3)^3 - x(2);
dx(2,1) = x(3);
dx(3,1) = u;

end

