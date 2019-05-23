function [dx] = dynamicsIP(t,x,u)
% Inverted pendulum dynamics
% x1 = position
% x2 = velocity
% x3 = theta
% x4 = omega

dx(1,1) = x(2);
dx(2,1) = 0.004300000000000637 * x(4) - 2.75 * x(3) + 1.9399999999986903 * u - 10.950000000011642 * x(2);
dx(3,1) = x(4);
dx(4,1) = 28.580000000016298 * x(3) - 0.04399999999998272 * x(4) - 4.440000000002328 * u + 24.919999999983702 * x(2);

end

