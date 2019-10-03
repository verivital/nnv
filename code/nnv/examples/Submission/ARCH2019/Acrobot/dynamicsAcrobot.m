function [dx, y] = dynamicsAcrobot(t,x,u)
% Parameters
l1 = 1; %length of link 1
l2 = 1; %length of link 2
m1 = 1; %mass of link 1
m2 = 1; %mass of link 2
lc1 = 0.5; %position of center of mass of link 1
lc2 = 0.5; %position of center of mass of link 2
moi = 1; %moments of inertia for both links
max_v1 = 4*pi; %Max velocity link 1
max_v2 = 9*pi; %Max velocity link 2
g = 9.8; %gravity constant

% Equations
d1 = (m1 * lc1^2 + m2* (l1^2 + lc2^2 + 2 * l1 * lc2 * cos(theta2)) + 2*moi);
d2 = (m2 * (lc2^2 + l1 * lc2 *cos(theta2)) + moi);
phi1 = (-m2 * l1 * lc2  *omega2^2 * sin(theta2) - 2 * m2  * l1 * lc2 * omega2 * omega1  *sin(theta2) + (m1 * lc1 + m2 * l1) * g * cos(theta1 - pi/2) + phi2);
phi2 = (m2 * lc2 * g * cos(theta1 + theta2 - pi/2));

% ODEs
dx(1,1) = x(3);
dx(2,1) = x(4);
dx(3,1) = -(d1^(-1)) * (d1 * x(4) + phi1);
dx(4,1) = (m2*lc1^2 + moi - (d2^2 / d1))^(-1) * (u + d2/d1 * phi1 -phi2);

y = [cos(x(1)); cos(x(2)); sin(x(1)); sin(x(2)); x(3); x(4)];
end

