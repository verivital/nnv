function [dx,y] = dynamicsMPCQ(t,x,u)
% dynamics of the MPC quadrotor
% u(1) = pitch
% u(2) = roll
% u(3) = thrust

% x1 = x-position of planner
% x2 = y-position of planner
% x3 = z-position of planner
% x4 = x-position of quadrotor
% x5 = y-position of quadrotor
% x6 = z-position of quadrotor
% x7 = x-velocity of quadrotor
% x8 = y-velocity of quadrotor
% x9 = z-velocity of quadrotor

dx(1,1) = 0.25;
dx(2,1) = -0.25;
dx(3,1) = 0.25;
dx(4,1) = x(7);
dx(5,1) = x(8);
dx(6,1) = x(9);
dx(7,1) = 9.81 * tan(u(1));
dx(8,1) = -9.81 * tan(u(2));
dx(9,1) = u(3) - 9.81;

y(1,1) = x(4) - x(1);
y(2,1) = x(5) - x(2);
y(3,1) = x(6) - x(3);
y(4,1) = x(7);
y(5,1) = x(8);
y(6,1) = x(9);
end

