function [dx] = dynamics10(t,x,u)
% Ex_car_model
  %vehicleODE Bicycle model of a vehicle with
  % states
  %       x(1), x(2): x,y positions
  %       x(3): Yaw angle (\psi)
  %       x(4): velocity
  % control inputs
  %       u(1): acceleration m/s^2
  %       u(2): steering angle of front wheel
  % disturbance input
  %       w: disturbance with a range (-10^-4,10^4)
  % Initial state range [9.5, 9.55] × [-4.5, -4.45] × [2.1, 2.11] × [1.5, 1.51]
  
dx(1,1) = x(4) * cos(x(3));
dx(2,1) = x(4) * sin(x(3));
dx(3,1) = u(2);
dx(4,1) = u(1); % No noise
% dx(4,1) = u(1) + w; % noise
end