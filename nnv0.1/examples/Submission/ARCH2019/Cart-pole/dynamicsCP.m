function [dx] = dynamicsCP(t,x,u)
% x1 = position of cart (x)
% x2 = velocity of cart (v)
% x3 = pendulum angle (theta)
% x4 = angular velocity of pendulum (omega)

dx(1,1) = x(2);
dx(2,1) =((u + 0.05*x(4)^2*sin(x(3)))/1.1)-0.05*(9.8*sin(x(3))-cos(x(3))*((u+0.05*x(4)^2*sin(x(3)))/1.1))/(0.5*(4/3-0.1*cos(x(3))*cos(x(3))/1.1))*cos(x(3))/1.1;
dx(3,1) = x(4);
dx(4,1) = (9.8*sin(x(3))-cos(x(3))*((u+0.05*x(4)^2*sin(x(3)))/1.1))/(0.5*(4/3-0.1*cos(x(3))*cos(x(3))/1.1));


end

