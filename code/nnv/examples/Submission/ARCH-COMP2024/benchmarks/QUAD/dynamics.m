function dx = dynamics(x,u)
% Parameters 
g = 9.81;    m = 1.4;     Jx = 0.054;
Jy = 0.054;  Jz = 0.104;  t = 0;
% Quadrotor with 12 states (x) and 3 inputs (u)
dx(1,1) = cos(x(8))*cos(x(9))*x(4) + (sin(x(7))*sin(x(8))*cos(x(9)) - cos(x(7))*sin(x(9)))*x(5) + (cos(x(7))*sin(x(8))*cos(x(9)) + sin(x(7))*sin(x(9)))*x(6) ;
dx(2,1) = cos(x(8))*sin(x(9))*x(4) + (sin(x(7))*sin(x(8))*sin(x(9)) - cos(x(7))*cos(x(9)))*x(5) + (cos(x(7))*sin(x(8))*sin(x(9)) + sin(x(7))*cos(x(9)))*x(6);
dx(3,1) = sin(x(8))*x(4) - sin(x(7))*cos(x(8))*x(5) - cos(x(7))*cos(x(8))*x(6);
dx(4,1) = x(12)*x(5) * x(11)*x(6) - g*sin(x(8));
dx(5,1) = x(10)*x(6) - x(11)*x(6) - g*sin(x(8));
dx(6,1) = x(11)*x(4) - x(10)*x(5) + g*cos(x(8))*cos(x(7)) - g - u(1)/m;
dx(7,1) = x(10) + sin(x(7))*tan(x(8))*x(11) + cos(x(7))*tan(x(8))*x(12);
dx(8,1) = cos(x(7))*x(11) - sin(x(7))*x(12);
dx(9,1) = sin(x(7))*x(11)/cos(x(8)) - cos(x(7))*x(12)/cos(x(8));
dx(10,1) = x(11)*x(12)*(Jy-Jz)/Jx + u(2)/Jx;
dx(11,1) = (Jz - Jx)*x(10)*x(12)/Jy + u(3)/Jy;
dx(12,1) = (Jx - Jy)*x(10)*x(11)/Jz + t/Jz;
end

