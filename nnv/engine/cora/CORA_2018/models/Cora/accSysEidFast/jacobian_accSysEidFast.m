function [A,B]=jacobian_accSysEidFast(x,u)

A=[0,1;...
0,-(1281*u(1))/(25*x(2)^2)];

B=[0;...
1281/(25*x(2))];

