function [A,B]=jacobian_vanderPolEq(x,u)

A=[0,1;...
- 2*x(1)*x(2) - 1,1 - x(1)^2];

B=[0;...
1];

