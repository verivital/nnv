function [dx] = dynamics9(t,x,u)
% Ex_tora
% Initial state range:
% [0.6, 0.7] × [-0.7, -0.6] × [-0.4, -0.3] × [0.5, 0.6]
  
dx(1,1) = x(2);
dx(2,1) = - x(1) + 0.1*sin(x(3));
dx(3,1) = x(4);
dx(4,1) = u;

end