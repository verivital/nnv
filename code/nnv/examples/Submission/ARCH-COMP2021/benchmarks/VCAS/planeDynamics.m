function dx = planeDynamics(x,u,T)
%PLANEDYNAMICS 
%   These are the equations used to describe the VCAS airplane dynamics
% 
% Dynamics: (T = 1)
%
% h(t+1) = h - h0_dot * T - 0.5* h0_dot_dot*T^2;
% h0_dot(t+1) = h0_dot + h0_dot_dot*T;
% tao(t+1) = tao - T;
% adv(t+1) = adv;
%
% adv = u(1); New advisory, the output of the neural network controllers
% h0_dot_dot = u(2); % Acceleration chosen on the basis of the new advisory
%
% Originally there are up to 3 possible acceleration values to choose from for each advisory.
% When proving safety, all possible values should be accounted for.
%
% The tools that cannot handle efficiently quantification over all possible
% acceleration can follow the following two strategies for picking up
% next acceleration:
%
% -- a worst-case scenario selection, where we choose the acceleration 
%    that will take the ownship closer to or less far apart from the intruder.
%
% -- always select the acceleration in the middle.

dx(1,1) = x(1) - x(2)*T - 0.5*u(2)*T^2;
dx(2,1) = x(2) + u(2)*T;
dx(3,1) = x(3) - T;
dx(4,1) = u(1);

end

