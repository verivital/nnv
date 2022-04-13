function ds = ACASXuVerifModel1Agent(s,u)
%
% Differential equation describing the dynamics of the physical component
% of the system. 
%
% Here the dynamics corresponds to 1 Dubins airplanes
% 
% INPUTS
%
% s: state vector of the physical component of the system
% u: control input
%
% OUTPUTS
%
% ds: derivative of the state vector.

x_own = s(1); 
y_own = s(2);
psi_own = s(3);
x_int = s(4);
y_int = s(5);
psi_int = s(6);
v_own = s(7); % ft/s
v_int = s(8); % ft/s


% Ownship
dx_own = v_own*cos(psi_own); % x (ft)
dy_own = v_own*sin(psi_own); % y (ft)
dpsi_own = u;           % heading (rad)
dv_own = 0;

% Intruder

dx_int = v_int*cos(psi_int); % x (ft)
dy_int = v_int*sin(psi_int); % y (ft)
dpsi_int = 0; % Constant heading (rad)
dv_int = 0;


ds = [dx_own;dy_own;dpsi_own;dx_int;dy_int;dpsi_int;dv_own;dv_int];
end