function ds = VCASVerifModel2Agent(s,u)

% Differential equation describing the dynamics of the physical component
% of the system
% 
% INPUTS
%
% s: state vector of the physical component of the system
% u: control input
%
% OUTPUTS
%
% ds: derivative of the state vector.

h = s(1); 
hown_dot = s(2);
hint_dot = s(3);
tau = s(4);

dh = hint_dot - hown_dot;
dhown_dot = u(1);
dhint_dot = u(2);
dtau = -1;

ds = [dh;dhown_dot;dhint_dot;dtau];