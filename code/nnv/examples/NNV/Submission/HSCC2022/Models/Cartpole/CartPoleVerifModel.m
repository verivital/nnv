function ds = CartPoleVerifModel(s,u)

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

m = 0.1;
l = 0.5;
g = 9.8;
mt = 1.1;

x = s(1); 
v = s(2);
theta = s(3);
omega = s(4);

dx = v;
dtheta = omega;
% Taken from Simulink Model --> theta = from 0 to +-2pi
dv =((u + m*l*omega*omega*sin(theta))/mt)-m*l*(g*sin(theta)-cos(theta)*...
    ((u+m*l*omega*omega*sin(theta))/mt))/...
    (l*(4/3-m*cos(theta)*cos(theta)/mt))*cos(theta)/mt;
domega = (g*sin(theta)-cos(theta)*((u+m*l*omega*omega*sin(theta))/mt))/...
    (l*(4/3-m*cos(theta)*cos(theta)/mt));

ds = [dx;dv;dtheta;domega];