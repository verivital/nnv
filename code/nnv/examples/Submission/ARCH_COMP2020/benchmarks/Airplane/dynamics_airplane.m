function [dstates] = dynamics(t,states,actions)
% Ex_airplane
% constants as per the documentation Ix=Iy=Iz = 1; Ixz = 0; m = 1; g= 1;
% description of the pararameters in these equations can be found here
% https://github.com/amaleki2/benchmark_closedloop_verification/blob/master/AINNC_benchmark.pdf
 
x = states(1);
y = states(2);
z = states(3);
u = states(4);
v = states(5);
w = states(6);
phi = states(7);
theta = states(8);
psi = states(9);
r = states(10);
p = states(11);
q = states(12);

Fx = actions(1);
Fy = actions(2);
Fz = actions(3);
Mx = actions(4);
My = actions(5);
Mz = actions(6);

T_psi = [cos(psi), -sin(psi), 0.;
         sin(psi),  cos(psi), 0.;
         0.,        0.,       1.];
        
T_theta = [cos(theta), 0., sin(theta);
           0.        , 1., 0.        ;
          -sin(theta), 0., cos(theta)];

T_phi = [1., 0.,        0.      ;
         0., cos(phi), -sin(phi);
         0., sin(phi),  cos(phi)];
     
mat_1 = T_psi * T_theta * T_phi;

mat_2 = [cos(theta), sin(theta) * sin(phi),  sin(theta) * cos(phi);
         0.,         cos(theta) * cos(phi), -cos(theta) * sin(phi);
         0.,         sin(phi),               cos(phi)];
mat_2 = 1 / cos(theta) * mat_2;

a1 = [u; v; w];
a2 = mat_1 * a1;

dx = a2(1);
dy = a2(2);
dz = a2(3);

a3 = [p; q; r];
a4 = mat_2 * a3;

dphi = a4(1);
dtheta = a4(2);
dpsi = a4(3);

du = -sin(theta)            + Fx - q * w + r * v;
dv =  cos(theta) * sin(phi) + Fy - r * u + p * w;
dw =  cos(theta) * cos(phi) + Fz - p * v + q * u;

dp = Mx;
dq = My;
dr = Mz;

dstates(1,1) = dx;
dstates(2,1) = dy;
dstates(3,1) = dz;
dstates(4,1) = du;
dstates(5,1) = dv;
dstates(6,1) = dw;
dstates(7,1) = dphi;
dstates(8,1) = dtheta;
dstates(9,1) = dpsi;
dstates(10,1) = dp;
dstates(11,1) = dq;
dstates(12,1) = dr;
