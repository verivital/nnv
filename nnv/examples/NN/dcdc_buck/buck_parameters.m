% This script loads the Buck converter parameters, simulates pid model
% generates data for training

clc; 
clear all; 

% initial conditions
i0 = 0;
v0 = 0;

% Define Buck converter parameters
T_sw = 1/10000; % switching time period
Vs = 30;%30
% Generate avriable ref sequence for training 
time1 = [0 .4 .8 1.2 1.6 2 2.4 2.8 3.2 3.6 4 4.4 4.8 5.2 5.6 6];
Vref = [10 10 10 11 12 13 11 10 9 8 6 7 8 9 10 10];
%
C = 2.2e-3;
L = 2.65e-3;
R = 8;
rs=100e-3;% switching loss 
rL=520e-3;% inductor loss

%SYSTEM Matrix
A = [-1*(rs+rL)/L, -1/L; 
    1/C, -1/(R*C)];
%Input matrices
u1 = [Vs/L 0]';
u2 = [0 0]';

%Stateflow Model Parameters

sys(1).Ac = A;
sys(1).Bc = u1;
sys(1).Bo = u2;
    
a00c = sys(1).Ac(1,1);
a01c = sys(1).Ac(1,2);
a10c = sys(1).Ac(2,1);
a11c = sys(1).Ac(2,2);
    
b0c = sys(1).Bc(1);
b1c = sys(1).Bc(2);
    
a00o = -rL/L;% sys(1).Ac(1,1) same as closed
a01o = sys(1).Ac(1,2);
a10o = sys(1).Ac(2,1);
a11o = sys(1).Ac(2,2);
    
b0o = sys(1).Bo(1);
b1o = sys(1).Bo(2);

sim ('buck_pid');
save 'buck_data';
save 'iodata'; % saves voltage data

input_data(1,:) =  (iodata.signals(1).values)';% error
input_data(2,:) =  (iodata.signals(3).values)';% Vref
input_data(3,:) =  (buck_data.signals(2).values)';% Vc
input_data(4,:) =  (buck_data.signals(1).values)';% iL
% input_cell =  input_data;
output_data =  (iodata.signals(2).values)'; %
% output_cell =  output_data; %