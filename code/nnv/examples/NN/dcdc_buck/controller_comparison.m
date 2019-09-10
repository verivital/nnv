% This script simulates pid and nn controlled models and generates comp
% plots
clc
clear
close all
% initial conditions
i0 = 0;
v0 = 0;

% Define Buck converter parameters
T_sw = 1/10000; % switching time period
Vs = 30;
Vref = 10;
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

sim ('test_buck_pid'); %pid model sim
sim ('test_buck_nn'); %nn-controlloed model sim

%plot comp figure
plot(nn_buck_data.time, nn_buck_data.signals(2).values)
hold on
plot(pid_buck_data.time, pid_buck_data.signals(2).values)
title('Output Voltage of the Buck Converter')
grid on
legend('NN contoller', 'PID controller');
xlabel('Time, Sec');%,'LineWidth',2)'FontWeight','bold',
ylabel('v_{C}, V');
hold off
figure
plot(nn_buck_data.time, nn_buck_data.signals(1).values)
hold on
plot(pid_buck_data.time, pid_buck_data.signals(1).values)
title('Inductor Current of the Buck Converter')
grid on
legend('NN contoller', 'PID controller');
hold off
xlabel('Time, Sec');%,'LineWidth',2)'FontWeight','bold',
ylabel('i_{L}, A');
figure
plot(nn_buck_data.signals(1).values,nn_buck_data.signals(2).values)
hold on
plot(pid_buck_data.signals(1).values,pid_buck_data.signals(2).values)
title('Phase Plane Plot of the Buck Converter')
grid on
legend('NN contoller', 'PID controller');
xlabel('i_{L}, A');%,'LineWidth',2)'FontWeight','bold',
ylabel('v_{C}, V');
hold off