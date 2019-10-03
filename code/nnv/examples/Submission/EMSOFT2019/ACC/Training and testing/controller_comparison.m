close all
mdl = 'test_mpcACCsystem';
open_system(mdl)

%Define the sample time, Ts, and simulation duration, T, in seconds.
Ts = 0.1;
T = 400;
seed = randi([1001,2000]);
%Specify the linear model for ego car.
G_ego = tf(1,[0.5,1,0]);

%Specify the initial position and velocity for the two vehicles.
x0_lead = 50;   % initial position for lead car (m)
v0_lead = 25;   % initial velocity for lead car (m/s)

x0_ego = 10;   % initial position for ego car (m)
v0_ego = 20;   % initial velocity for ego car (m/s)

t_gap = 1.4;
D_default = 10;

%Specify the driver-set velocity in m/s.
v_set = 30;

%The acceleration is constrained to the range [-3,2] (m/s^2).
amin_ego = -3;
amax_ego = 2;

%% Run the simulation with actual controller .
sim(mdl)

v_ego = logsout.get(4).Values.Data; % input3
d_rel = logsout.get(7).Values.Data; % input4
v_lead = logsout.get(5).Values.Data;
v_rel = v_lead - v_ego; % input5

a_ego = logsout.get(1).Values.Data; % output1

a_lead = logsout.get(6).Values.Data;

%% Run the simulation with NN controller .
sim('test_nncACCsystem')
times = logsout.get(1).Values.Time;
v_ego_nn = logsout.get(1).Values.Data; % input3
d_rel_nn = logsout.get(5).Values.Data; % input4
v_lead_nn = logsout.get(3).Values.Data;
v_rel_nn = v_lead_nn - v_ego_nn; % input5

a_ego_nn = logsout.get(2).Values.Data; % output1


%%
figure

subplot(4,1,1)
plot(times, a_lead)
title('Acceleration of the lead car')
grid on

subplot(4,1,2)
hold on
plot(times, a_ego_nn,'r')
plot(times, a_ego, 'b')
title('Acceleration of the ego car: NN controller gives a smoothier control singal')
legend('NN contoller', 'MPC controller')
grid on
hold off

subplot(4,1,3)
hold on
plot(times,v_ego_nn,'r')
plot(times,v_ego, 'b')
title('Velocity of the ego car: NN controller drives ego car with a smoothier velocity')
legend('NN contoller', 'MPC controller')
grid on
hold off

subplot(4,1,4)
hold on
plot(times, d_rel_nn,'r')
plot(times, d_rel, 'b')
safe_distance = D_default + t_gap*v_ego_nn;
plot(times, safe_distance, '-.k')
title('Relative distance vs. safe distance')
legend('NN contoller', 'MPC controller', 'Safe distance')
grid on
hold off


%% Remove example file folder from MATLAB path, and close Simulink model.
rmpath(fullfile(matlabroot,'examples','mpc','main'));