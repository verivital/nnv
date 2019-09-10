function [inputs_cell, output_cell] =generateOnetrace()
global mdl Ts T G_ego t_gap D_default v_set amin_ego amax_ego
global x0_lead v0_lead x0_ego v0_ego seed
%Uniform randome generator 
seed = randi([0,1000]);
%Define the sample time, Ts, and simulation duration, T, in seconds.
Ts = 0.1;
T = 1000;

%Specify the linear model for ego car.
G_ego = tf(1,[0.5,1,0]);

t_gap = 1.4;
D_default = 10;

%Specify the driver-set velocity in m/s.
v_set = 30;

%The acceleration is constrained to the range [-3,2] (m/s^2).
amin_ego = -3;
amax_ego = 2;

%Specify the initial position and velocity for the two vehicles.
x0_lead = 50 + 1.5 - 3*rand(1);   % initial position for lead car (m)
v0_lead = 25 + 1.5 - 3*rand(1);  % initial velocity for lead car (m/s)

x0_ego = 10 + 1.5 - 3*rand(1);   % initial position for ego car (m)
v0_ego = 20 + 1.5 - 3*rand(1);  % initial velocity for ego car (m/s)

%Run the simulation.
sim(mdl)

% input and output data for training
v_ego = logsout.get(4).Values.Data; % input3
d_rel = logsout.get(7).Values.Data; % input4
v_lead = logsout.get(5).Values.Data;
v_rel = v_lead - v_ego; % input5

v_sets = v_set*ones(length(v_ego),1); %input1
t_gaps = t_gap*ones(length(v_ego),1); %input2

inputs_mat = [v_sets, t_gaps, v_ego, d_rel, v_rel];
inputs_cell = num2cell(inputs_mat',1);

a_ego = logsout.get(1).Values.Data; % output1
output_mat = a_ego;
output_cell = num2cell(output_mat',1);
