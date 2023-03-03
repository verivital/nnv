%% Reachability analysis of the buck converter using NNV

% author: Diego Manzanas 
% date: December 6th 2020

%% 1) Setup
% Load hybrid automata
buck = HybridA(5,1,buck_v2,2);

% Setup reachability parameters
steps = 10; % Number of reachability steps
Ts = 1.6667e-05 / 10; % 1/10 switching frequency
buck.set_tFinal(steps*Ts); % final time
buck.set_timeStep(Ts); % set time step for reachability of the hybrid system

%% 2) Reachability and simulation
% Execute reachability analysis star set NNV
inp_set = Star; % Input set
lb = [0; 0; 0; 0.75; 0]; % lower bound initial state set
ub = [0.11; 0; 0; 0.75; 0.11]; % upper bound initial state set
init_set = Star(lb,ub); % initial state set as a Star
% Reachability execution
S = buck.stepReachStar(init_set,inp_set); % S corresponds to the Star set(s) at t=tFinal
Sall = buck.intermediate_reachSet; % get all reachable sets
% Simulation
x0sim = [0.01; 0; 0; 0.75; 0.01]; % initial state
[tT,yT,locT] = buck.evaluate(0,x0sim); % simulation trayectory

%% 3) Visualize results

figure;
Star.plotBoxes_2D_noFill(Sall,3,1,'b');
hold on;
for i=1:length(yT)
    plot(yT{i}(:,3),yT{i}(:,1),'r')
end
xlabel('x_3');
ylabel('x_1');

figure;
Star.plotBoxes_2D_noFill(Sall,3,5,'b');
for i=1:length(yT)
    plot(yT{i}(:,3),yT{i}(:,5),'r')
end
xlabel('x_3');
ylabel('x_5');

figure;
Star.plotBoxes_2D_noFill(Sall,1,5,'b');
for i=1:length(yT)
    plot(yT{i}(:,1),yT{i}(:,5),'r')
end
xlabel('x_1');
ylabel('x_5');

figure;
Star.plotBoxes_2D_noFill(Sall,3,2,'b');
for i=1:length(yT)
    plot(yT{i}(:,3),yT{i}(:,2),'r')
end
xlabel('x_3');
ylabel('x_2');

figure;
Star.plotBoxes_2D_noFill(Sall,4,5,'b');
for i=1:length(yT)
    plot(yT{i}(:,4),yT{i}(:,5),'r')
end
xlabel('x_4');
ylabel('x_5');
