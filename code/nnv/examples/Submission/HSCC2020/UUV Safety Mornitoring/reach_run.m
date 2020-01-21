% Create script to run the reachability anlysis for the LEC1 for entire
% experiment
clc;clear;close all;

% Load the first part of the NN controller (LEC)
Cont = load_nn(' ',' ',' ',' ',' ','model.ckpt-558138.mat');
% Load the second part of the NN controller (LEc)
Norm = load_nn(' ',' ',' ',' ',' ','LEC1_NormFcnRelu.mat');
% Load the SysID model of the UUV
load('blackbox_sys_6s_3'); % sys is loaded

% Load data from one of the systemId experiments
data = load('Exp10a_all.csv');
inputs = [data(:,34) data(:,11)];
outputs = [data(:,12:13) data(:,23)];
% Load the pipe segments info
pipe = load('Exp10_pipe.csv');
% Obstacle info
obstacle2 = [-1869.3 38.1835 -46.458];
% Create data object
Ts = 1; % time step
ida = iddata(outputs,inputs,Ts);

% Check if the SysID is continuous or discrete (time)
if sys.Ts == 0
    Plant = LinearODE(sys.A,sys.B,sys.C,sys.D);
else
    Plant = DLinearODE(sys.A,sys.B,sys.C,sys.D,sys.Ts);
end
ReachSystem = []; % Initialize trajectory reach set to plot
%% Compute the reachable set of the LEC1 combined
% Inputs to LEC1 are:
%     1. Pipeline Orientation           (rad?)
%     2. Pipe Range (Port Side)         (m?)
%     3. Pipe Range (Starboard Range)   (m?)
%     4. Time Since Last Detection      (s?)
%     5. 1" SAS Range (Port Side)       (?)
%     6. 1" SAS Range (Starboard Side)  (?)
%     7. Closest Point of Approach (CPA)(?)
%     8. Nearest Obstacle Range         (m?)

% Run the loop to generate the reach sets at each time step
for k=1:length(data)
    % Use current and future data to create reach sets
    iL = data(k,2:9);
    b = 0.02;
    lb = [iL(1)-b; iL(2)-b; iL(3)-b ; max(0,iL(4)-b); max(0,iL(5)-b); max(0,iL(6)-b) ;-60;-60]; 
    up = [iL(1)+b; iL(2)+b; iL(3)+b; iL(4); iL(5) ;iL(6); -60; -60];
    % Create start set for input to LEC1
    InitSetC = Star(lb,up); % input set to LEC1

    % Generate random numbers within the input set to simulate
    np = 100; %  number of points to simulate
    rsim = (up-lb).*rand(8,np) + lb;
    % compute reach set of controller
    OutC = Cont.reach(InitSetC,'exact-star',4); %Approximate or exact method using stars
    OutC_e = zeros(4,np); % memory allocation
    % Simulate Controller
    for i=1:np
        OutC_e(:,i) = Cont.evaluate(rsim(:,i));
    end

    % Compute reach set for the second part of the LEC
    OutLEC = Norm.reach(OutC,'exact-star',4);

    % System dynamics

    %estimate initial states for current data point
    x0 = findstates(sys,ida(k)); 
    % Create a set around the initial states of the system
    ab = 0.001;
    lb = zeros(length(x0),1);
    ub = zeros(length(x0),1);
    for i = 1:1:length(x0)
        lb(i) = x0(i) - ab; 
        ub(i) = x0(i) + ab;
    end
    % Create initial state set based on collected data
    init_set = Star(lb,ub);
    % Compute reach set of the states of the sysID
    n = length(OutLEC);
    OutPlant = [];
    if sys.Ts == 0 % continuous-time
        for i=1:n
            OutPlant = [OutPlant Plant.simReach('direct', init_set, OutLEC(i), 1, 1)];
        end
    else % discrete-time
        for i=1:n
            OutPlant = [OutPlant Plant.stepReachStar(init_set,OutLEC(i))];
        end
    end
    % Compute reach set of the outputs of the sysID
    n = length(OutLEC);
    OutSys = [];
    for i=1:n
        OutSys = [OutSys OutPlant(i).affineMap(sys.C,[])];
        OutSys(i) = OutSys(i).affineMap([1 0 0; 0 1 0],[]);
    end
    % Concatenate all output sets of the SysID to plot later
    ReachSystem = [ReachSystem OutSys];

%     clear OutSys OutC OutLEC OutPlant init_set initSetC;
end

% Begin the simulation of the system
time = length(data)-1; % simulation time
x = zeros(time,4); % memory allocation
x2 = zeros(time,2);
% Simulate the LEC
for i=1:time
    x(i,:) = Cont.evaluate(data(i,2:9)');
    x2(i,:) = Norm.evaluate(x(i,:)');
end
% Simulate the SysID
[y00,fit00,x_0] = compare(ida,sys);
x4 = lsim(sys,x2,linspace(0,time-1,time)',x_0);
%% Plot everything in same figure
%aa = figure('WindowState','maximized'); % doesn't work on codeocean
aa = figure();
aa1 = plot(x4(:,1),x4(:,2),'LineWidth',4);
hold on;
aa2 = plot(data(:,12),data(:,13),'LineWidth',4);
aa3 = plot(pipe(1:11,2),pipe(1:11,1),'LineWidth',4);
scatter(obstacle2(1), obstacle2(2),'s','LineWidth',3);
%Star.plotBoxes_2D_noFill(ReachSystem, 1, 2, 'red')
Star.plots(ReachSystem); % exact method, slower and problems with mpt
xlabel('x-pos (m)');
ylabel('y-pos (m)');
legend('SysID','UUVSim','pipe','obstacle','Reach Sets','FontSize',24,'Location','southeast');
title('ReachSet LEC1+SysID');
uistack(aa1,'top');
uistack(aa2,'top');
uistack(aa3,'top');
set(gca,'FontSize',18);
axis([-2000 -1625 30 200])
saveas(aa,'safety_monitoring','png');
% save('Exp10a_reachset_b001','ReachSystem');
