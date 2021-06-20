%% Compute the reachability analysis of the obstacle avoidance scenario for UUV
clc;clear;close all;

addpath(genpath(pwd))

% computation steps
% step 0: initial state of plant
% step 1: output set of the plant (input to controller)
% step 2: transform output of plant to input of Norm (NN)
% step 3: compute output of Norm (NN mdeling sensor info)
% step 4: compose output of Cont1 (NN)
% step 5: compute output of Cont2 (NN normalizing function of controller)
% step 6: compute reachable set of the plant with the new control input and initial state
% step 7: go back to step 1, ...

%% Load all elements 
% Discrete-time Plant
load('UUV_model.mat');
Plant = DLinearODE(sys.A,sys.B,sys.C,sys.D,sys.Ts); % Continuous
% Norm: sensor NN model
Norm = Load_nn('FNNsensor.mat'); 
% Controller
Cont = Load_nn('FNNcontroller.mat');
% Controller output normalization NN
Cont2 = Load_nn('FNNnorm.mat');
% Gazebo data
load('data_exp3.mat');% (28 variables, 1Hz)
% Obstacle locations
load('obstacle34');

%% Perform simulation of the system given initial conditions and obstacle location (Visualizing purpose)
% Find initial state of the plant based on recorded Gazebo data
ida = iddata([data(1,6:7) data(1,17)],[data(1,28) data(1,5)],1);
% x, y, yaw, heading change, speed
ipoint = data(1,6:7);
ihea = data(1,17);
fpoint = ipoint + [17*cos(0.3225) 17*sin(0.3225)];
x0 = findstates(sys,ida);
m = size(data,1); % length of data points
x = x0; % First step
t = [0; 1];
u = ida.u(1,:);

% Following same steps as reachability analysis algorithm
for i=1:m+6
    [y,temp,x] = lsim(sys,[u;u],[],x);
    y = y(1,:);x = x(2,:);
    out(i,:) = y;
    innorm = y - [obstacle(1,2) obstacle(1,1) 0];
    outnorm = Norm.evaluate(innorm');
    outlec = Cont.evaluate(outnorm);
    u = Cont2.evaluate(outlec);
    u = u';
end



%% Perform Reachability Analysis following the computation steps described

% Initial state set of the Plant
as = 0.0001;
lb = zeros(length(x0),1); % lower bound
ub = zeros(length(x0),1); % upper bound
for i = 1:1:length(x0)
    lb(i) = x0(i) - as; 
    ub(i) = x0(i) + as;
end
% Step 0
init_set = Star(lb,ub); % initial set of states (Star)

% Begin loop to compute the reachable sets of all the components
ReachSystem = []; % To store all output sets of the plant
% Uncertainty in obstacle position
unc = 1;
x_obsS = Star(obstacle(1,2)-unc,obstacle(1,2)+unc);
y_obsS = Star(obstacle(1,1)-unc,obstacle(1,1)+unc);
% Reachability methods for NNs
met = 'approx-star';


for i=1:m+6
    % Step 1.
    OutPlant = init_set.affineMap(sys.C,[]);
    ReachSystem = [ReachSystem OutPlant];
    % Step 2.
    In_Norm = uncertainObstacle(x_obsS,y_obsS,OutPlant);
    % Step 3.
    Out_Norm = Norm.reach(In_Norm,met,4);
    % Step 4.
    Out_Lec = Cont.reach(Out_Norm,met,4);
    % Step 5.
    In_Plant = Cont2.reach(Out_Lec,met,4);
    % Step 6.
    init_set = reachPlantStep(Plant,init_set,In_Plant);
    % Step 7. Start loop again
end

% Create real size obstacle (obstacle corners)
xo1 = obstacle(1,2) + 0.5;
xo2 = obstacle(1,2) - 0.5;
yo1 = obstacle(1,1) + 0.5;
yo2 = obstacle(1,1) - 0.5;

% Plot results
aa = figure('units','normalized','outerposition',[0 0 0.5 0.5]);
aa1 = plot(out(:,1),out(:,2),'b','LineWidth',3);
hold on;
aa2 = patch([xo1 ; xo1; xo2; xo2; xo1],[yo1 ; yo2; yo2; yo1; yo1],'g');
aa3 = viscircles([obstacle(1,2) obstacle(1,1)],2+unc,'Color','r');
aa4 = plot([ipoint(1) fpoint(1)],[ipoint(2) fpoint(2)], '--m','LineWidth',3);
Star.plotBoxes_2D_noFill(ReachSystem, 1, 2,'b');
xlabel('X Position (m)');
ylabel('Y Position (m)');
legend([aa1 aa2 aa3 aa4],{'Trajectory','Obstacle','Unsafe Region','N.O.P'},'Position',[0.3 0.7 0.1 0.1]);
set(gca,'FontSize',16);
set(gca,'DataAspectRatio',[1 1 1]);
title('Experiment 3');
grid;
saveas(aa,'figures/CaseStudy3_unc','png');


%%% Helper Functions
function X = uncertainObstacle(xobs,yobs,OutPlant)
    [xo(1),xo(2)] = xobs.getRanges;
    [yo(1),yo(2)] = yobs.getRanges;
    [m,M] = OutPlant.getRanges;
    xd = [];
    xd(1) = m(1)-xo(1);
    xd(2) = m(1)-xo(2);
    xd(3) = M(1)-xo(1);
    xd(4) = M(1)-xo(2);
    yd = [];
    yd(1) = m(2)-yo(1);
    yd(2) = m(2)-yo(2);
    yd(3) = M(2)-yo(1);
    yd(4) = M(2)-yo(2);
    m(1) = min(xd);
    M(1) = max(xd);
    m(2) = min(yd);
    M(2) = max(yd);
    X = Star(m,M);
end

function Y = reachPlantStep(plant,X0,controls)
    n = length(controls);
    Y = [];
    for i=1:n
        U = controls(i).orderReduction_box(length(X0.V));
        try
            yz =  stepReachPlant(Plant.A, Plant.B, X0, controls(i));
            Y = [Y yz];
        catch
            yz = plant.stepReachZono(X0.Z,controls(i).Z);
            Y = [Y yz.toStar];
        end
    end
    X = Star.get_hypercube_hull(Y);
%     X = Star.get_convex_hull(Y);
    X = X.toStar; 
end

% Discrete-time Systems
function X = stepReachPlant(A, B, X0, controls)
    
    n = length(controls);
    X = X0.affineMap(A, []);
    U = [];
    for i=1:n
        U = [U controls(i).affineMap(B, [])];
    end
    
    next_X = [];
    for i=1:n
        V = X.V + U(i).V;
        next_X = [next_X Star(V, U(i).C, U(i).d)];
    end
    
    X = Star.get_hypercube_hull(next_X);
    X = X.toStar;
    
end