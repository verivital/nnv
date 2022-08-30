%% Reachability analysis of Airplane Benchmark

%% load the controller
net = Load_nn('controller_airplane.mat');

% Specify the reach step, has to be smaller than the control period
reachStep = 0.05;
%% specify the control period as specified by the benchmark description
controlPeriod = 0.1;

% define the plant as specified by nnv
plant = NonLinearODE(12,6,@dynamics, reachStep, controlPeriod, eye(12));
plant.set_taylorTerms(20)
plant.set_zonotopeOrder(10);

%% Reachability analysis 
% Initial set
lb = [0; 0; 0; 0; 0.8; 0; 0; 0; 0; 0; 0; 0];
ub = [0; 0; 0; 0; 0.81; 0; 0; 0; 0; 0; 0; 0];
init_set = Star(lb,ub);
% Input set
lb = [0;0;0;0;0;0];
ub = [0;0;0;0;0;0];
input_set = Star(lb,ub);
% Store all reachable sets
reachAll = init_set;
% Execute reachabilty analysis
num_steps = 8;
t = tic;
for i=1:num_steps
    %  Compute controller output set
    input_set = net.reach(init_set,'approx-star');
    % Compute plant reachable set
    init_set = plantReach(plant,init_set,input_set,'lin');
    reachAll = [reachAll init_set];
end
t = toc(t);

path_out = [path_results(), filesep, 'Airplane', filesep];
mkdir(path_out)
save([path_out, 'reach.mat'],'reachAll','t','-v7.3');
%% Visualize results

f2 = figure;
rectangle('Position',[-0.5,-1,1,2],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
hold on;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,2,5,'b');
grid;
xlabel('x_2');
ylabel('x_5');
saveas(f2,[path_out, 'reach2vs5.pdf']);

f5 = figure;
rectangle('Position',[-1,-1,2,2],'FaceColor',[0 0.5 0 0.5],'EdgeColor','y', 'LineWidth',0.1)
hold on;
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,8,9,'b');
grid;
xlabel('x_8');
ylabel('x_9');
saveas(f5,[path_out, 'reach8vs9.pdf']);


%% Helper function
function init_set = plantReach(plant,init_set,input_set,algoC)
    nS = length(init_set);
    nL = length(input_set);
    ss = [];
    for k=1:nS
        for l=1:nL
            ss =[ss plant.stepReachStar(init_set(k), input_set(l),algoC)];
        end
    end
    init_set = ss;
end