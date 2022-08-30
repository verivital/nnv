%% Reachability analysis of Single Pendulum Benchmark

%% load the controller
net = Load_nn('controller_single_pendulum.mat');

% Specify the reach step, has to be smaller than the control period
reachStep = 0.01;
%% specify the control period as specified by the benchmark description
controlPeriod = 0.05;

% define the plant as specified by nnv
plant = NonLinearODE(3,1,@dynamics_sp, reachStep, controlPeriod, eye(3));

%% Reachability analysis
% Initial set
% lb = [1.0; 0.0; 0];
lb = [1.1999; 0.1999; 0];
ub = [1.2; 0.2; 0];
init_set = Star(lb,ub);
% Input set
lb = [0];
ub = [0];
input_set = Star(lb,ub);
% Store all reachable sets
reachAll = init_set;
% reachCell = {};
% reachCell{1} = init_set;
% Execute reachabilty analysis
% for i =1:steps
num_steps = 11;
t = tic;
for i=1:num_steps
    % Compute controller output set
    inp_nn = [];
    for i = 1:length(init_set)
        inp_nn = [inp_nn init_set(i).affineMap([1 0 0;0 1 0],[])];
    end
    input_set = net.reach(inp_nn,'approx-star');
    % Compute plant reachable set
    init_set = plantReach(plant,init_set, input_set,'lin');
    reachAll = [reachAll init_set];
%     reachCell{i+1} = init_set;
end
t = toc(t);

%% Random simulations
% Initial region
lb = [1.1999; 0.1999];
ub = [1.2; 0.2];
N = 20; % number of simulations
r1 = unifrnd(lb(1),ub(1),[N 1]);
r2 = unifrnd(lb(2),ub(2),[N 1]);
r = [r1 r2 zeros(N,1)];
sims = zeros(3,N,num_steps);
for k=1:N
    rT = r(k,:); % initial state
    sims(:,k,1) = rT(end,:);
    cA = 0; % control action
    for s = 2:num_steps+1
        cA = net.evaluate(rT(end,1:2)');
        [~,rT] = plant.evaluate(rT(end,:), cA);
        sims(:,k,s) = rT(end,:);
    end
end

path_out = [path_results(), filesep, 'SinglePendulum', filesep];
mkdir(path_out)
save([path_out, 'reach.mat'],'reachAll','t','-v7.3');

%% Visualize results

f = figure;
hold on;
rectangle('Position',[0.5,1,1,1],'FaceColor',[1 0 0 0.5],'EdgeColor','r', 'LineWidth',0.1)
Star.plotBoxes_2D_noFill(plant.intermediate_reachSet,3,1,'b');
for s = 1:N
    plot(squeeze(sims(3,s,:)), squeeze(sims(1,s,:)), 'Color', [0 0 1 0.2])
end
grid;
xlabel('Time (s)');
ylabel('\Theta');
xlim([0 0.6])
ylim([0.95 1.25])
saveas(f,[path_out, 'reach.pdf']);


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