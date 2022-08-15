%% Bug found during the reachability analysis of TORA (benchmark 9)
% Load components and set reachability parameters
net = Load_nn('controllerTora_nnv.mat');
reachStep = 0.05;
controlPeriod = 1;
plant = NonLinearODE(4,1,@dynamics9, reachStep, controlPeriod, eye(4));
time = 0:controlPeriod:20;
steps = length(time);
% Initial set
lb = [0.69; -0.7; -0.4; 0.5];
ub = [0.7; -0.69; -0.39; 0.51];
% init_set = Box(lb,ub);
offset = 10;
scale_factor = 1;
% plant.set_tensorOrder(3);
plant.options.intermediateOrder = 20;
plant.options.zonotopeOrder = 20;
plant.options.taylorTerms  = 8;
plant.set_maxError(0.005*ones(4,1))

%% Control step 1 - Fine
init_set = Star(lb,ub);
% Input set
lb = 0;
ub = 0;
input_set = Star(lb,ub);
% Store all reachable sets
reachLin = init_set;

% Execute reachabilty analysis
% Compute plant reachable set
init_set = plantReach(plant,init_set,input_set,'lin');
reachLin = [reachLin init_set];

% Compute controller output set
% Exact
inp_set_exact = net.reach(init_set,'exact-star');
input_set_exact = [];
for k = 1:length(inp_set_exact)
    input_set_exact = [input_set_exact inp_set_exact(k).affineMap(1,-offset)];
end
[me1,Me1] = inp_set_exact.getRanges()

% Approx
inp_set_approx= net.reach(init_set,'approx-star');
input_set_approx = [];
for k = 1:length(inp_set_approx)
    input_set_approx = [input_set_approx inp_set_approx(k).affineMap(1,-offset)];
end
[ma1,Ma1] = inp_set_approx.getRanges()

%% Control Step 2 - Bug 
% But taking a closer look at the second step, there is a mismatch on the output set of the controller
% Exact and approach seem to be the same, so we can use either for the next
% step
% Interval range of exact star > approx, while it should be the other way around. Exact star must be contained in approx star

% Execute reachabilty analysis
% Compute plant reachable set
init_set = plantReach(plant,init_set,input_set_approx,'lin');
reachLin = [reachLin init_set];
% Compute controller output set
% Exact
inp_set_exact = net.reach(init_set,'exact-star');
input_set_exact = [];
for k = 1:length(inp_set_exact)
    input_set_exact = [input_set_exact inp_set_exact(k).affineMap(1,-offset)];
end
me = zeros(length(inp_set_exact),1);
Me = zeros(length(inp_set_exact),1);
for s = 1:length(inp_set_exact)
    [me(s),Me(s)] = inp_set_exact(s).getRanges();
end
min(me)
max(Me)
% Approx
inp_set_approx= net.reach(init_set,'approx-star');
input_set_approx = [];
for k = 1:length(inp_set_approx)
    input_set_approx = [input_set_approx inp_set_approx(k).affineMap(1,-offset)];
end
[ma,Ma] = inp_set_approx.getRanges()

% If plotted, it seems like exact > approx, but when looking at the actual
% ranges, it's fine. We should probably disable the exact plots.
figure;
Star.plots(inp_set_approx)
hold on
Star.plots(inp_set_exact)
% vs
figure;
Star.plotBoxes_2D_noFill(inp_set_approx,1,1,'r');
hold on;
Star.plotBoxes_2D_noFill(inp_set_exact,1,1,'b');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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