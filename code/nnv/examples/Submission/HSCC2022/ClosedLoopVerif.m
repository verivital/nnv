function [allReach,verifResult,time_elapsed,time_dynamics_eval,time_nn_eval,max_n_branches] = ClosedLoopVerif(Param, init_set, init_labels, tf)

% Reachability analysis of the closedLoop system of the selected model
% 
% INPUTS
%
% Param: struct containing the information needed to execute the simulation.
% init_set: set of the possible initial states of the pysical component of the system. 
% init_labels: the initial labels, corresponding to the initial commands of
% the agents composing the system.
% tf: the time horizon for the simulation.
%
% OUTPUTS
%
% allReach: struct containing the results of the simulation.
% verifResult: a boolean indicating whether the property of interest holds
% or not.
% time_elaplsed: the time necessary for performing the reachability
% analysis.
% time_dynamics_eval: the time necessary for analyzing the physical
% components of the system.
% time_nn_eval: the time necessary for analyzing the neural networks.
% max_n_branches: the maximum number of branches created for analyzing the
% system

%% Measure time elapsed
t_start = tic;

%% Load components
numagents = Param.numagents;
reachStep = Param.reachStep;
controlPeriod = Param.controlPeriod;
outputMat = Param.outputMat; 
dynamics = Param.dynamics; 
Idx = init_labels; % scalar or array
reachMethod = Param.method;
plantdim = Param.plantDim; 
NNfiles = Param.NNfiles; % cell
NNformat = Param.NNformat;
NNreachMethod = Param.NNreachMethod;
listofcommands = Param.commands; 
preprocessing = Param.pre; 
mean = Param.scale_mean; % cell
range = Param.scale_range; % cell
lambda = Param.lambda; % cell
stopfcn = Param.stopfcn;
isSafefcn = Param.isSafefcn;

numinputs = numagents;

% Plant in Correct Format
plant = NonLinearODE(plantdim,numinputs,dynamics,reachStep,controlPeriod,outputMat);

for nn = 1:numagents
    
% NN in Correct Format 
NN{nn} = CreateNNStruct(NNfiles{nn},NNformat,NNreachMethod);

end

%% Reachability analysis (Multiple steps)
times = 0:controlPeriod:tf;

% Initialize variables
allReach.init_set = cell(1,length(times));
allReach.Ro = cell(numagents,length(times));
allReach.Unn = cell(numagents,length(times));
allReach.yNN = cell(numagents,length(times));
allReach.minIdx = cell(numagents,length(times));
allReach.Up = cell(numagents,length(times));
allReach.combos = cell(1,length(times));
allReach.init_set{1,1} = init_set;
verifResult = 'safe';
time_dynamics_eval = 0.0;
time_nn_eval = 0.0;
max_n_branches = 1;

cmbs = {Idx,1,Idx,1};  
step_sets = [init_set];
Ro = [];


% cmbs EXPLANATION 
%
% cmbs = [ Prev advisory, output set, current advisory, state set]
% cmbs(1): previous advisory. Used to select the current NN. Scalar or
%          vector (depending on agent number)
% cmbs(2): number of Ro in which current init_set is divided. Such
%          divisions might happen at PreProcessing or at Lambda.
%          At the end of one time step we refresh this column
%          with the value of column (4).         
% cmbs(3): current advisory. Used in Advisory to select what command send
%          to the plant. Scalar or vector (depending on agent number).
% cmbs(4): state set. Number of current init_set. The update of this 
%          value is done in PlantReach. 

stop = false;

% Start reachability loop
for j = 1:length(times)-1
        
    for i = 1:numagents
        
        % Adequates the values of combos when we have numagents > 1
        cmbs = preSetCombos(cmbs);
        
        % Output set --> adequates init_set to enter the network
        [Ro,cmbs] = preprocessing(init_set,outputMat,cmbs,i); % Ro_0_i = [Star(init_set_0)]
        
        % Normalize inputs
        Unn = Normalize(Ro,mean{i},range{i},cmbs,i);  % Unn_0_i = [Star(Ro_0_i)]
        
        % Compute NN outputs
        t_start_nn = tic; % measure time elapsed
        yNN = ReachNN(Unn,NN{i},reachMethod,cmbs,i); % yNN_0_i = {Star(Unn_0_i)}
        t_nn = toc(t_start_nn);
        
        % Compute advisory command
        [Idx,cmbs] = lambda{i}(yNN,cmbs,i); % Idx = {index or array of indexes of the agent i per branch}
        Up = PostProcessing(Idx,listofcommands); % Up = {chosen command per branch}
        
        % End cycle
        allReach.Ro{i,j} = Ro;
        allReach.Unn{i,j} = Unn;
        allReach.yNN{i,j} = yNN;
        allReach.minIdx{i,j} = Idx;
        allReach.Up{i,j} = Up;
        
        % update time_nn_eval variable
        time_nn_eval = time_nn_eval + t_nn;
       
    end
        
    allReach.combos{1,j} = cmbs; 
    
    % Reachability step plant
    t_start_dynamics = tic; % measure time elapsed
    [init_set,cmbs] = PlantReach(plant,init_set, cmbs, listofcommands); % init_set_1 = Star(init_set_0, commands (eq. to Up but with another srtuct)))
    t_dynamics = toc(t_start_dynamics);
    
    % End cycle
    allReach.init_set{1,j+1} = init_set;
    step_sets = [step_sets init_set];
    
    % update time_dynamics_eval variable
    time_dynamics_eval = time_dynamics_eval + t_dynamics;
    
    % update the max_n_branches variable
    if size(cmbs,1) > max_n_branches
        max_n_branches = size(cmbs,1);
    end
    
    % Checks the stopping criteria
    isSafe = isSafefcn(init_set); 
    stop = stopfcn(isSafe,init_set);
    allReach.isSafe{1,j} = isSafe;
    
    if ~isSafe
        verifResult = 'unknown';
    end
    
    if stop
        
        allReach.step_sets = step_sets;
        allReach.int_reachSet = plant.intermediate_reachSet;
        time_elapsed = toc(t_start);
        return;
    end
            
end
 
allReach.step_sets = step_sets;
allReach.int_reachSet = plant.intermediate_reachSet;

time_elapsed = toc(t_start);

end