function Param = SetUpCartPoleVerif

% Predefined Set Up for CartPole Verification model.
% With double relu neural network.
% 
% 
% SETTINGS
%
% numagents: number of agents in the system
% plantDim: dimension of the state of the physical component of the system
% geometry: class with the method used to represent the set of states
% dynamics: plant dynamics
% method: reach method used by the neural network
% reachStep: time step for the analysis of the dynamics
% controlPeriod: execution period of the controllers
% outputMat: observer matrix
% NNfiles: directory of NN path files
% NNformat: class that defining correct format for NN
% NNreachMethod: class that contains only reachability analysis method
% commands: list of possible outputs
% pre: preprocessing func (splits interval if needed)
% norm: normalization function
% scale_mean: mean values of the inputs of the NNs (for normalization purpose)
% scale_range: range values of the inputs of the NNs (for normalization
% purpose)
% lambda: function applied on the output of the NNs, for deciding which
% command applying
% stopfcn: function that returns True if we reach the stopping criteria
% isSafefcn: function that returns True if we reach a safe scenario
% plotResults: function that plots the results of the reachability analysis

if ismac
    % Code to run on Mac platform
elseif isunix
    addpath('Models/Cartpole/');
    addpath('Models/Cartpole/nnv_format');
elseif ispc
    addpath('Models\Cartpole\');
    addpath('Models\Cartpole\nnv_format');
else
    disp('Platform not supported')
end

Param.numagents = 1;
Param.plantDim = 4;
Param.geometry = @Star;
Param.dynamics = @CartPoleVerifModel;
Param.method = 'approx-star';
Param.reachStep = 10^-3;
Param.controlPeriod = 0.02;
Param.outputMat = eye(4);  % The observer --> selects the variables of interest
Param.NNfiles = {CartPoleVerifNeuralNetworks};
Param.NNformat = @FFNNS;
Param.NNreachMethod = @LayerS;
Param.commands = CartPoleCommands;
Param.pre = @CartPoleVerifPreProcessing;
Param.norm = @Normalize;
Param.scale_mean = {[-0.0157030,-0.0121955, -0.0002575, 0.0000999;
    -0.0314033,0.0074397,-0.0014457,0.0034358
    ]};
Param.scale_range = {[2.7572668, 2.8771458, 0.2425329, 0.9079171;
    2.7638733,2.8828793,0.2418418,0.9086830
    ]};
Param.lambda = {@ArgMaxVerif};
Param.stopfcn = @CartPoleVerifStoppingCriteria;
Param.isSafefcn = @CartPoleVerifSafetyness;
Param.plotResults = @CartpolePlot;

end