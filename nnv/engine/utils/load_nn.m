function [status,net_info] = load_nn(PyPath,nnmtPath,input,output,formatin)
%% This function loads a neural network and transforms it into a mat file using the nnmt tool
% % OUTPUTS
% status = 1 or 0, 1 if we run into an error, 0 if run was succesful
% net_info = object contains weights, bias matrices and other nn info
%
% % INPUTS
% PyPath = path to the python environment to use
% nnmtPath = path to the nnmt folder
% input = path to the neural network we want to transform (includes everything except for the 
% output = path to the folder we want to store the mat file at
% formatin = format of neural network to transform
%
% EXAMPLE 
%load_nn('...\envs\DL1','...\nnmt','...\nnmt\testing\neural_network_information_13','...\nnmt\examples','Sherlock');


% Check OS running on
osc = computer;
if contains(osc,'WIN')
    sh = '\'; % windows
else 
    sh = '/'; % mac and linux
end

% Generate the command to run in the terminal
commandPy = [PyPath sh 'python ' nnmtPath sh 'nnvmt.py -i ' input ' -o ' output ' -t ' formatin];
% Run the command in the terminal and convert the nn
[status,cmdout] = system(commandPy);

% Check if the nn was converted
if status == 0
    disp('Neural Netowrk converted succesfully');
else
    error(cmdout);
end

name = split(input,sh);
name = name(end);
name = [sh name{1}];
net_info = load([output name '.mat']);

end

