function [Controller] = load_nn(PyPath,nnmtPath,input,output,formatin,opt)
%% This function loads a neural network and transforms it into a mat file using the nnmt tool
% % OUTPUTS
% status = 1 or 0, 1 if we run into an error, 0 if run was succesful
% net_info = object contains weights, bias matrices and other nn info
%
% % INPUTS
% PyPath = path to the python environment to use
% nnmtPath = path to the nnmt folder
% input = path to the neural network we want to transform 
% output = path to the folder we want to store the mat file at
% formatin = format of neural network to transform
% opt = corresponds to the path for the extra inputs necessary for the
%       nnmt tool such as the (optional) json file (Keras), the path to
%       checkpoint file for tensorflow, or if we have a mat file
%       already, we can input it here and load it without using nnmt
% For more info about the nnmt tool go to:
%       https://github.com/verivital/nnmt
%
% EXAMPLE 
%load_nn('...\envs\DL1','...\nnmt','...\nnmt\testing\neural_network_information_13','...\nnmt\examples','Sherlock',[]);
%cont = load_nn('...\envs\DL1','...\nnmt','...\model.ckpt-558138.meta','...\MATLAB','Tensorflow','...\UUV-LEC1-558138')

% If the is a mat file, no need to parse it through nnmt  
if ~isempty(opt) && contains(opt,'.mat')
    net_info = load(opt);
    disp('Neural network loaded');
else
    % Check OS running on
    osc = computer;
    if contains(osc,'WIN')
        sh = '\'; % windows
    else 
        sh = '/'; % mac and linux
    end
    if ~isempty(opt)
        % Generate the command to run in the terminal
        commandPy = strcat(PyPath,sh,'python ',[' ' nnmtPath],sh,'nnvmt.py -i ',[' ' input],' -o ',[' ' output],' -t ',[' ' formatin],' -j ',[' ' opt]);
    else
        commandPy = strcat(PyPath,sh,'python ',[' ' nnmtPath],sh,'nnvmt.py -i ',[' ' input],' -o ',[' ' output],' -t ',[' ' formatin]);
    end
    % Run the command in the terminal and convert the nn
    %disp(commandPy);
    [status,cmdout] = system(commandPy);

    % Check if the nn was converted
    if status == 0
        disp('Neural Network converted succesfully');
    else
        error(cmdout);
    end
    if contains(formatin,'ensorflow')
        a = importdata(input);
        a = a(1);
        a = a{1};
        a = split(string(a),'"');
        name = a(2);
    else
        name = split(string(input),sh);
        name = name(end);
        name = split(name,'.');
        if length(name) ~= 1
            name = name(1:end-1);
        end
        name = strjoin(name,'.');
    end
    %ouch = strcat(output,sh,name,'.mat');
    net_info = load(strcat(output,sh,name,'.mat'));
    disp('Neural network loaded');
end

% Create the NN object for nnv
W = net_info.W; % weights
b = net_info.b; % bias
n = length(b); % number of layers
acf = net_info.act_fcns;

Layers = [];
for i=1:n 
    L = LayerS(W{i}, b{i}, ActFunction(acf(i,:)));
    Layers = [Layers L];
end


Controller = FFNNS(Layers); % feedforward neural network controller
disp('NN Controller created');

end

