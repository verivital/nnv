function [Controller] = Load_nn(varargin)
%% This function loads a neural network and transforms it into a mat file using the nnmt tool
%
% ----  OUTPUTS  ----
% net = NNV Feedforward Neural Network object contains weights, bias matrices and other nn info
%
% ---- INPUTS  ----
% input = path to the neural network we want to transform
% output (optional) = path to the folder we want to store the mat file at
% opt (optional) = corresponds to the path for the extra inputs necessary for the
%       nnmt tool such as the (optional) json file (Keras) 
%
% For more info about the nnmt tool go to:
%       https://github.com/verivital/nnmt
%
% ----  EXAMPLES  ----
% net = Load_nn('...\nnmt\testing\neural_network_information_13');
% net = Load_nn('...\model.ckpt-558138.meta');
% net = Load_nn('acasxu_1_1.mat');

if isempty(varargin)
    error('Number of inputs must be between 1 and 3');
end
% If there is a mat file, no need to parse it through nnmt  
if nargin == 1 && contains(varargin{1},'.mat')
    NNpath = varargin{1};
    net_info = load(NNpath);
    disp('Neural network loaded');
else
    % Check for python and nnmt paths to load the networks
    sh = filesep;
    [~,e,~] = pyversion;
    allpaths = path;
    allpaths = split(string(allpaths),':');
    matches = find(contains(allpaths,'nnmt'));
    nnmtPath = char(allpaths(matches(1)));
    % Default opt value
    opt = [];
    % Check the input format
    NNpath = varargin{1};
    NNname = split(string(NNpath),sh);
    NNname = NNname(end);
    NNformat = split(NNname,'.');
    if NNformat(end) == "h5" || NNformat(end) == "hdf5"
        formatin = 'keras';
    elseif NNformat(end) == "checkpoint"
        formatin = 'tensorflow';
    elseif NNformat(end) == "nnet"
        formatin = 'reluplex';
    elseif NNformat(end) == "onnx"
        formatin = 'onnx';
    elseif NNformat(end) == "txt" || length(NNformat) == 1
        formatin = 'sherlock';
    else
        error('Wrong input format');
    end
    % Set defauls output path to be current working directory
    if nargin > 1 && strcmp(formatin,'keras')
        if nargin > 2
            outPath = char(varargin{2});
            opt = char(varargin{3});
        else
            temp = split(string(varargin{2}),sh);
            if contains(temp,'json')
                opt = char(varargin{2});
                outPath = pwd;
            else
                outPath = char(varargin{2});
            end
        end
    elseif nargin > 1
        outPath = char(varargin{2});
    else
        outPath = pwd;
    end
    input = char(NNpath);
    if ~isempty(opt)
        % Generate the command to run in the terminal
        commandPy = strcat(e ,[' ' nnmtPath],sh,'nnvmt.py -i ',[' ' input],' -o ',[' ' outPath],' -t ',[' ' formatin],' -j ',[' ' opt]);
    else
        commandPy = strcat(e ,[' ' nnmtPath],sh,'nnvmt.py -i ',[' ' input],' -o ',[' ' outPath],' -t ',[' ' formatin]);
    end
    % Run the command in the terminal and convert the nn
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
    net_info = load(strcat(outPath,sh,name,'.mat'));
    disp('Neural network loaded');
end

% Create the NN object for nnv
W = net_info.W; % weights
b = net_info.b; % bias
n = length(b); % number of layers
acf = net_info.act_fcns;

Layers = [];
for i=1:n 
    if W{i} == 0
        Layers(end).f = ActFunction(acf(i,:));
    else
        L = LayerS(W{i}, b{i}, ActFunction(acf(i,:)));
        Layers = [Layers L];
    end
end


Controller = FFNNS(Layers); % feedforward neural network controller
disp('NN Controller created');

end

