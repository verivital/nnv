function cnn = load_cnn(tool,Input,varargin)
%LOAD_CNN Loads a CNN from multiple formats into Matlab's default structure
%   For the moment, nnv only supports FullyConnected networks with 1-to-1
%   connections. This means, that we will only return the Layers in order,
%   without specifying the connections, as we assume that layer1 is
%   connected to layer2, layer2 to layer3,..., and layerN-1 to layerN.
%
%   The output returns a CNN object used by nnv for its robustness analysis.
%
%   Currently we support 4 different formats:
%       1) Matlab --> tool = matlab, Input = 'name_of_cnn.mat' (Defined as lgraph, DAGNetwork or SeriesNetwork)      
%       2) ONNX --> tool = onnx, Input = 'name_of_cnn.onnx', varargin = outputType
%       3) Keras --> tool = keras, Input = 'name_of_cnn.h5'
%       4) Caffe --> tool = caffe, Input = 'name_of_cnn.prototxt', varargin = 'name_of_cnn.caffemodel'
%       For more info about these tools and their integration with MATLAB, please see:
%       - help SeriesNetwork/lgraph/DAGNetwork
%       - help importONNXNetwork
%       - help importKerasNetwork
%       - help importCaffeNetwork

% Load the networks
if contains('matlab',tool)
    nn = load(Input);
elseif contains('onnx',tool)
    if length(varargin) == 1
        nn = importONNXNetwork(Input,'OutputLayerType',varargin{1});
    elseif length(varargin) == 2
        nn = importONNXNetwork(Input,'OutputLayerType',varargin{1},'Classes',varargin{2});
    else
        error('Incorrect input arguments. Please, specify the name of the cnn,'...
        +' the OutputLayerType and (optionally), add a list of output classes');
    end
elseif contains('keras',tool)
    if isempty(varargin)
        nn = importKerasNetwork(Input);
    elseif length(varargin) == 1
        nn = importKerasNetwork(Input,'Classes',varargin{1});
    else
        error('Incorrect input arguments. Please, specify the name of the cnn, and (optionally), add a list of output classes');
    end
elseif contains('caffe',tool)
    nn = importCaffeNetwork(Input,varargin{1});
else
    error('The format of the cnn does not correspond to any of the current supported formats. Please, use our transformation tool (NNVMT) available at github.com/verivital/nnvmt to convert it to one the of supported formats');
end

% Parse the networks into nnv
if string(class(nn)) == "struct"
    if isempty(nn.net.Layers(1,1).AverageImage)
        temp = nn.net.saveobj;
        layer1 = imageInputLayer([nn.net.Layers(1,1).InputSize]);
        layer1.AverageImage = zeros(layer1.InputSize);
        temp.Layers(1,1) = layer1;
%         nn.net = nn.net.loadobj(temp);
        nn = SeriesNetwork(temp.Layers);
    end
    cnn = CNN.parse(nn.net,'CNN_c');
else
    if isempty(nn.Layers(1,1).AverageImage)
        temp = nn.saveobj;
        layer1 = imageInputLayer([nn.Layers(1,1).InputSize]);
        layer1.AverageImage = zeros(layer1.InputSize);
        temp.Layers(1,1) = layer1;
%         nn = nn.loadobj(temp);
        nn = SeriesNetwork(temp.Layers);
    end
    cnn = CNN.parse(nn,'CNN_c');
end

end

