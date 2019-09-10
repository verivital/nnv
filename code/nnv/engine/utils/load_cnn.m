function cnn = load_cnn(Input,tool)
%LOAD_CNN Summary of this function goes here
%   Detailed explanation goes here
if contains('atlab',tool)
    nn = load(Input);
elseif contains('nnx',tool)
    nn = importONNXNetwork(Input);
else
    error('Need to convert the cnn using the nnvmt tool first);
end

cnn = CNN.parse(nn,'CNN_c');


end

