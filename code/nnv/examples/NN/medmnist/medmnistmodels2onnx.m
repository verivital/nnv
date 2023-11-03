% Transform all models into onnx

models = dir("models/*.mat");

if ~isfolder('onnx')
    mkdir('onnx')
end

for i=1:length(models)
    filename = string(models(i).name);
    load("models/" + filename);
    onnxfile = split(filename, '.');
    exportONNXNetwork(net,['onnx/', onnxfile{1}, '.onnx']);
end
