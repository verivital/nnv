 function nnv_net = parse_onnx_to_nnvnet(onnxfile)
 %% loads an onnx model and parses it to nnv network
 %  as of now works for vnncomp2021 benchmarks:
 %                                      verivital, 
 %                                      eran, 
 %                                      cifar2020(except pgd),
 %                                      oval21 
 % need to check for others
 
  lgraph = importONNXLayers(onnxfile,'OutputLayerType','classification','ImportWeights',true);  
  [placeholderLayers, indices] = findPlaceholderLayers(lgraph);
  
  if isempty(placeholderLayers) 
     matlabnet = assembleNetwork(lgraph);
  elseif ~isempty(placeholderLayers)
     numPH = length(indices)+1;
     idx = 1;
     for i=1:length(lgraph.Layers)
         if ~ismember(i,indices) && idx < numPH && i< indices(idx) 
            layers(i-idx+1,1) = lgraph.Layers(i,1);
         elseif ismember(i,indices) && idx < numPH && i == indices(idx)
             idx = idx + 1;
         elseif ~ismember(i,indices) && i> indices(numPH-1)
            layers(i-idx+1,1) = lgraph.Layers(i,1);
         end
     end
     
     
     %if class(layers(1,1)) ~= 'nnet.cnn.layer.ImageInputLayer'
         
     modif_lgraph = layerGraph(layers);
     matlabnet = assembleNetwork(modif_lgraph);
  end
  nnv_net = CNN.parse(matlabnet);
 end